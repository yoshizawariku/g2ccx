"""CalculiX input generation, supervised execution and history reduction."""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

import numpy as np
from fea_geometry import CORNERS, FEATURE_SETS
from g2ccx_paths import resolve_solver



def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def write_deck(grid, config, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    diagnostic = grid.diagnostics()
    save_json(output/"geometry.json", diagnostic)
    if diagnostic["face_components"] != 1:
        raise ValueError(f"Mesh has {diagnostic['face_components']} face-connected components; inspect geometry.json. No auto-repair performed.")
    nodes, elements, cells, masks = grid.mesh(config["max_elements"])
    lo, hi = nodes.min(axis=0), nodes.max(axis=0)
    if not 0 < config["displacement_mm"] < hi[2]-lo[2]:
        raise ValueError("displacement_mm must be positive and smaller than the ROI thickness")
    ne, nn = len(elements), len(nodes)
    # Two finite-thickness rigid bricks. Independent control nodes are not body members.
    margin = max(2.0, (hi-lo)[:2].max()*0.1)
    plate_nodes = []
    for z in (lo[2]-1, hi[2]):
        p0 = np.array([lo[0]-margin, lo[1]-margin, z])
        size = np.array([hi[0]-lo[0]+2*margin, hi[1]-lo[1]+2*margin, 1.0])
        plate_nodes.extend(p0+CORNERS*size)
    refs = {"FLOOR_REF": nn+17, "FLOOR_ROT": nn+18,
            "LOAD_REF": nn+19, "LOAD_ROT": nn+20}
    bottom_nodes = np.flatnonzero(np.isclose(nodes[:,2],lo[2],atol=1e-8))+1
    centre = (lo+hi)/2
    anchor = int(bottom_nodes[np.argmin(np.linalg.norm(nodes[bottom_nodes-1,:2]-centre[:2],axis=1))])
    same_row = bottom_nodes[np.isclose(nodes[bottom_nodes-1,1], nodes[anchor-1,1])]
    rotation_anchor = int(same_row[np.argmax(np.abs(nodes[same_row-1,0]-nodes[anchor-1,0]))])
    if rotation_anchor == anchor:
        raise ValueError("ROI too narrow for three independent in-plane stabilization constraints")
    pad = np.pad(grid.masks != 0, 1)
    lower = ~pad[cells[:,0]+1, cells[:,1]+1, cells[:,2]]
    upper = ~pad[cells[:,0]+1, cells[:,1]+1, cells[:,2]+2]
    feature_sets = {FEATURE_SETS[name]: (np.flatnonzero(masks & (1 << i))+1).tolist()
                    for i,name in enumerate(grid.features)}
    path = output/"model.inp"
    with path.open("w", encoding="ascii", newline="\n") as f:
        def line(s):
            f.write(s+"\n")
        def ids(seq):
            for i in range(0,len(seq),16):
                line(", ".join(str(int(n)) for n in seq[i:i+16]))
        line("*HEADING\nG-code voxel TPU compression; N mm MPa")
        line("*NODE")
        for i,p in enumerate(nodes,1):
            line(f"{i}, {p[0]:.10g}, {p[1]:.10g}, {p[2]:.10g}")
        for i,p in enumerate(plate_nodes,nn+1):
            line(f"{i}, {p[0]:.10g}, {p[1]:.10g}, {p[2]:.10g}")
        for name,n in refs.items():
            z = hi[2] if name.startswith("LOAD") else lo[2]
            line(f"{n}, {centre[0]:.10g}, {centre[1]:.10g}, {z:.10g}")
        line(f"*ELEMENT, TYPE={config['element_type']}, ELSET=TPU")
        for i,e in enumerate(elements,1):
            line(f"{i}, "+", ".join(map(str,e)))
        for k,name in enumerate(("FLOOR", "PLATEN")):
            line(f"*ELEMENT, TYPE=C3D8, ELSET={name}")
            line(f"{ne+k+1}, "+", ".join(map(str,range(nn+8*k+1,nn+8*k+9))))
            line(f"*NSET, NSET={name}_BODY")
            ids(list(range(nn+8*k+1,nn+8*k+9)))
        for name,seq in feature_sets.items():
            if seq:
                line(f"*ELSET, ELSET={name}")
                ids(seq)
        for name,n in refs.items():
            line(f"*NSET, NSET={name}\n{n}")
        line(f"*NSET, NSET=STABILIZATION\n{anchor}, {rotation_anchor}")
        line("*NSET, NSET=TPU_NODES, GENERATE")
        line(f"1, {nn}, 1")
        for name,selector,face in (("TPU_DOWN",lower,"S1"),("TPU_UP",upper,"S2")):
            line(f"*SURFACE, NAME={name}, TYPE=ELEMENT")
            for eid in np.flatnonzero(selector)+1:
                line(f"{eid}, {face}")
        line(f"*SURFACE, NAME=FLOOR_TOP, TYPE=ELEMENT\n{ne+1}, S2")
        line(f"*SURFACE, NAME=PLATEN_BOTTOM, TYPE=ELEMENT\n{ne+2}, S1")
        line("*MATERIAL, NAME=TPU90A\n*HYPERELASTIC, MOONEY-RIVLIN")
        m = config["material"]
        line(f"{m['c10']}, {m['c01']}, {m['d1']}")
        line("*SOLID SECTION, ELSET=TPU, MATERIAL=TPU90A")
        line("*MATERIAL, NAME=PLATE_MAT\n*ELASTIC\n210000., 0.3")
        for name in ("FLOOR", "PLATEN"):
            line(f"*SOLID SECTION, ELSET={name}, MATERIAL=PLATE_MAT")
        line(f"*RIGID BODY, NSET=FLOOR_BODY, REF NODE={refs['FLOOR_REF']}, ROT NODE={refs['FLOOR_ROT']}")
        line(f"*RIGID BODY, NSET=PLATEN_BODY, REF NODE={refs['LOAD_REF']}, ROT NODE={refs['LOAD_ROT']}")
        line("*SURFACE INTERACTION, NAME=CONTACT\n*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=LINEAR")
        line(f"{config['contact_stiffness_mpa_per_mm']}, 1.e-6")
        if config['friction'] > 0:
            line(f"*FRICTION\n{config['friction']}, {config['contact_stiffness_mpa_per_mm']/10}")
        line("*CONTACT PAIR, INTERACTION=CONTACT, TYPE=SURFACE TO SURFACE\nTPU_DOWN, FLOOR_TOP")
        line("*CONTACT PAIR, INTERACTION=CONTACT, TYPE=SURFACE TO SURFACE\nTPU_UP, PLATEN_BOTTOM")
        line("*BOUNDARY\nFLOOR_REF, 1, 3\nFLOOR_ROT, 1, 3\nLOAD_REF, 1, 2\nLOAD_ROT, 1, 3")
        line(f"{anchor}, 1, 2\n{rotation_anchor}, 2, 2")
        line(f"*STEP, NLGEOM, INC=1000\n*STATIC, SOLVER={config['linear_solver']}")
        line(f"{config['initial_increment']}, 1., {config['minimum_increment']}, {config['maximum_increment']}")
        line(f"*BOUNDARY\nLOAD_REF, 3, 3, {-config['displacement_mm']}")
        freq = config['output_frequency']
        line(f"*NODE FILE, FREQUENCY={freq}\nU, RF\n*EL FILE, FREQUENCY={freq}\nS, E")
        line(f"*CONTACT FILE, FREQUENCY={freq}\nCDIS, CSTR")
        for name in ("LOAD_REF", "FLOOR_REF", "STABILIZATION"):
            line(f"*NODE PRINT, NSET={name}, FREQUENCY={freq}\nU, RF")
        line(f"*EL PRINT, ELSET=TPU, FREQUENCY={freq}\nS, E")
        line("*END STEP")
    diagnostic.update({"nodes": nn, "references": refs, "stabilization_nodes": [anchor,rotation_anchor],
                       "bounds_mm": [lo.tolist(),hi.tolist()], "feature_element_counts": {k:len(v) for k,v in feature_sets.items()}})
    save_json(output/"geometry.json", diagnostic)
    # Compact membership arrays also support per-feature integration-point summaries.
    np.savez_compressed(output/"features.npz", masks=masks, names=np.array([FEATURE_SETS[x] for x in grid.features]))
    save_json(output/"config.resolved.json", config)
    save_json(output/"build.json", {"input_sha256": digest(config["input"]) if config.get("input") else None,
                                    "deck_sha256": digest(path), "config": config})
    return diagnostic


def solve(config, output):
    output = Path(output).resolve()
    exe = Path(resolve_solver(config.get("solver")))
    geometry = json.loads((output/"geometry.json").read_text(encoding="utf-8"))
    dofs = 3*geometry["nodes"]
    if dofs > config.get("max_solver_dofs",500_000):
        raise ValueError(f"{dofs:,} approximate TPU DOFs exceed max_solver_dofs={config.get('max_solver_dofs',500_000):,}. Reduce ROI or explicitly raise the limit on a suitably sized machine; sparse factorization memory is not predictable from voxel count alone.")
    build = json.loads((output/"build.json").read_text(encoding="utf-8"))
    if digest(output/"model.inp") != build["deck_sha256"]:
        raise ValueError("model.inp changed after build; rebuild to keep results traceable")
    # Refuse stale-result reuse; a new run directory is required.
    if any((output/f"model.{ext}").exists() for ext in ("frd", "dat", "sta")):
        raise FileExistsError("Solver results already exist; choose a new output directory")
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(config["threads"])
    start = time.monotonic()
    status = {"status": "running", "solver": str(exe), "deck_sha256": build["deck_sha256"],
              "runtime":{k:config[k] for k in ('threads','timeout_seconds','max_solver_dofs')}}
    save_json(output/"run.json", status)
    try:
        with (output/"solver.log").open("w", encoding="utf-8") as log:
            proc = subprocess.Popen([str(exe), "-i", "model"], cwd=output, env=env,
                                    stdout=log, stderr=subprocess.STDOUT,
                                    creationflags=getattr(subprocess,"CREATE_NO_WINDOW",0))
            status['pid'] = proc.pid
            save_json(output/'run.json',status)
            try:
                deadline = start+config["timeout_seconds"]
                next_report = start+30
                while proc.poll() is None:
                    now = time.monotonic()
                    if now >= deadline:
                        raise subprocess.TimeoutExpired(proc.args,config["timeout_seconds"])
                    if now >= next_report:
                        print(f"CalculiX running: {now-start:.0f} s; progress: {output/'model.sta'}",flush=True)
                        with (output/'solver.log').open('rb') as current:
                            current.seek(max(0,current.seek(0,2)-262144))
                            tail = current.read().decode('utf-8',errors='replace')
                        if re.search(r'\bnan\b|\binfinity\b',tail,re.I):
                            status['numerical_failure'] = 'Non-finite solver state detected'
                            proc.kill()
                            proc.wait()
                            break
                        next_report = now+30
                    time.sleep(0.5)
                code = proc.returncode
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait()
                code = proc.returncode
                status['timed_out'] = True
            except KeyboardInterrupt:
                proc.kill()
                proc.wait()
                raise
        log_text = (output/"solver.log").read_text(encoding="utf-8",errors="replace")
        version = re.search(r'CalculiX Version (\S+)',log_text)
        status.update(returncode=code, solver_finished="Job finished" in log_text,
                      solver_version=version[1].rstrip(',') if version else None,
                      solver_errors=bool(re.search(r'\bnan\b|\binfinity\b',log_text,re.I)) or any(x in log_text.upper() for x in ('*ERROR','COULD NOT CONVERGE')))
        status["status"] = "solver_finished" if code == 0 and status["solver_finished"] and not status["solver_errors"] else "failed"
    except BaseException as exc:
        status.update(status="failed", error=str(exc) or type(exc).__name__)
        raise
    finally:
        status["elapsed_seconds"] = time.monotonic()-start
        save_json(output/"run.json",status)
    return status


NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
HEADER = re.compile(r"^\s*(displacements|forces|stresses|strains)\s+.*?for set\s+(\S+)\s+and time\s+("+NUMBER+r")", re.I)


def dat_blocks(path):
    """Stream one output block at a time; avoid retaining the full history field."""
    header, rows = None, []
    with Path(path).open(encoding="utf-8", errors="replace") as f:
        for line in f:
            match = HEADER.match(line)
            if match:
                if header:
                    yield (*header, np.asarray(rows,dtype=float))
                header = (match[1].lower(),match[2],float(match[3].replace("D","E")))
                rows = []
            elif header:
                tokens = line.split()
                if tokens and all(re.fullmatch(NUMBER,t) for t in tokens):
                    rows.append([float(t.replace("D","E")) for t in tokens])
        if header:
            yield (*header,np.asarray(rows,dtype=float))


def tensor_stats(values):
    # CalculiX DAT order: xx yy zz xy xz yz (tensor, not engineering shear).
    a = np.zeros((len(values),3,3))
    a[:,0,0],a[:,1,1],a[:,2,2] = values[:,0],values[:,1],values[:,2]
    a[:,0,1]=a[:,1,0]=values[:,3]
    a[:,0,2]=a[:,2,0]=values[:,4]
    a[:,1,2]=a[:,2,1]=values[:,5]
    eigen = np.linalg.eigvalsh(a)
    vm = np.sqrt(((eigen[:,0]-eigen[:,1])**2+(eigen[:,1]-eigen[:,2])**2+(eigen[:,2]-eigen[:,0])**2)/2)
    return {"principal_min": float(eigen.min()), "principal_max": float(eigen.max()), "von_mises_measure_max": float(vm.max())}


def postprocess(output):
    output = Path(output)
    config = json.loads((output/"config.resolved.json").read_text(encoding="utf-8"))
    run = json.loads((output/"run.json").read_text(encoding="utf-8"))
    if digest(output/'model.inp') != run['deck_sha256']:
        raise ValueError('Deck differs from the executed input; refusing to associate stale results')
    membership = np.load(output/"features.npz")
    masks,names = membership["masks"],membership["names"]
    history, summaries = {}, []
    for kind, name, t, rows in dat_blocks(output/"model.dat"):
        if rows.size == 0:
            continue
        if name in {"LOAD_REF", "FLOOR_REF", "STABILIZATION"}:
            history.setdefault(t,{})[(name,kind)] = rows[:,1:4].sum(axis=0)
        elif name == "TPU" and kind in {"stresses", "strains"}:
            eid = rows[:,0].astype(int)-1
            for i, feature in enumerate(names):
                selection = (masks[eid] & (1 << i)) != 0
                if selection.any():
                    stats = tensor_stats(rows[selection,2:8])
                    measure = stats.pop('von_mises_measure_max')
                    stats['stress_mises_max_mpa'] = measure if kind == 'stresses' else ''
                    stats['strain_deviatoric_equivalent_max'] = measure*2/3 if kind == 'strains' else ''
                    extra = {'stretch_min':'','stretch_max':'','log_strain_min':'','log_strain_max':''}
                    if kind == 'strains' and 1+2*stats['principal_min'] > 0:
                        for suffix in ('min','max'):
                            squared = 1+2*stats['principal_'+suffix]
                            extra['stretch_'+suffix] = float(np.sqrt(squared))
                            extra['log_strain_'+suffix] = float(np.log(squared)/2)
                    summaries.append({"time":t,"feature":str(feature),"field":kind,**stats,**extra})
    records = []
    for t,data in sorted(history.items()):
        required = [("LOAD_REF","displacements"),("LOAD_REF","forces"),("FLOOR_REF","forces")]
        if not all(k in data for k in required):
            continue
        u,load,floor = (data[k] for k in required)
        mismatch = abs(load[2]+floor[2])/max(abs(load[2]),abs(floor[2]),1e-12)
        stab = data.get(("STABILIZATION","forces"), np.zeros(3))
        records.append({"time":t,"indentation_mm":-u[2],"compression_force_n":-load[2],
                        "floor_force_n":floor[2],"balance_relative":mismatch,
                        "stabilization_force_n":float(np.linalg.norm(stab))})
    if records:
        with (output/"force_displacement.csv").open("w",newline="",encoding="utf-8") as f:
            writer = csv.DictWriter(f,fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot([0]+[r["indentation_mm"] for r in records],[0]+[r["compression_force_n"] for r in records],".-")
        ax.set(xlabel="Indentation [mm]", ylabel="Compression force [N]")
        ax.grid(True)
        fig.tight_layout()
        fig.savefig(output/"force_displacement.png",dpi=150)
        plt.close(fig)
    if summaries:
        with (output/"feature_fields.csv").open("w",newline="",encoding="utf-8") as f:
            writer = csv.DictWriter(f,fieldnames=list(summaries[0]))
            writer.writeheader()
            writer.writerows(summaries)
    last = records[-1] if records else {}
    complete = (run.get("solver_finished",False) and not run.get("solver_errors",True)
                and run.get("returncode") == 0 and bool(records)
                and abs(last["indentation_mm"]-config["displacement_mm"]) <= max(1e-6,config["displacement_mm"]*1e-5)
                and abs(last["time"]-1) < 1e-6)
    result = {"status":"complete" if complete else "incomplete", "saved_increments":len(records), "last":last,
              "force_balance_pass": bool(records) and bool(max(r["balance_relative"] for r in records) < 0.01),
              "strain_measure":"CalculiX E: Lagrangian strain, tensor shear; not logarithmic strain",
              "material_note":"Surrogate TPU90A. D1=0 invokes CalculiX default compressibility."}
    save_json(output/"summary.json",result)
    run["status"] = result["status"]
    save_json(output/"run.json",run)
    return result
