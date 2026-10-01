"""Whole-insole homogenised model: tile calibration, column mesh, gait load cases, CalculiX decks."""
from __future__ import annotations

import csv
import json
import re
import time
from pathlib import Path

import numpy as np

import homog
from fea_solver import digest, save_json, solve
from gait_load import pressure_field
from insole_loads import load_setup

DEFAULTS = {
    "gcode": None, "sensors": None, "gait": None, "side": "left", "body_mass_kg": 53.0,
    "dx_mm": 1.2, "layers": 6, "classes": 5, "tile_mm": 10.0, "tile_voxel_mm": 0.2,
    "tile_strain": 0.25, "tile_element_type": "C3D8", "tile_timeout_seconds": 10800,
    "calibration_max_stress_mpa": 0.25, "sensor_rotation_deg": 0.0, "medial_sign": None, "peak_factor": 1.1, "standing_factor": 0.5,
    "threads": 16, "tile_threads": 6, "timeout_seconds": 7200, "max_solver_dofs": 1_500_000,
    "linear_solver": "PARDISO", "initial_increment": 0.1, "maximum_increment": 0.25, "minimum_increment": 1e-5,
    "solver": None,
}


def load_config(path, create_output=False):
    from g2ccx_paths import resolve_output
    path = Path(path).resolve()
    raw = json.loads(path.read_text(encoding="utf-8"))
    unknown = set(raw) - set(DEFAULTS) - {"output"}
    if unknown:
        raise ValueError(f"Unknown config keys: {sorted(unknown)}")
    c = {**DEFAULTS, **raw}
    for key in ("gcode", "sensors", "gait", "output", "solver"):
        if key == "solver" and not c.get(key):
            continue  # looked up when the solver starts
        if not c.get(key):
            raise ValueError(f"Config requires {key}")
        p = Path(c[key])
        c[key] = str(p if p.is_absolute() else (path.parent / p).resolve())
    c["output"] = resolve_output(c["output"], create_output)
    return c


def prepare(config):
    """Column table, classes and tile input files. Returns the plan (also saved as plan.json)."""
    out = Path(config["output"])
    table = homog.column_table(config["gcode"], config["dx_mm"], out)
    labels, centres = homog.classify(table, config["classes"])
    tiles = homog.representative_tiles(table, labels, centres, config["tile_mm"])
    from bambu_gcode3mf_to_fea_stl import DEFAULT_FEATURES, parse_segments, read_gcode
    from fea_geometry import voxelize
    segments, _ = parse_segments(read_gcode(Path(config["gcode"])), DEFAULT_FEATURES)
    plan = {"classes": [], "dx_mm": config["dx_mm"]}
    (out / "tiles").mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "labels.npz", labels=labels)
    def single_piece(class_height, tile):
        # The solver refuses disconnected meshes, so screen candidates the same way it will.
        cx, cy = tile["centre_mm"]
        h = tile["tile_mm"] / 2
        crop = [cx - h, cx + h, cy - h, cy + h, 0, round(float(np.ceil((class_height + 0.6) * 5) / 5), 1)]
        try:
            return voxelize(segments, config["tile_voxel_mm"], crop, 50_000_000).diagnostics()["face_components"] == 1
        except ValueError:
            return False

    for c, (centre, candidates) in enumerate(zip(centres, tiles)):
        tile = next((t for t in candidates if single_piece(float(centre[0]), t)), None)
        entry = {"class": c, "columns": int((labels == c).sum()), "height_mm": float(centre[0]),
                 "phi": float(centre[1]), "phi_bottom": float(centre[2]), "phi_top": float(centre[3]), "tile": tile}
        if tile:
            cx, cy = tile["centre_mm"]
            h = tile["tile_mm"] / 2
            zmax = float(centre[0]) + 0.6
            name = f"tile_c{c}"
            save_json(out / "tiles" / f"{name}.json", {
                "input": str(Path(config["gcode"])), "output": name,
                "crop_mm": [cx - h, cx + h, cy - h, cy + h, 0, round(float(np.ceil(zmax * 5) / 5), 1)],
                "voxel_mm": config["tile_voxel_mm"], "element_type": config["tile_element_type"],
                "displacement_mm": round(config["tile_strain"] * float(centre[0]), 3),
                "initial_increment": 0.01, "maximum_increment": 0.04, "output_frequency": 1,
                "timeout_seconds": config["tile_timeout_seconds"],
                **({"solver": config["solver"]} if config["solver"] else {}),
                "threads": config["tile_threads"]})
            entry["config"] = str(out / "tiles" / f"{name}.json")
        plan["classes"].append(entry)
    save_json(out / "plan.json", plan)
    return plan


def run_tiles(config):
    """Run all tile compressions in parallel (each a normal bambu_fea run)."""
    import subprocess
    import sys
    out = Path(config["output"])
    plan = json.loads((out / "plan.json").read_text(encoding="utf-8"))
    script = Path(__file__).with_name("bambu_fea.py")
    procs = []
    import shutil
    for entry in plan["classes"]:
        if "config" in entry:
            cfg_path = out / "tiles" / Path(entry["config"]).name
            run = out / "tiles" / json.loads(cfg_path.read_text(encoding="utf-8"))["output"]
            done = run / "summary.json"
            if done.exists() and json.loads(done.read_text(encoding="utf-8")).get("status") == "complete":
                print(f"tile class {entry['class']}: already complete, skipped", flush=True)
                continue
            shutil.rmtree(run, ignore_errors=True)  # discard a partial run; the solver refuses stale results
            log = open(cfg_path.with_suffix(".log"), "w")
            if procs:
                time.sleep(60)  # stagger: each run parses the whole G-code at start-up
            procs.append((entry["class"], subprocess.Popen([sys.executable, str(script), "run", "--config", str(cfg_path)],
                                                            stdout=log, stderr=subprocess.STDOUT)))
    for c, proc in procs:
        print(f"tile class {c}: exit code {proc.wait()}", flush=True)


def calibrate(config):
    """Fit compressible Ogden parameters per class from finished tile runs."""
    out = Path(config["output"])
    plan = json.loads((out / "plan.json").read_text(encoding="utf-8"))
    result = {}
    for entry in plan["classes"]:
        if "config" not in entry:
            continue
        tile_cfg = json.loads((out / "tiles" / Path(entry["config"]).name).read_text(encoding="utf-8"))
        run = out / "tiles" / tile_cfg["output"]
        csvp = run / "force_displacement.csv"
        if not csvp.exists():
            continue
        geo = json.loads((run / "geometry.json").read_text(encoding="utf-8"))
        lo, hi = np.array(geo["bounds_mm"][0]), np.array(geo["bounds_mm"][1])
        crop = tile_cfg["crop_mm"]
        area = (crop[1] - crop[0]) * (crop[3] - crop[2])
        height = hi[2] - lo[2]
        with csvp.open(encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        strain = np.r_[0, np.array([float(r["indentation_mm"]) for r in rows]) / height]
        stress = np.r_[0, np.array([float(r["compression_force_n"]) for r in rows]) / area]  # MPa, compression +
        # A diverging run ends with increments whose force falls again (buckling / lost convergence): keep
        # only the monotone, converged part of the curve.
        peak = np.maximum.accumulate(stress)
        bad = np.flatnonzero(stress < 0.98 * peak)
        keep = int(bad[0]) if len(bad) else len(stress)
        strain, stress = strain[:keep], stress[:keep]
        # Foot pressures stay below ~0.15 MPa: fit the load range rather than a buckling plateau beyond it.
        within = int(np.searchsorted(stress, config["calibration_max_stress_mpa"], side="right"))
        if within >= 6:
            keep = min(keep, within)
        strain, stress = strain[:keep], stress[:keep]
        if keep < 4:
            raise ValueError(f"class {entry['class']}: fewer than 3 usable increments in {run}")
        fit = homog.fit_ogden(strain, stress)
        status = json.loads((run / "summary.json").read_text(encoding="utf-8")).get("status")
        fit.update(tile_height_mm=float(height), tile_area_mm2=float(area), max_strain=float(strain[-1]),
                   max_stress_mpa=float(stress[-1]), tile_status=status,
                   strain=strain.tolist(), stress_mpa=stress.tolist())
        fit["valid_max_stress_mpa"] = float(stress[-1])
        result[str(entry["class"])] = fit
    save_json(out / "materials.json", result)
    return result


def plot_fits(config):
    """Tile FEA curves with the fitted Ogden curves -> <output>/figures/tile_fits.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = Path(config["output"])
    materials = json.loads((out / "materials.json").read_text(encoding="utf-8"))
    fig, axes = plt.subplots(1, len(materials), figsize=(4 * len(materials), 3.8))
    for ax, (c, m) in zip(np.atleast_1d(axes), materials.items()):
        e, s = np.array(m["strain"]), np.array(m["stress_mpa"]) * 1000
        ee = np.linspace(0, max(e.max(), 1e-3), 40)
        fit = [-homog.ogden_uniaxial(1 - x, m["mu_mpa"], m["alpha"], m["d1_per_mpa"])[0] * 1000 for x in ee]
        ax.plot(e, s, "o", label="tile FEA")
        ax.plot(ee, fit, "-", label="Ogden fit")
        ax.set(title=f"class {c}  (H={m['tile_height_mm']:g} mm)", xlabel="nominal strain", ylabel="nominal stress [kPa]")
        ax.grid(True)
    np.atleast_1d(axes)[0].legend()
    fig.tight_layout()
    (out / "figures").mkdir(exist_ok=True)
    fig.savefig(out / "figures" / "tile_fits.png", dpi=110)
    plt.close(fig)


def build_cases(config, cases_wanted=None):
    """Write one CalculiX deck per gait case in <output>/cases/<name>."""
    out = Path(config["output"])
    table = homog.column_table(config["gcode"], config["dx_mm"], out)
    labels = np.load(out / "labels.npz")["labels"]
    materials = json.loads((out / "materials.json").read_text(encoding="utf-8"))
    used = {int(c) for c in np.unique(labels[labels >= 0])}
    missing = sorted(used - {int(k) for k in materials})
    if missing:
        raise ValueError(f"No calibrated material for classes {missing}; run calibrate after all tiles finish")
    _, _, frame, sensor_xy, _, cases = load_setup(config["gcode"], config["sensors"], config["gait"], config["side"],
                                                  config["body_mass_kg"], config["medial_sign"],
                                                  sensor_rotation_deg=config["sensor_rotation_deg"])
    dx, nz = config["dx_mm"], config["layers"]
    ox, oy = table["origin"]
    valid = labels >= 0
    nx, ny = valid.shape
    # Node heights: mean over the valid columns sharing the node -> conforming mesh.
    zz = np.stack([table["z0"], table["z0"] + table["H"]], axis=-1) * valid[..., None]
    zsum = np.zeros((nx + 1, ny + 1, 2))
    cnt = np.zeros((nx + 1, ny + 1))
    for di in (0, 1):
        for dj in (0, 1):
            sl = (slice(di, nx + di), slice(dj, ny + dj))
            cnt[sl] += valid
            zsum[sl] += zz
    node_col = cnt > 0
    zbot, ztop = zsum[..., 0] / np.maximum(cnt, 1), zsum[..., 1] / np.maximum(cnt, 1)
    ij = np.argwhere(node_col)
    node_index = -np.ones((nx + 1, ny + 1), int)
    node_index[node_col] = np.arange(len(ij))
    nn = len(ij) * (nz + 1)
    nid = lambda i, j, k: node_index[i, j] * (nz + 1) + k + 1
    if 3 * nn > config["max_solver_dofs"]:
        raise ValueError(f"{3*nn:,} DOFs exceed max_solver_dofs; increase dx_mm or reduce layers")
    xyz = np.empty((nn, 3))
    for m, (i, j) in enumerate(ij):
        z = zbot[i, j] + (ztop[i, j] - zbot[i, j]) * np.arange(nz + 1) / nz
        xyz[m * (nz + 1):(m + 1) * (nz + 1)] = np.c_[np.full(nz + 1, ox + i * dx), np.full(nz + 1, oy + j * dx), z]
    cols = np.argwhere(valid)
    elem, elem_class, top_elem, col_centre = [], [], [], []
    for (i, j) in cols:
        for k in range(nz):
            elem.append([nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                         nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)])
            elem_class.append(labels[i, j])
        top_elem.append(len(elem))  # 1-based id of the column's top element
        col_centre.append((ox + (i + 0.5) * dx, oy + (j + 0.5) * dx))
    elem, elem_class, top_elem = np.array(elem), np.array(elem_class), np.array(top_elem)
    col_centre = np.array(col_centre)
    ne = len(elem)
    bottom = np.flatnonzero(np.arange(nn) % (nz + 1) == 0) + 1
    xb = xyz[bottom - 1]
    a = bottom[np.argmin(np.linalg.norm(xb[:, :2] - col_centre.mean(axis=0), axis=1))]
    row = bottom[np.isclose(xyz[bottom - 1, 1], xyz[a - 1, 1])]
    b = row[np.argmax(np.abs(xyz[row - 1, 0] - xyz[a - 1, 0]))]
    stats = {"nodes": int(nn), "elements": int(ne), "columns": int(len(cols)), "dofs": int(3 * nn),
             "dx_mm": dx, "layers": nz, "area_mm2": float(len(cols) * dx * dx)}
    for name, case in cases.items():
        if cases_wanted and "all" not in cases_wanted and name not in cases_wanted:
            continue
        d = out / "cases" / name
        d.mkdir(parents=True, exist_ok=True)
        p = pressure_field((col_centre, dx * dx), sensor_xy, case["pressure"], case["force_n"])
        with (d / "model.inp").open("w", encoding="ascii", newline="\n") as f:
            w = lambda s: f.write(s + "\n")
            w(f"*HEADING\nWhole-insole homogenised model, case {name}; N mm MPa\n*NODE")
            for n, q in enumerate(xyz, 1):
                w(f"{n}, {q[0]:.6g}, {q[1]:.6g}, {q[2]:.6g}")
            w("*ELEMENT, TYPE=C3D8, ELSET=INSOLE")
            for n, e in enumerate(elem, 1):
                w(f"{n}, " + ", ".join(map(str, e)))
            for c in sorted(set(elem_class.tolist())):
                w(f"*ELSET, ELSET=CL{c}")
                ids = np.flatnonzero(elem_class == c) + 1
                for s in range(0, len(ids), 16):
                    w(", ".join(map(str, ids[s:s + 16])))
                m = materials[str(c)]
                w(f"*MATERIAL, NAME=OGDEN{c}\n*HYPERELASTIC, OGDEN, N=1\n{m['mu_mpa']:.8g}, {m['alpha']:.8g}, {m['d1_per_mpa']:.8g}")
                w(f"*SOLID SECTION, ELSET=CL{c}, MATERIAL=OGDEN{c}")
            w("*NSET, NSET=BOTTOM")
            for s in range(0, len(bottom), 16):
                w(", ".join(map(str, bottom[s:s + 16])))
            w(f"*BOUNDARY\nBOTTOM, 3, 3\n{a}, 1, 2\n{b}, 2, 2")
            # CalculiX scales the residual tolerance by the mean nodal force, which is tiny here (N over 90k nodes), so the default stalls although the
            # structure is in equilibrium. Tolerances are relaxed; the force balance is checked afterwards (summary.json).
            w("*STEP, NLGEOM, INC=1000\n*CONTROLS, PARAMETERS=FIELD\n20.,5.,,,2.,2.,,,,")
            w(f"*STATIC, SOLVER={config['linear_solver']}")
            w(f"{config['initial_increment']}, 1., {config['minimum_increment']}, {config['maximum_increment']}")
            # Fixed-direction nodal forces (p * column area shared by the four top nodes): the follower-pressure
            # load stiffness made Newton stall on this stiff, step-rich top surface, and the strains are small.
            w("*CLOAD")
            fz = np.zeros(nn + 1)
            for eid, pv in zip(top_elem, p):
                for node in elem[eid - 1][4:]:
                    fz[node] += pv * dx * dx / 4
            for node in np.flatnonzero(fz):
                w(f"{node}, 3, {-fz[node]:.6g}")
            w("*NODE FILE, FREQUENCY=1\nU, RF\n*EL FILE, FREQUENCY=1\nS, E")
            w("*NODE PRINT, NSET=BOTTOM, TOTALS=ONLY\nRF")
            w("*END STEP")
        save_json(d / "geometry.json", {**stats, "case": name, "origin_mm": [float(ox), float(oy), 0.0],
                                       "voxel_mm": dx, "grid_shape": [int(nx), int(ny), int(nz)],
                                       "column_mesh": True})
        beyond = {}
        for c in sorted(used):
            pc = p[labels[valid] == c]
            lim = materials[str(c)]["valid_max_stress_mpa"]
            if len(pc) and pc.max() > lim:
                beyond[str(c)] = {"max_pressure_kpa": float(pc.max() * 1000), "calibrated_to_kpa": float(lim * 1000),
                                  "columns_beyond": int((pc > lim).sum())}
        save_json(d / "build.json", {"deck_sha256": digest(d / "model.inp"), "case": name, "force_n": case["force_n"],
                                    "time_s": case["time"], "pressure_max_kpa": float(p.max() * 1000),
                                    "classes_beyond_calibration": beyond})
        save_json(d / "config.resolved.json", config)
    return stats


def solve_case(config, name):
    d = Path(config["output"]) / "cases" / name
    status = solve(config, d)
    summarize(d)
    return status


def summarize(d):
    """Reaction balance from the DAT file of one case."""
    d = Path(d)
    build = json.loads((d / "build.json").read_text(encoding="utf-8"))
    text = (d / "model.dat").read_text(errors="replace") if (d / "model.dat").exists() else ""
    m = re.findall(r"total force[^\n]*\n\s*([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)", text, re.I)
    rz = float(m[-1][2]) if m else None
    result = {"case": build["case"], "applied_force_n": build["force_n"], "reaction_z_n": rz,
              "balance_relative": abs(rz - build["force_n"]) / build["force_n"] if rz is not None else None,
              "pressure_max_kpa": build["pressure_max_kpa"]}
    save_json(d / "summary.json", result)
    return result


def main(args):
    config = load_config(args.config, create_output=args.action == "prepare")
    if args.action == "prepare":
        plan = prepare(config)
        print(json.dumps(plan, indent=2))
    elif args.action == "tiles":
        run_tiles(config)
    elif args.action == "calibrate":
        for k, v in calibrate(config).items():
            print(k, {n: (round(x, 4) if isinstance(x, float) else x) for n, x in v.items() if n not in ("strain", "stress_mpa")})
    elif args.action == "fits":
        plot_fits(config)
    elif args.action == "build":
        print(json.dumps(build_cases(config, args.case or ["all"]), indent=2))
    elif args.action == "solve":
        for name in args.case or sorted(p.name for p in (Path(config["output"]) / "cases").iterdir()):
            print(name, json.dumps(solve_case(config, name)), flush=True)
    return 0
