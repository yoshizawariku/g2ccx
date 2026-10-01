"""Reproducible Bambu G-code -> voxel -> CalculiX pipeline."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

import numpy as np
from bambu_gcode3mf_to_fea_stl import DEFAULT_FEATURES, arc_points, parse_segments, read_gcode
from fea_geometry import bounds, voxelize
from fea_solver import postprocess, save_json, solve, write_deck
from g2ccx_paths import resolve_output

DEFAULTS = {
    "voxel_mm":0.2, "roi_size_mm":[60.0,60.0], "crop_mm":None,
    "displacement_mm":2.0, "material":{"c10":-0.44,"c01":4.99,"d1":0.0},
    "friction":0.0, "contact_stiffness_mpa_per_mm":1000.0,
    "initial_increment":0.01, "minimum_increment":1e-6, "maximum_increment":0.05,
    "max_elements":2_000_000, "max_grid_cells":50_000_000, "max_solver_dofs":500_000,
    "threads":min(4,os.cpu_count() or 1), "timeout_seconds":7200,
    "solver":None, "linear_solver":"PARDISO", "element_type":"C3D8R", "output_frequency":5, "export_stl":False,
}


def load_config(path, create_output=False):
    path = Path(path).resolve()
    supplied = json.loads(path.read_text(encoding="utf-8"))
    unknown = set(supplied)-set(DEFAULTS)-{"input","output"}
    if unknown:
        raise ValueError(f"Unknown config keys: {sorted(unknown)}")
    config = {**DEFAULTS, **supplied}
    config["material"] = {**DEFAULTS["material"], **supplied.get("material",{})}
    for key in ("input","output","solver"):
        if key not in config:
            raise ValueError(f"Config requires {key}")
        if config[key] is None:  # solver: looked up when the solver starts (config > $G2CCX_SOLVER > local file > PATH)
            continue
        p = Path(config[key])
        config[key] = str((path.parent/p).resolve()) if not p.is_absolute() else str(p.resolve())
        if key == "output":
            config[key] = resolve_output(config[key],create_output)
    for key in ("voxel_mm","displacement_mm","contact_stiffness_mpa_per_mm","initial_increment",
                "minimum_increment","maximum_increment","timeout_seconds"):
        if not isinstance(config[key],(int,float)) or not np.isfinite(config[key]) or config[key] <= 0:
            raise ValueError(f"{key} must be finite and positive")
    for key in ("threads","max_elements","max_grid_cells","max_solver_dofs","output_frequency"):
        if not isinstance(config[key],int) or config[key] <= 0:
            raise ValueError(f"{key} must be a positive integer")
    if not 0 < config["minimum_increment"] <= config["initial_increment"] <= config["maximum_increment"] <= 1:
        raise ValueError("Require 0 < minimum <= initial <= maximum increment <= 1")
    if config['linear_solver'] not in {'PARDISO','PASTIX','SPOOLES'}:
        raise ValueError('linear_solver must be PARDISO, PASTIX or SPOOLES')
    if config['element_type'] not in {'C3D8R','C3D8'}:
        raise ValueError('element_type must be C3D8R or C3D8')
    if not np.isfinite(config["friction"]) or config["friction"] < 0:
        raise ValueError("friction must be finite and nonnegative")
    m = config["material"]
    if set(m) != {"c10","c01","d1"} or not all(np.isfinite(v) for v in m.values()) or m["c10"]+m["c01"] <= 0 or m["d1"] < 0:
        raise ValueError("Invalid Mooney-Rivlin material parameters")
    if len(config["roi_size_mm"]) != 2 or not all(np.isfinite(x) and x>0 for x in config["roi_size_mm"]):
        raise ValueError("roi_size_mm needs two positive lengths")
    return config


def inspect_input(path):
    segments,report = parse_segments(read_gcode(Path(path)), DEFAULT_FEATURES)
    lo,hi = bounds(segments)
    report.update(bounds_mm=[lo.tolist(),hi.tolist()], layer_count=len({s.layer for s in segments}),
                  arc_segments=sum(s.cmd in {"G2","G3"} for s in segments))
    return segments, report


def projection(segments, ax):
    from matplotlib.collections import LineCollection
    # Top/bottom contours and sparse infill retain a legible overview.
    layers = sorted({s.layer for s in segments})
    selected = {layers[0],layers[len(layers)//2],layers[-1]}
    curves = [arc_points(s,0.5)[:,:2] for s in segments
              if s.layer in selected and s.feature in {"Outer wall","Sparse infill"}]
    ax.add_collection(LineCollection(curves,linewidths=0.35,colors="#277da8",alpha=0.7))
    ax.autoscale()
    ax.set_aspect("equal")
    ax.set(xlabel="G-code X [mm]", ylabel="G-code Y [mm]")


def select_roi(config, path):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.widgets import Button
    segments,report = inspect_input(config["input"])
    lo,hi = np.array(report["bounds_mm"])
    fig,ax = plt.subplots(figsize=(9,8))
    fig.subplots_adjust(bottom=0.15)
    projection(segments,ax)
    w,h = config["roi_size_mm"]
    rect = Rectangle((0,0),w,h,fill=False,edgecolor="crimson",linewidth=2,visible=False)
    ax.add_patch(rect)
    selected = []
    title = ax.set_title("Click heel centre, then Save ROI. Close to cancel.")
    def click(event):
        toolbar = getattr(fig.canvas.manager,'toolbar',None)
        if event.inaxes is ax and event.button == 1 and not (toolbar and toolbar.mode):
            x,y = event.xdata,event.ydata
            selected[:] = [x-w/2,x+w/2,y-h/2,y+h/2,float(lo[2]),float(hi[2])]
            rect.set_xy((x-w/2,y-h/2))
            rect.set_visible(True)
            title.set_text(f"Centre ({x:.2f}, {y:.2f}); {w:g} x {h:g} mm. Save ROI to confirm.")
            fig.canvas.draw_idle()
    def save(_):
        if not selected:
            title.set_text("Select a centre first")
            fig.canvas.draw_idle()
            return
        config["crop_mm"] = selected.copy()
        # Edit only crop_mm in the user's file; do not rewrite its relative paths / {timestamp} placeholders.
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        raw["crop_mm"] = selected.copy()
        save_json(path,raw)
        fig.savefig(Path(path).with_suffix(".roi.png"),dpi=150)
        plt.close(fig)
        print(f"Saved ROI to {path}")
    button = Button(fig.add_axes((0.4,0.025,0.2,0.06)),"Save ROI")
    button.on_clicked(save)
    fig.canvas.mpl_connect("button_press_event",click)
    plt.show()


def build(config):
    if config["crop_mm"] is None:
        raise ValueError("ROI not selected. Run select-roi --config ... first, or set crop_mm explicitly.")
    out = Path(config["output"])
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {out}. Choose a new directory.")
    out.mkdir(parents=True,exist_ok=True)
    save_json(out/"config.resolved.json",config)
    segments,report = inspect_input(config["input"])
    save_json(out/"input.json",report)
    grid = voxelize(segments,config["voxel_mm"],config["crop_mm"],config["max_grid_cells"])
    diagnostic = grid.diagnostics()
    save_json(out/"geometry.json",diagnostic)
    print(json.dumps(diagnostic,indent=2),flush=True)
    if config["export_stl"]:
        grid.stl(out/"geometry.stl")
    return write_deck(grid,config,out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="command",required=True)
    parser = sub.add_parser("inspect",help="Inspect sliced input without creating a mesh")
    parser.add_argument("input",type=Path)
    parser.add_argument("--preview",type=Path,help="Save an XY preview PNG")
    parser = sub.add_parser("insole",help="Whole-insole homogenised gait-load analysis")
    parser.add_argument("action",choices=["prepare","tiles","calibrate","build","solve"])
    parser.add_argument("--config",type=Path,required=True)
    parser.add_argument("--case",nargs="*",help="Gait case names (default: all)")
    parser = sub.add_parser("view",help="Open the heat-map slice viewer for a finished analysis")
    parser.add_argument("--context",type=float,nargs="?",const=0.4,metavar="VOXEL_MM",help="With --3d, also draw the whole insole (coarse voxel, default 0.4 mm) around the analysed ROI")
    parser.add_argument("--3d",dest="three_d",action="store_true",help="Open the 3-D viewer instead")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--config",type=Path)
    group.add_argument("--output",type=Path)
    parser.add_argument("--save",type=Path,help="Render the last increment to a PNG instead of opening a window")
    for name in ("select-roi","build","solve","postprocess","run"):
        parser = sub.add_parser(name)
        parser.add_argument("--config",type=Path,required=True)
    args = ap.parse_args(argv)
    try:
        if args.command == "inspect":
            segments,report = inspect_input(args.input)
            print(json.dumps(report,indent=2,ensure_ascii=False))
            if args.preview:
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.pyplot as plt
                fig,ax = plt.subplots(figsize=(8,8))
                projection(segments,ax)
                fig.tight_layout()
                args.preview.parent.mkdir(parents=True,exist_ok=True)
                fig.savefig(args.preview,dpi=150)
                plt.close(fig)
            return 0
        if args.command == "insole":
            import insole_global
            return insole_global.main(args)
        if args.command == "view":
            from fea_viewer import view, view3d
            target = args.output or load_config(args.config)["output"]
            if args.three_d:
                view3d(target,args.save,args.context)
            else:
                view(target,args.save)
            return 0
        config = load_config(args.config,create_output=args.command in {"build","run","select-roi"})
        if args.command == "select-roi":
            select_roi(config,args.config)
            return 0
        if args.command in {"build","run"}:
            build(config)
        if args.command in {"solve","postprocess"}:
            resolved = Path(config["output"])/"config.resolved.json"
            original = load_config(resolved)
            # Analysis parameters must correspond to the built deck.
            runtime = {'solver','threads','timeout_seconds','max_solver_dofs'}
            if any(config[k] != original[k] for k in config if k not in runtime):
                raise ValueError("Configuration differs from the built model; use config.resolved.json or rebuild")
        if args.command in {"solve","run"}:
            status = solve(config,config["output"])
            print(json.dumps(status,indent=2))
        if args.command in {"postprocess","solve","run"}:
            if not (Path(config["output"])/"model.dat").exists():
                raise RuntimeError("No DAT output; inspect solver.log")
            summary = postprocess(config["output"])
            print(json.dumps(summary,indent=2))
            return 0 if summary["status"] == "complete" and summary["force_balance_pass"] else 2
        return 0
    except (ValueError,RuntimeError,OSError) as exc:
        print(f"ERROR: {exc}",file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
