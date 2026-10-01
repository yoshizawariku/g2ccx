#!/usr/bin/env python3
"""
Convert Bambu Studio .gcode.3mf sliced toolpaths into a voxel-unioned STL
suitable for meshing in PrePoMax / Gmsh / other FEA preprocessors.

Designed for Bambu G-code containing:
  ; FEATURE: ...
  ; LINE_WIDTH: ...
  G1 linear extrusion moves
  G2/G3 arc extrusion moves using I/J center offsets

The generated solid approximates each deposited bead as an elliptic/rounded
swept volume with width = slicer LINE_WIDTH and height = layer height.
This preserves wall / skin / sparse-infill topology (e.g. gyroid) far better
than rebuilding from the original CAD body.

Dependencies:
    pip install numpy scikit-image

Example:
    python bambu_gcode3mf_to_fea_stl.py Insole-L2.gcode.3mf insole_fdm.stl --voxel 0.20

Local submodel (recommended for stress convergence studies):
    python bambu_gcode3mf_to_fea_stl.py Insole-L2.gcode.3mf heel_roi.stl \
        --voxel 0.10 --crop 40 100 40 110 0 8
"""
from __future__ import annotations

import argparse
import json
import math
import re
import struct
import zipfile
from dataclasses import dataclass
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np
from skimage.measure import marching_cubes

NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)"
PARAM_RE = re.compile(r"([XYZIJER])(" + NUM + r")")
MOVE_RE = re.compile(r"^(G[0123])\b")

DEFAULT_FEATURES = {
    "Outer wall",
    "Inner wall",
    "Sparse infill",
    "Internal solid infill",
    "Top surface",
    "Bottom surface",
    "Gap infill",
    "Bridge",
    "Floating vertical shell",
    "Overhang wall",
}

@dataclass
class Segment:
    p0: np.ndarray
    p1: np.ndarray
    cmd: str
    ij: tuple[float, float] | None
    width: float
    height: float
    feature: str
    layer: int


def read_gcode(path: Path) -> str:
    if path.suffix.lower() == ".3mf" or path.name.lower().endswith(".gcode.3mf"):
        with zipfile.ZipFile(path, "r") as zf:
            candidates = [n for n in zf.namelist() if n.lower().endswith(".gcode")]
            if not candidates:
                raise RuntimeError("No .gcode member found inside 3MF")
            if len(candidates) != 1:
                raise ValueError("Multiple G-code plates found; export one plate or extract the desired .gcode: " + ", ".join(candidates))
            preferred = next((n for n in candidates if n.endswith("plate_1.gcode")), candidates[0])
            return zf.read(preferred).decode("utf-8", errors="ignore")
    return path.read_text(encoding="utf-8", errors="ignore")


def config_value(gcode: str, key: str, default=None):
    m = re.search(r"^; " + re.escape(key) + r" = (.*)$", gcode, flags=re.MULTILINE)
    return m.group(1).strip() if m else default


def parse_float_first(s, default):
    if s is None:
        return default
    try:
        return float(str(s).split(",")[0].strip())
    except Exception:
        return default


def parse_segments(gcode: str, included_features: set[str]) -> tuple[list[Segment], dict]:
    layer_h = parse_float_first(config_value(gcode, "layer_height"), 0.2)
    default_w = parse_float_first(config_value(gcode, "line_width"), 0.42)
    nozzle = parse_float_first(config_value(gcode, "nozzle_diameter"), None)
    infill = config_value(gcode, "sparse_infill_pattern", None)
    infill_density = config_value(gcode, "sparse_infill_density", None)

    pos = np.array([0.0, 0.0, 0.0], dtype=float)
    coordinate_offset = np.zeros(3)
    feature = "Custom"
    width = default_w
    layer = 0
    relative_e = True
    relative_xyz = False
    plane = "G17"
    e_abs = 0.0
    segments: list[Segment] = []
    counts = Counter()
    widths = defaultdict(list)

    for raw in gcode.splitlines():
        s = raw.strip()
        if not s:
            continue
        if s.startswith("; LAYER_HEIGHT:"):
            layer_h = float(s.split(":", 1)[1])
            if layer_h <= 0:
                raise ValueError("Layer height must be positive")
            continue
        if not s.startswith(";"):
            s = s.split(";", 1)[0].strip()
        if s in {"G90", "G91"}:
            relative_xyz = s == "G91"
            continue
        if s == "G90.1":
            raise ValueError("Absolute arc centres are unsupported; use relative I/J centres")
        if s in {"G17", "G18", "G19"}:
            plane = s
            continue
        if s == "G20":
            raise ValueError("Inch G-code is not supported; export in mm")
        if s == "M82":
            relative_e = False
            continue
        if s == "M83":
            relative_e = True
            continue
        if re.match(r"^G92\b", s):
            vals = dict((k, float(v)) for k, v in PARAM_RE.findall(s))
            if "E" in vals:
                e_abs = vals["E"]
            for i, k in enumerate("XYZ"):
                if k in vals:
                    coordinate_offset[i] = pos[i] - vals[k]
            continue
        if s.startswith("; FEATURE:"):
            feature = s.split(":", 1)[1].strip()
            continue
        if s.startswith("; LINE_WIDTH:"):
            try:
                width = float(s.split(":", 1)[1].strip())
                widths[feature].append(width)
            except ValueError:
                pass
            continue
        if s.startswith("; layer num/total_layer_count:"):
            try:
                layer = int(s.split(":", 1)[1].split("/", 1)[0].strip())
            except Exception:
                pass
            continue

        mm = MOVE_RE.match(s)
        if not mm:
            if re.match(r"^G\d+", s) and "E" in dict(PARAM_RE.findall(s)) and feature in included_features:
                raise ValueError(f"Unsupported extrusion command: {s}")
            continue
        cmd = mm.group(1)
        vals = {k: float(v) for k, v in PARAM_RE.findall(s)}
        new = pos.copy()
        for i, k in enumerate("XYZ"):
            if k in vals:
                new[i] = pos[i] + vals[k] if relative_xyz else vals[k] + coordinate_offset[i]

        extruding = False
        if "E" in vals:
            if relative_e:
                extruding = vals["E"] > 0.0
            else:
                extruding = vals["E"] > e_abs + 1e-12
                e_abs = vals["E"]

        xy_motion = ((abs(new[0] - pos[0]) + abs(new[1] - pos[1])) > 1e-12
                     or (cmd in {"G2", "G3"} and ("I" in vals or "J" in vals)))
        if cmd == "G0" and extruding and xy_motion and feature in included_features:
            raise ValueError(f"Rapid extrusion is unsupported: {s}")
        if cmd in {"G1", "G2", "G3"} and extruding and xy_motion and feature in included_features:
            ij = None
            if cmd in {"G2", "G3"}:
                if plane != "G17" or "R" in vals or not ("I" in vals or "J" in vals):
                    raise ValueError(f"Only XY arcs with I/J offsets are supported: {s}")
                ij = (vals.get("I", 0.0), vals.get("J", 0.0))
            if width <= 0 or layer_h <= 0:
                raise ValueError("Bead width and height must be positive")
            segments.append(Segment(pos.copy(), new.copy(), cmd, ij, width, layer_h, feature, layer))
            counts[feature] += 1
        pos = new

    report = {
        "layer_height_mm": layer_h,
        "nozzle_diameter_mm": nozzle,
        "default_line_width_mm": default_w,
        "sparse_infill_pattern": infill,
        "sparse_infill_density": infill_density,
        "included_features": sorted(included_features),
        "segment_counts": dict(counts),
        "line_width_median_mm": {
            k: float(np.median(v)) for k, v in widths.items() if v
        },
        "segment_total": len(segments),
    }
    return segments, report


def arc_points(seg: Segment, max_step: float) -> np.ndarray:
    p0, p1 = seg.p0, seg.p1
    if seg.cmd == "G1" or seg.ij is None:
        length = float(np.linalg.norm(p1 - p0))
        n = max(1, int(math.ceil(length / max_step)))
        return np.linspace(p0, p1, n + 1)

    cx = p0[0] + seg.ij[0]
    cy = p0[1] + seg.ij[1]
    r0x, r0y = p0[0] - cx, p0[1] - cy
    r1x, r1y = p1[0] - cx, p1[1] - cy
    r = math.hypot(r0x, r0y)
    if r < 1e-10:
        return np.vstack([p0, p1])
    a0 = math.atan2(r0y, r0x)
    a1 = math.atan2(r1y, r1x)
    if seg.cmd == "G3":  # CCW
        da = (a1 - a0) % (2.0 * math.pi)
        if da < 1e-12 and np.linalg.norm(p1[:2] - p0[:2]) < 1e-9:
            da = 2.0 * math.pi
    else:  # G2 CW
        da = -((a0 - a1) % (2.0 * math.pi))
        if abs(da) < 1e-12 and np.linalg.norm(p1[:2] - p0[:2]) < 1e-9:
            da = -2.0 * math.pi
    length = abs(da) * r
    n = max(1, int(math.ceil(length / max_step)))
    t = np.linspace(0.0, 1.0, n + 1)
    a = a0 + da * t
    z = p0[2] + (p1[2] - p0[2]) * t
    return np.column_stack([cx + r * np.cos(a), cy + r * np.sin(a), z])


def crop_segment_points(points: np.ndarray, crop):
    if crop is None:
        return points
    xmin, xmax, ymin, ymax, zmin, zmax = crop
    m = ((points[:,0] >= xmin) & (points[:,0] <= xmax) &
         (points[:,1] >= ymin) & (points[:,1] <= ymax) &
         (points[:,2] >= zmin) & (points[:,2] <= zmax))
    return points[m]


def binary_stl_write(path: Path, verts: np.ndarray, faces: np.ndarray):
    header = b"Bambu G-code reconstructed FDM solid".ljust(80, b" ")[:80]
    with path.open("wb") as f:
        f.write(header)
        f.write(struct.pack("<I", len(faces)))
        for tri in faces:
            a, b, c = verts[tri[0]], verts[tri[1]], verts[tri[2]]
            n = np.cross(b - a, c - a)
            ln = np.linalg.norm(n)
            if ln > 0:
                n = n / ln
            else:
                n = np.zeros(3)
            f.write(struct.pack("<12fH",
                                float(n[0]), float(n[1]), float(n[2]),
                                float(a[0]), float(a[1]), float(a[2]),
                                float(b[0]), float(b[1]), float(b[2]),
                                float(c[0]), float(c[1]), float(c[2]), 0))


def reconstruct(segments: list[Segment], voxel: float, crop=None):
    from fea_geometry import voxelize
    grid = voxelize(segments, voxel, crop)
    occupied = np.pad(grid.masks != 0, 1)
    verts, faces, _, _ = marching_cubes(occupied, 0.5, spacing=(voxel,)*3)
    verts += grid.origin - voxel/2
    hi = grid.origin + np.array(grid.masks.shape)*voxel
    return verts.astype(np.float32), faces.astype(np.int32), grid.origin, hi, np.array(grid.masks.shape)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input", type=Path, help="Bambu .gcode.3mf or plain .gcode")
    ap.add_argument("output", type=Path, help="output STL")
    ap.add_argument("--voxel", type=float, default=0.20,
                    help="voxel pitch in mm (default 0.20; use 0.10-0.15 for local submodels)")
    ap.add_argument("--crop", nargs=6, type=float, metavar=("XMIN","XMAX","YMIN","YMAX","ZMIN","ZMAX"),
                    help="optional ROI crop in mm")
    ap.add_argument("--features", nargs="*", default=None,
                    help="feature names to include; default = all structural printed features")
    ap.add_argument("--report", type=Path, default=None, help="write JSON parse report")
    args = ap.parse_args()

    features = set(args.features) if args.features else set(DEFAULT_FEATURES)
    gcode = read_gcode(args.input)
    segments, report = parse_segments(gcode, features)
    print(json.dumps(report, indent=2, ensure_ascii=False))

    verts, faces, lo, hi, shape = reconstruct(segments, args.voxel, args.crop)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    binary_stl_write(args.output, verts, faces)
    report.update({
        "voxel_pitch_mm": args.voxel,
        "crop_mm": args.crop,
        "grid_shape": [int(x) for x in shape],
        "bbox_min_mm": [float(x) for x in lo],
        "bbox_max_mm": [float(x) for x in hi],
        "stl_vertices": int(len(verts)),
        "stl_triangles": int(len(faces)),
        "output_stl": str(args.output),
    })
    rp = args.report or args.output.with_suffix(args.output.suffix + ".json")
    rp.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Wrote {args.output} ({len(faces):,} triangles)")
    print(f"Wrote {rp}")

if __name__ == "__main__":
    main()
