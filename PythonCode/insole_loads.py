"""Check figures for the foot frame, sensor placement and static gait load cases."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from bambu_gcode3mf_to_fea_stl import DEFAULT_FEATURES, parse_segments, read_gcode
from gait_load import Gait, pressure_field
from insole_frame import FootFrame, footprint, read_sensors


def rotate_clockwise(xy, centre, deg):
    """Rotate points clockwise (as drawn with X right, Y up) about `centre` by `deg` degrees."""
    a = -np.radians(deg)
    r = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    return (np.asarray(xy) - centre) @ r.T + centre


def load_setup(gcode, sensors, gait, side="left", body_mass_kg=53.0, medial_sign=None, pitch=1.0,
               sensor_rotation_deg=0.0):
    segments, _ = parse_segments(read_gcode(Path(gcode)), DEFAULT_FEATURES)
    mask, origin, p = footprint(segments, 0.5)
    frame = FootFrame(mask, origin, p, side, medial_sign=medial_sign)
    xy, names, ref_length = read_sensors(sensors, side)
    if ref_length:  # scale the reference insole (e.g. Moticon size 5) onto the printed one
        xy = xy * (frame.length / ref_length)
    sensor_xy = frame.to_xy(xy[:, 0], xy[:, 1])
    if sensor_rotation_deg:  # reference insole differs from the printed one: turn the layout about the insole centre
        sensor_xy = rotate_clockwise(sensor_xy, frame.centre, sensor_rotation_deg)
    cases = Gait(gait, side).cases(body_mass_kg)
    return segments, (mask, origin, p), frame, sensor_xy, names, cases


def field_on_mask(mask, origin, p, stride=2):
    """Cell centres (subsampled by `stride`) inside the footprint and their cell area."""
    ij = np.argwhere(mask)[:: 1]
    ij = ij[(ij[:, 0] % stride == 0) & (ij[:, 1] % stride == 0)]
    return origin + (ij + 0.5) * p, (p * stride) ** 2


def plot(gcode, sensors, gait, out, **kw):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    segments, (mask, origin, p), frame, sensor_xy, names, cases = load_setup(gcode, sensors, gait, **kw)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(mask.T, origin="lower", extent=(origin[0], origin[0] + mask.shape[0] * p, origin[1], origin[1] + mask.shape[1] * p),
              cmap="Greys", alpha=0.25)
    ax.plot(*sensor_xy.T, "ro")
    for i, (q, nm) in enumerate(zip(sensor_xy, names), 1):
        ax.annotate(f"{i}", q, fontsize=8, xytext=(3, 3), textcoords="offset points")
    tip = frame.centre + (frame.s0 + frame.length) * frame.t
    heel = frame.centre + frame.s0 * frame.t
    ax.annotate("heel", heel, color="blue"); ax.annotate("toe", tip, color="blue")
    med = frame.centre + 10 * frame.medial_sign * frame.n
    ax.annotate("medial side ->", frame.centre + frame.medial_sign * frame.n * 45, color="green")
    ax.set_aspect("equal"); ax.set(xlabel="G-code X [mm]", ylabel="G-code Y [mm]", title=f"Sensors on insole (auto frame, rotated {kw.get('sensor_rotation_deg', 0):g} deg clockwise)")
    fig.savefig(out, dpi=130); plt.close(fig)

    pts = field_on_mask(mask, origin, p, 2)
    names_c = [k for k in cases if not k.startswith("walk_")]
    fig, axes = plt.subplots(1, len(names_c), figsize=(5 * len(names_c), 6))
    for ax, k in zip(np.atleast_1d(axes), names_c):
        c = cases[k]
        pf = pressure_field(pts, sensor_xy, c["pressure"], c["force_n"])
        sc = ax.scatter(*pts[0].T, c=pf * 1000, s=2, cmap="turbo", vmin=0)
        ax.plot(*sensor_xy.T, "k.", ms=3)
        ax.set_aspect("equal"); ax.set_title(f"{k}  F={c['force_n']:.0f} N  t={c['time']:.2f} s")
        fig.colorbar(sc, ax=ax, label="pressure [kPa]", shrink=0.7)
    fig.tight_layout(); fig.savefig(Path(out).with_name(Path(out).stem + "_cases.png"), dpi=120); plt.close(fig)
    return frame, cases
