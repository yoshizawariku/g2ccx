"""Foot frame of a sliced insole and mapping of insole-sensor coordinates onto G-code XY. Units: mm."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy import ndimage

from bambu_gcode3mf_to_fea_stl import arc_points


def footprint(segments, pitch=0.5):
    """Filled XY outline: (mask[nx, ny], origin_xy, pitch)."""
    pts = np.vstack([arc_points(s, pitch)[:, :2] for s in segments if s.feature in {"Outer wall", "Inner wall"}])
    lo = pts.min(axis=0) - 3 * pitch
    n = np.ceil((pts.max(axis=0) + 3 * pitch - lo) / pitch).astype(int) + 1
    mask = np.zeros(n, bool)
    idx = np.floor((pts - lo) / pitch).astype(int)
    mask[idx[:, 0], idx[:, 1]] = True
    mask = ndimage.binary_closing(mask, iterations=2)
    return ndimage.binary_fill_holes(mask), lo, pitch


class FootFrame:
    """s: heel->toe distance, v: lateral distance from the local mid-line (positive = medial)."""

    def __init__(self, mask, origin, pitch, side="left", toe="min_xy", medial_sign=None):
        self.side = side
        ij = np.argwhere(mask)
        xy = origin + (ij + 0.5) * pitch
        centre = xy.mean(axis=0)
        w, vec = np.linalg.eigh(np.cov((xy - centre).T))
        t = vec[:, 1]
        toward_small = -np.array([1.0, 1.0])
        if toe == "min_xy":
            t = t if t @ toward_small > 0 else -t
        elif toe == "max_xy":
            t = t if t @ toward_small < 0 else -t
        self.t = t
        self.n = np.array([-t[1], t[0]])  # +90 deg from the toe direction
        s = (xy - centre) @ self.t
        v = (xy - centre) @ self.n
        self.centre, self.s0 = centre, s.min()
        self.length = float(s.max() - s.min())
        # Lateral profile: mid-line and both edges as functions of s (1 mm bins).
        bins = np.arange(0, self.length + 1, 1.0)
        k = np.clip(((s - self.s0)).astype(int), 0, len(bins) - 1)
        lo = np.full(len(bins), np.nan)
        hi = np.full(len(bins), np.nan)
        for b in range(len(bins)):
            sel = v[k == b]
            if len(sel):
                lo[b], hi[b] = sel.min(), sel.max()
        ok = np.isfinite(lo)
        lo, hi = np.interp(bins, bins[ok], lo[ok]), np.interp(bins, bins[ok], hi[ok])
        self.bins, self.edge_lo, self.edge_hi = bins, lo, hi
        smooth = lambda a: ndimage.gaussian_filter1d(a, 20.0, mode="nearest")
        self.mid = smooth((lo + hi) / 2)
        if medial_sign is None:
            # The arch is the deeper waist: compare each edge with the chord between 15 % and 85 % of the length.
            i0, i1 = int(0.15 * self.length), int(0.85 * self.length)
            def depth(e):
                chord = np.interp(np.arange(i0, i1), [i0, i1], [e[i0], e[i1]])
                return (e[i0:i1] - chord)
            # For the +n edge an inward dent is negative deviation; for the -n edge it is positive.
            dent_hi = -depth(hi).min()
            dent_lo = depth(lo).max()
            medial_sign = 1.0 if dent_hi > dent_lo else -1.0
        self.medial_sign = float(medial_sign)

    def to_xy(self, x_mm, y_mm):
        """Sensor (x medial-positive, y from heel edge) -> G-code XY."""
        s = self.s0 + np.asarray(y_mm, float)
        v = np.interp(np.asarray(y_mm, float), self.bins, self.mid) + self.medial_sign * np.asarray(x_mm, float)
        return self.centre + s[..., None] * self.t + v[..., None] * self.n


def read_sensors(path, side, heel_edge_from_imu_mm=None):
    """Sensor table -> (x medial-positive, y from heel edge) in the *reference insole's* millimetres, names, ref length.

    Two layouts are understood: `side,channel,anatomical_region,x_mm,y_mm` (already heel-referenced), and the
    Moticon IMU-referenced table `sensor_no,X_mm_from_IMU,Y_mm_from_IMU,region,coordinate_note` (X toward the toe,
    Y medial) whose note gives the reference insole size. There the heel edge is assumed to lie half of the
    reference length behind the IMU unless `heel_edge_from_imu_mm` is given (negative)."""
    import re
    with Path(path).open(encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if "sensor_no" in rows[0]:
        rows.sort(key=lambda r: int(r["sensor_no"]))
        m = re.search(r"([\d.]+)\s*x\s*([\d.]+)\s*mm", rows[0].get("coordinate_note", ""))
        length = float(m[2]) if m else None
        edge = heel_edge_from_imu_mm if heel_edge_from_imu_mm is not None else -(length or 261.1) / 2
        xy = np.array([[float(r["Y_mm_from_IMU"]), float(r["X_mm_from_IMU"]) - edge] for r in rows])
        return xy, [r["region"] for r in rows], length
    rows = sorted((r for r in rows if r["side"] == side), key=lambda r: int(r["channel"]))
    return (np.array([[float(r["x_mm"]), float(r["y_mm"])] for r in rows]),
            [r["anatomical_region"] for r in rows], None)
