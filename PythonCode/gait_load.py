"""Moticon OpenGo gait records -> plantar pressure maps for static load cases. Units: N, mm, MPa."""
from __future__ import annotations

import numpy as np

G = 9.80665
N_PER_CM2_TO_MPA = 0.01
STANCE_MIN_N = 50.0
REGIONS = {"heel": range(0, 3), "midfoot": range(3, 5), "forefoot": range(5, 12), "toes": range(12, 16)}


class Gait:
    def __init__(self, path, side="left"):
        lines = open(path, encoding="utf-8").read().splitlines()
        header = [l for l in lines if l.startswith("# time")][0][2:].split("\t")
        rows = [l.split("\t") for l in lines if l and not l.startswith("#")]
        data = np.array([[np.nan if x == "" else float(x) for x in r] for r in rows])
        self.t = data[:, 0]
        cols = [header.index(f"{side} pressure {i}[N/cm²]") for i in range(1, 17)]
        self.pressure = np.nan_to_num(data[:, cols])  # N/cm^2, channel order
        self.force = np.nan_to_num(data[:, header.index(f"{side} total force[N]")])
        self.side = side

    def stances(self):
        """List of (start, end) index pairs of complete stance phases (loaded, within the record)."""
        on = self.force > STANCE_MIN_N
        edges = np.flatnonzero(np.diff(on.astype(int)))
        starts, ends = edges[on[edges + 1]] + 1, edges[~on[edges + 1]] + 1
        out = []
        for a in starts:
            b = ends[ends > a]
            if len(b) and a > 0:
                out.append((int(a), int(b[0])))
        return out

    def region_sum(self, index, region):
        return self.pressure[index][..., list(REGIONS[region])].sum(axis=-1)

    def representative_stance(self):
        stances = self.stances()
        peaks = np.array([self.force[a:b].max() for a, b in stances])
        return stances[int(np.argsort(peaks)[len(peaks) // 2])], peaks

    def cases(self, body_mass_kg, peak_factor=1.1, standing_factor=0.5, series=8):
        """Static cases: {name: {'pressure': 16 x N/cm^2 (shape), 'force_n': total, 'time': s}}."""
        (a, b), peaks = self.representative_stance()
        bw = body_mass_kg * G
        norm = peaks.mean() / (peak_factor * bw)  # recorded N per (N of physical load)
        sl = slice(a, b)
        idx = np.arange(a, b)
        n = len(idx)
        heel = self.region_sum(sl, "heel")
        fore = self.region_sum(sl, "forefoot") + self.region_sum(sl, "toes")
        heel_i = a + int(np.argmax(heel[: max(2, n // 2)]))
        push_i = a + n // 2 + int(np.argmax(fore[n // 2:]))
        mid = slice(a + int(0.3 * n), a + int(0.7 * n) + 1)
        out = {}
        flat = self.pressure[mid].mean(axis=0)
        out["standing"] = {"pressure": flat, "force_n": standing_factor * bw, "time": float(self.t[mid].mean())}
        for name, i in (("heel_strike", heel_i), ("push_off", push_i)):
            out[name] = {"pressure": self.pressure[i], "force_n": float(self.force[i] / norm), "time": float(self.t[i])}
        for k, i in enumerate(np.linspace(a, b - 1, series).round().astype(int)):
            out[f"walk_{k + 1:02d}"] = {"pressure": self.pressure[i], "force_n": float(self.force[i] / norm),
                                          "time": float(self.t[i])}
        return out


def pressure_field(points_xy, sensors_xy, sensor_p, force_n, sigma=15.0, floor_weight=0.05):
    """Pressure [MPa] at points_xy from 16 sensors, scaled to a total force (N) over the points' area.

    `points_xy` are cell centres of equal area `area`, supplied as (xy, area)."""
    xy, area = points_xy
    d2 = ((xy[:, None, :] - sensors_xy[None, :, :]) ** 2).sum(axis=2)
    w = np.exp(-d2 / (2 * sigma ** 2))
    p = (w * sensor_p[None, :]).sum(axis=1) / (w.sum(axis=1) + floor_weight)
    total = p.sum() * area
    return p * (force_n / total) if total > 0 else p
