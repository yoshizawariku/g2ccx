"""Column-wise homogenisation of a sliced insole: geometry table, classes and tile calibration."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.cluster.vq import kmeans2
from scipy.optimize import least_squares

from bambu_gcode3mf_to_fea_stl import DEFAULT_FEATURES, parse_segments, read_gcode
from fea_geometry import voxelize

FINE = 0.2


def column_table(gcode, dx, cache_dir):
    """Per column (dx x dx): bottom z0, height H, volume fraction phi (whole column), skin fractions.

    The bead-resolved 0.2 mm occupancy is reduced by *counting* (not any-pooling), so phi is unbiased."""
    cache_dir = Path(cache_dir)
    cache = cache_dir / f"columns_{dx:g}.npz"
    if cache.exists():
        return dict(np.load(cache))
    k = round(dx / FINE)
    if abs(k * FINE - dx) > 1e-9:
        raise ValueError("dx must be a multiple of 0.2 mm")
    fine = cache_dir / "fine_occupancy.npz"
    if fine.exists():
        z = np.load(fine)
        occ, origin = z["occ"], z["origin"]
    else:
        print(f"Voxelizing {gcode} at {FINE:g} mm (one-off)...", flush=True)
        segments, _ = parse_segments(read_gcode(Path(gcode)), DEFAULT_FEATURES)
        grid = voxelize(segments, FINE, None, 200_000_000)
        occ, origin = grid.masks != 0, np.asarray(grid.origin, float)
        cache_dir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(fine, occ=occ, origin=origin)
    nx, ny, nz = occ.shape[0] // k, occ.shape[1] // k, occ.shape[2]
    occ = occ[: nx * k, : ny * k].reshape(nx, k, ny, k, nz)
    count = occ.sum(axis=(1, 3)).astype(np.float32) / (k * k)  # (nx, ny, nz): fill of each 0.2-mm slab
    has = count > 0.02
    any_ = has.any(axis=2)
    z0 = np.where(any_, has.argmax(axis=2), 0) * FINE
    top = np.where(any_, nz - has[:, :, ::-1].argmax(axis=2), 0) * FINE
    H = top - z0
    cum = np.concatenate([np.zeros((nx, ny, 1), np.float32), count.cumsum(axis=2)], axis=2)
    layers = np.clip(np.round(H / FINE).astype(int), 0, nz)
    i0 = np.round(z0 / FINE).astype(int)
    def frac(a, b):
        ia = np.clip(i0 + a, 0, nz)
        ib = np.clip(i0 + b, 0, nz)
        vol = np.take_along_axis(cum, ib[..., None], 2)[..., 0] - np.take_along_axis(cum, ia[..., None], 2)[..., 0]
        return vol / np.maximum(ib - ia, 1)
    table = {"origin": origin[:2],
             "dx": np.array(dx), "z0": z0, "H": H, "valid": any_ & (H >= 1.0),
             "phi": frac(0, layers), "phi_bot": frac(0, np.minimum(layers, 5)),
             "phi_top": frac(np.maximum(layers - 5, 0), layers), "layers": layers}
    cache_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, **table)
    return table


def classify(table, k=5, seed=1):
    """K-means classes over (H, phi, phi_bot, phi_top). Returns labels (-1 = void) and class centres."""
    v = table["valid"]
    feat = np.stack([table[n][v] for n in ("H", "phi", "phi_bot", "phi_top")], axis=1)
    scale = feat.std(axis=0) + 1e-9
    rng = np.random.default_rng(seed)
    centres, lab = kmeans2(feat / scale, k, minit="++", seed=rng, missing="raise")
    order = np.argsort(centres[:, 1])  # class 0 = lowest phi
    remap = np.empty(k, int)
    remap[order] = np.arange(k)
    labels = np.full(v.shape, -1, int)
    labels[v] = remap[lab]
    return labels, (centres[order] * scale)


def representative_tiles(table, labels, centres, tile_mm=12.0, per_class=8):
    """Candidate tile centres per class, best (highest class purity) first and at least one tile apart."""
    dx = float(table["dx"])
    half = int(round(tile_mm / dx / 2))
    out = []
    for c in range(len(centres)):
        found = []
        for i, j in np.argwhere(labels == c)[:: 2]:
            win = labels[max(i - half, 0): i + half, max(j - half, 0): j + half]
            if win.shape != (2 * half, 2 * half) or (win < 0).any():
                continue
            found.append(((win == c).mean(), i, j))
        found.sort(reverse=True)
        chosen = []
        for purity, i, j in found:
            if all(abs(i - a) > 2 * half or abs(j - b) > 2 * half for _, a, b in chosen):
                chosen.append((purity, i, j))
            if len(chosen) == per_class:
                break
        out.append([{"class": c, "purity": float(p),
                     "centre_mm": [float(v) for v in table["origin"] + (np.array([i, j]) + 0.5) * dx],
                     "tile_mm": tile_mm} for p, i, j in chosen])
    return out


# ---- Compressible one-term Ogden (CalculiX *HYPERELASTIC,OGDEN,N=1) under uniaxial compression -----
# W = 2 mu/alpha^2 (sum lam_bar^alpha - 3) + (J-1)^2 / D ; lateral faces free (zero lateral stress).
# Verified against a one-element CalculiX run (sigma_z and lateral stretch agree to 6 digits).
def _nominal(lam, ll, mu, alpha, d):
    li = np.array([lam, ll, ll])
    j = lam * ll * ll
    lb = j ** (-1 / 3) * li
    return (2 * mu / alpha) / li * (lb ** alpha - (lb ** alpha).sum() / 3) + (2 / d) * (j - 1) * j / li


def ogden_uniaxial(lam, mu, alpha, d):
    """Axial nominal stress (tension +) and lateral stretch for axial stretch `lam`."""
    from scipy.optimize import brentq
    if abs(lam - 1) < 1e-12:
        return 0.0, 1.0
    f = lambda ll: _nominal(lam, ll, mu, alpha, d)[1]
    ll = brentq(f, 0.2, 5.0)
    return float(_nominal(lam, ll, mu, alpha, d)[0]), float(ll)


def fit_ogden(strain, stress_mpa, nu0=0.15):
    """Fit (mu, alpha) to compressive nominal stress (positive) vs nominal strain (positive).

    Uniaxial data cannot identify the lateral behaviour, so D1 follows from a fixed initial Poisson ratio `nu0`
    (K = 2 mu (1+nu)/(3 (1-2 nu)), D1 = 2/K); cellular structures are only weakly lateral-coupled (nu ~ 0.1-0.2)."""
    d_of = lambda mu: 3 * (1 - 2 * nu0) / (mu * (1 + nu0))
    lam = 1 - np.asarray(strain, float)
    target = -np.asarray(stress_mpa, float)
    scale = np.abs(target).max() + 1e-12
    def model(q):
        return np.array([ogden_uniaxial(l, np.exp(q[0]), q[1], d_of(np.exp(q[0])))[0] for l in lam])
    best = None
    for a0 in (-2, 2, 4):
        q0 = [np.log(max(scale, 1e-3)), a0]
        try:
            r = least_squares(lambda q: (model(q) - target) / scale, q0,
                              bounds=([-12, -4], [6, 6]), x_scale="jac")
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    q = best.x
    rms = float(np.sqrt(np.mean((model(q) - target) ** 2)) / scale)
    mu, alpha = float(np.exp(q[0])), float(q[1])
    d = float(d_of(mu))
    return {"mu_mpa": mu, "alpha": alpha, "d1_per_mpa": d, "rel_rms_error": rms,
            "poisson_lateral_at_max_strain": float(ogden_uniaxial(lam[-1], mu, alpha, d)[1] - 1) / max(strain[-1], 1e-9)}
