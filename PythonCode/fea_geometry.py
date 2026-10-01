"""Cell-centred, globally aligned voxel geometry. Units: mm."""
from __future__ import annotations

from dataclasses import dataclass
import math
from functools import lru_cache
import numpy as np
from scipy import ndimage
from skimage.measure import marching_cubes
from bambu_gcode3mf_to_fea_stl import arc_points, binary_stl_write

CORNERS = np.array([[0,0,0], [1,0,0], [1,1,0], [0,1,0],
                    [0,0,1], [1,0,1], [1,1,1], [0,1,1]], dtype=np.int32)
FEATURE_SETS = {
    "Outer wall": "OUTER_WALL", "Inner wall": "INNER_WALL",
    "Sparse infill": "SPARSE_INFILL", "Internal solid infill": "SOLID_INFILL",
    "Top surface": "TOP_SURFACE", "Bottom surface": "BOTTOM_SURFACE",
    "Gap infill": "GAP_INFILL", "Bridge": "BRIDGE",
    "Floating vertical shell": "FLOATING_VERTICAL_SHELL", "Overhang wall": "OVERHANG_WALL",
}


@lru_cache(maxsize=1)
def local_connectivity_lookup():
    result = np.zeros(256, dtype=np.uint8)
    for pattern in range(256):
        cube = np.array([(pattern >> i) & 1 for i in range(8)]).reshape(2,2,2)
        result[pattern] = ndimage.label(cube)[1]
    return result


def bounds(segments):
    lo, hi = np.full(3, np.inf), np.full(3, -np.inf)
    for s in segments:
        p = arc_points(s, 0.1).copy()
        p[:, 2] -= s.height / 2
        radius = np.array([s.width/2, s.width/2, s.height/2])
        lo = np.minimum(lo, p.min(axis=0)-radius)
        hi = np.maximum(hi, p.max(axis=0)+radius)
    if not np.isfinite(lo).all():
        raise ValueError("No structural extrusion found")
    return lo, hi


@dataclass
class Grid:
    masks: np.ndarray
    origin: np.ndarray
    pitch: float
    features: list[str]

    def diagnostics(self):
        occupied = self.masks != 0
        labels, count = ndimage.label(occupied)
        sizes = np.bincount(labels.ravel())[1:]
        _, touch_count = ndimage.label(occupied, structure=np.ones((3,3,3)))
        # Detect local edge/vertex joints even when a longer face-connected path exists.
        node_patterns = np.zeros(np.array(occupied.shape)+1,dtype=np.uint8)
        for bit,offset in enumerate(np.ndindex(2,2,2)):
            slices = tuple(slice(d,d+n) for d,n in zip(offset,occupied.shape))
            node_patterns[slices] |= occupied.astype(np.uint8) << bit
        weak_nodes = int(np.count_nonzero(local_connectivity_lookup()[node_patterns] > 1))
        n = int(occupied.sum())
        # Conservative assembly estimate, not a solver factorization guarantee.
        return {"grid_shape": list(self.masks.shape), "elements": n,
                "volume_mm3": n*self.pitch**3, "face_components": count,
                "touch_components": touch_count,
                "component_sizes": sorted(sizes.tolist(), reverse=True),
                "edge_or_vertex_only_connections": count != touch_count,
                "local_nonmanifold_vertices": weak_nodes,
                "assembly_estimate_mib": math.ceil(n*2500/2**20),
                "origin_mm": self.origin.tolist(), "voxel_mm": self.pitch}

    def mesh(self, max_elements=2_000_000):
        cells = np.argwhere(self.masks != 0).astype(np.int32)
        if len(cells) == 0 or len(cells) > max_elements:
            raise ValueError(f"Occupied elements {len(cells):,}; allowed 1..{max_elements:,}. Adjust ROI/voxel/limit.")
        shape = np.array(self.masks.shape) + 1
        corners = cells[:, None, :] + CORNERS
        keys = np.ravel_multi_index(corners.reshape(-1, 3).T, shape)
        unique, inverse = np.unique(keys, return_inverse=True)
        nodes = np.column_stack(np.unravel_index(unique, shape))*self.pitch+self.origin
        elements = inverse.reshape(-1,8)+1
        masks = self.masks[tuple(cells.T)]
        return nodes, elements, cells, masks

    def stl(self, path):
        occupied = np.pad(self.masks != 0, 1)
        v, f, _, _ = marching_cubes(occupied, 0.5, spacing=(self.pitch,)*3)
        v += self.origin-self.pitch/2
        binary_stl_write(path, v, f)
        return v, f


def voxelize(segments, pitch=0.2, crop=None, max_grid_cells=50_000_000):
    if not np.isfinite(pitch) or pitch <= 0:
        raise ValueError("voxel must be finite and positive")
    if crop is None:
        lo, hi = bounds(segments)
    else:
        box = np.asarray(crop, dtype=float)
        if box.shape != (6,) or not np.isfinite(box).all() or np.any(box[1::2] <= box[::2]):
            raise ValueError("crop must be xmin xmax ymin ymax zmin zmax, with min < max")
        lo, hi = box[::2], box[1::2]
    # Include cells whose centres lie in the requested ROI. Faces snap to global lattice.
    first = np.ceil(lo/pitch-0.5-1e-9).astype(np.int64)
    end = np.ceil(hi/pitch-0.5-1e-9).astype(np.int64)
    shape = end-first
    ngrid = math.prod(int(x) for x in shape)
    if np.any(shape <= 0) or ngrid > max_grid_cells:
        raise ValueError(f"Grid has {ngrid:,} cells; limit {max_grid_cells:,}. Reduce ROI or increase voxel.")
    features = sorted({s.feature for s in segments})
    if len(features) > 16:
        raise ValueError("At most 16 structural features supported")
    bits = {f: np.uint16(1 << i) for i, f in enumerate(features)}
    masks = np.zeros(tuple(shape), dtype=np.uint16)
    origin = first*pitch
    cache = {}
    for i, s in enumerate(segments):
        p = arc_points(s, pitch*0.3).copy()
        p[:, 2] -= s.height/2
        radius = np.array([s.width/2, s.width/2, s.height/2])
        p = p[np.all((p >= origin-radius-pitch) & (p <= end*pitch+radius+pitch), axis=1)]
        if not len(p):
            continue
        key = tuple(np.ceil(radius/pitch).astype(int)+1)
        if key not in cache:
            cache[key] = np.stack(np.meshgrid(*[np.arange(-n,n+1) for n in key], indexing="ij"), axis=-1).reshape(-1,3)
        offsets = cache[key]
        # Bound temporary arrays even for very long straight moves.
        for start in range(0, len(p), 512):
            points = p[start:start+512]
            centres = np.floor((points-origin)/pitch).astype(np.int32)
            ids = centres[:,None,:]+offsets
            delta = (origin+(ids+0.5)*pitch-points[:,None,:])/radius
            good = (np.sum(delta**2, axis=2) <= 1+1e-9) & np.all((ids >= 0) & (ids < shape),axis=2)
            q = ids[good]
            if len(q):
                masks[tuple(q.T)] |= bits[s.feature]
        if i and i % 10000 == 0:
            print(f"Voxelized {i:,}/{len(segments):,}", flush=True)
    if not masks.any():
        raise ValueError("ROI has no occupied cells; check ROI and voxel size")
    return Grid(masks, origin, pitch, features)
