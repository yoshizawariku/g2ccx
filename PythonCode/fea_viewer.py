"""Interactive slice/heat-map viewer for CalculiX results written by this pipeline."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

FIELDS = {
    "von Mises stress [MPa]": "mises",
    "Max principal stress [MPa]": "smax",
    "Min principal stress [MPa]": "smin",
    "Equivalent strain (E, deviatoric)": "eqstrain",
    "Displacement |U| [mm]": "umag",
    "Displacement Uz [mm]": "uz",
}


class Frd:
    """Random-access reader for the ASCII .frd nodes, hex8 elements and nodal result blocks."""

    def __init__(self, path):
        self.buf = np.memmap(path, dtype=np.uint8, mode="r")
        nl = np.flatnonzero(self.buf == 10)
        self.starts = np.concatenate(([0], nl[:-1] + 1))
        self.ends = nl  # exclusive; a trailing '\r' is trimmed by the fixed-width slices
        self._index()

    def _line(self, i):
        return bytes(self.buf[self.starts[i]:self.ends[i]]).decode("ascii", "replace").rstrip("\r")

    def _rows(self, first, last, width):
        """Fixed-width rows [first, last) as a (n, width) uint8 matrix."""
        idx = self.starts[first:last, None] + np.arange(width)
        return np.asarray(self.buf[idx])

    def _index(self):
        s = self.starts
        c = [np.asarray(self.buf[s + k]) for k in range(6)]
        tag = lambda text: np.flatnonzero(np.all([c[k] == ord(ch) for k, ch in enumerate(text)], axis=0))
        self.node_block = int(tag("    2C")[0])
        self.elem_block = int(tag("    3C")[0])
        end_marks = np.flatnonzero((c[0] == 32) & (c[1] == ord("-")) & (c[2] == ord("3")))
        self.results = {}
        for h in tag("  100C"):
            name = self._line(h + 1)[5:17].strip()
            comps = 0
            j = h + 2
            while self._line(j).startswith(" -5"):
                comps += 1
                j += 1
            end = int(end_marks[np.searchsorted(end_marks, j)])
            time = float(self._line(h)[12:24])
            self.results.setdefault(name, []).append((time, j, end))
        self._end_marks = end_marks
        self.node_ids, self.node_xyz = self._nodes()
        self.elem_ids, self.elem_nodes = self._elements()

    def _fixed(self, first, last, start, widths, count):
        rows = self._rows(first, last, start + widths * count)
        cols = rows[:, start:start + widths * count].copy()
        return cols.view(f"S{widths}").reshape(len(rows), count).astype(float)

    def _nodes(self):
        first = self.node_block + 1
        last = int(self._end_marks[np.searchsorted(self._end_marks, first)])
        ids = self._fixed(first, last, 3, 10, 1)[:, 0].astype(int)
        xyz = self._fixed(first, last, 13, 12, 3)
        return ids, xyz

    def _elements(self):
        first = self.elem_block + 1
        last = int(self._end_marks[np.searchsorted(self._end_marks, first)])
        s = self.starts
        head = np.arange(first, last)
        is2 = np.asarray(self.buf[s[head] + 2]) == ord("2")
        hdr, body = head[~is2], head[is2]
        if len(hdr) != len(body):
            raise ValueError("Only single-line 8-node hexahedra are supported in the FRD")
        rows = self._rows_at(hdr, 13)
        ids = rows[:, 3:13].copy().view("S10").ravel().astype(float).astype(int)
        nrows = self._rows_at(body, 3 + 80)
        nodes = nrows[:, 3:83].copy().view("S10").reshape(len(body), 8).astype(float).astype(int)
        return ids, nodes

    def _rows_at(self, lines, width):
        idx = self.starts[lines, None] + np.arange(width)
        return np.asarray(self.buf[idx])

    def times(self, name="STRESS"):
        return [t for t, _, _ in self.results.get(name, [])]

    def nodal(self, name, step):
        """Nodal component matrix (n_nodes_in_frd_order, ncomp) for the step-th block of `name`."""
        _, first, end = self.results[name][step]
        ncomp = {"STRESS": 6, "TOSTRAIN": 6, "DISP": 3}[name]
        return self._fixed(first, end, 13, 12, ncomp)


class Result:
    """Per-element scalar fields on the voxel grid of one output directory."""

    def __init__(self, output):
        self.output = Path(output)
        geo = json.loads((self.output / "geometry.json").read_text(encoding="utf-8"))
        self.origin = np.array(geo["origin_mm"], float)
        self.voxel = float(geo["voxel_mm"])
        self.shape = tuple(geo["grid_shape"])
        self.n_tpu = int(geo["elements"])
        frd = self.frd = Frd(self.output / "model.frd")
        keep = frd.elem_ids <= self.n_tpu
        self.elem_nodes = frd.elem_nodes[keep]
        lookup = np.full(frd.node_ids.max() + 1, -1)
        lookup[frd.node_ids] = np.arange(len(frd.node_ids))
        self.rows = lookup[self.elem_nodes]
        centre = frd.node_xyz[self.rows].mean(axis=1)
        self.column_mesh = bool(geo.get("column_mesh"))  # variable-height columns: 3-D view only
        self.cell = np.floor((centre - self.origin) / self.voxel + 1e-6).astype(int)
        if self.column_mesh:
            pass
        elif (self.cell < 0).any() or (self.cell >= np.array(self.shape)).any():
            raise ValueError("FRD elements fall outside geometry.json grid; results and geometry differ")
        self.times = frd.times("STRESS")
        self.summary = self._read_json("summary.json")
        self.force = self._force_curve()
        self._cache = {}

    def _read_json(self, name):
        p = self.output / name
        return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}

    def _force_curve(self):
        p = self.output / "force_displacement.csv"
        if not p.exists():
            return None
        with p.open(encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        return (np.array([float(r["time"]) for r in rows]),
                np.array([float(r["indentation_mm"]) for r in rows]),
                np.array([float(r["compression_force_n"]) for r in rows]))

    def indentation(self, step):
        if self.force is None:
            return None
        k = int(np.argmin(np.abs(self.force[0] - self.times[step])))
        return float(self.force[1][k]), float(self.force[2][k])

    def _element_tensor(self, name, step):
        v = self.frd.nodal(name, step)[self.rows].mean(axis=1)
        a = np.zeros((len(v), 3, 3))
        a[:, 0, 0], a[:, 1, 1], a[:, 2, 2] = v[:, 0], v[:, 1], v[:, 2]
        a[:, 0, 1] = a[:, 1, 0] = v[:, 3]
        a[:, 1, 2] = a[:, 2, 1] = v[:, 4]
        a[:, 0, 2] = a[:, 2, 0] = v[:, 5]
        return a

    def elements(self, key, step):
        """Scalar per TPU element (mean of its 8 nodal values)."""
        if (key, step) in self._cache:
            return self._cache[key, step]
        if key in {"mises", "smax", "smin"}:
            eig = np.linalg.eigvalsh(self._element_tensor("STRESS", step))
            out = {"smax": eig[:, 2], "smin": eig[:, 0],
                   "mises": np.sqrt(((eig[:, 0]-eig[:, 1])**2 + (eig[:, 1]-eig[:, 2])**2 + (eig[:, 2]-eig[:, 0])**2)/2)}[key]
        elif key == "eqstrain":
            e = self._element_tensor("TOSTRAIN", step)
            dev = e - np.trace(e, axis1=1, axis2=2)[:, None, None] / 3 * np.eye(3)
            out = np.sqrt(2/3 * np.einsum("nij,nij->n", dev, dev))
        else:
            u = self.frd.nodal("DISP", step)[self.rows].mean(axis=1)
            out = np.linalg.norm(u, axis=1) if key == "umag" else u[:, 2]
        self._cache[key, step] = out
        return out

    def volume(self, key, step):
        if self.column_mesh:
            raise RuntimeError("Slice view needs a regular voxel model; use view --3d for the whole-insole case")
        vol = np.full(self.shape, np.nan)
        vol[tuple(self.cell.T)] = self.elements(key, step)
        return vol


def _draw(res, key, step, pos, axes, scale, cmap="turbo"):
    vol = res.volume(key, step)
    lo, hi = scale
    n = np.array(res.shape)
    org, h = res.origin, res.voxel
    edges = [org[a] + h * np.arange(n[a] + 1) for a in range(3)]
    views = (("XY", 2, 0, 1), ("XZ", 1, 0, 2), ("YZ", 0, 1, 2))
    mesh = []
    for ax, (label, fixed, u, v) in zip(axes, views):
        sl = np.take(vol, pos[fixed], axis=fixed)
        ax.clear()
        m = ax.pcolormesh(edges[u], edges[v], np.ma.masked_invalid(sl).T, cmap=cmap, vmin=lo, vmax=hi, shading="flat")
        ax.set_aspect("equal")
        ax.set_facecolor("#dddddd")
        ax.set_title(f"{label}  ({'XYZ'[fixed]} = {org[fixed] + h * (pos[fixed] + 0.5):.2f} mm)", fontsize=9)
        ax.set_xlabel(f"{'XYZ'[u]} [mm]")
        ax.set_ylabel(f"{'XYZ'[v]} [mm]")
        mesh.append(m)
    return mesh[0]


def _range(res, key, step, per_step):
    vals = res.elements(key, step if per_step else len(res.times) - 1)
    lo, hi = float(np.min(vals)), float(np.max(vals))
    return (lo, hi) if hi > lo else (lo - 1e-9, hi + 1e-9)


def view(output, save=None):
    """Open the viewer; with `save`, render the last increment to a PNG without a window."""
    import matplotlib
    if save:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.widgets import CheckButtons, RadioButtons, Slider

    res = Result(output)
    if not res.times:
        raise RuntimeError("No stress results in model.frd; nothing to show")
    state = {"key": "mises", "step": len(res.times) - 1, "per_step": False,
             "pos": [s // 2 for s in res.shape]}
    fig = plt.figure(figsize=(14, 8))
    fig.canvas.manager.set_window_title(f"g2ccx viewer - {res.output.name}")
    gs = fig.add_gridspec(2, 3, left=0.3, right=0.93, top=0.93, bottom=0.27, hspace=0.35, wspace=0.35)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])]
    curve = fig.add_subplot(gs[1, :])
    cax = fig.add_axes((0.94, 0.55, 0.012, 0.38))
    holder = {}

    def title():
        ind = res.indentation(state["step"])
        text = f"increment {state['step'] + 1}/{len(res.times)}  t={res.times[state['step']]:.4f}"
        if ind:
            text += f"   indentation {ind[0]:.3f} mm   force {ind[1]:.3f} N"
        fig.suptitle(text + f"   [{res.summary.get('status', '?')}]", fontsize=10)

    def redraw(_=None):
        key = FIELDS[state["field"]] if "field" in state else state["key"]
        state["key"] = key
        scale = _range(res, key, state["step"], state["per_step"])
        mesh = _draw(res, key, state["step"], state["pos"], axes, scale)
        cax.clear()
        fig.colorbar(mesh, cax=cax).set_label(next(k for k, v in FIELDS.items() if v == key), fontsize=8)
        curve.clear()
        if res.force is not None:
            curve.plot(res.force[1], res.force[2], ".-", color="#277da8")
            ind = res.indentation(state["step"])
            curve.plot([ind[0]], [ind[1]], "o", color="crimson", ms=9)
        curve.set(xlabel="Indentation [mm]", ylabel="Compression force [N]")
        curve.grid(True)
        title()
        fig.canvas.draw_idle()

    holders = []
    for i, (label, n) in enumerate(zip("XYZ", res.shape)):
        s = Slider(fig.add_axes((0.34, 0.17 - 0.04 * i, 0.5, 0.025)), f"{label} slice", 0, max(n - 1, 1),
                   valinit=state["pos"][i], valstep=1)
        s.on_changed(lambda v, i=i: (state["pos"].__setitem__(i, int(v)), redraw()))
        holders.append(s)
    inc = Slider(fig.add_axes((0.34, 0.04, 0.5, 0.025)), "Increment", 1, max(len(res.times), 2),
                 valinit=len(res.times), valstep=1)
    inc.on_changed(lambda v: (state.__setitem__("step", int(v) - 1), redraw()))
    holders.append(inc)
    names = list(FIELDS)
    state["field"] = names[0]
    radio = RadioButtons(fig.add_axes((0.01, 0.35, 0.2, 0.3)), names)
    radio.on_clicked(lambda label: (state.__setitem__("field", label), redraw()))
    check = CheckButtons(fig.add_axes((0.01, 0.2, 0.2, 0.08)), ["Scale per increment"], [False])
    check.on_clicked(lambda _: (state.__setitem__("per_step", not state["per_step"]), redraw()))
    holders += [radio, check]

    def key_press(event):
        if event.key in {"left", "right"}:
            inc.set_val(min(max(inc.val + (1 if event.key == "right" else -1), 1), len(res.times)))
    fig.canvas.mpl_connect("key_press_event", key_press)
    fig._g2ccx_widgets = holders  # keep widgets alive
    redraw()
    if save:
        fig.savefig(save, dpi=130)
        plt.close(fig)
        print(f"Saved {save}")
    else:
        plt.show()


def _context_surface(res, pitch):
    """Whole-insole surface from the G-code named in config.resolved.json, ROI removed.

    Beads are thinner than a coarse voxel, so the grid is built at the layer height and only
    XY is pooled (`pitch` mm, a multiple of the layer voxel) to keep the display light.
    """
    import pyvista as pv
    from skimage.measure import marching_cubes
    from bambu_gcode3mf_to_fea_stl import DEFAULT_FEATURES, parse_segments, read_gcode
    from fea_geometry import voxelize

    fine = 0.2
    k = max(1, round(pitch / fine))
    config = json.loads((res.output / "config.resolved.json").read_text(encoding="utf-8"))
    cache = res.output / f"context_{k}x{fine:g}.npz"
    if cache.exists():
        z = np.load(cache)
        occ, origin = z["occ"], z["origin"]
    else:
        print(f"Voxelizing the whole insole at {fine:g} mm, pooled x{k} in XY (cached afterwards)...", flush=True)
        segments, _ = parse_segments(read_gcode(Path(config["input"])), DEFAULT_FEATURES)
        grid = voxelize(segments, fine, None, 200_000_000)
        occ, origin = grid.masks != 0, grid.origin
        nx, ny = (occ.shape[0] // k) * k, (occ.shape[1] // k) * k
        occ = occ[:nx, :ny].reshape(nx // k, k, ny // k, k, -1).any(axis=(1, 3))
        np.savez_compressed(cache, occ=occ, origin=origin)
    spacing = (fine * k, fine * k, fine)
    if config.get("crop_mm"):
        c = np.asarray(config["crop_mm"], float)
        inside = [(origin[a] + spacing[a] * (np.arange(occ.shape[a]) + 0.5) >= c[2*a]) &
                  (origin[a] + spacing[a] * (np.arange(occ.shape[a]) + 0.5) <= c[2*a+1]) for a in range(3)]
        occ = occ.copy()
        occ[np.ix_(*inside)] = False
    v, f, _, _ = marching_cubes(np.pad(occ, 1), 0.5, spacing=spacing)
    v += origin - np.array(spacing)
    surface = pv.PolyData(v, np.hstack([np.full((len(f), 1), 3), f]).ravel())
    return surface, config.get("crop_mm")


def view3d(output, save=None, context=None):
    """3-D view of the TPU hexahedra coloured by an element field, optionally on the deformed shape."""
    import pyvista as pv

    res = Result(output)
    if not res.times:
        raise RuntimeError("No stress results in model.frd; nothing to show")
    n = len(res.rows)
    cells = np.hstack([np.full((n, 1), 8), res.rows]).ravel()
    grid = pv.UnstructuredGrid(cells, np.full(n, pv.CellType.HEXAHEDRON), res.frd.node_xyz.copy())
    names = list(FIELDS)
    state = {"field": names[0], "step": len(res.times) - 1, "deform": True, "scale": 1.0, "clip": False,
             "per_step": False, "cut": 1.0}
    pl = pv.Plotter(off_screen=bool(save), window_size=(1400, 850))
    pl.set_background("white")

    cy = res.frd.node_xyz[res.rows][:, :, 1].mean(axis=1)

    def update(*_):
        key = FIELDS[state["field"]]
        step = state["step"]
        disp = res.frd.nodal("DISP", step)
        grid.points = res.frd.node_xyz + (state["scale"] * disp if state["deform"] else 0)
        grid.cell_data["value"] = res.elements(key, step)
        lo, hi = _range(res, key, step, state["per_step"])
        keep = cy <= cy.min() + state["cut"] * (cy.max() - cy.min()) + 1e-9
        shown = grid if keep.all() else grid.extract_cells(np.flatnonzero(keep))
        pl.add_mesh(shown, scalars="value", cmap="turbo", clim=(lo, hi), show_edges=False, name="tpu",
                    scalar_bar_args={"title": state["field"], "color": "black", "vertical": True,
                                     "position_x": 0.90, "position_y": 0.2, "height": 0.6, "width": 0.06})
        ind = res.indentation(step)
        text = f"increment {step + 1}/{len(res.times)}  t={res.times[step]:.4f}"
        if ind:
            text += f"\nindentation {ind[0]:.3f} mm  force {ind[1]:.3f} N"
        pl.add_text(text + f"\ndeformation x{state['scale']:g}" if state["deform"] else text, name="info",
                    font_size=10, color="black")

    def set_field(index):
        state["field"] = names[index]
        update()

    def step_by(delta):
        state["step"] = min(max(state["step"] + delta, 0), len(res.times) - 1)
        update()

    platen = None
    if context:
        surface, _ = _context_surface(res, context)
        pl.add_mesh(surface, color="#b9b9b9", opacity=1.0, smooth_shading=False, name="insole")
        xyz = res.frd.node_xyz[res.rows].reshape(-1, 3)
        lo, hi = xyz.min(axis=0), xyz.max(axis=0)
        margin = max(2.0, float((hi - lo)[:2].max()) * 0.1)  # same as the solver's rigid platen
        platen = pv.Plane(center=((lo[0]+hi[0])/2, (lo[1]+hi[1])/2, 0), direction=(0, 0, 1),
                          i_size=hi[0]-lo[0]+2*margin, j_size=hi[1]-lo[1]+2*margin)
        ztop = float(hi[2])

    def move_platen():
        if platen is not None:
            ind = res.indentation(state["step"])
            z = ztop - (ind[0] if ind else 0)
            pl.add_mesh(platen.copy().translate((0, 0, z), inplace=False), color="#3b6ea5", opacity=0.25,
                        name="platen")

    _update = update
    def update(*a):
        _update(*a)
        move_platen()

    update()
    pl.add_axes(color="black")
    pl.add_slider_widget(lambda v: (state.__setitem__("step", int(round(v)) - 1), update()),
                         [1, max(len(res.times), 2)], value=len(res.times), title="Increment",
                         fmt="%.0f", pointa=(0.05, 0.06), pointb=(0.3, 0.06), style="modern")
    pl.add_slider_widget(lambda v: (state.__setitem__("scale", float(v)), update()), [0, 5], value=1.0,
                         title="Deformation scale", fmt="%.1f", pointa=(0.37, 0.06), pointb=(0.62, 0.06),
                         style="modern")
    pl.add_slider_widget(lambda v: (state.__setitem__("cut", float(v)), update()), [0, 1], value=1.0,
                         title="Cut (Y)", fmt="%.2f", pointa=(0.69, 0.06), pointb=(0.94, 0.06), style="modern")
    for i in range(len(names)):
        pl.add_key_event(str(i + 1), lambda i=i: set_field(i))
    pl.add_key_event("Right", lambda: step_by(1))
    pl.add_key_event("Left", lambda: step_by(-1))
    pl.add_key_event("d", lambda: (state.__setitem__("deform", not state["deform"]), update()))
    pl.add_key_event("s", lambda: (state.__setitem__("per_step", not state["per_step"]), update()))
    pl.add_text("1-6 field / Left,Right increment\nd deformed / s scale per increment",
                position="upper_right", font_size=7, color="black", name="help")
    pl.camera_position = "iso"
    if save:
        pl.screenshot(str(save))
        pl.close()
        print(f"Saved {save}")
    else:
        pl.show(title=f"g2ccx 3D viewer - {res.output.name}")
