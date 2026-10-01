import numpy as np
import pytest
from bambu_gcode3mf_to_fea_stl import parse_segments, DEFAULT_FEATURES, arc_points
from fea_geometry import Grid, voxelize


def parse(body):
    return parse_segments('; FEATURE: Outer wall\n; LINE_WIDTH: 0.45\nG90\nM83\n'+body,DEFAULT_FEATURES)[0]


def test_modes_and_comments():
    s = parse('G0 X1 Y2 Z0.2\nG1 X2 E1 ; X900\nG91\nG1 X1 E1\nG90\nM82 ; comment\nG92 E0\nG1 X4 E1\nG1 X5 E0.5\nG92 X10\nG1 X11 E2')
    assert len(s) == 4
    assert [x.p1[0] for x in s] == [2,3,4,6]
    assert s[-1].p0[0] == 5


def test_full_circle_and_layer_height():
    s = parse('G0 X1 Y0 Z0.2\nG3 I-1 J0 E2\n; LAYER_HEIGHT: 0.3\nG1 X2 E1')
    assert len(s) == 2
    p = arc_points(s[0],0.05)
    assert p[:,1].max() > 0.99 and p[:,1].min() < -0.99
    assert s[1].height == 0.3


@pytest.mark.parametrize('command',['G2 X1 R2 E1','G5 X1 E1','G20','G0 X1 E1','G90.1'])
def test_unsupported(command):
    with pytest.raises(ValueError):
        parse(command)


def test_cell_alignment_and_crop_crossing():
    s = parse('G0 X-1 Y0.1 Z0.2\nG1 X2 E1')
    g = voxelize(s,0.2,[0,1,0,0.4,0,0.2])
    assert g.masks.shape == (5,2,1)
    assert g.masks[:,0,0].all()
    assert g.origin.tolist() == [0,0,0]
    assert g.diagnostics()['face_components'] == 1


def test_shared_nodes_and_jacobian():
    g = Grid(np.ones((2,1,1),dtype=np.uint16),np.zeros(3),0.2,['Outer wall'])
    n,e,_,_ = g.mesh()
    assert len(n) == 12
    assert len(set(e[0]) & set(e[1])) == 4
    p = n[e[0]-1]
    assert np.linalg.det(np.stack([p[1]-p[0],p[3]-p[0],p[4]-p[0]])) == pytest.approx(0.2**3)


def test_feature_overlap_and_disconnected():
    s = parse('G0 X0 Y0.1 Z0.2\nG1 X1 E1\n; FEATURE: Sparse infill\nG1 X0 E1')
    g = voxelize(s,0.2,[0,1,0,0.2,0,0.2])
    assert np.all(g.masks == 3)
    masks = np.zeros((2,2,2),dtype=np.uint16)
    masks[0,0,0] = masks[1,1,1] = 1
    d = Grid(masks,np.zeros(3),1,['Outer wall']).diagnostics()
    assert d['face_components'] == 2 and d['edge_or_vertex_only_connections']


def test_invalid_grid():
    s = parse('G0 Z0.2\nG1 X1 E1')
    for pitch,crop in [(0,None),(0.2,[1,0,0,1,0,1]),(0.2,[5,6,5,6,0,1])]:
        with pytest.raises(ValueError):
            voxelize(s,pitch,crop)
