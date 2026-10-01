import json
import sys
from pathlib import Path
import numpy as np
import pytest
from bambu_fea import DEFAULTS
from fea_geometry import Grid
from fea_solver import write_deck, dat_blocks, tensor_stats, postprocess, save_json, digest, solve


def test_deck(tmp_path):
    g = Grid(np.ones((4,4,4),dtype=np.uint16),np.zeros(3),0.5,['Outer wall'])
    c = {**DEFAULTS,'input':None,'displacement_mm':0.2}
    d = write_deck(g,c,tmp_path)
    s = (tmp_path/'model.inp').read_text()
    assert '*HYPERELASTIC, MOONEY-RIVLIN' in s
    assert '*CONTACT PAIR' in s and 'TPU_DOWN, FLOOR_TOP' in s
    assert '*FRICTION' not in s  # CalculiX rejects a zero friction coefficient card.
    assert '*STATIC, SOLVER=PARDISO' in s
    assert d['references']['LOAD_REF'] > d['nodes']+16
    assert d['elements'] == 64


def test_dat_and_tensor(tmp_path):
    p = tmp_path/'model.dat'
    p.write_text(' displacements (vx,vy,vz) for set LOAD_REF and time 0.1000000E+01\n\n 100 0 0 -0.2\n\n forces (fx,fy,fz) for set LOAD_REF and time 0.1000000E+01\n\n 100 0 0 -12\n')
    blocks = list(dat_blocks(p))
    assert len(blocks) == 2 and blocks[0][2] == 1
    assert blocks[1][3][0,3] == -12
    stats = tensor_stats(np.array([[0,0,-10,0,0,0]]))
    assert stats['von_mises_measure_max'] == pytest.approx(10)


def test_disconnected_deck_rejected(tmp_path):
    a = np.zeros((3,3,3),dtype=np.uint16)
    a[0,0,0] = a[2,2,2] = 1
    with pytest.raises(ValueError,match='face-connected'):
        write_deck(Grid(a,np.zeros(3),0.2,['Outer wall']),DEFAULTS,tmp_path)


@pytest.mark.parametrize('finished,target,status',[(True,0.2,'complete'),(False,0.2,'incomplete'),(True,0.5,'incomplete')])
def test_postprocess_completion_and_stale_results(tmp_path, finished, target, status):
    g = Grid(np.ones((2,2,2),dtype=np.uint16),np.zeros(3),0.5,['Outer wall'])
    write_deck(g,{**DEFAULTS,'input':None,'displacement_mm':target},tmp_path)
    save_json(tmp_path/'run.json',{'solver_finished':finished,'solver_errors':False,'returncode':0,
                                 'deck_sha256':digest(tmp_path/'model.inp')})
    text = ''
    for field,name,vector in [('displacements','LOAD_REF','0 0 -0.2'),('forces','LOAD_REF','0 0 -10'),('forces','FLOOR_REF','0 0 10')]:
        text += f' {field} (x,y,z) for set {name} and time 1.0\n\n 100 {vector}\n'
    (tmp_path/'model.dat').write_text(text)
    result = postprocess(tmp_path)
    assert result['status'] == status and result['force_balance_pass']
    (tmp_path/'model.inp').write_text('changed')
    with pytest.raises(ValueError,match='stale'):
        postprocess(tmp_path)


def test_real_sample_parse():
    from bambu_fea import inspect_input
    path = Path(__file__).resolve().parents[1]/'sample/Insole-L2.gcode.3mf'
    if not path.exists():
        pytest.skip('Local sample not present')
    segments, report = inspect_input(path)
    assert len(segments) == 90233
    assert report['segment_counts']['Sparse infill'] == 52514
    assert report['layer_count'] == 40


@pytest.mark.parametrize('failure',['nan','timeout'])
def test_solver_stops_on_invalid_state(tmp_path,monkeypatch,failure):
    import fea_solver
    c = {**DEFAULTS,'input':None,'displacement_mm':0.2,'solver':sys.executable,
         'timeout_seconds':10 if failure == 'timeout' else 300}
    write_deck(Grid(np.ones((2,2,2),dtype=np.uint16),np.zeros(3),0.5,['Outer wall']),c,tmp_path)
    class Process:
        pid = 99999
        args = ['fake']
        returncode = None
        def __init__(self,*args,**kwargs):
            kwargs['stdout'].write('CalculiX Version 2.22,\nline search factor=nan\n')
            kwargs['stdout'].flush()
        def poll(self):
            return self.returncode
        def kill(self):
            self.returncode = -9
        def wait(self):
            return self.returncode
    monkeypatch.setattr(fea_solver.subprocess,'Popen',Process)
    times = iter([0,31,32])
    monkeypatch.setattr(fea_solver.time,'monotonic',lambda: next(times))
    status = solve(c,tmp_path)
    assert status['status'] == 'failed' and status['returncode'] == -9
    assert status.get('timed_out') if failure == 'timeout' else status.get('numerical_failure')


def test_resource_limit(tmp_path):
    c = {**DEFAULTS,'input':None,'displacement_mm':0.2,'solver':sys.executable,'max_solver_dofs':1}
    write_deck(Grid(np.ones((2,2,2),dtype=np.uint16),np.zeros(3),0.5,['Outer wall']),c,tmp_path)
    with pytest.raises(ValueError,match='DOFs'):
        solve(c,tmp_path)
