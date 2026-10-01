import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseEvent
import numpy as np
import bambu_fea
from bambu_gcode3mf_to_fea_stl import Segment


def test_roi_click_save_and_cancel(tmp_path,monkeypatch):
    segments = [Segment(np.array([0.,0.,.2]),np.array([100.,100.,.2]),'G1',None,.45,.2,'Outer wall',1)]
    monkeypatch.setattr(bambu_fea,'inspect_input',lambda _: (segments,{'bounds_mm':[[0,0,0],[100,100,8]]}))
    config = {**bambu_fea.DEFAULTS,'input':'fixture','output':'output'}
    path = tmp_path/'roi.json'
    path.write_text(json.dumps(config))
    def choose():
        fig = plt.gcf()
        fig.canvas.draw()
        ax, button = fig.axes
        px,py = ax.transData.transform((50,50))
        event = MouseEvent('button_press_event',fig.canvas,px,py,button=1)
        fig.canvas.callbacks.process('button_press_event',event)
        bx,by = button.transAxes.transform((.5,.5))
        for name in ('button_press_event','button_release_event'):
            fig.canvas.callbacks.process(name,MouseEvent(name,fig.canvas,bx,by,button=1))
    monkeypatch.setattr(plt,'show',choose)
    bambu_fea.select_roi(config,path)
    saved = json.loads(path.read_text())
    assert np.allclose(saved['crop_mm'],[20,80,20,80,0,8])
    assert path.with_suffix('.roi.png').is_file()
    previous = path.read_bytes()
    monkeypatch.setattr(plt,'show',lambda: plt.close('all'))
    bambu_fea.select_roi(config,path)
    assert path.read_bytes() == previous
