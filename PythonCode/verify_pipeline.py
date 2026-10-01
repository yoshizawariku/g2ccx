"""Real-solver checks. Writes fresh output directories; never overwrites results."""
from __future__ import annotations
import argparse
from pathlib import Path
import json
import numpy as np
from bambu_fea import DEFAULTS
from fea_geometry import Grid
from fea_solver import write_deck, solve, postprocess


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--pitch",type=float,default=0.2)
    args = ap.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        ap.error("output must be empty")
    n = round(2/args.pitch)
    if n < 2 or not np.isclose(n*args.pitch,2):
        ap.error("pitch must divide the 2 mm block")
    grid = Grid(np.ones((n,n,n),dtype=np.uint16),np.zeros(3),args.pitch,['Outer wall'])
    config = {**DEFAULTS,'input':None,'output':str(args.output.resolve()),
              'voxel_mm':args.pitch,'displacement_mm':0.2,'timeout_seconds':300}
    write_deck(grid,config,args.output)
    solve(config,args.output)
    result = postprocess(args.output)
    print(json.dumps(result,indent=2))
    if result['status'] != 'complete' or not result['force_balance_pass'] or result['last']['compression_force_n'] <= 0:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
