#!/usr/bin/env python3
"""The16-short/three-long direction creates the94th coherent wall."""
import argparse,json
from fractions import Fraction as F
from fresh_h20_edge_family_completion_20261013_exact import ROOT,read,run
PREVIOUS='data/enumerations/fresh_h20_edge_family_completion_20261013_exact.json'
OUT='data/enumerations/fresh_h20_edge_family_completion_20261013_four_exact.json'
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run([('four_angles',38,F(1))],PREVIOUS)
    if a.verify:assert out==read(OUT)
    else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
