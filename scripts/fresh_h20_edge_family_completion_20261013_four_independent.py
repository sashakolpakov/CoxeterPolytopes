#!/usr/bin/env python3
"""Separate arithmetic replay of the94th H20 wall, using the prior auditor."""
import argparse,json
from copy import deepcopy
from fractions import Fraction as Q
import fresh_h20_edge_family_completion_20261013_independent as E
SOURCE='data/enumerations/fresh_h20_edge_family_completion_20261013_four_exact.json'
VOLUME='data/enumerations/fresh_h20_edge_family_completion_20261013_four_volume.json'
OUT='results/fresh_h20_edge_family_completion_20261013_four_independent.json'
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    E.SOURCE=SOURCE;E.VOLUME=VOLUME;data=E.read(SOURCE);out=E.check(data)
    for field in('normal','norm','strict_target'):
        bad=deepcopy(data);r=bad['new_caps'][0]
        if field=='normal':r[field][0][0][0]=str(Q(r[field][0][0][0])+1)
        elif field=='norm':r[field][0][0]='0'
        else:r[field]['new_pairing'][0][0]=str(Q(r[field]['new_pairing'][0][0])+1)
        try:E.check(bad)
        except AssertionError:pass
        else:raise AssertionError('forged '+field+' accepted')
    out['forged_normal_norm_cut_rejected']=True
    if a.verify:assert E.read(OUT)==out
    else:(E.ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_pair_checks'],out['finite_vertices_checked'])
