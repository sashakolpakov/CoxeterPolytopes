#!/usr/bin/env python3
"""Literal spherical-deletion census for higher simplex incidence in the saved triangle dual."""
import argparse,hashlib,json
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank13_triangle_simplex_flags_20261017'
INPUT='fresh_sqrt5_rank13_triangle_discriminant_gate_20261017'
def run():
    data=B.read(INPUT);finite={B.C.code(g):name for name,g in B.C.models() if len(g)<=13}
    records=[];counts=Counter();mass=Counter()
    for r in data['records']:
        g={i:{} for i in range(14)}
        for i,j,m in r['edges']:g[i][j]=g[j][i]=m
        spherical=[];rejected=[]
        for omitted in range(14):
            h={i:{j:m for j,m in row.items() if j!=omitted} for i,row in g.items() if i!=omitted}
            pieces=B.C.pieces(h);codes=[];names=[]
            for p in pieces:
                if sum(map(len,p.values()))!=2*len(p)-2:break
                code=B.C.code(p)
                if code not in finite:break
                codes.append(code);names.append(finite[code])
            else:
                spherical.append(dict(omitted=omitted,forest=sorted(names),unmarked_forest_code=sorted(codes)))
                continue
            rejected.append(omitted)
        S=[s['omitted'] for s in spherical];assert set(r['black'])<=set(S)
        m=len(S);counts[m]+=1;mass[m]+=F(r['multiplier'])
        records.append(dict(face=r['face'],multiplier=r['multiplier'],edges=r['edges'],original_black=r['black'],
            simplex_dimension=m-1,spherical_deletions=spherical,nonspherical_deletions=rejected,
            triangle_root_triples=[list(t) for t in combinations(S,3)],aggregate_required_payment=2*len(list(combinations(S,3)))))
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{B.C.TAG}.py',f'results/{INPUT}.json']
    return dict(status='PASS_EXACT_TRIANGLE_SIMPLEX_FLAG_CENSUS',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        counts_by_simplex_dimension={str(m-1):n for m,n in sorted(counts.items())},
        dual_multiplier_mass_by_simplex_dimension={str(m-1):str(a) for m,a in sorted(mass.items())},records=records,
        scope='The census identifies all spherical one-root deletions in each saved14-root triangle star. Its aggregate-payment proposal requires the written acute simplex-extension proof and a new exact optimization; this census alone neither repairs the weighting rule nor proves a dimension bound.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()))
    out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts_by_simplex_dimension'],d['dual_multiplier_mass_by_simplex_dimension'],flush=True)
