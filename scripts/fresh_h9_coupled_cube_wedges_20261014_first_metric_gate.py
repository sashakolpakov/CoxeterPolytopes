#!/usr/bin/env python3
"""Exact obstruction to all metric completions of the first label survivor.

Nine actual-vertex affine minors force two Lanner triangles to be orthogonal.
Unassigned meeting angles and dotted distances are not fixed numerically.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_coupled_cube_wedges_20261014_first_metric_gate'
GEOMETRY='results/fresh_h9_coupled_cube_wedges_20261014.json'
LABELS='results/fresh_h9_coupled_cube_wedges_20261014_csp.json'
WITNESSES=[(0,9,53,[0,1,4,7,9,10,14,16,17]),
 (0,10,21,[0,9,10,16,17]),(0,11,23,[0,9,11,16,17]),
 (1,9,53,[0,1,9,10,14,16,17]),(1,10,21,[1,9,10,16,17]),
 (1,11,23,[1,9,11,16,17]),(2,9,149,[0,2,9,10,14,16,17]),
 (2,10,117,[2,9,10,16,17]),(2,11,119,[2,9,11,16,17])]

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    geometry,labels=[json.loads((ROOT/p).read_text()) for p in [GEOMETRY,LABELS]]
    for doc in [geometry,labels]:
        for p,h in doc['sources'].items():assert sha(p)==h,p
    g=geometry['geometry'];vs=g['vertices']
    fixed={(i,j):m for i,j,m in labels['solution']['pair_labels']}
    assert all(m in [3,4] for m in fixed.values())
    dotted={tuple(m) for m in g['minimal_nonfaces'] if len(m)==2}
    records=[]
    for i,j,vi,subset in WITNESSES:
        assert i in subset and j in subset and set(subset)<=set(vs[vi]['incident'])
        assert (i,j) not in fixed and (i,j) not in dotted
        lab=dict(fixed);lab[i,j]=3
        # These witnesses use only simple bonds, so their kernels are rational.
        assert all(lab.get(tuple(sorted((a,b))),2) in [2,3] for a in subset for b in subset if a!=b)
        M=S.Matrix([[S.Integer(1) if a==b else -S.Rational(1,2) if lab.get(tuple(sorted((a,b))),2)==3 else S.Integer(0) for b in subset] for a in subset])
        ns=M.nullspace();assert ns
        v=ns[0]
        if all(x<0 for x in v):v=-v
        assert all(x>0 and x.is_Rational for x in v) and M*v==S.zeros(len(subset),1)
        assert M.det()==0
        # Any nonright label has cos(pi/m)>=1/2. Increasing this or
        # any other omitted contact only decreases v^T G v from zero.
        records.append(dict(pair=[i,j],actual_vertex_index=vi,actual_vertex=vs[vi]['incident'],
                            principal_subset=subset,gram=[[str(x) for x in row] for row in M.tolist()],
                            positive_null_vector=[str(x) for x in v],forced_label=2))
    triangles=[[0,1,2],[9,10,11]]
    for t in triangles:
        assert t in g['minimal_nonfaces']
        assert sorted(fixed[tuple(sorted((i,j)))] for k,i in enumerate(t) for j in t[k+1:])==[3,3,4]
    assert {tuple(r['pair']) for r in records}=={(i,j) for i in triangles[0] for j in triangles[1]}
    q=S.sqrt(2)/2;G=S.Matrix([[1,-S.Rational(1,2),-S.Rational(1,2)],[-S.Rational(1,2),1,-q],[-S.Rational(1,2),-q,1]])
    assert S.simplify(G.det()+S.sqrt(2)/4)==0
    assert G[:2,:2].det()==S.Rational(3,4)
    paths=[GEOMETRY,LABELS,f'scripts/{TAG}.py']
    return dict(status='PASS_FIRST_MANDATORY_ASSIGNMENT_HAS_NO_LORENTZIAN_COMPLETION',
        sources={p:sha(p) for p in paths},forced_right_witnesses=records,
        orthogonal_Lanner_blocks=triangles,each_block_determinant='-sqrt(2)/4',negative_index_lower_bound=2,
        scope='Excludes every completion of the one archived mandatory assignment, for all omitted finite labels and all nine dotted distances. Nine positive null-vector tests at actual vertices force the blocks orthogonal. Other mandatory assignments and the whole incidence remain open.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(path.read_text())
    else:path.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
