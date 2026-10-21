#!/usr/bin/env python3
"""One paid quad and two literal marked forests force T>=6042/169."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import fresh_sqrt5_h3_a4_actual_count_split_20261016 as E
M=E.M;ROOT=M.ROOT;TAG='fresh_sqrt5_h3_a3_multiplicity_barrier_20261016'
QUAD='fresh_sqrt5_adjacent_discount_frontier_20261016'

def run():
    E.install()
    data=json.loads((ROOT/f'results/{QUAD}.json').read_text())
    row=next(r['row']for r in data['immediate_sealed_weight_obligations']if r['class_index']==7668)
    white=row['white'];es=[(a,b,m)for(a,b),m in zip(combinations(range(row['nodes']),2),row['labels'])]
    compiled=[M.compile_corner(white+[a,b],es,[a,b],[4])for a,b in product((0,2),(1,3))]
    face,tables=M.add_bonus(compiled,(1,0,4,0,0,0,0))
    for table in tables:face=M.plus(face,table[4])
    bonus=next(j for j in face if j>=2)
    assert face=={-1:F(1),1:F(-1),bonus:F(1,15)}
    fronts=[dict(context=[1,0,0,0,0,0,0,0],marked_hub='A1',marked_component='A1',facet_root=0,residual_components=['D20']),
            dict(context=[3,0,1,0,0,0,0,0],marked_hub='A1',marked_component='A1',facet_root=0,residual_components=['A3']*5)]
    v,w=[M.literal_profile(f)for f in fronts]
    assert v[-1]==38 and v[0]==v[1]==-19 and v.get(bonus,0)==0
    assert w[-1]==18 and w.get(0,0)==w.get(1,0)==0 and w[bonus]==10
    assert all(z>=0 for p in(v,w)for j,z in p.items()if j>=2)
    columns={M.N:F(-1)}
    for p,scale in[(v,F(150,169)),(w,F(19,169))]:
        for j,z in p.items():
            if j>=0:columns[j]=columns.get(j,F(0))+scale*z
    columns[0]+=F(2850,169);columns[1]+=F(2850,169);columns[bonus]-=F(190,169)
    assert columns[0]==columns[1]==columns[bonus]==0
    assert all(z>=0 for j,z in columns.items()if j!=M.N)
    lower=(150*F(38)+19*F(18))/169;assert lower==F(6042,169)>35
    prior=json.loads((ROOT/'results/fresh_sqrt5_h3_a4_actual_count_split_20261016.json').read_text())
    paths=[f'scripts/{TAG}.py',f'scripts/{E.TAG}.py','scripts/fresh_sqrt5_marked_h3_center_cross_20261016.py',f'results/{QUAD}.json',
           'results/fresh_sqrt5_h3_a4_actual_count_split_20261016.json']
    return dict(status='PASS_EXACT_H3_MULTIPLICITY_MODEL_BARRIER',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        lower_bound=str(lower),lower_decimal=float(lower),bonus_index=bonus,bonus_definition=prior['variable_order'][bonus],
        face=dict(index=7668,row=row,external_counts=[1,0,4,0,0,0,0],facet_root=4,profile=[[j,str(z)]for j,z in sorted(face.items())]),
        vertices=[dict(metadata=f,profile=[[j,str(z)]for j,z in sorted(p.items())])for f,p in zip(fronts,[v,w])],
        exact_combination=dict(vertex_multipliers=['150/169','19/169'],quad_multiplier='2850/169',distance2_zero_multiplier='2850/169',
            remaining_columns=[[j,str(z)]for j,z in sorted(columns.items())if z],rhs=str(-lower)),
        scope='Exact obstruction in the current presence-only H3 cross-payment family with distance-two discount0. The face has actual H3 count5 while the second vertex has count1, but both use the same coefficient. It is not a geometric nonexistence theorem. Splitting actual H3 count removes this algebraic coupling; it does not establish successor feasibility or full paid-face coverage.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['lower_bound'],flush=True)
