#!/usr/bin/env python3
"""Small exact obstruction in the current uniform-discount payment family."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import fresh_sqrt5_marked_h3_center_cross_20261016 as M
ROOT=M.ROOT;TAG='fresh_sqrt5_uniform_adjacent_discount_barrier_20261016'
QUAD='fresh_sqrt5_adjacent_discount_frontier_20261016'

def run():
    doc=json.loads((ROOT/f'results/{QUAD}.json').read_text())
    bonus_free=[];killers=[]
    for rec in doc['immediate_sealed_weight_obligations']:
        row=rec['row'];white=row['white'];es=[(a,b,m)for(a,b),m in zip(combinations(range(row['nodes']),2),row['labels'])]
        compiled=[M.compile_corner(white+[a,b],es,[a,b],white)for a,b in product((0,2),(1,3))]
        common,tables=M.add_bonus(compiled,(0,)*7)
        corners=[]
        for a,b in product((0,2),(1,3)):
            pcs=M.R.components(M.R.graph(white+[a,b],es));g=next(g for g in pcs if a in g)
            if b in g and len(g)>=6 and g[a].get(b)==3:
                corners.append(dict(pair=[a,b],component=M.R.type_of(g),rank=len(g),
                    root_degrees=[len(g[a]),len(g[b])],pair_orbit=M.R.tree_key(g,[a,b]),
                    component_edges=[[u,v,m]for u in g for v,m in g[u].items()if u<v]))
        for mark in white:
            p=common
            for tab in tables:p=M.plus(p,tab[mark])
            if set(p)<=set([-1,0,1]):
                bonus_free.append(dict(index=rec['class_index'],facet_root=mark,row=row,
                    external_counts=[0]*7,profile=[[j,str(z)]for j,z in sorted(p.items())],
                    discounted_corners=corners))
                if p=={-1:F(1),1:F(-1)}:killers.append(bonus_free[-1])
    assert {r['index']for r in killers}=={4266,7746}
    assert len({c['pair_orbit']for r in killers for c in r['discounted_corners']})==1
    front=dict(context=[1,0,0,0,0,0,0,0],marked_hub='A1',marked_component='A1',facet_root=0,residual_components=['D20'])
    p=M.literal_profile(front)
    assert p[-1]==38 and p[0]==p[1]==-19 and all(z>=0 for j,z in p.items()if j>=2)
    # Vertex profile - T <= 0; add19 times each zero-discount bound.
    columns={j:z for j,z in p.items()if j>=0};columns[M.N]=F(-1)
    columns[0]+=19;columns[1]+=19
    assert columns[0]==columns[1]==0 and all(z>=0 for j,z in columns.items()if j!=M.N)
    bound=F(38)
    paths=[f'scripts/{TAG}.py',f'scripts/{M.TAG}.py',f'scripts/{QUAD}.py',f'results/{QUAD}.json']
    return dict(status='PASS_EXACT_UNIFORM_ADJACENT_DISCOUNT_MODEL_BARRIER',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        bonus_free_rows=bonus_free,zero_discount_witnesses=killers,
        exact_bound=str(bound),marked_vertex=front,marked_vertex_profile=[[j,str(z)]for j,z in sorted(p.items())],
        combined_columns=[[j,str(z)]for j,z in sorted(columns.items())if z],
        exact_combination=dict(vertex_multiplier='1',distance2_zero_bound_multiplier='19',quad4266_multiplier='19',rhs='-38'),
        explanation='Either quad4266 or7746 imposes adjacent_discount<=0 in this family. Nonnegativity forces it0. With distance2_discount fixed0, the marked A1+D20 vertex has baseline38 and only nonnegative bonuses, so every feasible target is at least38.',
        scope='Exact obstruction in the stated current payment family with the two literal quad rows imposed. This does not exclude geometric polytopes. The two rows share one E6 edge orbit; exempting that orbit or adding a payment at another corner can break this small algebraic obstruction. The finite37-row scan and zero-reserved-spectator contexts are bounded diagnostics, not complete paid-quad coverage.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['exact_bound'],len(result['zero_discount_witnesses']),flush=True)
