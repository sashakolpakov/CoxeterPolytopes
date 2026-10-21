#!/usr/bin/env python3
"""Literal maximizing forests and two exact failed faces of the selected rule."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import fresh_sqrt5_marked_h3_center_cross_20261016 as M
ROOT=M.ROOT;TAG='fresh_sqrt5_h3_center_selected_diagnostic_20261016'
RULE='fresh_sqrt5_marked_h3_center_sparse_recovery_20261016'
TRI='fresh_sqrt5_triangle_extension_frontier_20261016_rank9'
QUAD='fresh_sqrt5_adjacent_discount_frontier_20261016'

def run():
    doc=json.loads((ROOT/f'results/{RULE}.json').read_text());x=list(map(F,doc['certificate']['primal']))
    assert doc['certificate']['kind']=='exact_feasible_restricted_rule'and F(doc['verification']['maximum'])==35
    assert x[0]==0 and 0<x[1]<F(2,3)
    active=[]
    for front in doc['verification']['fronts']:
        if F(front['maximum'])!=35:continue
        profile=M.literal_profile(front);assert M.value(profile,x)==35
        counts=Counter(dict(zip(M.RESERVED,front['context'][:7])))
        name=front['marked_component'];mark=front['facet_root']
        if name in M.RESERVED:counts[name]-=1
        counts.update(front['residual_components']);counts[name]+=1
        assert sum(len(M.MODELS[n])*c for n,c in counts.items())==21
        active.append(dict(**front,forest={n:c for n,c in sorted(counts.items())if c},
            marked_root_key=M.R.tree_key(M.MODELS[name],[mark]),
            marked_component_edges=[[u,v,m]for u in M.MODELS[name]for v,m in M.MODELS[name][u].items()if u<v],
            profile=[[j,str(z)]for j,z in sorted(profile.items())]))
    assert len(active)==306
    tri=json.loads((ROOT/f'results/{TRI}.json').read_text())
    quad=json.loads((ROOT/f'results/{QUAD}.json').read_text())
    rows=[('triangle',136,tri['records'][136]['row'],[2,1,1,0,0,0,0],3),
          ('quad',4202,next(r['row']for r in quad['immediate_sealed_weight_obligations']if r['class_index']==4202),[0,2,0,0,0,0,0],4)]
    failures=[]
    for kind,index,row,ext,mark in rows:
        n=row['nodes'];white=row['white'];es=[(a,b,m)for(a,b),m in zip(combinations(range(n),2),row['labels'])]
        pairs=[[b for b in row['black']if b!=o]for o in row['black']]if kind=='triangle'else list(product((0,2),(1,3)))
        compiled=[M.compile_corner(white+list(pair),es,pair,[mark])for pair in pairs]
        common,tables=M.add_bonus(compiled,ext);profile=common
        for table in tables:profile=M.plus(profile,table[mark])
        target=F(2 if kind=='triangle'else 1);paid=M.value(profile,x)
        assert paid==target-x[1]<target
        failures.append(dict(kind=kind,index=index,row=row,external_counts=ext,facet_root=mark,
            target=str(target),payment=str(paid),payment_decimal=float(paid),shortfall=str(x[1]),
            profile=[[j,str(z)]for j,z in sorted(profile.items())],
            variables=[dict(index=j,definition=doc['variable_order'][j],weight=str(x[j]),coefficient=str(z))for j,z in sorted(profile.items())if j>=0]))
    paths=[f'scripts/{TAG}.py',f'scripts/{M.TAG}.py',f'scripts/{RULE}.py',f'results/{RULE}.json',f'results/{TRI}.json',f'results/{QUAD}.json']
    return dict(status='PASS_EXACT_SELECTED_RULE_MAXIMIZERS_AND_FAILED_FACES',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        dimension=21,facet_dimension=20,selected_model_maximum='35',discounts=list(map(str,x[:2])),
        maximizing_marked_forests=len(active),marked_component_distribution=dict(sorted(Counter(r['marked_component']for r in active).items())),
        maximizing_witnesses=active,exact_failed_faces=failures,
        scope='306 chosen maximizing DP witnesses, one per active context/mode front; not an enumeration of every tied maximizer. Literal root-pair recounts verify all306 payments exactly35. Two new literal faces fail by the exact adjacent discount; local arithmetic realizability is inherited or checked separately, not proved here. The selected-domain rule is therefore not a full compact21 exclusion. No issue is implied for the older compact<=21 certificate, which has no discounts.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['maximizing_marked_forests'],len(result['exact_failed_faces']),flush=True)
