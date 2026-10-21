#!/usr/bin/env python3
"""Exact payment diagnostic for a preserved decimal numerical vector.

This does not certify its vertex budget or exhaust the full face domain.
Every exported literal face is tested over all shared spectator contexts and
all possible facet-root categories using the frozen proposed payment rule.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import fresh_sqrt5_marked_h3_center_cross_20261016 as M

ROOT=M.ROOT
TAG='fresh_sqrt5_h3_center_new_face_screen_20261016'
INPUTS=[
    'scripts/fresh_sqrt5_marked_h3_center_cross_20261016.py',
    'data/enumerations/fresh_sqrt5_marked_h3_center_secondary_20261016_checkpoint.json',
    'results/fresh_sqrt5_triangle_extension_frontier_20261016_rank9.json',
    'results/fresh_sqrt5_adjacent_discount_frontier_20261016.json',
]

def run():
    checkpoint=json.loads((ROOT/INPUTS[1]).read_text())
    x=list(map(lambda z:F(str(z)),checkpoint['latest_numeric_weights'][:M.N]))
    assert x[0]==0 and 0<=x[1]<=F(2,3)
    tri=json.loads((ROOT/INPUTS[2]).read_text())
    quad=json.loads((ROOT/INPUTS[3]).read_text())
    todo=[('triangle',r['index'],r['row'],r['status'])for r in tri['records']
          if r['status']in('EXPANDED','OPEN_RANK_BOUNDARY')]
    todo += [('quad',r['class_index'],r['row'],'SEALED_SINGULAR')
             for r in quad['immediate_sealed_weight_obligations']]
    todo += [('quad',r['class_index'],r['row'],'SINGULAR_EXTENSIBLE_ZERO_TOUCH')
             for r in quad['records']if r['class_index']in(4446,4448)]
    rows=[];total=0
    for kind,index,row,status in todo:
        n=row['nodes'];white=row['white'];remaining=(22 if kind=='triangle'else 23)-n
        es=[(a,b,m)for(a,b),m in zip(combinations(range(n),2),row['labels'])]
        pairs=[[b for b in row['black']if b!=o]for o in row['black']]if kind=='triangle'else list(product((0,2),(1,3)))
        marks=white+([-1,-2,-3]if remaining else[])
        compiled=[M.compile_corner(white+list(pair),es,pair,marks)for pair in pairs]
        target=F(2 if kind=='triangle'else 1);minimum=None;violations=0;cases=0
        for ext in M.contexts(remaining):
            common,tables=M.add_bonus(compiled,ext)
            valid=white+([-1]if remaining-4*ext[3]-ext[0]>0 else[])+([-2]if ext[3]>0 else[])+([-3]if ext[0]>0 else[])
            for mark in valid:
                p=common
                for table in tables:p=M.plus(p,table[mark])
                paid=M.value(p,x);cases+=1
                if paid<target:violations+=1
                if minimum is None or paid<minimum[0]:minimum=(paid,ext,mark,p)
        assert minimum is not None
        paid,ext,mark,p=minimum;total+=cases
        rows.append(dict(kind=kind,index=index,producer_status=status,row=row,
            target=str(target),minimum=str(paid),minimum_decimal=float(paid),
            checked_context_marks=cases,violating_context_marks=violations,
            witness=dict(external_counts=list(ext),facet_root=mark,
                         profile=[[j,str(z)]for j,z in sorted(p.items())])))
        print(kind,index,'minimum',float(paid),'violating',violations,'cases',cases,flush=True)
    return dict(status='PASS_EXACT_DECIMAL_VECTOR_FACE_DIAGNOSTIC',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in INPUTS+[f'scripts/{TAG}.py']},
        vector_definition='Each preserved numerical coefficient is interpreted as the exact rational given by its JSON decimal representation. This is not an exact feasible primal certificate.',
        discounts=list(map(str,x[:2])),checked_rows=len(rows),checked_context_marks=total,
        violating_rows=sum(r['violating_context_marks']>0 for r in rows),records=rows,
        scope='All exported expanded/open rank9 triangle rows and all35 sealed paid quads plus two singular-extensible stopping rows; all common spectator counts and facet-root categories. This bounded diagnostic is not complete face coverage and does not certify the vertex budget.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',result['checked_rows'],result['checked_context_marks'],result['violating_rows'],flush=True)
