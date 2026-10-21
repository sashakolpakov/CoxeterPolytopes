#!/usr/bin/env python3
"""Repair the H3-center proposal with the complete affected-triangle domain.

The 37 added paid quadrilaterals are still only a partial quad domain.
The old overcomplete selected triangles are replaced, not retained.
"""
import argparse
import json
import os
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_marked_h3_center_sparse_recovery_20261016 as S

M=S.M;ROOT=M.ROOT
TAG='fresh_sqrt5_h3_center_complete_triangle_repair_20261016'
BASE_TAG='fresh_sqrt5_marked_h3_center_cross_20261016'
TRI='fresh_sqrt5_triangle_extension_complete_20261016'
QUAD='fresh_sqrt5_adjacent_discount_frontier_20261016'
INITIAL=f'data/enumerations/{TAG}_initial.json'
CHECKPOINT=Path('/private/tmp/'+TAG+'_checkpoint.json')
OLD_FACES=M.faces

def faces():
    inherited,oldcounts=OLD_FACES()
    obligations={(M.key(p),target):meta for p,target,meta in inherited if meta['kind']=='quad'}
    tri=json.loads((ROOT/f'results/{TRI}_rank22.json').read_text())
    quad=json.loads((ROOT/f'results/{QUAD}.json').read_text())
    assert len(tri['records'])==996 and tri['complete_through_bound'] and not tri['pending']
    todo=[('triangle',r['index'],r['row'],r['status'])for r in tri['records']]
    todo += [('quad',r['class_index'],r['row'],'SEALED_SINGULAR')for r in quad['immediate_sealed_weight_obligations']]
    todo += [('quad',r['class_index'],r['row'],'SINGULAR_EXTENSIBLE_STOP')for r in quad['records']if r['class_index']in(4446,4448)]
    assert sum(kind=='quad'for kind,_,_,_ in todo)==37
    counts=dict(inherited_baseline_zero_quad_rows=oldcounts['quad_rows'],
                inherited_quad_common_contexts=oldcounts['quad_common_contexts'],
                triangle_rows=996,paid_quad_stopping_rows=37,
                triangle_marked_contexts=0,paid_quad_marked_contexts=0)
    for pos,(kind,index,row,status)in enumerate(todo):
        n=row['nodes'];white=row['white'];remaining=(22 if kind=='triangle'else 23)-n
        edges=[(a,b,m)for(a,b),m in zip(combinations(range(n),2),row['labels'])]
        pairs=[[b for b in row['black']if b!=o]for o in row['black']]if kind=='triangle'else list(product((0,2),(1,3)))
        marks=white+([-1,-2,-3]if remaining else[])
        compiled=[M.compile_corner(white+list(pair),edges,pair,marks)for pair in pairs]
        for ext in M.contexts(remaining):
            common,tables=M.add_bonus(compiled,ext)
            valid=white+([-1]if remaining-4*ext[3]-ext[0]>0 else[])+([-2]if ext[3]>0 else[])+([-3]if ext[0]>0 else[])
            for mark in valid:
                p=common
                for table in tables:p=M.plus(p,table[mark])
                obligations.setdefault((M.key(p),F(2 if kind=='triangle'else 1)),dict(
                    kind=kind,family='complete_triangles'if kind=='triangle'else'paid_quad_stopping',
                    index=index,producer_status=status,external_counts=list(ext),facet_root=mark))
                counts['triangle_marked_contexts'if kind=='triangle'else'paid_quad_marked_contexts']+=1
        if pos%100==0:print('FACE_REPAIR',pos,'profiles',len(obligations),flush=True)
    # Explicit rows make these modeling restrictions part of exact weak duality.
    # They also justify inherited hereditary-safe pruning in the full catalogue.
    obligations[M.key({0:F(-1)}),F(0)]=dict(kind='model_bound',name='distance2_discount_zero')
    obligations[M.key({-1:F(2,3),1:F(-1)}),F(0)]=dict(kind='model_bound',name='adjacent_discount_at_most_two_thirds')
    return [(dict(p),target,meta)for(p,target),meta in obligations.items()],counts

def checkpoint(vertices,steps,weights):
    paths=[f'scripts/{BASE_TAG}.py',f'scripts/{TAG}.py',
           'scripts/fresh_sqrt5_marked_h3_center_sparse_recovery_20261016.py']
    doc=dict(status='PARTIAL_CHECKPOINT_NOT_A_CERTIFICATE',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
             vertex_cuts=vertices,steps=steps,latest_numeric_weights=list(map(float,weights))if weights is not None else None)
    temporary=CHECKPOINT.with_suffix('.tmp');temporary.write_text(json.dumps(doc,separators=(',',':'))+'\n');os.replace(temporary,CHECKPOINT)

def run(saved=None):
    S.TAG=TAG;S.CHECKPOINT=CHECKPOINT;S.RESUME=True;S.checkpoint=checkpoint
    M.TAG=TAG;M.faces=faces;M.solve=S.solve
    result=M.run(saved)
    paths=[f'scripts/{BASE_TAG}.py','scripts/fresh_sqrt5_marked_h3_center_sparse_recovery_20261016.py',
           f'scripts/{TRI}.py',f'results/{TRI}_rank22.json',f'scripts/{QUAD}.py',f'results/{QUAD}.json',INITIAL]
    result['sources'].update({p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths})
    result['status']='PASS_EXACT_COMPLETE_TRIANGLE_PARTIAL_PAID_QUAD_PROPOSAL'
    result['scope']='Same nonnegative bonus family, distance-two discount fixed0 and adjacent discount at most2/3. All996 states of the complete affected-triangle catalogue replace the old196 selected triangles, including removal of their conjugate-inadmissible D5-root3 row. All inherited baseline-zero quad profiles and37 newly affected paid-quad stopping rows are enforced. Complete marked finite-forest budgets are retained. Further paid-quadrilateral extension coverage remains open. No dimension-bound promotion.'
    result['optimizer_scope']='T fixed35, minimize sum of all nonnegative bonus coefficients. If infeasible, recover an exact restricted-model lower-bound dual from the same literal constraints. The objective is bounded below0; every bonus is bounded on each feasible objective sublevel. Inherited vertex cuts are exact valid marked-forest rows, independent of the face-domain change.'
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text())if args.verify else None
    if not args.verify and not CHECKPOINT.exists():CHECKPOINT.write_bytes((ROOT/INITIAL).read_bytes())
    result=json.loads(json.dumps(run(saved),default=str))
    if args.verify:assert result==saved
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    proof=result['verification'];number=proof.get('maximum',proof.get('lower_bound'))
    print(result['status'],result['certificate']['kind'],float(F(number)),flush=True)
