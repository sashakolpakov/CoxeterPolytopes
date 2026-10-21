#!/usr/bin/env python3
"""Exact finite first-extension singular-stop screen; no deeper closure claim."""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_global_discount_20261016 as G
import fresh_sqrt5_nonsingular_rank29_fast_20261015 as A
import fresh_sqrt5_sealed_subcore_closure_20261015_review as J
I=A.I;P=A.P;R=A.R;ROOT=G.ROOT
TAG='fresh_sqrt5_global_b1_finite_first_20261016'
INPUT='results/fresh_sqrt5_adjacent_discount_frontier_20261016.json'
RULE='results/fresh_sqrt5_global_discount_20261016.json'

def classify(row):
    assert 0 not in row['labels'];n=row['nodes'];edges=P.edges(row)
    assert all(I.D.finite(len(ids),tuple(edges[tuple(sorted((i,j)))]for i,j in combinations(ids,2)))
               for a,b in product((0,2),(1,3))for ids in[row['white']+[a,b]])
    matrix=I.finite_matrix(row)
    if matrix is None:return dict(status='REJECT_CYCLE_FIELD')
    physical,conjugate,kernel=A.analyze(matrix,row['white'])
    out=dict(physical=physical,conjugate=conjugate,kernel=[[R.enc(x)for x in v]for v in kernel])
    if physical[1]!=1:return dict(**out,status='REJECT_PHYSICAL_SIGNATURE')
    if conjugate[1]:return dict(**out,status='REJECT_CONJUGATE_SIGNATURE')
    if physical[2]>1:return dict(**out,status='REJECT_MULTIPLE_RELATIONS')
    assert conjugate==[n-physical[2],0,physical[2]]
    cyclic=P.cycles(matrix)
    if not all(ok for value,ok in cyclic.values()):return dict(**out,status='REJECT_CYCLIC_INTEGRALITY')
    out['integral_cycles']=len(cyclic)
    if kernel and(all(R.sign(x)>=0 for x in kernel[0])or all(R.sign(x)<=0 for x in kernel[0])):
        return dict(**out,status='REJECT_NONNEGATIVE_ROOT_RELATION')
    if n==23 and not kernel:return dict(**out,status='REJECT_RANK_EXCEEDS_DIMENSION21')
    return dict(**out,status='ACCEPT_SINGULAR_STOP'if kernel else'OPEN_NONSINGULAR')

def payment(row,weights,cache):
    X=G.X;M=G.M;n=row['nodes']
    es=[(a,b,m)for(a,b),m in zip(combinations(range(n),2),row['labels'])]
    compiled=[X.compile_corner(row['white']+[a,b],es,[a,b],[None])for a,b in product((0,2),(1,3))]
    best=None;worst=None;profile=None;cases=0;failed=0
    for ext in X.external_contexts(23-n):
        p,tables=X.add_bonus(compiled,ext)
        for t in tables:p=M.plus(p,t[None])
        key=M.key(p)
        if key not in cache:cache[key]=M.value(p,weights)
        value=cache[key];cases+=1;failed+=value<1
        if best is None or value<best:best=value;worst=ext;profile=p
    return dict(minimum=str(best),common_contexts=cases,failed_contexts=failed,
                worst_external_counts=list(worst),worst_profile=[[j,str(z)]for j,z in sorted(profile.items())])

def run():
    G.install();source=json.loads((ROOT/INPUT).read_text());rule=json.loads((ROOT/RULE).read_text())
    weights=list(map(F,rule['certificate']['primal']));assert weights[0]==0 and weights[1]<F(1,2)
    inputs=[r for r in source['outgoing_roots']if r['analysis']['baseline']==1 and r['parameter_type']=='finite']
    assert len(inputs)==2464
    buckets=defaultdict(list);rows=[];maps=[]
    for k,rec in enumerate(inputs):
        row=rec['row'];mat=J.matrix(row);invariant=J.invariant(mat);found=None
        for j,target in buckets[invariant]:
            mp=J.isomorphism(mat,target)
            if mp is not None:J.check_map(mat,target,mp);found=(j,mp);break
        if found is None:
            j=len(rows);rows.append(dict(index=j,row=row,analysis=rec['analysis'],classification=classify(row)))
            buckets[invariant].append((j,mat));found=(j,list(range(row['nodes'])))
        maps.append(dict(literal_index=k,parent_class=rec['parent_class'],contacts=rec['contacts'],target=found[0],map=found[1]))
        if k%200==0:print('FINITE_FIRST_CLASSIFY',k,'classes',len(rows),flush=True)
    counts=Counter(r['classification']['status']for r in rows);print('FINITE_FIRST_COUNTS',dict(counts),flush=True)
    tested=[];failure=None;cache={}
    ordered=sorted((r for r in rows if r['classification']['status']=='ACCEPT_SINGULAR_STOP'),key=lambda r:(r['row']['nodes'],r['index']))
    for pos,rec in enumerate(ordered):
        pay=payment(rec['row'],weights,cache);tested.append(dict(index=rec['index'],**pay))
        if F(pay['minimum'])<1:
            M=I.finite_matrix(rec['row'])
            failure=dict(**rec,payment=pay,scaled_gram=[[R.enc(x)for x in row]for row in M],
                         origins=[m for m in maps if m['target']==rec['index']])
            print('UNPAID_ADMISSIBLE_SINGULAR_FIRST',rec['index'],rec['row']['nodes'],pay['minimum'],flush=True)
            break
        if pos%100==0:print('FINITE_FIRST_PAID',pos,len(ordered),flush=True)
    paths=[f'scripts/{TAG}.py',INPUT,RULE,f'scripts/{G.TAG}.py',f'scripts/{A.TAG}.py',
           f'scripts/{I.TAG}.py',f'scripts/{J.TAG}.py']
    return dict(status='PASS_EXACT_FINITE_FIRST_EXTENSION_CLASSIFICATION',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        literal_first_extensions=len(inputs),marked_isomorphism_classes=len(rows),counts=dict(counts),
        seed_maps=maps,records=rows,singular_payment_cases=tested,payment_profiles=len(cache),
        all_singular_first_stops_paid=failure is None,first_unpaid_singular=failure,
        scope='All2464finiteB1first extensions are classified and deduplicated with complete opposite-pair-preserving marked maps. Exact separate tuple-Fraction Gram signs, cyclic integrality and root relations are checked. Singular stopping payments are tested until the first failure, if any. An admissible singular Gram is not by itself an independently verified compact quadrilateral or an ambient compact polytope. No second-or-later white-attachment closure and no dimension claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['marked_isomorphism_classes'],result['all_singular_first_stops_paid'],flush=True)
