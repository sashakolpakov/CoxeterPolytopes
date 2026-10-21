#!/usr/bin/env python3
"""Exact B4 folded mixed-H5 seed and its extreme-ray polar closure."""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_alternative_mixed_couplings_20261014_b4_fold'
SOURCE='results/fresh_h9_mixed_d4_coupling_20261014.json'

def pair(x,y):return A.pair(x,y,Q)
def add(x,y):return [a+b for a,b in zip(x,y)]
def scale(a,x):return [a*b for b in x]
def gram(N):return [[pair(x,y)for y in N]for x in N]
def reflect(v,b):return add(v,scale(-2*pair(v,b)/pair(b,b),b))

def census(N,h):
    rows=[A.mv(Q,n)for n in N]
    dual=[sum((n[j]for n in N),A.Z)for j in range(10)]
    assert pair(dual,dual).sign()<0 and pair(dual,h).sign()<0
    found={};drop=zero=0
    for ids in combinations(range(len(N)),9):
        v=A.kernel_one([rows[i]for i in ids])
        if v is None:drop+=1;continue
        s=pair(v,dual).sign()
        if not s:zero+=1;continue
        if s>0:v=scale(-1,v)
        vals=[A.dot(row,v)for row in rows]
        if any(x.sign()>0 for x in vals):continue
        v=A.normkey(v);q=pair(v,v)
        active=[i for i,row in enumerate(rows)if not A.dot(row,v)]
        assert pair(v,dual).sign()<0
        found[v]=dict(vector=list(v),norm=q,incident=active,
            kind={-1:'finite',0:'ideal',1:'spacelike'}[q.sign()])
    rays=[found[k]for k in sorted(found,key=lambda v:str(A.enc(v)))]
    counts=dict(Counter(r['kind']for r in rays))
    print('census',len(N),counts,flush=True)
    return dict(tested_supports=sum(1 for _ in combinations(range(len(N)),9)),
        rank_deficient_supports=drop,zero_dual_nullspaces=zero,
        dual_sum=dual,dual_sum_norm=pair(dual,dual),counts=counts,rays=rays)

def polar(N,h,cs):
    out=[];passing=[]
    norms=[pair(n,n)for n in N]
    for rid,row in enumerate(cs['rays']):
        if row['kind']!='spacelike':continue
        n,q=row['vector'],row['norm'];rels=[];bad=[]
        assert pair(n,h).sign()<0
        for j,v in enumerate(N):
            g=pair(n,v);t,m=A.relation(g,q,norms[j])
            rels.append(dict(wall=j,product=g,squared_cosine=g*g/(q*norms[j]),type=t,label=m))
            if t not in('right','finite','parallel','dotted'):bad.append(j)
        if not bad:passing.append(len(out))
        out.append(dict(ray=rid,normal=n,norm=q,relations=rels,bad_old_walls=bad,
            passes=not bad,center_pairing=pair(n,h)))
    mutual=[]
    for i,j in combinations(passing,2):
        a,b=out[i],out[j];g=pair(a['normal'],b['normal'])
        t,m=A.relation(g,a['norm'],b['norm'])
        mutual.append(dict(candidates=[i,j],product=g,type=t,label=m))
    print('polars',len(out),'accepted',passing,'mutual',[(x['type'],x['label'])for x in mutual],flush=True)
    return dict(wall_count=len(N),candidates=out,passing_indices=passing,mutual_pairs=mutual)

def run():
    global Q
    d=json.loads((ROOT/SOURCE).read_text())
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    Q=A.dm(d['form']);old=A.dm(d['seed_normals']);h0=A.dv(d['center'])
    assert A.inertia(Q)==[9,1]
    beta=add(old[10],scale(-1,old[9]));assert pair(beta,beta)==4
    assert reflect(old[9],beta)==old[10] and reflect(old[10],beta)==old[9]
    assert all(reflect(n,beta)==n for n in old[:9])
    N=old[:7]+[old[8],old[7],old[9],beta]
    seed=N[:];h=add(h0,scale(A.E(F(-1,1000)),beta))
    assert pair(h,h).sign()<0 and all(pair(h,n).sign()<0 for n in N)
    support=gram(N[7:]);assert support==A.dm([
        [A.E(x).enc()for x in r]for r in [[2,-1,0,0],[-1,2,-1,0],[0,-1,2,-2],[0,0,-2,4]]])
    # Nine independent spatial walls plus a sixth source wall span the full form.
    A.inverse([A.mv(Q,N[i])for i in list(range(6))+list(range(7,11))])
    stages=[];cs=census(N,h);seedcs=cs
    assert cs['counts']=={'spacelike':20,'finite':10}
    # At most two polar rounds are tested. The second is diagnostic, not presumed complete.
    for k in range(2):
        pr=polar(N,h,cs);stages.append(pr)
        assert all(x['type']in('right','finite','parallel','dotted')for x in pr['mutual_pairs'])
        if not pr['passing_indices']:break
        for i in pr['passing_indices']:N.append(pr['candidates'][i]['normal'])
        cs=census(N,h)
    pairs=[];counts=Counter();facets=[]
    norms=[pair(n,n)for n in N]
    for i,n in enumerate(N):
        f=add(h,scale(-pair(h,n)/norms[i],n))
        assert pair(f,f).sign()<0 and pair(f,h).sign()<0
        fp=[pair(f,v)for v in N]
        assert all(x.sign()==(0 if i==j else -1)for j,x in enumerate(fp))
        facets.append(dict(wall=i,point=f,norm=pair(f,f),all_pairings=fp))
        for j in range(i):
            g=pair(n,N[j]);t,m=A.relation(g,norms[i],norms[j])
            assert t in('right','finite','parallel','dotted')
            counts[t+(str(m)if m else'')]+=1
            pairs.append(dict(pair=[j,i],product=g,squared_cosine=g*g/(norms[i]*norms[j]),type=t,label=m))
    finite_old=[]
    for row in seedcs['rays']:
        if row['kind']!='finite':continue
        vals=[pair(row['vector'],v)for v in N[len(seed):]]
        assert all(x.sign()<0 for x in vals)
        finite_old.append(dict(vector=row['vector'],new_pairings=vals))
    v=next(r['vector']for r in cs['rays']if r['kind']=='spacelike')
    for k in range(1,40):
        e=add(v,scale(A.E(F(1,2**k)),h));ev=[pair(e,n)for n in N]
        if pair(e,e).sign()>0 and pair(e,h).sign()<0 and all(x.sign()<0 for x in ev):break
    else:raise AssertionError('no strict escape perturbation')
    # The mixed-source cyclic-field obstruction is unchanged under folding.
    G=A.dm(d['nonquasiarithmetic_obstruction']['fixed_principal_gram'])
    cycle=pair(N[5],N[6]);cycle=cycle*cycle
    assert cycle.b and cycle.sigma()!=cycle and A.inertia(G)==[4,1]
    assert all(x.sigma()==x for row in G for x in row)
    assert len(N)==12 and A.normkey(N[11])==A.normkey(A.dv(d['final_normals'][11]))
    assert not pair(N[11],beta)
    files=[SOURCE,'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',f'scripts/{TAG}.py']
    return A.enc(dict(status='PASS_EXACT_B4_FOLD_POLAR_EXPERIMENT_INFINITE_VOLUME',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in files},form=Q,
        center=h,center_norm=pair(h,h),reflection_normal=beta,support_gram=support,
        seed_normals=seed,seed_complete_ray_census=seedcs,polar_rounds=stages,
        final_normals=N,final_complete_ray_census=cs,final_pairs=pairs,final_pair_counts=dict(counts),
        actual_facet_witnesses=facets,old_finite_vertex_preservation=finite_old,
        strict_escape=dict(vector=e,norm=pair(e,e),center_pairing=pair(e,h),all_pairings=ev,perturbation_denominator=2**k),
        conditional_nonquasiarithmetic=dict(cycle=cycle,sigma_cycle=cycle.sigma(),fixed_indefinite_block=G),
        inherited_D4_polar=True,polar_bisector_pairing=pair(N[11],beta),
        exhausted_polar_closure=not stages[-1]['passing_indices'],
        scope='Exact B4 half of the fixed mixed D4 seed and at most two simultaneous coherent extreme-ray polar rounds. Infinite volume certified by strict spacelike escape; no arbitrary-cap exclusion.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text());print('SAVED_REPLAY_PASS')
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],len(out['final_normals']),out['final_complete_ray_census']['counts'])
