#!/usr/bin/env python3
"""Target-rank continuous near-envelope with E7/E8 center omissions.

Finite rank coverage is legitimate only conditional on the inherited n<=163
candidate: this scout covers every affine rank through162, not unbounded rank.
"""
import argparse,json
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_nikulin_two_budget_frontier_20261015 as P

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nikulin_center_omission_scout_20261016'
BOUND=162
B=F(20,19)
CAP=F(2)

def classical(typ,r):
    if typ=='C':
        ds=[F(1)]+[F(2)]*(r-1)+[F(1)]
        return ds,lambda i,j:F(abs(i-j),2)
    ds=([F(1),F(1)]+[F(4)]*(r-2)+[F(2)]) if typ=='B' else ([F(1),F(1)]+[F(4)]*(r-3)+[F(1),F(1)])
    def R(i,j):
        if i==j:return F(0)
        if typ=='D' and (i>=r-1 or j>=r-1):
            if i>=r-1 and j>=r-1:return F(1)
            k=min(i,j)
            return F(r if k<=1 else r-k,4)
        if i>j:i,j=j,i
        if j==1:return F(1)
        return F(j if i<=1 else j-i,4)
    return ds,R

def models():
    out=[]
    for r in range(1,BOUND+1):
        m=r+1
        def R(i,j,m=m):
            d=(i-j)%m
            return F(d*(m-d),m)
        out.append(dict(name='A',rank=r,modified=False,om=(0,m//2),ds=[F(1)]*m,R=R))
    for typ in ('B','C','D'):
        for r in range(4 if typ=='D' else 2,BOUND+1):
            ds,R=classical(typ,r)
            out.append(dict(name=typ,rank=r,modified=False,om=(0,r),ds=ds,R=R))
            if typ=='C' or typ=='B' and r>=3:
                out.append(dict(name=typ,rank=r,modified=True,om=(1,r) if typ=='C' else (r-1,0),ds=ds,R=R))
    for model in P.models():
        if model['name'] not in ('E6','E7','E8','F4','G2','B2modified'):continue
        model=dict(model)
        if model['name'] in ('E7','E8') and model.get('forbidden') not in (0,1):
            model['om']=(0,1)
        out.append(model)
    return out

def lower_height(ds,U,V):
    a,b=[sum((P.sqrt_lower(ds[i]) for i in ids),F(0)) for ids in (U,V)]
    return max(min(a,b)**2,max(a,b)**2/B**2)

def domain(model):
    ds,R,om=[model[k] for k in ('ds','R','om')]
    answer=[]
    for o in om:
        cheap=[k for k in range(len(ds)) if k!=o and ds[k]*R(o,k)<=CAP]
        row=[]
        for ids in P.H.subsets(cheap):
            cost=P.support_cost(ids,o,ds,R)
            if cost>CAP:continue
            exact=sum((ds[i]*R(o,i) for i in ids),F(0))
            for i,j in combinations(ids,2):
                cross=R(o,i)+R(o,j)-R(i,j)
                if cross:exact+=P.H.sqrtq(ds[i]*ds[j])*cross
            assert exact==cost
            row.append((ids,cost))
        answer.append(row)
    return answer

def envelope():
    records=[];positive=[];models_seen=0;case_count=0;far=[]
    for model in models():
        name,r,modified,om,ds,R=[model[k] for k in ('name','rank','modified','om','ds','R')]
        T=F(2*(r-int(modified)),9);dom=domain(model);cases=[]
        for (left,p),(right,q) in product(*dom):
            for a,b in product((False,True),repeat=2):
                U=left+((om[0],) if a else ());V=right+((om[1],) if b else ())
                if not U or not V:continue
                e,dual=P.closest(U,V,R);h=lower_height(ds,U,V);loss=T-h*e
                row=dict(left=U,right=V,costs=[str(p),str(q)],minimum_energy=str(e),height_lower=str(h),loss_upper=str(loss),dual=dual)
                cases.append(row);case_count+=1
                if loss>0:
                    assert p+q>=F(1,2)
                    positive.append(dict(model=len(records),case=len(cases)-1,name=name,rank=r,modified=modified,p=p,q=q,loss=loss))
                if name in ('E7','E8') and om==(0,1) and p+q<F(3,2):
                    a,b=[sum((P.H.sqrtq(ds[i]) for i in ids),F(0)) for ids in (U,V)]
                    margin=max(a,b)**2*e+F(4,3)*(p+q)-2
                    far.append(dict(model=len(records),case=len(cases)-1,max_height_squared=str(max(a,b)**2),margin=str(margin)))
        records.append(dict(name=name,rank=r,modified=modified,forbidden=model.get('forbidden'),omissions=om,kernel_squares=list(map(str,ds)),target=str(T),cases=cases))
        models_seen+=1
        if models_seen%100==0:print('models',models_seen,'cases',case_count,flush=True)
    unique={}
    for z in positive:unique.setdefault((z['modified'],z['p'],z['q'],z['loss']),z)
    cs=list(unique.values());table=[]
    for z in cs:
        if any(y is not z and y['modified']==z['modified'] and y['p']<=z['p'] and y['q']<=z['q'] and y['loss']>=z['loss'] and (y['p']<z['p'] or y['q']<z['q'] or y['loss']>z['loss']) for y in cs):continue
        table.append(z)
    for z in positive:
        indices=[i for i,y in enumerate(table) if y['modified']==z['modified'] and y['p']<=z['p'] and y['q']<=z['q'] and y['loss']>=z['loss']]
        assert indices
        records[z['model']]['cases'][z['case']]['dominator']=indices[0]
    @lru_cache(None)
    def best(start,p,q,modified):
        value=F(0);witness=()
        for i in range(start,len(table)):
            z=table[i]
            if z['p']>p or z['q']>q or int(z['modified'])>modified:continue
            v,w=best(i,p-z['p'],q-z['q'],modified-int(z['modified']));v+=z['loss']
            if v>value:value,witness=v,(i,)+w
        return value,witness
    maximum,witness=best(0,F(2),F(2),1)
    print('near maximum',maximum,'far minimum',min(F(z['margin']) for z in far),'table',len(table),flush=True)
    return dict(models=records,counts=dict(models=len(records),domains=case_count,positive=len(positive),frontier=len(table)),pareto_table=table,
        near_budget_packing=dict(maximum=str(maximum),witness=witness,states=best.cache_info().currsize),
        new_far_lemma_cases=far,new_far_lemma_minimum_margin=str(min(F(z['margin']) for z in far)))

def run():
    data=envelope()
    paths=[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py','results/fresh_nikulin_two_budget_frontier_20261015_dimension163.json']
    return dict(status='EXACT_TARGET_RANK_CENTER_OMISSION_SCOUT',sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        affine_rank_limit=BOUND,near_ratio_cap=str(B),endpoint_cap=str(CAP),**data,
        scope='Continuous convex-support bounds for all affine types/ranks up to162 under changed E7/E8 omissions. This is an investigation, not a global dimension theorem. Any general conclusion must bind the preceding163 ceiling, justify all omitted-node choices and both height ranges, and independently review all arithmetic and scope.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=json.loads(json.dumps(run(),default=str));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
