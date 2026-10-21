#!/usr/bin/env python3
"""Exact finite budgets conditional on the new rank16/17 and D14-role rule."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_e6_d11_quarter_coupon_forest_budget_20261017 as Q

B=Q.B;ROOT=Q.ROOT;TAG='fresh_d16_d17_d14_coupon_forest_budget_20261017'


def base(t):
    kind,r=t
    if r==16:return {'A':232,'B':208,'D':236}[kind]
    if r==17:return {'A':228,'B':210,'D':250}[kind]
    return B.base(t)


def cross(a,b):
    out=Q.cross(a,b)
    for small,large in ((a,b),(b,a)):
        if small in (('A',1),('A',2)) and large==('D',14):out+=small[1]
    return out


def rep(t):return Q.rep(t)
def rtype(r):return rep(('A' if r<=9 else 'D',r))
def actual(ts):return sum(map(base,ts))+sum(cross(a,b) for i,a in enumerate(ts) for b in ts[i+1:])


def dominance():
    types=([('A',r) for r in range(1,76)]+[('B',r) for r in range(2,76)]+
           [('D',r) for r in range(4,76)]+[('E',6),('E',7),('E',8),('F',4),('H',3),('H',4),('I2*',2)])
    rows=[]
    for t in types:
        a=rep(t);assert a in types and a[1]==t[1] and base(a)>=base(t)
        deltas=[]
        for u in types:
            d=cross(a,u)-cross(t,u);assert d>=0,(t,a,u)
            if d:deltas.append([B.label(u),d])
        rows.append(dict(original=B.label(t),representative=B.label(a),doubled_base_gain=base(a)-base(t),positive_cross_gains=deltas))
    return dict(types=len(types),partner_checks=len(types)**2,rows=rows)


def domain():
    out=[]
    for old in B.count_domain():
        c=old['component_counts'];large=Q.large_types(old)
        key=old['interaction_key']+[c[4]]
        out.append(dict(component_counts=c,used_rank=old['used_rank'],interaction_key=key,doubled_base_cost=actual(large)))
    return out


def costs(key,small):
    large_ranks=B.RESIDUAL if not small else tuple(range(1,10))+(17,)
    broad,upgrade,a4,d11,d14=key
    out={}
    for r in large_ranks:
        out[r]=base(rtype(r))
        if r<=5:out[r]+=r*broad
        if r==1:out[r]+=2*upgrade+d14
        if r==2:out[r]+=4*d11+2*d14
        if r==4:out[r]+=2*a4
        if r==6:out[r]+=2*d11
    return out


def table(key,small):
    c=costs(key,small);v=[0]+[None]*75
    for r,w in c.items():
        for n in range(r,76):
            if v[n-r] is not None:
                z=v[n-r]+w
                if v[n] is None or z>v[n]:v[n]=z
    assert all(x is not None for x in v)
    return dict(key=list(key),small_components_only=small,values=v,
        last=[0]+[max(r for r,w in c.items() if r<=n and v[n]==v[n-r]+w) for n in range(1,76)])


def unpack(t,n):
    out=[]
    while n:
        r=t['last'][n];assert 0<r<=n;out.append(rtype(r));n-=r
    return out


def run(verify=False):
    old=json.loads((ROOT/f'results/{Q.TAG}.json').read_text())
    for p,h in old['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    ds=domain();keys=sorted({tuple(r['interaction_key']) for r in ds});path=ROOT/f'results/{TAG}.json'
    saved=json.loads(path.read_text()) if verify else None
    all_tables=[];results=[];checks=dict(interaction_formulas=0,bellman_upper=0,attaining_predecessors=0,global_alternatives=0,typed_global_attainers=0)
    for small in (False,True):
        if verify:
            tt=[t for t in saved['tables'] if t['small_components_only']==small]
            tables={tuple(t['key']):t for t in tt};assert len(tables)==len(tt) and set(tables)==set(keys)
        else:tables={k:table(k,small) for k in keys}
        for r in ds:
            large=Q.large_types(r);k=tuple(r['interaction_key'])
            assert actual(large)==r['doubled_base_cost']
            for rank,w in costs(k,small).items():
                t=rtype(rank);assert w==base(t)+sum(cross(t,u) for u in large);checks['interaction_formulas']+=1
        for k,t in tables.items():
            c=costs(k,small);v=t['values'];last=t['last'];assert len(v)==len(last)==76 and v[0]==last[0]==0
            for n in range(1,76):
                assert isinstance(v[n],int)
                for r,w in c.items():
                    if r<=n:assert v[n]>=v[n-r]+w;checks['bellman_upper']+=1
                r=last[n];assert r in c and r<=n and v[n]==v[n-r]+c[r];checks['attaining_predecessors']+=1
        answers=[]
        for n in range(76):
            alternatives=[]
            for r in ds:
                if r['used_rank']>n:continue
                t=tables[tuple(r['interaction_key'])];left=n-r['used_rank']
                value=r['doubled_base_cost']+t['values'][left]
                forest=B.canonical(Q.large_types(r)+unpack(t,left))
                assert sum(r for _,r in forest)==n and actual(forest)==value
                if small:assert max((r for _,r in forest),default=0)<=17
                alternatives.append((value,tuple(map(B.label,forest))))
            value,forest=max(alternatives);checks['global_alternatives']+=len(alternatives);checks['typed_global_attainers']+=1
            if not small:assert value<=old['dimension_maxima'][n]['doubled_cost']
            answers.append(dict(dimension=n,doubled_cost=value,forest=forest,curvature=str(B.curvature(n)),slack=str(B.curvature(n)-F(value,2))))
        results.append(dict(small_components_only=small,dimension_maxima=answers))
        all_tables.extend(tables[k] for k in keys)
    assert results[0]['dimension_maxima'][65]['doubled_cost']==974
    assert results[1]['dimension_maxima'][65]['doubled_cost']==958
    assert results[1]['dimension_maxima'][65]['forest']==('D16','D16','D16','D17')
    paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}' for s in(Q.TAG,B.TAG,'fresh_finite_d17_long_pair_support_20261017')
        for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    out=dict(status='PASS_D16_D17_D14_ROLE_COUPON_EXACT_FINITE_BUDGET',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        conditional_costs=dict(A16=116,B16=104,D16=118,A17=114,B17=105,D17=125),
        new_coupon='Weight1/2 between every mark of whole A1 or A2 and the unique long-arm neighbor of the fork in separate whole D14; added totals1/2 and1 per component pair respectively.',
        geometric_coverage_proved_here=False,type_domination=dominance(),large_domain=ds,tables=all_tables,results=results,checks=checks,
        large_vectors=len(ds),distinct_interaction_keys=len(keys),maximum_rank=75,
        max65_all='487',max65_components_at_most17='479',curvature65='480',
        scope='Complete exact conditional typed finite-forest counts through75 with all inherited coupons plus the stated D14-role coupon and rank16/17 costs. Two domains: arbitrary components and every component of rank<=17. Generation uses item-first DP; replay checks all last-item upper inequalities, equality predecessors and typed global attainers. Geometric justification and ideal composition are separate; no ambient theorem asserted here.')
    out=json.loads(json.dumps(out))
    if verify:assert out==saved
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['large_vectors'],out['distinct_interaction_keys'],out['checks'],out['max65_all'],out['max65_components_at_most17'],flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();run(args.verify)
