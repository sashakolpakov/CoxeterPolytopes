#!/usr/bin/env python3
"""Exact rank64 forest budget with enlarged E6/D11 and A1/B10 coupons.

Standard-library arithmetic only. Generation uses item-first dynamic programming;
verification checks saved all-last-item upper bounds and attaining predecessors.
The enlarged payment rule is explicit; no ambient dimension claim is made.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_e6_b10_half_coupon_forest_budget_20261018'
BASE='fresh_ambient64_finite_closure_20261018'
FRONT='fresh_ambient64_payment_frontier_20261018'
LOCAL='fresh_order17_e6_half_coupon_transfer_20261018'
N=64

def representative(r):return ('E',6)if r==6 else('A'if r<=9 else'D',r)

def base(t):
    k,r=t
    if r==16:return {'A':232,'B':208,'D':232}[k]
    if r==17:return {'A':228,'B':210,'D':254}[k]
    if r<=15:return r*(r-1)
    old=16*r-(66 if k=='D'and r>=20 else 64 if k=='D'else 68)
    deep=max(0,r-(22+{'A':8,'B':9,'D':10}[k]))
    return old-2*deep

def cross(a,b):
    result=0
    for s,l in((a,b),(b,a)):
        if s[1]<=5:
            if l==('D',10):result+=2*s[1]
            elif l[0]=='D'and 11<=l[1]<=16:result+=3*s[1]
            elif l in(('A',13),('A',15)):result+=s[1]
            elif l==('A',17):result+=2*s[1]
        if(s,l)in((('A',1),('D',11)),(('A',1),('D',12)),(('A',2),('B',10)),
            (('A',4),('D',13)),(('A',1),('A',11)),(('A',2),('A',11))):result+=2
        if(s,l)==(('A',2),('D',11)):result+=4
        if(s,l)==(('E',6),('D',11)):result+=4
        if(s,l)==(('A',1),('B',10)):result+=1
        if s in(('A',1),('A',2))and l==('D',14):result+=s[1]
    return result

def domain():
    rows=[]
    for count in range(N//10+1):
        for rs in combinations_with_replacement(range(10,18),count):
            if sum(rs)>N:continue
            c=Counter(rs);key=(2*c[10]+3*sum(c[r]for r in range(11,17))+2*c[17],c[11]+c[12],c[13],c[11],c[14])
            rows.append(dict(counts=[c[r]for r in range(10,18)],used_rank=sum(rs),key=list(key)))
    return sorted(rows,key=lambda r:r['counts'])

def items(key):
    broad,extra,a4,d11,d14=key;out={}
    for r in tuple(range(1,10))+tuple(range(17,N+1)):
        w=base(representative(r))+(r*broad if r<=5 else 0)
        w+=2*extra+d14 if r==1 else 4*d11+2*d14 if r==2 else 2*a4 if r==4 else 4*d11 if r==6 else 0
        out[r]=w
    return out

def generate_table(key):
    weights=items(key);v=[0]+[None]*N
    for r,w in weights.items():
        for n in range(r,N+1):
            if v[n-r]is not None:
                z=v[n-r]+w
                if v[n]is None or z>v[n]:v[n]=z
    last=[0]+[max(r for r,w in weights.items()if r<=n and v[n]==v[n-r]+w)for n in range(1,N+1)]
    return dict(key=list(key),values=v,last=last)

def verify_table(table,counts):
    key=tuple(table['key']);w=items(key);v,last=table['values'],table['last']
    assert len(v)==len(last)==N+1 and v[0]==last[0]==0
    assert all(type(x)is int for x in v+last)
    for n in range(1,N+1):
        for r,c in w.items():
            if r<=n:assert v[n]>=v[n-r]+c;counts['bellman_upper']+=1
        r=last[n];assert r in w and r<=n and v[n]==v[n-r]+w[r];counts['attaining_predecessors']+=1

def actual_cost(forest):return sum(base(t)for t in forest)+sum(cross(a,b)for i,a in enumerate(forest)for b in forest[i+1:])

def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in(BASE,FRONT,LOCAL)}
    for tag,d in docs.items():
        assert d['status'].startswith('PASS')
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    assert docs[BASE]['finite_maximum']=='464'and docs[BASE]['finite_budget_hypotheses_discharged']
    assert docs[BASE]['simultaneous_rules']['D16']==116 and docs[BASE]['simultaneous_rules']['D17']==127
    assert docs[LOCAL]['new_cost']=='118'and docs[LOCAL]['finite_coupon']['new_total']=='2'
    types=([('A',r)for r in range(1,N+1)]+[('B',r)for r in range(2,N+1)]+[('D',r)for r in range(4,N+1)]+
        [('E',6),('E',7),('E',8),('F',4),('H',3),('H',4),('I2*',2)])
    domination=[]
    for t in types:
        rep=t if t==('A',17)else representative(t[1]);assert base(rep)>=base(t)
        gaps=[cross(rep,u)-cross(t,u)for u in types];assert min(gaps)>=0
        domination.append(dict(type=list(t),representative=list(rep),base_gap=base(rep)-base(t),minimum_partner_gap=min(gaps)))
    assert len(types)==195 and cross(('D',10),('A',1))==2>=cross(('B',10),('A',1))==1
    assert cross(('D',10),('A',2))==4>=cross(('B',10),('A',2))==2
    dom=domain();assert len(dom)==753
    keys=sorted({tuple(r['key'])for r in dom});assert len(keys)==348
    if saved is None:tables={k:generate_table(k)for k in keys}
    else:
        assert [r['key']for r in saved['tables']]==list(map(list,keys))
        tables={tuple(r['key']):r for r in saved['tables']}
    counts=Counter()
    for table in tables.values():verify_table(table,counts)
    alternatives=[]
    for row in dom:
        large=[('A'if r==17 else'D',r)for r,c in zip(range(10,18),row['counts'])for _ in range(c)]
        key=tuple(row['key']);table=tables[key];weights=items(key)
        for r,w in weights.items():
            assert base(representative(r))+sum(cross(representative(r),u)for u in large)==w;counts['literal_interactions']+=1
        rem=N-row['used_rank'];value=sum(base(t)for t in large)+table['values'][rem];forest=list(large)
        while rem:r=table['last'][rem];forest.append(representative(r));rem-=r
        assert sum(r for k,r in forest)==N and actual_cost(forest)==value
        alternatives.append(dict(counts=row['counts'],doubled_upper=value,attaining_forest=sorted([k+str(r)for k,r in forest])))
        counts['global_alternatives']+=1
    maximum=max(r['doubled_upper']for r in alternatives);attainers=[r for r in alternatives if r['doubled_upper']==maximum]
    assert maximum==928 and any(r['attaining_forest']==['D16']*4 for r in attainers)
    old=next(r for r in docs[FRONT]['finite_conditional_scenarios']if(r['feature_depth'],r['deep_pair_weight'],r['D16_drop'],r['D17_drop'])==(11,'0',2,-2))
    assert old['exact_maximum']=='464'and old['checks']==dict(counts)
    paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}'for s in(BASE,FRONT,LOCAL)for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_RANK64_E6_B10_HALF_COUPON_BUDGET',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        ambient_dimension=64,finite_maximum='464',finite_curvature='464',finite_slack='0',attaining_forest=['D16']*4,
        added_coupons=[dict(types=['E6','D11'],old_component_pair_cost='1',new_component_pair_cost='2',increment='1',role_source=LOCAL),
            dict(types=['A1','B10'],old_component_pair_cost='0',new_component_pair_cost='1/2',increment='1/2',marked_B10_role='doubled terminal')],
        domination=domination,type_partner_comparisons=len(types)**2,large_domain=dom,tables=[tables[k]for k in keys],
        alternatives=alternatives,checks=dict(counts),old_maximum_preserved=True,
        representative_type_scope='All finite Coxeter types under the inherited classification; rank2 arbitrary labels represented by I2*. E6 and A17 retained distinctly.',
        new_dimension_bound=None,
        scope='Exact all-forest rank64 maximum for the explicit enlarged finite rule. Item-first generation; saved replay checks all-last-item inequalities, predecessors, complete component-count domain, literal coupon interactions and all global attainers. Original type classification and same-component geometric rules inherited. A1/B10 coupon is included in the counting rule without asserting its new ideal application here.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();path=ROOT/f'results/{TAG}.json'
    saved=json.loads(path.read_text())if a.verify else None;out=json.loads(json.dumps(run(saved)))
    if saved is not None:assert out==saved
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['finite_maximum'],out['checks'],flush=True)
