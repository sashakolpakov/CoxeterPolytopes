#!/usr/bin/env python3
"""All-support, all-label common escape for symmetric finite chain arcs.

Standalone Fraction Q(sqrt2) arithmetic; exact support-budget reduction.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_symmetric_arc_escape_gate_20261016'
BASE='results/fresh_h20_s61_multicap_20261014_exact.json'
CAP='results/fresh_h20_dotted_arc_cap_20261016.json'
ABEL='results/fresh_h20_chain_abelian_gate_20261016.json'
REVIEW='results/fresh_h20_chain_abelian_gate_20261016_review.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def add(x,y):return (x[0]+y[0],x[1]+y[1])
def mul(x,y):return (x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def scale(x,a):return (x[0]*a,x[1]*a)
def sign(x):
    a,b=x
    if not b:return (a>0)-(a<0)
    if not a:return (b>0)-(b<0)
    if a*b>0:return (a>0)-(a<0)
    d=a*a-2*b*b;assert d
    return ((a>0)-(a<0))*((d>0)-(d<0))
def enc(x):return list(map(str,x))

def run():
    base,cap,abel,rev=[json.loads((ROOT/p).read_text())for p in [BASE,CAP,ABEL,REVIEW]]
    for doc in [base,cap,abel,rev]:
        for p,h in doc['sources'].items():assert sha(p)==h,p
    assert rev['sources'][ABEL]==sha(ABEL)
    S=base['support'];J=[[F(q)for q in row]for row in base['support_inverse']]
    alo,ahi=map(F,base['parameter']['isolating_interval']);assert alo*alo>280
    assert len(S)==18 and all(q>=0 for row in J for q in row)
    inside=[S.index(i)for i in range(1,11)];outside=[i for i in range(18)if i not in inside]
    assert [[J[i][j]for j in inside]for i in inside]==[
        [F(min(i,j)*(11-max(i,j)),11)for j in range(1,11)]for i in range(1,11)]
    assert all(J[i][j]==0 for i in inside for j in outside)
    assert min(J[i][i]for i in outside)==F(1,2)
    assert sum(J[i][j]for i in outside for j in outside)==13
    total=F(281,70);small=total-2;assert small==F(141,70)
    assert 4+4/(alo*alo)<total
    assert total-F(40,11)<F(1,2)
    assert 4*F(18,11)>total
    assert F(24,11)>small
    assert F(10,11)*F(5,2)>small # label>=5: squared magnitude>5/2
    assert F(18,11)*2>small     # second end: label>=4
    assert 3*F(10,11)>small
    eb=list(map(F,cap['strict_escape']['dual_coefficients']))
    assert eb==[F(-1)if i==5 else F(-1,100000)for i in S]
    zz=[-sum(J[i][j]*eb[j]for j in range(18))for i in range(18)]
    K=F(9999,10000);eps=F(1,100000);allow=F(1,1000)
    assert eps*eps*13*total<allow*allow
    # A dotted A10 contact can only be a single endpoint with all other
    # support contacts zero. Exact reverse-Lorentz-CS margin on its whole ray.
    endpoint=[]
    for wall in [1,10]:
        j=S.index(wall);u=4*J[j][j]-2
        lead=J[j][j]*K-zz[j]*zz[j]
        margin=u*lead-2*zz[j]*zz[j]
        assert lead>0 and margin>0
        endpoint.append(dict(wall=wall,u_lower=str(u),linear_coefficient=str(lead),strict_margin=str(margin)))
    # Otherwise some dotted outside contact costs>2. The A10 cost is<141/70.
    # Only endpoints1/10 at labels3/4 and second endpoints2/9 atlabel3 remain.
    zero=(F(0),F(0));one=(F(1),F(0));sqrt2=(F(0),F(1))
    slots=[1,2,9,10];domains=[[zero,one,sqrt2],[zero,one],[zero,one],[zero,one,sqrt2]]
    kept=[];rejected=[]
    for co in product(*domains):
        d=[zero]*18
        for i,x in zip(slots,co):d[S.index(i)]=x
        cost=zero;pay=zero
        for i in inside:
            pay=add(pay,scale(d[i],zz[i]))
            for j in inside:cost=add(cost,scale(mul(d[i],d[j]),J[i][j]))
        key=[0 if x==zero else 3 if x==one else 4 for x in co]
        if sign(add(cost,(-small,F(0))))>=0:
            rejected.append(dict(labels=key,cost=enc(cost)));continue
        assert cost[1]==0
        u=max(cost[0],1/(ahi*ahi))
        rhs=add(pay,(allow,F(0)))
        margin=add((K*u,F(0)),scale(mul(rhs,rhs),-1))
        assert sign(margin)>0
        kept.append(dict(labels=key,A10_cost=str(cost[0]),escape_A10_pair_magnitude=enc(pay),
                         u_lower=str(u),outside_pair_upper=str(allow),strict_squared_margin=enc(margin)))
    assert len(kept)==8 and len(rejected)==28
    expected=[[0,0,0,0],[3,0,0,0],[4,0,0,0],[0,3,0,0],[0,0,3,0],[0,0,0,3],[0,0,0,4],[3,0,0,3]]
    assert sorted(r['labels']for r in kept)==sorted(expected)
    paths=[BASE,CAP,ABEL,REVIEW,f'scripts/{TAG}.py']
    return dict(status='PASS_ALL_SUPPORT_ALL_LABEL_SYMMETRIC_ARC_COMMON_ESCAPE',
        sources={p:sha(p)for p in paths},centers=list(range(108,132)),support=S,
        all_label_budget=dict(total_support_cost_upper=str(total),A10_cost_upper_if_dotted_outside=str(small),
                              outside_inverse_sum='13',outside_escape_pair_upper=str(allow)),
        dotted_A10_endpoint_cases=endpoint,A10_slots=slots,A10_profiles=kept,A10_excluded_profiles=rejected,
        proof_reductions=[
            'B<=A would put all26 chain contacts in[-2,0], contradicting the reviewed abelian gate. Thus B>A>=0 and B>=1.',
            'If all support contacts were finite or parallel, their rational-Gram cost would lie in Q^ab and force C=B^2/(A^2+2cost-4) in Q^ab. Hence at least one support contact is dotted.',
            'Dotted A10 contact forces a single endpoint and no other nonzero support contact by the total budget.',
            'Otherwise an outside dotted coordinate costs>2. Complete A10 finite-label reduction gives the eight recorded profiles.',
            'Reverse Lorentz CS, outside positive-definite CS and the exact ten margins prove r.E<0 for every allowed cap.'],
        escape_dependency=CAP,
        conclusion='The saved simple strict spacelike escape survives every coherent symmetric finite/parallel three-chain-contact cap with arbitrary18 support contacts. Therefore any mutually compatible finite collection of such caps remains infinite volume.',
        scope='All24 consecutive triple centers and all integer finite labels or parallel contacts, equal outer contacts only. All18 support contacts may independently be right, finite, parallel or dotted; arbitrary real coefficients. The other133-wall pair conditions are only relaxed, never replaced by sampling. Unequal outer contacts, fewer than three finite/parallel chain contacts or changed retained walls remain open.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['dotted_A10_endpoint_cases']),len(d['A10_profiles']),flush=True)
