#!/usr/bin/env python3
"""Independent exact local H17 half-coupon transfer; no global bound."""
import argparse
import json
from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h17_half_coupon_transfer_20261017_review'
OLD='fresh_nonincident_carrier_incidence_20261017'
N=set(range(19))
EDGES=[(0,1),(0,2),(0,4),(2,3),(4,5),(5,6),(6,7),(7,8),(8,9)]
EDGES += [(a+10 if a<9 else 9,b+10 if b<9 else 9)for a,b in EDGES.copy()]
ADJ={i:set()for i in N}
for a,b in EDGES:ADJ[a].add(b);ADJ[b].add(a)
G=[[F(2 if i==j else -int(j in ADJ[i]))for j in range(19)]for i in range(19)]
Z=[-6,-3,-4,-2,-5,-4,-3,-2,-1,0,6,3,4,2,5,4,3,2,1]


def inertia(a):
    a=[r.copy()for r in a];out=[0,0,0]
    while a:
        n=len(a);i=next((i for i in range(n)if a[i][i]),None)
        if i is not None:
            q=a[i][i];out[0 if q>0 else 1]+=1;rest=[j for j in range(n)if j!=i]
            a=[[a[j][k]-a[j][i]*a[i][k]/q for k in rest]for j in rest]
        else:
            e=next(((i,j)for i in range(n)for j in range(i+1,n)if a[i][j]),None)
            if e is None:out[2]+=n;break
            i,j=e;q=a[i][j];out[0]+=1;out[1]+=1;rest=[k for k in range(n)if k not in e]
            a=[[a[k][l]-(a[k][i]*a[j][l]+a[k][j]*a[i][l])/q for l in rest]for k in rest]
    return out


def sub(nodes):return [[G[i][j]for j in nodes]for i in nodes]
def component(nodes,a):
    nodes=set(nodes);ds={a:0};q=deque([a])
    while q:
        x=q.popleft()
        for y in sorted((ADJ[x]&nodes)-ds.keys()):ds[y]=ds[x]+1;q.append(y)
    return ds


def name_and_role(nodes,mark):
    nodes=set(nodes);deg={i:len(ADJ[i]&nodes)for i in nodes};fork=[i for i in nodes if deg[i]==3]
    assert sum(deg.values())==2*(len(nodes)-1)
    if not fork:
        assert max(deg.values())<=2
        return 'A'+str(len(nodes)),None
    assert len(fork)==1
    p=fork[0];short=[i for i in ADJ[p]&nodes if deg[i]==1]
    assert len(short)==2
    role='fork_center'if mark==p else'short_fork_leaf'if mark in short else'other'
    return 'D'+str(len(nodes)),role


def run():
    assert inertia(G)==[17,1,1]
    assert all(sum(G[i][j]*Z[j]for j in N)==0 for i in N)
    affine=sorted(N-{3,13});assert inertia(sub(affine))==[16,0,1]
    # A literal positive null vector of the17-node affineD16 block.
    K={i:(1 if len(ADJ[i]&set(affine))==1 else 2)for i in affine}
    assert all(sum(G[i][j]*K[j]for j in affine)==0 for i in affine)
    old=json.loads((ROOT/f'results/{OLD}.json').read_text())
    saved={tuple(f['side']):f for f in old['face_ledger']['two_faces']if f['polygon']==4}
    rows=[];coupon_needed=[]
    for a,b in product((0,1,2),(10,11,12)):
        black={a,3,b,13};W=sorted(N-black)
        assert inertia(sub(W))==[15,0,0]
        face=saved[tuple(sorted(black))]
        assert sorted(face['white'])==W and len(face['vertices'])==4
        core_corner=[]
        for left,right in product((a,3),(b,13)):
            marks=[left,right];ids=sorted(set(W)|set(marks));sig=inertia(sub(ids))
            if marks==[a,b]:
                assert sig==[16,0,1] and ids==affine
                distance=component(ids,left)[right]
                assert distance>=13 or marks==[0,10]
                core_corner.append(dict(marks=marks,kind='affine_D16',distance=distance,
                    dropped_ideal_payment=True,base_payment='0',coupon='0'))
                continue
            assert sig==[17,0,0]
            c=component(ids,left);distance=c.get(right)
            baseline=int(distance is not None and(distance<=9 or len(c)<=17))
            if distance is not None:
                name,role=name_and_role(c,left);parts=[dict(nodes=sorted(c),type=name,mark=left,role=role)]
                coupon=F(0)
            else:
                d=component(ids,right)
                parts=[]
                for block,mark in ((c,left),(d,right)):
                    name,role=name_and_role(block,mark);parts.append(dict(nodes=sorted(block),type=name,mark=mark,role=role))
                types=[p['type']for p in parts]
                coupon=F(0)
                for x,y in ((0,1),(1,0)):
                    if (types[x],types[y])in(('A1','D16'),('A2','D14')) and parts[y]['role']in('fork_center','short_fork_leaf'):
                        coupon=F(1,2)
            core_corner.append(dict(marks=marks,kind='finite',signature=sig,components=parts,
                distance=distance,base_payment=str(baseline),coupon=str(coupon)))
        base=sum(F(c['base_payment'])for c in core_corner)
        coupons=sum(F(c['coupon'])for c in core_corner)
        assert base+coupons>=1
        if not base:
            assert coupons==1 and sum(F(c['coupon'])==F(1,2)for c in core_corner)==2
            coupon_needed.append([a,b])
        rows.append(dict(dropped_pair=[a,b],white=W,black=sorted(black),corners=core_corner,
            base_total=str(base),coupon_total=str(coupons),face_total=str(base+coupons)))
    assert coupon_needed==[[0,10],[0,12],[2,10],[2,12]]
    paths=[Path(__file__)]
    for stem in(OLD,'fresh_h17_d16_payment_bottleneck_20261017'):
        paths += [ROOT/f'scripts/{stem}.py',ROOT/f'results/{stem}.json',ROOT/f'results/{stem}.md']
    return dict(status='PASS_INDEPENDENT_H17_HALF_COUPON_LOCAL_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        gram=[[str(x)for x in r]for r in G],core_signature=[17,1,1],radical=Z,
        affine_nodes=affine,affine_positive_kernel=[K[i]for i in affine],
        remaining_radical_coefficients={'3':-2,'13':2},cases=rows,
        coupon_required_cases=coupon_needed,old_finite_alternative_cases=5,
        scope='Nine actual-star local cases and a general affine-kernel/radical sealing proof in the companion report. No complete old108 embedding census, finite global budget, ideal pair count, or dimension composition is claimed.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],'9 cases / 27 finite corners / 4 two-half transfers',flush=True)
