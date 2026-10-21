#!/usr/bin/env python3
"""Two literal rational polar facets on the frozen98 chamber.

Exact new-pair geometry only; no numerical scout is a proof dependency.
The separate finite-vertex census is not implemented in this file.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from functools import lru_cache
from math import isqrt
from collections import Counter
import fresh_h20_e7_causal_residuals_20261014_bridge_review as B

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_polar_pair_20261014'
BASE='results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json'
BRIDGE='results/fresh_h20_e7_bridge_face_20261014_exact.json'
PRE='results/fresh_h20_bridge_edge_closure_20261014.json'
SUM='results/fresh_h20_bridge_faces_20261014_review.json'
RAD=[1,2,79,158,11,22,869,1738,3,6,237,474,33,66,2607,5214]
RECIPES=[
 [7,14,21,18,15,12,9,6,3,1,-1,-3,-5,-6,-7,-8,-9,-10,10,0,0],
 [14,27,40,35,30,25,20,15,10,5,0,-5,-10,-12,-14,-16,-18,-20,20,0,0]]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def parse(x,k):return tuple(map(F,x))if k==0 else tuple(parse(y,k-1)for y in x)
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(tuple,list)):return[enc(y)for y in x]
    if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
    return x
def run():
    old,bridge,pre,summed=[json.loads((ROOT/p).read_text())for p in[BASE,BRIDGE,PRE,SUM]]
    for doc in[old,bridge,pre,summed]:
        for p,h in doc['sources'].items():assert sha(p)==h,p
    assert old['field']['base_square_root_basis']==RAD
    zeros=[B.Z]
    for _ in range(3):zeros.append((zeros[-1],zeros[-1]))
    radicals=[parse(old['field']['discriminant'],0),parse(bridge['field']['new_positive_sqrt_radical'],1),parse(pre['field']['new_positive_sqrt_radical'],2)]
    def const(q,k=3):
        x=(F(q),)+(F(0),)*15
        for i in range(k):x=(x,zeros[i])
        return x
    def lift(x,k,target=3):
        for j in range(k,target):x=(x,zeros[j])
        return x
    def add(x,y,k=3):return B.add(x,y)if k==0 else tuple(add(a,b,k-1)for a,b in zip(x,y))
    def scale(x,q,k=3):return B.scale(x,q)if k==0 else tuple(scale(a,q,k-1)for a in x)
    def mul(x,y,k=3):
        if x==zeros[k]or y==zeros[k]:return zeros[k]
        if k==0:return B.mul(x,y)
        a,b=x;c,d=y;z=zeros[k-1]
        if b==z:return(mul(a,c,k-1),mul(a,d,k-1))
        if d==z:return(mul(a,c,k-1),mul(b,c,k-1))
        return(add(mul(a,c,k-1),mul(radicals[k-1],mul(b,d,k-1),k-1),k-1),add(mul(a,d,k-1),mul(b,c,k-1),k-1))
    def total(xs,k=3):
        z=zeros[k]
        for x in xs:z=add(z,x,k)
        return z
    @lru_cache(None)
    def rb(bits):
        den=1<<bits;out=[]
        for d in RAD:
            a=isqrt(d*den*den);out.append((F(a,den),F(a+int(a*a!=d*den*den),den)))
        return out
    def interval(x,k,bits):
        if k==0:
            lo=hi=F(0)
            for q,(l,h)in zip(x,rb(bits)):
                lo+=q*(l if q>=0 else h);hi+=q*(h if q>=0 else l)
            return lo,hi
        l,h=interval(x[0],k-1,bits)
        if x[1]==zeros[k-1]:return l,h
        a,b=interval(x[1],k-1,bits);c,d=sqrtbounds(k,bits);v=[a*c,a*d,b*c,b*d]
        return l+min(v),h+max(v)
    @lru_cache(None)
    def sqrtbounds(k,bits):
        l,h=interval(radicals[k-1],k-1,bits);assert l>0
        den=1<<bits;a=isqrt(l.numerator*den*den//l.denominator);b=isqrt(h.numerator*den*den//h.denominator)+1
        assert F(a*a,den*den)<=l and F(b*b,den*den)>h
        return F(a,den),F(b,den)
    maxbits=0
    def sign(x,k=3):
        nonlocal maxbits
        if x==zeros[k]:return 0
        for bits in[192,384,768,1536,3072]:
            try:l,h=interval(x,k,bits)
            except AssertionError:continue
            if l>0 or h<0:maxbits=max(maxbits,bits);return 1 if l>0 else -1
        raise AssertionError(('unresolved sign',k))
    Q=[[F(q)for q in row]for row in old['form']]
    assert Q==[[F(q)for q in row]for row in pre['form']]
    def qv(v):return[sum(q*x for q,x in zip(row,v))for row in Q]
    def rdot(v,w):return sum(x*y for x,y in zip(v,qv(w)))
    def linear(r,v,k=3):return total((scale(x,q,k)for x,q in zip(v,qv(r))if q),k)
    N=[tuple(lift(parse(x,1),1)for x in v)for v in old['normals']]
    N+=[tuple(lift(parse(x,2),2)for x in bridge['new_normal'])]
    N +=[tuple(parse(x,3)for x in v)for v in pre['new_normals']]
    assert len(N)==98
    norms=[const(F(q))for q in old['norms']]+[const(2)]+[parse(pre['new_normal_norm'],3)]*2
    caps=[tuple(map(F,r))for r in RECIPES];newq=[rdot(r,r)for r in caps]
    assert newq==[F(2),F(4)]and rdot(*caps)==0
    pairs=[[linear(r,n)for n in N]for r in caps]
    def relation(g,prod):
        if g==zeros[3]:return'right'
        assert sign(g)<0
        sq=mul(g,g)
        for c,name in[(F(1,4),'pi3'),(F(1,2),'pi4'),(F(3,4),'pi6'),(F(1),'parallel')]:
            if sq==scale(prod,c):return name
        assert sign(add(sq,scale(prod,-1)))>0,'Unexpected finite angle: inspect explicitly.'
        return'dotted'
    rels=[[relation(g,scale(q,newq[j]))for g,q in zip(row,norms)]for j,row in enumerate(pairs)]
    counts=Counter(x for row in rels for x in row);counts['right']+=1;assert sum(counts.values())==197
    # T is a strict future timelike positive sum certified independently at95.
    # Reconstruct its physical coordinates and inherit the hash-bound old signs.
    T=tuple(total((parse(old['normals'][i][j],1)for i in range(95)),1)for j in range(21))
    qT=parse(summed['positive_wall_sum_norm'],1);assert sign(qT,1)<0
    tp=[lift(parse(x,1),1)for x in summed['positive_wall_sum_all_pairings']]
    assert len(tp)==95 and all(sign(x)<0 for x in tp)
    assert [x['wall']for x in bridge['new_pairs']]==list(range(95))
    tp.append(lift(total((parse(row['pairing'],2)for row in bridge['new_pairs']),2),2))
    tp +=[total(parse(g,3)for g in row[:95])for row in pre['new_pairings']]
    newtp=[linear(r,T,1)for r in caps];tp +=[lift(x,1)for x in newtp]
    assert len(tp)==100 and all(sign(x)<0 for x in tp)
    # Exact facets via the acute-center projection lemma; save both new points.
    facets=[]
    for k,(r,q,tr)in enumerate(zip(caps,newq,newtp)):
        alpha=scale(tr,F(1,q),1)
        P=tuple(add(t,scale(alpha,-c,1),1)for t,c in zip(T,r))
        norm=add(qT,scale(mul(tr,tr,1),-F(1,q),1),1);assert sign(norm,1)<0
        full=pairs[k]+[const(rdot(r,s))for s in caps]
        signs=[add(t,scale(mul(lift(alpha,1),g),-1))for t,g in zip(tp,full)]
        assert all(sign(x)==0 if i==98+k else sign(x)<0 for i,x in enumerate(signs))
        facets.append(dict(wall=98+k,point=P,norm=norm,all100_pairings=signs,point_formula='T-(T.r/r²)r'))
    cusps=[]
    points=[old['completed_cusp']['point']]+[c['point']for c in old['preserved_old_bounded_cusps']]
    for v in points:
        nu=tuple(parse(x,0)for x in v);g=[linear(r,nu,0)for r in caps]
        assert all(sign(x,0)<0 for x in g)
        cusps.append(dict(point=nu,new_pairings=g))
    E=tuple(parse(x,2)for x in pre['strict_escape']['point'])
    Er=[]
    for x in E:
        assert x[1]==zeros[1]and x[0][1]==zeros[0]and not any(x[0][0][1:]);Er.append(x[0][0][0])
    qE=rdot(Er,Er);eg=[rdot(r,Er)for r in caps]
    assert qE==F(20181174579467302741,10**18)and all(g<0 for g in eg)
    assert sign(linear(Er,T,1),1)<0
    # The two spacelike normal vectors themselves are feasible old rays.
    assert all(sign(g)<=0 for row in pairs for g in row)
    finite_neighbors=[[i for i,x in enumerate(row)if x in['right','pi3','pi4','pi6']]+[99-k]for k,row in enumerate(rels)]
    paths=[BASE,BRIDGE,PRE,SUM,f'scripts/{TAG}.py','scripts/fresh_h20_e7_causal_residuals_20261014_bridge_review.py','results/acute_extension_preservation_20261014.md']
    return enc(dict(status='PASS_100_WALL_RATIONAL_POLAR_PAIR_STILL_INFINITE',sources={p:sha(p)for p in paths},form=Q,field=pre['field'],recipes=RECIPES,new_normals=caps,new_norms=newq,mutual_pairing=F(0),new_pairings=pairs,new_relations=rels,new_pair_counts=dict(counts),complete_finite_neighbors=finite_neighbors,old98_zero_sets=[[i for i,g in enumerate(row)if g==zeros[3]]for row in pairs],strict_center=dict(point=T,norm=qT,all100_pairings=tp,source='sum of the first95 physical normals; old95 signs and norm independently certified in the bound whole-face review'),new_facet_points=facets,preserved_bounded_cusps=cusps,known_finite_vertex_count_at_least=118,strict_escape=dict(rational_point=Er,norm=qE,new_pairings=eg,all_old98_signs_inherited=True),maximum_interval_bits=maxbits,scope='Two literal rational polar walls, exact197 new Coxeter pairs, common strict future center, all100 facets actual, three old bounded cusps and at least118 old finite vertices preserved. The strict rational escape remains. New finite vertices are not enumerated here; a separate census may strengthen the vertex count. No finite-volume claim.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(p.read_text())==out
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_pair_counts'],flush=True)
