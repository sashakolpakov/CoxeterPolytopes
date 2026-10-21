#!/usr/bin/env python3
"""Two nonright caps closing the specified open edge of the frozen96 chamber.

Exact finite construction only. All frozen predecessor files are read-only.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from functools import lru_cache
from math import isqrt
from collections import Counter
import fresh_h20_e7_causal_residuals_20261014_bridge_review as B
from fresh_h20_e8_interior_caps_20261014_exact_review import pd

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_bridge_edge_closure_20261014'
BASE='results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json'
PRE='results/fresh_h20_e7_bridge_face_20261014_exact.json'
FACE='results/fresh_h20_bridge_faces_20261014_review.json'
RAD=[1,2,79,158,11,22,869,1738,3,6,237,474,33,66,2607,5214]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def parse(x,k):return tuple(map(F,x))if k==0 else tuple(parse(y,k-1)for y in x)
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(tuple,list)):return[enc(y)for y in x]
    if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
    return x
def run():
    old,pre,face=[json.loads((ROOT/p).read_text())for p in[BASE,PRE,FACE]]
    for doc in[old,pre,face]:
        for p,h in doc['sources'].items():assert sha(p)==h,p
    assert old['field']['base_square_root_basis']==RAD
    zeros=[B.Z]
    for _ in range(3):zeros.append((zeros[-1],zeros[-1]))
    radicals=[parse(old['field']['discriminant'],0),parse(pre['field']['new_positive_sqrt_radical'],1)]
    def const(q,k=2):
        x=(F(q),)+(F(0),)*15
        for i in range(k):x=(x,zeros[i])
        return x
    def add(x,y,k=2):return B.add(x,y)if k==0 else tuple(add(a,b,k-1)for a,b in zip(x,y))
    def scale(x,q,k=2):return B.scale(x,q)if k==0 else tuple(scale(a,q,k-1)for a in x)
    def mul(x,y,k=2):
        if x==zeros[k]or y==zeros[k]:return zeros[k]
        if k==0:return B.mul(x,y)
        a,b=x;c,d=y;z=zeros[k-1]
        if b==z:return(mul(a,c,k-1),mul(a,d,k-1))
        if d==z:return(mul(a,c,k-1),mul(b,c,k-1))
        return(add(mul(a,c,k-1),mul(radicals[k-1],mul(b,d,k-1),k-1),k-1),add(mul(a,d,k-1),mul(b,c,k-1),k-1))
    def total(xs,k=2):
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
        l,h=interval(x[0],k-1,bits);a,b=interval(x[1],k-1,bits)
        if x[1]==zeros[k-1]:return l,h
        c,d=sqrtbounds(k,bits);v=[a*c,a*d,b*c,b*d]
        return l+min(v),h+max(v)
    @lru_cache(None)
    def sqrtbounds(k,bits):
        l,h=interval(radicals[k-1],k-1,bits);assert l>0
        den=1<<bits;a=isqrt(l.numerator*den*den//l.denominator);b=isqrt(h.numerator*den*den//h.denominator)+1
        assert F(a*a,den*den)<=l and F(b*b,den*den)>h
        return F(a,den),F(b,den)
    maxbits=0
    def sign(x,k=2):
        nonlocal maxbits
        if x==zeros[k]:return 0
        for bits in[192,384,768,1536,3072]:
            try:l,h=interval(x,k,bits)
            except AssertionError:continue
            if l>0 or h<0:maxbits=max(maxbits,bits);return 1 if l>0 else -1
        raise AssertionError(('unresolved sign',k))
    Q=[[F(q)for q in row]for row in old['form']];rows=[[(j,q)for j,q in enumerate(row)if q]for row in Q]
    def vparse(v,k):return tuple(parse(x,k)for x in v)
    def lift(v):return tuple((x,zeros[1])for x in vparse(v,1))
    N=[lift(v)for v in old['normals']]+[vparse(pre['new_normal'],2)]
    H=vparse(pre['strict_center'],2);E=vparse(pre['strict_escape']['point'],2)
    cache={}
    def dot(v,w):
        if w not in cache:cache[w]=tuple(total(scale(w[j],q)for j,q in row)for row in rows)
        return total(mul(x,y)for x,y in zip(v,cache[w]))
    def lin(vs,cs):return tuple(total(mul(c,v[j])for c,v in zip(cs,vs))for j in range(21))
    def rational(x):
        for k in[2,1]:assert x[1]==zeros[k-1];x=x[0]
        assert not any(x[1:]);return x[0]
    one=const(1);S=pre['recipe']['right_support'];assert 61 in S
    G=[[rational(dot(N[i],N[j]))for j in S]for i in S];assert pd(G);GI=B.inv(G)
    omega=lin([N[i]for i in S],[const(q)for q in GI[S.index(61)]]);assert dot(omega,omega)==const(F(2,3))
    assert all(dot(omega,N[i])==const(i==61)for i in S)
    r=N[95];assert dot(r,r)==const(2)
    X=vparse(pre['actual_edge']['point'],2);qX=dot(X,X);A=scale(qX,-1);assert sign(A)>0
    assert all(dot(X,N[i])==zeros[2]for i in S+[95])
    def project(v):
        a=[dot(v,N[i])for i in S];coef=[total(scale(x,q)for x,q in zip(a,row))for row in GI]
        return lin([v]+[N[i]for i in S],[one]+[scale(x,-1)for x in coef])
    A0=project(N[0]);aX=dot(A0,X)
    W=lin([A0,r,X],[qX,scale(qX,F(1,2)),scale(aX,-1)])
    w2=dot(W,W);assert sign(w2)>0
    assert dot(W,X)==zeros[2]and all(dot(W,N[i])==zeros[2]for i in S+[95])
    p=lin([omega,r],[const(-1),const(F(-1,2))]);assert dot(p,p)==const(F(7,6))
    assert dot(p,X)==zeros[2]and dot(p,W)==zeros[2]
    norms=list(map(F,old['norms']))+[F(2)]
    xN=[dot(X,n)for n in N];wN=[dot(W,n)for n in N];pN=[dot(p,n)for n in N]
    xH,wH,pH=[dot(v,H)for v in[X,W,p]]
    xE,wE,pE=[dot(v,E)for v in[X,W,p]]
    # Each iteration is an exact test. No floating point chooses acceptance.
    m=F(1);found=None
    for power in range(80):
        mw2=scale(w2,m*m)
        if sign(add(mw2,const(-2)))<=0:m*=2;continue
        Delta=mul(A,add(mw2,const(F(-5,6))))
        assert sign(Delta)>0
        if len(radicals)==2:radicals.append(Delta)
        else:radicals[2]=Delta
        sqrtbounds.cache_clear()
        pairvec=[];centers=[]
        for side in[1,-1]:
            pairs=[(mul(A,add(pn,scale(wn,side*m))),xn)for pn,wn,xn in zip(pN,wN,xN)]
            ch=(mul(A,add(pH,scale(wH,side*m))),xH)
            ok=sign(ch,3)<0
            for i,g in enumerate(pairs):
                if i in S+[95]:
                    expected=scale(A,-1)if i in[61,95]else zeros[2]
                    assert g==(expected,zeros[2])
                else:
                    # All remaining old walls must be strictly dotted.
                    target=(scale(mul(A,A),2*norms[i]),zeros[2])
                    if sign(g,3)>=0 or sign(add(mul(g,g,3),scale(target,-1,3),3),3)<=0:ok=False;break
            if not ok:break
            pairvec.append(pairs);centers.append(ch)
        if len(pairvec)==2:found=(pairvec,centers);break
        m*=2
    assert found is not None
    pairvec,centerpairs=found;Delta=radicals[2];A2=mul(A,A);newnorm=scale(A2,2)
    capnorm=(newnorm,zeros[2]);mutual=(mul(A2,add(const(2),scale(w2,-2*m*m))),zeros[2])
    assert sign(mutual,3)<0 and sign(add(mul(mutual,mutual,3),scale(mul(capnorm,capnorm,3),-1,3),3),3)>0
    caps=[];endpoints=[]
    for side in[1,-1]:
        caps.append(tuple((mul(A,add(v,scale(w,side*m))),x)for v,w,x in zip(p,W,X)))
        endpoints.append(tuple((scale(mul(w2,x),m),scale(w,side))for x,w in zip(X,W)))
    # Formula-derived exact geometry; no inverse in the new radical field.
    pointnorm=scale(mul(A,w2),F(-5,6));assert sign(pointnorm)<0
    mate=(zeros[2],scale(mul(A,w2),-2*m));assert sign(mate,3)<0
    pointrows=[]
    for k,side in enumerate([1,-1]):
        pairs=[(scale(mul(w2,xn),m),scale(wn,side))for xn,wn in zip(xN,wN)]
        assert all(sign(g,3)==0 if i in S+[95]else sign(g,3)<0 for i,g in enumerate(pairs))
        future=(scale(mul(w2,xH),m),scale(wH,side));assert sign(future,3)<0
        pairs+= [(zeros[2],zeros[2]),mate]if k==0 else[mate,(zeros[2],zeros[2])]
        # Normalizing the new normal by A gives the exact rational20 Gram.
        ids=S+[95];M=[[rational(dot(N[i],N[j]))for j in ids]for i in ids]
        c=[F(-1)if i in[61,95]else F(0)for i in ids]
        gram=[row+[ci]for row,ci in zip(M,c)]+[c+[F(2)]];assert pd(gram)
        assert F(2)-F(2,3)-F(1,2)==F(5,6)
        pointrows.append(dict(wall=96+k,point=endpoints[k],norm=(pointnorm,zeros[2]),future_pairing=future,incident=ids+[96+k],all98_pairings=pairs,normalized_support_gram=gram,final_Schur_pivot=F(5,6)))
    # The inherited center is strict at both new caps, so every listed facet
    # and all old finite vertices/bounded cusps are retained by projection.
    assert sign(dot(H,H))<0 and all(sign(dot(H,n))<0 for n in N)
    cusp_points=[tuple((tuple(map(F,a)),B.Z)for a in old['completed_cusp']['point'])]+[tuple((tuple(map(F,a)),B.Z)for a in c['point'])for c in old['preserved_old_bounded_cusps']]
    cusp_rows=[]
    for nu0 in cusp_points:
        nu=tuple((x,zeros[1])for x in nu0);xx,ww,pp=[dot(v,nu)for v in[X,W,p]]
        pair=[(mul(A,add(pp,scale(ww,side*m))),xx)for side in[1,-1]];assert all(sign(g,3)<0 for g in pair)
        cusp_rows.append(dict(point=nu,new_pairings=pair))
    escape_pairs=[(mul(A,add(pE,scale(wE,side*m))),xE)for side in[1,-1]]
    escape_survives=all(sign(g,3)<0 for g in escape_pairs)
    assert escape_survives,'Need a separately certified remaining escape; do not claim finite volume.'
    assert sign(dot(E,E))>0 and sign(dot(E,H))<0 and all(sign(dot(E,n))<0 for n in N)
    counts=Counter()
    for pairs in pairvec:
        for i,g in enumerate(pairs):counts['pi3'if i in[61,95]else'right'if i in S else'dotted']+=1
    counts['dotted']+=1
    assert sum(counts.values())==193
    paths=[BASE,PRE,FACE,f'scripts/{TAG}.py','scripts/fresh_h20_e7_causal_residuals_20261014_bridge_review.py','scripts/fresh_h20_e8_interior_caps_20261014_exact_review.py']
    return enc(dict(status='PASS_98_WALL_FINITE_EDGE_CLOSURE_STILL_INFINITE',sources={f:sha(f)for f in paths},form=Q,field=dict(base=pre['field'],new_positive_sqrt_radical=Delta,encoding='[old64 a,old64 b] = a+b sqrt(Delta); old64 follows the bound96 certificate'),recipe=dict(support=S,contact_wall=61,dual_weight=omega,dual_weight_norm=F(2,3),X=X,A=A,W=W,W_norm=w2,particular=p,particular_norm=F(7,6),boost_integer=m,Delta=Delta,formula='R±=A(-omega-r95/2±mW)+sqrt(Delta)X; Delta=A(m²W²-5/6); R±²=2A²'),new_normals=caps,new_normal_norm=capnorm,new_pairings=pairvec,mutual_pairing=mutual,new_pair_counts=dict(counts),strict_center=H,new_center_pairings=centerpairs,finite_endpoints=pointrows,finite_edge=dict(shared_support=S+[95],endpoint_walls=[96,97],proof='Both future timelike endpoints satisfy all98 inequalities, have precisely the stated20 active roots, and are strict at the opposite cap. Their joining feasible segment has exactly the sharedPD19 support active.'),preserved_bounded_cusps=cusp_rows,old_finite_vertex_count_at_least=116,known_finite_vertex_count_at_least=118,strict_escape=dict(point=E,norm=dot(E,E),new_pairings=escape_pairs,all_old96_signs_rechecked=True),maximum_interval_bits=maxbits,scope='Exactly two coherent nonright endpoint caps, retaining96 and closing its specified rank19 edge into a finite edge. All193 new pairs, two actual distinct PD20 vertices, common strict center, three cusp signs and a surviving strict escape are exact. No complete new vertex census or finite-volume claim.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(p.read_text())==out
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_pair_counts'],'m=',out['recipe']['boost_integer'],flush=True)
