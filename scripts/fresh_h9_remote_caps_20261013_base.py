#!/usr/bin/env python3
"""Exact rank-five remote-face frame and coupled all-right cap formulas."""
from scripts.fresh_h9_split_second_frontier_20261013_base import *

PREFIX='data/enumerations/fresh_h9_remote_caps_20261013'
REMOTE='results/fresh_h9_split_multiple_edge_escape_20261013.json'
S_NAMES=('seed9','seed10','seed11','seed14','seed15')
SIGNS=((1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1))

def source_frame():
    G,I,h,old,geometry=selected();point=tuple(map(parse,read(REMOTE)['point']))
    center=tuple(K.neg(x)for x in matvec(K,I,h));w=tuple(K.neg(x)for x in matvec(K,I,point))
    def proj(v):
        for n in S_NAMES:
            r=old[n];assert r['field']is K
            v=tuple(sub(x,mul(div(dot(K,v,r['p']),scale(O,2)),y))for x,y in zip(v,r['c'],strict=True))
        return v
    t=proj(center);wp=proj(w);a=K.neg(pairing(K,t,G,t));q=pairing(K,wp,G,wp);pt=pairing(K,wp,G,t)
    assert K.sign(a)>0 and K.sign(q)>0
    v=tuple(add(x,mul(div(pt,a),y))for x,y in zip(wp,t,strict=True));b=pairing(K,v,G,v)
    basis=[t,v];norms=[K.neg(a),b];chosen=[]
    for name in('seed8','seed12','seed13','seed0','seed1','seed2'):
        z=proj(old[name]['c'])
        for u,nu in zip(basis,norms,strict=True):
            c=div(pairing(K,z,G,u),nu);z=tuple(sub(x,mul(c,y))for x,y in zip(z,u,strict=True))
        qz=pairing(K,z,G,z)
        if K.sign(qz)==0:continue
        assert K.sign(qz)>0;basis.append(z);norms.append(qz);chosen.append(name)
        if len(basis)==5:break
    assert len(basis)==5
    assert all(pairing(K,basis[i],G,basis[j])==(norms[i]if i==j else Z)for i in range(5)for j in range(5))
    assert all(dot(K,u,old[n]['p'])==Z for u in basis for n in S_NAMES)
    return G,I,h,old,geometry,point,basis,norms,chosen,wp

def make_caps(G,basis,norms,k):
    assert k>1
    ct=Q(k*k+1,2*k);st=Q(k*k-1,2*k);assert ct*ct-st*st==1
    a=K.neg(norms[0]);b,c,d,e=norms[1:]
    radicands=[div(scale(O,Q(1,2)),q)for q in(a,b,c,d,e)]
    F=K;radicals=[]
    for r in radicands:
        assert K.sign(r)>0
        new=Quadratic(F,lift_K(F,r));radicals=[new.lift(x)for x in radicals]+[(F.zero,F.one)];F=new
    vectors=[tuple(F.mul(radicals[i],lift_K(F,x))for x in basis[i])for i in range(5)]
    common=tuple(F.add(F.mul(F.rational(st),x),F.mul(F.rational(ct),y))for x,y in zip(vectors[0],vectors[1],strict=True))
    coefficients=[]
    for signs in SIGNS:
        coefficients.append(tuple(F.add(x,add_many(F,[F.mul(F.rational(s),vectors[j+2][i])for j,s in enumerate(signs)]))for i,x in enumerate(common)))
    point=tuple(F.add(F.mul(F.rational(ct),x),F.mul(F.rational(st),y))for x,y in zip(vectors[0],vectors[1],strict=True))
    GF=tuple(tuple(lift_K(F,x)for x in row)for row in G)
    rows=[]
    for i,cvec in enumerate(coefficients):
        name=f'remote{i}';pvec=matvec(F,GF,cvec)
        rows.append({'name':name,'field':F,'c':cvec,'p':pvec,'q':scale(O,2)})
    return rows,point,radicands,(ct,st)

def add_many(F,values):
    out=F.zero
    for x in values:out=F.add(out,x)
    return out
