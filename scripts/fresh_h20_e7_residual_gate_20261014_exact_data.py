#!/usr/bin/env python3
"""Exact physical coefficients for the E7 all-residual optimization model."""
import json
from decimal import Decimal, localcontext
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h20_e8_interior_caps_20261014_exact_review as K

ROOT=Path(__file__).resolve().parents[1]
BASE='results/fresh_h20_e7_intrinsic_caps_20261014_cusps.json'
PARENT='results/fresh_h20_e7_intrinsic_caps_20261014.json'
FACE='results/fresh_h20_e7_face_20261014.json'
RADICANDS=[1,2,79,158,11,22,869,1738]


def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()


def value(x):
    """Exact zeros stay zero; only final field values are rounded to float."""
    if not any(x):return 0.0
    with localcontext() as ctx:
        ctx.prec=90
        v=sum(Decimal(a.numerator)/Decimal(a.denominator)*Decimal(r).sqrt() for a,r in zip(x,RADICANDS))
        assert v!=0
        assert (1 if v>0 else -1)==K.sign(x)
        return float(v)


def load():
    base,parent,face=[json.loads((ROOT/p).read_text())for p in [BASE,PARENT,FACE]]
    for doc in [base,parent,face]:
        for p,h in doc['sources'].items():assert digest(p)==h,p
    assert base['normals'][:85]==parent['normals'] and base['form']==parent['form']
    Q=[[F(a)for a in row]for row in base['form']]
    sparse=[[(j,a)for j,a in enumerate(row)if a]for row in Q];cache={}
    def dot(x,y):
        if y not in cache:
            cache[y]=tuple(tuple(sum(a*y[j][k]for j,a in row)for k in range(8))for row in sparse)
        out=K.Z
        for a,b in zip(x,cache[y]):
            if any(a)and any(b):out=K.add(out,K.multiply(a,b))
        return out
    N=list(map(K.vector,base['normals']));E=K.vector(base['strict_escape']['point']);H=K.vector(base['strict_center'])
    T=K.vector(face['intrinsic_normalization']['timelike_basis'])
    S=base['E7_support'];A=parent['face_facet_walls'][:13]
    assert len(N)==87 and len(S)==7 and len(A)==13 and not set(S)&set(A)
    I=K.inverse([[K.rational(dot(N[i],N[j]))for j in S]for i in S])
    J=K.inverse([[K.rational(dot(N[i],N[j]))for j in A]for i in A])
    assert K.pd(I)and K.pd(J)
    def combination(ids,row):
        return tuple(tuple(sum(row[j]*N[i][a][k]for j,i in enumerate(ids))for k in range(8))for a in range(21))
    W=[combination(S,row)for row in I];Lam=[combination(A,row)for row in J]
    assert all(dot(W[i],N[j])==K.scalar(int(i==k))for i in range(7)for k,j in enumerate(S))
    assert all(dot(Lam[i],N[j])==K.scalar(int(i==k))for i in range(13)for k,j in enumerate(A))
    assert all(dot(w,l)==K.Z for w in W for l in Lam)
    assert dot(T,T)==K.scalar(-2)and K.sign(dot(T,H))<0
    assert all(dot(T,N[i])==K.Z for i in S+A)
    coords=[[dot(l,n)for n in N]for l in Lam]
    times=[dot(T,n)for n in N];wp=[[dot(w,n)for n in N]for w in W]
    a=[K.neg(dot(l,E))for l in Lam];B=K.neg(dot(T,E));weightE=[K.neg(dot(w,E))for w in W]
    assert K.sign(B)>0
    assert all(coords[i][s]==K.Z for i in range(13)for s in S)
    assert all(times[s]==K.Z for s in S+A)
    assert all(wp[i][j]==K.scalar(int(i==k))for i in range(7)for k,j in enumerate(S))
    assert all(wp[i][j]==K.Z for i in range(7)for j in A)
    return dict(base=base,parent=parent,face=face,Q=Q,N=N,E=E,H=H,T=T,S=S,A=A,I=I,J=J,W=W,Lam=Lam,
                coords=coords,times=times,wp=wp,a=a,B=B,weightE=weightE,dot=dot)
