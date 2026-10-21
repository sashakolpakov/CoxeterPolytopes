#!/usr/bin/env python3
"""Separate physical-vector verification of the 144 Lorentz dual witnesses.

Reuses the existing exact field arithmetic, but not the primary coefficient
builder, optimization model or dual checker.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h20_e8_interior_caps_20261014_exact_review as K

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_e7_residual_gate_20261014_physical_review'
PRIMARY='results/fresh_h20_e7_residual_gate_20261014_dual.json'
BASE='results/fresh_h20_e7_intrinsic_caps_20261014_cusps.json'
PARENT='results/fresh_h20_e7_intrinsic_caps_20261014.json'

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def scale(a,q):return tuple(x*q for x in a)
def scalar(q):return K.scalar(q)+K.Z
def mul(a,b):
    return K.add(K.multiply(a[:8],b[:8]),K.scale(K.multiply(a[8:],b[8:]),3))+K.add(K.multiply(a[:8],b[8:]),K.multiply(a[8:],b[:8]))
def sign(a):
    sa,sb=K.sign(a[:8]),K.sign(a[8:])
    if not sb:return sa
    if not sa or sa==sb:return sb
    return sa*K.sign(K.sub(K.multiply(a[:8],a[:8]),K.scale(K.multiply(a[8:],a[8:]),3)))
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    return x

def run():
    primary,b,parent=[json.loads((ROOT/p).read_text()) for p in [PRIMARY,BASE,PARENT]]
    for doc in [primary,b,parent]:
        for p,h in doc['sources'].items():assert sha(p)==h,p
    Q=[[F(x) for x in row] for row in b['form']]
    sparse=[[(j,x) for j,x in enumerate(row) if x] for row in Q];cache={}
    def dot(x,y):
        if y not in cache:cache[y]=tuple(tuple(sum(a*y[j][k] for j,a in row) for k in range(8)) for row in sparse)
        v=K.Z
        for a,bb in zip(x,cache[y]):
            if any(a) and any(bb):v=K.add(v,K.multiply(a,bb))
        return v
    def linear(vs,cs):
        out=(K.Z,)*21
        for v,c in zip(vs,cs):out=K.vadd(out,K.vscale(v,c))
        return out
    N=list(map(K.vector,b['normals']));E=K.vector(b['strict_escape']['point']);H=K.vector(b['strict_center'])
    S=b['E7_support'];A=parent['face_facet_walls'][:13]
    I=K.inverse([[K.rational(dot(N[i],N[j])) for j in S] for i in S])
    W=[linear([N[i] for i in S],list(map(K.scalar,row))) for row in I]
    sq2=(F(0),F(1))+(F(0),)*6
    T=K.vscale(K.vector(parent['face_basis'][-1]),K.scale(sq2,F(1,2)))
    assert dot(T,T)==K.scalar(-2) and K.sign(dot(T,H))<0
    assert all(dot(T,w)==K.Z for w in W)
    records=[]
    for rec in primary['records']:
        X=E
        for key,val in rec['old_wall_multipliers'].items():
            i,q=int(key),F(val);assert 0<=i<87 and q>0
            X=K.vadd(X,K.vscale(N[i],K.scalar(-q)))
        for key,val in rec['nonnegative_contact_multipliers'].items():
            i,q=int(key),F(val);assert 0<=i<13 and q>0
            X=K.vadd(X,K.vscale(N[A[i]],K.scalar(-q)))
        U=K.vadd(X,linear(W,[K.neg(dot(N[i],X)) for i in S]))
        assert all(dot(U,N[i])==K.Z for i in S)
        delta=K.neg(dot(U,U));h=K.neg(dot(T,U))
        assert K.sign(delta)>0 and K.sign(h)>0
        assert delta==K.field(rec['Lorentz_dual_square'])
        assert h==K.field(rec['time_coefficient'])
        c=[scalar(0) for _ in range(7)]
        label_coeff={3:scalar(1),4:sq2+K.Z,6:K.Z+K.scalar(1)}
        for j,label,coeff in zip(rec['support'],rec['labels'],rec['contact_coefficients']):
            c[j]=label_coeff[label]
            assert c[j]==tuple(map(F,coeff))
        mass=scalar(-2);constant=scalar(0)
        for i in range(7):
            constant=add(constant,mul(c[i],K.neg(dot(W[i],X))+K.Z))
            for j in range(7):mass=add(mass,scale(mul(c[i],c[j]),I[i][j]))
        assert mass==tuple(map(F,rec['mass_square'])) and sign(mass)>=0
        assert constant==tuple(map(F,rec['upper_bound_constant']))
        shifted=add(constant,scalar(F(1,100000)))
        if sign(shifted)>=0:
            gap=add(mul(mass,delta+K.Z),scale(mul(shifted,shifted),-1))
            assert sign(gap)>0
        else:gap=None
        records.append(dict(support=rec['support'],labels=rec['labels'],dual_projection=U,
                            dual_square=K.neg(delta),constant=constant,mass=mass,squared_margin=gap))
    assert len(records)==144
    paths=[PRIMARY,BASE,PARENT,f'scripts/{TAG}.py','scripts/fresh_h20_e8_interior_caps_20261014_exact_review.py']
    return enc(dict(status='PASS_144_DIRECT_PHYSICAL_VECTOR_LORENTZ_DUAL_REVIEWS',sources={p:sha(p) for p in paths},
        records=records,scope='Independent physical-vector reconstruction of all144 dual witnesses, using the same existing exact field engine. The primary coefficient builder, solver and dual checker are not imported. Completeness and causality of the H13 face cone remain inherited exact dependencies.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(path.read_text())
    else:path.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
