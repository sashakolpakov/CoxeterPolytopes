#!/usr/bin/env python3
"""Separate sparse-field replay of the 96th bridge wall and its geometry.

Imports only the independent whole-face review's Fraction sparse arithmetic;
no construction arithmetic or normal-generation routine is imported.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import isqrt
from pathlib import Path
import fresh_h20_bridge_faces_20261014_review as A
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_e7_bridge_face_20261014_exact_review'
BASE='results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json'
INPUT='results/fresh_h20_e7_bridge_face_20261014_exact.json'
FACE='results/fresh_h20_bridge_faces_20261014_review.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def parse(x,level):
    if not level:return tuple(map(F,x))
    return tuple(parse(y,level-1)for y in x)
def run():
    old,dat,face=[json.loads((ROOT/p).read_text())for p in[BASE,INPUT,FACE]]
    for doc in[old,dat,face]:
        for p,h in doc['sources'].items():assert sha(p)==h,p
    assert old['field']['base_square_root_basis']==A.RAD
    assert dat['field']['base']==old['field']
    D=parse(old['field']['discriminant'],0)
    K=parse(dat['field']['new_positive_sqrt_radical'],1)
    zeros=[A.Z,(A.Z,A.Z),((A.Z,A.Z),(A.Z,A.Z))]
    def const(q,level=2):
        v=A.c(q)
        for k in range(level):v=(v,zeros[k])
        return v
    def add(x,y,level=2):
        if not level:return A.ba(x,y)
        return tuple(add(a,b,level-1)for a,b in zip(x,y))
    def scale(x,q,level=2):
        if not level:return A.bs(x,q)
        return tuple(scale(a,q,level-1)for a in x)
    def product(x,y,level=2):
        if x==zeros[level]or y==zeros[level]:return zeros[level]
        if not level:return A.bm(x,y)
        a,b=x;c,d=y;z=zeros[level-1]
        if b==z:return(product(a,c,level-1),product(a,d,level-1))
        if d==z:return(product(a,c,level-1),product(b,c,level-1))
        rad=D if level==1 else K
        ac=product(a,c,level-1);bd=product(b,d,level-1)
        return(add(ac,product(rad,bd,level-1),level-1),
               add(product(a,d,level-1),product(b,c,level-1),level-1))
    def total(xs,level=2):
        v=zeros[level]
        for x in xs:v=add(v,x,level)
        return v
    @lru_cache(None)
    def basebounds(bits):
        den=1<<bits;out=[]
        for q in A.RAD:
            n=isqrt(q*den*den);out.append((F(n,den),F(n+int(n*n!=q*den*den),den)))
        return out
    def interval(x,level,bits):
        if not level:
            lo=hi=F(0)
            for q,(l,h)in zip(x,basebounds(bits)):
                lo+=q*(l if q>=0 else h);hi+=q*(h if q>=0 else l)
            return lo,hi
        l,h=interval(x[0],level-1,bits);a,b=interval(x[1],level-1,bits)
        c,d=radicalbounds(level,bits);v=[a*c,a*d,b*c,b*d]
        return l+min(v),h+max(v)
    @lru_cache(None)
    def radicalbounds(level,bits):
        l,h=interval(D if level==1 else K,level-1,bits)
        assert l>0
        den=1<<bits;p=isqrt(l.numerator*den*den//l.denominator)
        q=isqrt(h.numerator*den*den//h.denominator)+1
        assert F(p*p,den*den)<=l and F(q*q,den*den)>h
        return F(p,den),F(q,den)
    used=0
    def sign(x,level=2):
        nonlocal used
        if x==zeros[level]:return 0
        for bits in[192,384,768,1536]:
            l,h=interval(x,level,bits)
            if l>0 or h<0:used=max(used,bits);return 1 if l>0 else -1
        raise AssertionError(('unresolved interval',level,x))
    assert sign(D,0)>0 and sign(K,1)>0
    Q=[list(map(F,row))for row in old['form']]
    assert Q==[list(map(F,row))for row in dat['form']]
    assert A.inertia(Q)==[20,1,0]
    rows=[[(j,q)for j,q in enumerate(row)if q]for row in Q]
    cache={}
    def dot(x,y):
        if y not in cache:cache[y]=tuple(total(scale(y[j],q)for j,q in row)for row in rows)
        return total(product(a,b)for a,b in zip(x,cache[y]))
    def vparse(v,level=2):return tuple(parse(x,level)for x in v)
    def lift(v):return tuple((x,zeros[1])for x in vparse(v,1))
    def lin(vs,cs):return tuple(total(product(c,v[k])for c,v in zip(cs,vs))for k in range(21))
    N=[lift(v)for v in old['normals']];H=lift(old['strict_center']);E=lift(old['strict_escape']['point'])
    r=vparse(dat['new_normal']);N.append(r)
    def expect(x,raw):assert x==parse(raw,2)
    norms=list(map(F,old['norms']))+[F(2)]
    assert all(dot(n,n)==const(q)for n,q in zip(N,norms))
    assert sign(dot(H,H))<0 and sign(dot(E,E))>0
    rec=dat['recipe'];S=rec['right_support']
    assert S==list(range(1,11))+[12,13,14,15,17,31,61,87]
    def rational(x):
        assert x[1]==zeros[1]and x[0][1]==A.Z and not any(x[0][0][1:])
        return x[0][0][0]
    G=[[rational(dot(N[i],N[j]))for j in S]for i in S]
    assert G==[list(map(F,row))for row in rec['support_gram']]
    piv=A.pd(G);J=A.inv(G);assert min(q for row in J for q in row)>=0
    part=lift(rec['particular']);z=lift(rec['complement']);branch=rec['branch'];assert branch in[-1,1]
    assert r==tuple((p[0],scale(zz[0],branch,1))for p,zz in zip(part,z))
    assert all(dot(part,N[i])==zeros[2]and dot(z,N[i])==zeros[2]for i in S)
    sq3=(tuple(F(i==8)for i in range(16)),A.Z)
    assert dot(part,N[0])==const(-1)and dot(part,N[28])==(scale(sq3,-1,1),zeros[1])
    assert dot(z,N[0])==dot(z,N[28])==dot(part,z)==zeros[2]
    p2=dot(part,part);z2=dot(z,z)
    assert p2==(parse(rec['particular_norm'],1),zeros[1])
    assert z2==(parse(rec['complement_norm'],1),zeros[1])and sign(z2)>0
    assert add(p2,product((K,zeros[1]),z2))==const(2)
    counts=Counter();pairs=[]
    for i,n in enumerate(N[:-1]):
        g=dot(r,n);assert sign(g)<=0
        expect(g,dat['new_pairs'][i]['pairing']);assert dat['new_pairs'][i]['wall']==i
        sq=product(g,g);prod=2*norms[i]
        if g==zeros[2]:kind='right'
        elif sq==const(prod/4):kind='pi3'
        elif sq==const(prod/2):kind='pi4'
        elif sq==const(3*prod/4):kind='pi6'
        elif sq==const(prod):kind='parallel'
        else:assert sign(add(sq,const(-prod)))>0;kind='dotted'
        assert kind==dat['new_pairs'][i]['relation'];counts[kind]+=1;pairs.append(g)
    assert dict(counts)==dat['new_pair_counts']=={'pi3':1,'right':18,'dotted':75,'pi6':1}
    assert [i for i,g in enumerate(pairs)if g==zeros[2]]==S
    expect(dot(E,r),dat['old_escape_pairing']);assert sign(dot(E,r))>0
    C=vparse(dat['strict_center']);P=vparse(dat['new_facet_point']);V=vparse(dat['actual_edge']['point'])
    X=vparse(dat['strict_escape']['point'])
    assert all(q[1]==zeros[1]and q[0][1]==A.Z and not any(q[0][0][1:])for q in X)
    for point,key in[(C,'strict_center_norm'),(V,None)]:
        assert sign(dot(point,point))<0 and sign(dot(point,H))<0
        if key:expect(dot(point,point),dat[key])
    assert all(sign(dot(C,n))<0 for n in N)
    assert P==lin([C,r],[const(1),scale(dot(C,r),F(-1,2))])
    assert sign(dot(P,P))<0 and sign(dot(P,H))<0
    assert dot(P,r)==zeros[2]and all(sign(dot(P,n))<0 for n in N[:-1])
    support=S+[95];assert dat['actual_edge']['support']==support
    vp=[dot(V,n)for n in N]
    assert all(sign(g)==(0 if i in support else -1)for i,g in enumerate(vp))
    expect(dot(V,V),dat['actual_edge']['norm'])
    assert vp==[parse(g,2)for g in dat['actual_edge']['all96_pairings']]
    # G is PD18 and r is orthogonal with norm2, so the support is PD19.
    assert all(dot(r,N[i])==zeros[2]for i in S)and dot(r,r)==const(2)
    assert sign(dot(X,X))>0 and sign(dot(X,H))<0
    xp=[dot(X,n)for n in N];assert all(sign(g)<0 for g in xp)
    expect(dot(X,X),dat['strict_escape']['norm'])
    assert xp==[parse(g,2)for g in dat['strict_escape']['all96_pairings']]
    oldcusps=[old['completed_cusp']['point']]+[c['point']for c in old['preserved_old_bounded_cusps']]
    cusp_pairs=[]
    for raw in oldcusps:
        nu=tuple(((tuple(map(F,a)),A.Z),zeros[1])for a in raw)
        assert dot(nu,nu)==zeros[2]and sign(dot(nu,H))<0
        g=dot(nu,r);assert sign(g)<0;cusp_pairs.append(g)
    assert len(cusp_pairs)==3
    assert cusp_pairs==[parse(g,2)for g in dat['preserved_bounded_cusp_pairings']]
    print('PASS_95_NEW_PAIRS_AND_96_WALL_GEOMETRY',dict(counts),flush=True)
    paths=[BASE,INPUT,FACE,f'scripts/{TAG}.py',
           'scripts/fresh_h20_e7_bridge_face_20261014_exact.py',
           'scripts/fresh_h20_bridge_faces_20261014_review.py']
    return A.enc(dict(status='PASS_INDEPENDENT_96_WALL_BRIDGE_REVIEW',sources={p:sha(p)for p in paths},
        arithmetic='Independent sparse Fraction radical multiplication and outward rational intervals; no construction imports',
        maximum_sign_enclosure_bits=used,all_norms_checked=96,new_pairs_checked=95,new_pair_counts=dict(counts),
        actual_edge_support=support,PD18_pivots=piv,actual_edge_rank=19,
        center_strict_pairings=96,new_facet_strict_old_pairings=95,old_escape_cut=True,
        rational_new_escape_norm=dot(X,X),new_escape_strict_pairings=96,preserved_bounded_cusp_pairings=cusp_pairs,
        finite_vertex_preservation='All old finite vertices are strictly retained by the reviewed acute-extension lemma. The inherited lower bound116 is not a fresh vertex census.',
        scope='Reconstructs saved cap recipe identities, all95 new Coxeter gates,96 norms, exact center/facet/PD19 edge, three cusp signs, old escape cut and new rational strict escape. Old95 pair/corner/link censuses are inherited. No new finite vertex or finite-volume claim.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run()
    p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],flush=True)
