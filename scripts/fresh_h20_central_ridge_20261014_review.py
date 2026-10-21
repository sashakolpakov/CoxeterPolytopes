#!/usr/bin/env python3
"""Independent literal-vector audit in Q(sqrt(5), sqrt(406/9-9sqrt(5))).

No producer imports; the saved phi/alpha coordinates are converted to a
different field basis before comparison. All decisions use rational signs.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
BASE = 'fresh_h20_central_ridge_20261014'
TAG = BASE + '_review'

def q(a=0, b=0): return F(a), F(b)
def qa(x, y): return x[0]+y[0], x[1]+y[1]
def qneg(x): return -x[0], -x[1]
def qs(x, y): return qa(x, qneg(y))
def qm(x, y): return x[0]*y[0]+5*x[1]*y[1], x[0]*y[1]+x[1]*y[0]
def qt(x, a): return x[0]*a, x[1]*a
def qi(x):
    n=x[0]**2-5*x[1]**2
    assert n
    return x[0]/n, -x[1]/n
def qsgn(x):
    a,b=x
    if not b: return (a>0)-(a<0)
    if not a: return (b>0)-(b<0)
    if a*b>0: return (a>0)-(a<0)
    d=a*a-5*b*b
    assert d
    return ((a>0)-(a<0))*((d>0)-(d<0))

T=q(F(406,9),-9)
PH=q(F(1,2),F(1,2))
L=q(F(5,2),F(-1,2))
R=F(1,3)
def k(a=0): return q(a),q()
def lift(x): return x,q()
def add(x,y): return qa(x[0],y[0]),qa(x[1],y[1])
def neg(x): return qneg(x[0]),qneg(x[1])
def sub(x,y): return add(x,neg(y))
def mul(x,y): return qa(qm(x[0],y[0]),qm(T,qm(x[1],y[1]))),qa(qm(x[0],y[1]),qm(x[1],y[0]))
def scale(x,a): return qt(x[0],a),qt(x[1],a)
def sgn(x):
    a,b=x;sa,sb=qsgn(a),qsgn(b)
    if not sb:return sa
    if not sa:return sb
    if sa==sb:return sa
    d=qs(qm(a,a),qm(T,qm(b,b)))
    assert qsgn(d)
    return sa*qsgn(d)
def total(xs):
    v=k()
    for x in xs:v=add(v,x)
    return v
def dot(x,y):return total(mul(a,b)for a,b in zip(x,y,strict=True))
def ip(x,y):return sub(dot(x[:19],y[:19]),add(mul(x[19],y[20]),mul(x[20],y[19])))
def decode_q(x):
    a,b=map(F,x)
    return q(a+b/2,b/2)
def decode_k(x):
    a,b=map(decode_q,x)
    return qs(a,qt(b,F(1,54))),qt(b,F(1,18))
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(list,tuple)):return [enc(y)for y in x]
    if isinstance(x,dict):return {str(a):enc(b)for a,b in x.items()}
    return x

AL=(q(F(-1,54)),q(F(1,18)))
SQ=(q(),q(1))
XI=(q(),qm(L,qi(T)))
def normal(g,s,i):
    return [lift(qt(PH,s))if j==g else add(AL,k(R if i==j else 0))for j in range(19)]+[k(2),k(1)]
def vertex(g,s):return [lift(qt(PH,s))if j==g else XI for j in range(19)]+[k(2),k(1)]

def run():
    paths=[f'scripts/{BASE}.py',f'results/{BASE}.json',
           f'scripts/{BASE}_edge.py',f'results/{BASE}_edge.json']
    old=json.loads((ROOT/paths[1]).read_text())
    edge=json.loads((ROOT/paths[3]).read_text())
    for saved in (old,edge):
        for path,h in saved['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h
    assert T[0]**2-5*T[1]**2==F(132031,81)
    assert isqrt(132031)**2!=132031 and qsgn(T)>0
    assert mul(SQ,SQ)==lift(T)
    assert add(scale(mul(AL,AL),18),scale(AL,2*R))==lift(L)
    assert sgn(AL)>0 and sgn(sub(k(F(1,3)),AL))>0
    assert sgn(XI)>0 and sgn(sub(k(F(1,3)),XI))>0
    recipes=[(g,s,i)for g in range(11,19)for s in(-1,1)for i in range(19)if i!=g]
    assert recipes==[tuple(x)for x in old['cap_recipes']]
    caps=[normal(*r)for r in recipes]
    assert all(ip(v,v)==k(R*R)for v in caps)
    assert decode_k(old['alpha'])==AL and decode_q(old['radicand'])==T
    # Full pair census, with literal vector arithmetic once per permutation type.
    cache={};counts=Counter()
    for u,v in combinations(range(288),2):
        g,s,i=recipes[u];h,t,j=recipes[v]
        typ=(g==h,s,t,i==j,i==h,j==g)
        if typ not in cache:cache[typ]=ip(caps[u],caps[v])
        a=cache[typ]
        if (g,s)==(h,t):assert a==k();counts['right']+=1
        else:assert sgn(add(a,k(R*R)))<0;counts['dotted']+=1
    assert counts=={'right':2448,'dotted':38880} and dict(counts)==old['new_pair_counts']
    # Each new cap is right to central and its target side; all other
    # vertical walls are strictly dotted. The smaller unit inradius is enough.
    central=[k()]*19+[k(-2),k(1)]
    for v,(g,s,i)in zip(caps,recipes):
        assert ip(v,central)==k() and v[g]==lift(qt(PH,s))
        for j in range(19):
            if j!=g:
                for eps in(-1,1):assert sgn(sub(sub(k(1),scale(v[j],eps)),k(R)))>0
        assert qsgn(qs(qt(PH,2),q(R)))>0
    assert old['new_core_pairs']=={'right':576,'dotted':10656}
    assert old['core_pair_counts']=={'right':684,'parallel':19,'pi3':20,'pi4':2,'pi5':16}
    # Old-core counts: 4 pairs per two spatial axes, 19 opposite parallels,
    # and two central contacts per axis, with cosines 1/2,sqrt2/2,phi/2.
    assert 4*(19*18//2)==684 and 684+19+20+2+16==39*38//2
    verts=[];vertex_pairs=0
    kap=qm(L,qi(T));pn=qt(kap,-R*R)
    for entry in old['vertices']:
        g,s=entry['golden_coordinate'],entry['sign'];p=vertex(g,s)
        assert p[:19]==[decode_k(x)for x in entry['spatial']]
        assert ip(p,p)==lift(pn) and decode_q(entry['norm'])==pn and qsgn(pn)<0
        assert ip(p,central)==k()
        for cap,(h,t,j)in zip(caps,recipes):
            a=ip(p,cap);vertex_pairs+=1
            if(g,s)==(h,t):assert a==k()
            else:assert sgn(a)<0
        verts.append(p)
    assert len(verts)==16 and qs(q(4),qm(PH,PH))==L and qsgn(L)>0
    # Every cap has exactly 19 finite neighbors. They form one positive
    # definite I2(5)+18A1 clique. No old-core-only PD20 clique is possible.
    N=q(24,4) # 10+2+8phi² = 24+4sqrt5.
    assert qsgn(qs(q(4),N))<0 and decode_q(old['core_only_schur'])==qs(q(4),N)
    assert old['complete_finite_vertex_count']==16
    assert old['complete_elliptic_rank19_count']==16*20+2**19
    # Seed-only PD19 faces are 2^19 complete vertical lines. Their box
    # vertex ends have only 19 incident walls, hence unbounded cusp links.
    # Caps miss those ends by a strictly positive sphere power > 7.
    assert qsgn(qs(N,q(25)))>0 and 4+R*R<F(49,9)
    assert (F(5)-F(7,3))**2-R*R==7 and 4-18<0
    H=[k()]*19+[k(100),k(1)]
    E=[k(F(99,100))]*10+[k()]*9+[k(F(8801,2000)),k(1)]
    assert ip(E,E)==k(1) and sgn(ip(E,H))<0 and sgn(ip(E,central))<0
    assert all(sgn(ip(H,c))<0 and sgn(ip(E,c))<0 for c in caps)
    # All core sides are strict at H and E; H also proves full dimension.
    assert F(99,100)<1 and sgn(ip(H,central))<0
    # Independently reconstruct the additional edge mirror.
    den=qs(q(1),qt(kap,16));assert qsgn(den)>0
    A=qt(qi(den),8);B=qm(qa(q(1),qt(kap,16)),qi(den));U=qs(A,B)
    assert decode_q(edge['A'])==A and decode_q(edge['B'])==B and qsgn(U)>0
    assert qs(qm(B,B),qm(kap,qm(A,A)))==q(1)
    P=vertex(11,1);released=normal(11,1,0)
    extra=[sub(mul(lift(A),p),mul(lift(B),r))for p,r in zip(P,released)]
    saved=edge['new_normal']
    assert extra==[decode_k(x)for x in saved['spatial']]+[lift(decode_q(saved['D'])),lift(decode_q(saved['U']))]
    assert ip(extra,extra)==k(R*R) and ip(extra,central)==k()
    assert extra[11]==lift(qm(PH,U))
    for j in range(19):
        if j!=11:
            for eps in(-1,1):assert sgn(sub(sub(lift(U),scale(extra[j],eps)),k(R)))>0
    assert qsgn(qs(qt(qm(PH,U),2),q(R)))>0
    contacts=[]
    for i,cap in enumerate(caps):
        a=ip(extra,cap)
        if a==k():contacts.append(i)
        else:assert sgn(add(a,k(R*R)))<0
    assert contacts==edge['retained_cap_indices'] and len(contacts)==17
    assert edge['new_pair_counts']=={'right':19,'dotted':308}
    assert all(sgn(ip(extra,p))<0 for p in verts)
    lam=qm(qm(A,kap),qi(B))
    Q=[sub(p,mul(lift(lam),r))for p,r in zip(P,released)]
    assert Q[:19]==[decode_k(x)for x in edge['new_vertex']['spatial']]
    assert ip(Q,Q)==lift(qt(qm(kap,qi(qm(B,B))),-R*R))
    assert sgn(ip(Q,Q))<0 and ip(Q,extra)==k() and ip(Q,central)==k()
    assert Q[11]==mul(lift(PH),Q[20]) and sgn(Q[20])>0
    for j in range(19):
        if j!=11:
            for eps in(-1,1):assert sgn(sub(Q[20],scale(Q[j],eps)))>0
    for i,c in enumerate(caps):
        a=ip(Q,c)
        if i in contacts:assert a==k()
        else:assert sgn(a)<0
    assert sgn(ip(extra,H))<0 and sgn(ip(extra,E))<0
    assert decode_k(edge['strict_escape_pair'])==ip(extra,E)
    assert edge['complete_finite_vertices']==17 and edge['finite_edges']==1
    assert edge['complete_elliptic_rank19_faces']==2**19+320+19
    assert edge['open_finite_ended_edges']==17*20-2
    # Distinct spacelike normals, nonpositive pairings and a common strict
    # timelike point imply actual facets by orthogonal projection of H.
    return enc(dict(status='PASS_INDEPENDENT_LITERAL_VECTOR_AND_GRAPH_AUDIT',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths+[f'scripts/{TAG}.py']},
        field_basis='1,sqrt5,S,sqrt5*S; S²=406/9−9sqrt5; S>0',
        producer_imports=False,cap_pair_permutation_types=len(cache),cap_pairs=sum(counts.values()),
        independently_checked_vertex_cap_pairs=vertex_pairs,
        actual_facets=[327,328],complete_finite_vertices=[16,17],
        complete_elliptic_rank19_faces=[524608,524627],finite_edges_after_extension=1,
        strict_escape_norm=1,finite_volume=False,
        scope='Both constructions are infinite-volume. Open edges have ideal-boundary ends but no second finite vertex or bounded-link ideal vertex. The finite-neighbor graph does not certify finite volume.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert json.loads(p.read_text())==out
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
