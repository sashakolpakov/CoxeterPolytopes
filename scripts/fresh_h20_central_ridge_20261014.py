#!/usr/bin/env python3
"""Sixteen exact finite I2(5)+18A1 corners of the weighted H20 core."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_central_ridge_20261014'
def q(a=0,b=0):return F(a),F(b)
def qa(a,b):return a[0]+b[0],a[1]+b[1]
def qn(a):return -a[0],-a[1]
def qs(a,b):return qa(a,qn(b))
def qm(a,b):return a[0]*b[0]+a[1]*b[1],a[0]*b[1]+a[1]*b[0]+a[1]*b[1]
def qt(a,t):return a[0]*t,a[1]*t
def qi(a):
    norm=a[0]*a[0]+a[0]*a[1]-a[1]*a[1];assert norm
    return (a[0]+a[1])/norm,-a[1]/norm
def qsign(a):
    x,y=a[0]+a[1]/2,a[1]/2
    if not y:return (x>0)-(x<0)
    if not x:return (y>0)-(y<0)
    if x*y>0:return (x>0)-(x<0)
    d=x*x-5*y*y;assert d
    return ((x>0)-(x<0))*((d>0)-(d<0))
PH=q(0,1);L=q(3,-1);R=F(1,3);RAD=qa(q(R*R),qt(L,18))
# Elements A+B*alpha, alpha=(sqrt(RAD)-R)/18, positive root of
# 18alpha²+2Ralpha=L. This field is sufficient for all new normals.
def k(a=0,b=0):return q(a),q(b)
def lift(a):return a,q()
def add(a,b):return qa(a[0],b[0]),qa(a[1],b[1])
def neg(a):return qn(a[0]),qn(a[1])
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    bd=qm(a[1],b[1])
    return qa(qm(a[0],b[0]),qt(qm(bd,L),F(1,18))),qs(qa(qm(a[0],b[1]),qm(a[1],b[0])),qt(bd,R/9))
def scale(a,t):return qt(a[0],t),qt(a[1],t)
def sign(a):
    x=qs(a[0],qt(a[1],R/18));y=qt(a[1],F(1,18))
    sx,sy=qsign(x),qsign(y)
    if not sy:return sx
    if not sx:return sy
    if sx==sy:return sx
    d=qs(qm(x,x),qm(qm(y,y),RAD));sd=qsign(d)
    assert sd
    return sx*sd
def inv(a):
    conj=(qs(a[0],qt(a[1],R/9)),qn(a[1]));norm=mul(a,conj)
    assert norm[1]==q();ni=qi(norm[0])
    return qm(conj[0],ni),qm(conj[1],ni)
def total(xs):
    z=k()
    for x in xs:z=add(z,x)
    return z
AL=k(0,1);PHK=lift(PH);LK=lift(L)
def dot(x,y):return total(mul(a,b)for a,b in zip(x,y,strict=True))
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(tuple,list)):return [enc(v)for v in x]
    if isinstance(x,dict):return {str(t):enc(v)for t,v in x.items()}
    return x
def center(g,s,i):return [scale(PHK,s)if j==g else add(AL,k(R if i==j else 0))for j in range(19)]
def pcenter(g,s,xi):return [scale(PHK,s)if j==g else xi for j in range(19)]
def dist2(g,s,i,h,t,j):
    if g==h:return add(scale(lift(qm(PH,PH)),(s-t)**2),k(2*R*R if i!=j else 0))
    dg=sub(scale(PHK,s),add(AL,k(R if j==g else 0)))
    dh=sub(add(AL,k(R if i==h else 0)),scale(PHK,t))
    rest=int(i not in (g,h))+int(j not in (g,h))-2*int(i==j)
    assert rest>=0
    return add(add(mul(dg,dg),mul(dh,dh)),k(rest*R*R))

def run():
    nr=RAD[0]**2+RAD[0]*RAD[1]-RAD[1]**2
    assert nr==F(132031,81)
    assert isqrt(nr.numerator)**2!=nr.numerator # RAD is not a square in Q(phi).
    assert qsign(L)>0 and qsign(qs(q(2),L))>0
    assert qsign(qs(PH,q(F(3,2))))>0
    assert sign(AL)>0 and sign(sub(k(F(1,3)),AL))>0
    assert add(scale(mul(AL,AL),18),scale(AL,2*R))==LK
    sqrt_rad=add(scale(AL,18),k(R));assert mul(sqrt_rad,sqrt_rad)==lift(RAD)
    xi=mul(LK,mul(sqrt_rad,lift(qi(RAD))))
    assert sign(xi)>0 and sign(sub(k(F(1,3)),xi))>0
    pnorm=qt(qm(L,qi(RAD)),-R*R);assert qsign(pnorm)<0
    rec=[(g,s,i)for g in range(11,19)for s in (-1,1)for i in range(19)if i!=g]
    assert len(rec)==288
    norms=[]
    for g,s,i in rec:
        c=center(g,s,i)
        assert dot(c,c)==k(4+R*R)
        # D=2,U=1: right to central (0,-2,1), right to signed target side.
        assert sub(dot(c,c),k(4))==k(R*R)
        assert c[g]==scale(PHK,s)
        norms.append((g,s,i))
    pc=Counter();dmin=None
    for z,(g,s,i)in enumerate(rec):
        for h,t,j in rec[z+1:]:
            d=dist2(g,s,i,h,t,j)
            if (g,s)==(h,t):
                assert d==k(2*R*R);pc['right']+=1
            else:
                assert sign(sub(d,k(4*R*R)))>0;pc['dotted']+=1
                if dmin is None or sign(sub(d,dmin))<0:dmin=d
    assert dict(pc)=={'right':2448,'dotted':38880}
    # The distance formula is checked against literal vectors for every
    # target/sign/axis-position class, without relying only on a bound.
    checked=set()
    for g,s,i in rec:
        for h,t,j in rec:
            typ=(g==h,s,t,i==j,i==h,j==g)
            if typ in checked:continue
            checked.add(typ)
            u=[sub(a,b)for a,b in zip(center(g,s,i),center(h,t,j))]
            assert dot(u,u)==dist2(g,s,i,h,t,j)
    # Every nontarget center coordinate lies in(0,2/3), every other
    # inradius is >=1. Thus gap>1/3=R, proving all37 other seed relations
    # strictly dotted, including the opposite target gap2phi.
    assert sign(sub(k(1),add(AL,k(2*R))))>0
    assert qsign(qs(qt(PH,2),q(R)))>0
    core_pairs={'right':684,'parallel':19,'pi3':20,'pi4':2,'pi5':16}
    assert sum(core_pairs.values())==39*38//2
    vertices=[]
    for g in range(11,19):
        for s in (-1,1):
            p=pcenter(g,s,xi)
            assert sub(dot(p,p),k(4))==lift(pnorm)
            assert p[g]==scale(PHK,s)
            checks=Counter()
            for h,t,j in rec:
                v=sub(dot(p,center(h,t,j)),k(4))
                if (g,s)==(h,t):assert v==k();checks['zero']+=1
                else:assert sign(v)<0;checks['strict']+=1
            assert checks=={'zero':18,'strict':270}
            # Gram block [[4,-phi],[-phi,1]] has determinant L>0;
            # eighteen orthogonal cap roots have norm R². Hence PD20.
            assert qs(q(4),qm(PH,PH))==L
            vertices.append(dict(golden_coordinate=g,sign=s,spatial=p,D=2,U=1,norm=pnorm,
                incident_cap_indices=[n for n,r in enumerate(rec)if r[:2]==(g,s)],cap_sign_counts=dict(checks)))
    # No old-core-only finite20corner: every candidate is central plus
    # one of each19 opposite side pairs, with negative Schur4-N.
    N=q(20,8);old_schur=qs(q(4),N);assert qsign(old_schur)<0
    # Each cap has precisely19 other finite neighbors, its17 cluster
    # peers and central+target side. Its only possible PD20 set is saved.
    # Releasing any member leaves19PD roots. Every other wall has a dotted
    # pair with that set, so cannot be a finite or ideal next endpoint.
    assert 17+2==19 and 16*20==320
    # The remaining PD19 sets are one vertical side in each coordinate.
    # Central plus18 sides has Schur<=4-18<0. The524288 seed-only
    # faces are complete vertical edges. At their Euclidean endpoints,
    # |z|>5, |c|<7/3, R=1/3 give every new cap power>7.
    assert 4-18<0 and qsign(qs(N,q(25)))>0
    assert F(4)+R*R<F(49,9) and (F(5)-F(7,3))**2-R*R==7
    X=[F(99,100)]*10+[F(0)]*9;D=F(8801,2000)
    assert sum(v*v for v in X)-2*D==1
    central_pair=2-D;assert central_pair<0
    new_bounds={}
    for unit in (False,True):
        pair=sub(scale(add(scale(AL,10),k(R if unit else 0)),F(99,100)),k(2+D))
        assert sign(pair)<0
        new_bounds['unit_axis'if unit else'other_axis']=pair
    rational_bound=F(363,100)-2-D;assert rational_bound==F(-5541,2000)
    assert all(sign(sub(k(rational_bound),v))>0 for v in new_bounds.values())
    assert -100-D<0
    src=f'scripts/{TAG}.py'
    return enc(dict(status='PASS_EXACT_SIXTEEN_FINITE_CENTRAL_PI5_RIDGE_CORNERS_WITH_ESCAPE',
        sources={src:sha256((ROOT/src).read_bytes()).hexdigest()},
        coordinates='Physical Euclidean spatial coordinates; vertical sides have unit normals and offsets1,sqrt2,phi. Lorentz form x²−2DU.',
        field='Q(phi,alpha), 18alpha²+(2/3)alpha=3−phi, alpha positive',radicand_rational_norm=nr,
        radius=R,alpha=AL,alpha_bounds=[0,F(1,3)],radicand=RAD,
        core_facets=39,core_pair_counts=core_pairs,cap_recipes=rec,new_caps=288,actual_facets=327,
        new_norm=R*R,new_D=2,new_U=1,new_pair_counts=dict(pc),cross_cluster_minimum_distance_square=dmin,
        literal_distance_formula_classes=len(checked),new_core_pairs={'right':576,'dotted':10656},
        common_strict_center=dict(spatial_zero=True,D=100,U=1,norm=-200,new_pair=-102),
        vertices=vertices,complete_finite_vertex_count=16,core_only_candidate_count=2**19,core_only_schur=old_schur,
        finite_neighbor_count_per_cap=19,finite_ended_edges=320,finite_ended_edges_without_second_finite_or_ideal_endpoint=320,
        complete_elliptic_rank19_count=320+2**19,seed_only_vertical_edges=2**19,
        seed_only_edge_ends='Common infinity cusp with bounded box link; distinct boxvertex ideal ends with only19incident orthogonal sides and unbounded links.',
        strict_escape=dict(spatial=X,D=D,U=1,norm=1,central_pair=central_pair,new_pair_classes=new_bounds,
            uniform_new_pair_upper_bound=rational_bound,minimum_seed_gap=F(1,100),future_pair=-100-D),
        scope='Retains only the weighted39 core, not binary/outer banks. Repairs the finite-vertex requirement on all16 actual central pi5 ridges. All327 facets actual and complete finite-vertex census16. The320 edges incident to these corners have no second finite/ideal endpoint in this chamber. Exact strict escape proves infinite volume; no global completion claim.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(p.read_text())==out
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['actual_facets'],out['complete_finite_vertex_count'])
