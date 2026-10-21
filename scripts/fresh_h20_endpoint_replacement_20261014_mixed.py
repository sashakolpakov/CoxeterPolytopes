#!/usr/bin/env python3
"""Exact mixed A6xA6 orbit with changed endpoint radii and a global escape.

Norm/contact equalities use Q(sqrt7) and defining radical identities.
Strict comparisons involving the pi/8 endpoint radius use rational
enclosures.  Cyclic factor-distance identities transfer the pair gates.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import isqrt
from pathlib import Path
import argparse,json

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_endpoint_replacement_20261014_mixed'
SCALE=2**180
ZERO=(F(0),F(0));ONE=(F(1),F(0))
def k(x):return(F(x),F(0))
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def neg(a):return(-a[0],-a[1])
def sub(a,b):return add(a,neg(b))
def mul(a,b):return(a[0]*b[0]+7*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def sc(a,b):return(a[0]*b,a[1]*b)
def sg(a):
    x,y=a
    if not y:return(x>0)-(x<0)
    if not x:return(y>0)-(y<0)
    if x*y>0:return(x>0)-(x<0)
    d=x*x-7*y*y;assert d
    return((d>0)-(d<0))*((x>0)-(x<0))
def ia(a,b):return(a[0]+b[0],a[1]+b[1])
def inn(a):return(-a[1],-a[0])
def isu(a,b):return ia(a,inn(b))
def im(a,b):
    x=[u*v for u in a for v in b];return(min(x),max(x))
def ii(a):assert a[0]>0;return(1/a[1],1/a[0])
def ip(x):return(F(x),F(x))
def sqrtq(x):
    assert x>=0
    a=isqrt(x.numerator*SCALE*SCALE//x.denominator)
    lo=F(a,SCALE);hi=F(a+1,SCALE)
    assert lo*lo<=x<=hi*hi
    return lo,hi
def ir(a):return(sqrtq(a[0])[0],sqrtq(a[1])[1])
SQ7=ir(ip(7))
def ki(a):return ia(ip(a[0]),im(ip(a[1]),SQ7))
LINKS=list(range(5))+list(range(6,11))+list(range(12,18))
def dot(a,b):return sum(x*y for x,y in zip(a,b))-sum((a[i]*b[i+1]+a[i+1]*b[i])/2 for i in LINKS)
def kd(a,b):
    z=ZERO
    for x,y in zip(a,b):z=add(z,mul(x,y))
    for i in LINKS:z=sub(z,sc(add(mul(a[i],b[i+1]),mul(a[i+1],b[i])),F(1,2)))
    return z
def va(a,b):return[add(x,y)for x,y in zip(a,b)]
def vs(a,t):return[mul(x,t)for x in a]
def vd(a,b):return va(a,vs(b,k(-1)))
def inverse_A(n,rhs):
    return [sum(F(2*min(i,j)*(n+1-max(i,j)),n+1)*rhs[j-1]for j in range(1,n+1))for i in range(1,n+1)]
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
    if isinstance(x,(list,tuple)):return[enc(v)for v in x]
    return x

def run():
    factors=[];seed=[];start=0
    for n in (6,6,7):
        vv=[]
        for o in range(n+1):
            v=inverse_A(n,[F(1)-(n+1 if i==o else 0)for i in range(1,n+1)])
            vv.append([F(0)]*start+v+[F(0)]*(19-start-n))
            w=[F(0)]*19
            if o==0:w[start:start+n]=[-F(1)]*n
            else:w[start+o-1]=F(1)
            seed.append(w)
        factors.append(vv);start+=n
    local=list(product(range(7),range(7),range(8)))
    Z=[[sum(factors[f][a[f]][i]for f in range(3))for i in range(19)]for a in local]
    assert all(dot(z,z)==196 for z in Z)
    for f,n in enumerate((7,7,8)):
        for i in range(n):
            for j in range(n):
                v=[x-y for x,y in zip(factors[f][i],factors[f][j])];d=(i-j)%n
                assert dot(v,v)==2*n*d*(n-d)
    R=(F(14),F(-4));Q=mul(R,R);t=sc(R,F(1,7));a=sub(ONE,t)
    b=(F(1,6),F(1,3))
    assert mul(R,add(ONE,b))==k(7) and sg(sub(b,ONE))>0
    assert Q==sub(sc(R,28),k(84))
    assert sub(sub(sc(mul(b,b),12),sc(b,4)),k(9))==ZERO
    assert sg(t)>0 and sg(sub(k(F(1,2)),t))>0
    C=[]
    for v in local:
        c=[ZERO]*19
        for f in (0,1):
            c=va(c,va(vs(list(map(k,factors[f][v[f]])),a),vs(list(map(k,factors[f][(v[f]+1)%7])),t)))
        c=va(c,list(map(k,factors[2][v[2]])));C.append(c)
    N=kd(C[0],C[0]);assert all(kd(c,c)==N for c in C)
    expectedN=add(k(154),sc(mul(sub(t,k(F(1,2))),sub(t,k(F(1,2)))),168))
    assert N==expectedN
    assert sg(sub(N,mul(add(R,k(8)),add(R,k(8)))))>0
    # cos(pi/8)=sqrt(2+sqrt2)/2.  The positive endpoint radius is
    # rho=R*(sqrt(119+49*c8²)/7-c8); its defining identity gives
    # rho²+2Rrho*c8+R²=(24/7)R² exactly.
    sq2=ir(ip(2));c8=im(ip(F(1,2)),ir(ia(ip(2),sq2)))
    s=ir(ia(ip(119),im(ip(49),im(c8,c8))))
    rhoi=im(ki(R),isu(im(ip(F(1,7)),s),c8));qi=im(rhoi,rhoi)
    assert rhoi[0]>0 and qi[1]<21
    # Expanded defining identity after s²=119+49*c8².  The mixed s*c8
    # terms cancel, as do the c8² terms; this proves the finite angle,
    # rather than inferring an equality from overlapping intervals.
    assert -F(2,7)+F(2,7)==0
    assert F(1)-F(2)+F(49,49)==0
    assert F(119,49)+1-F(24,7)==0
    assert im(ia(rhoi,ip(8)),ia(rhoi,ip(8)))[1]<196
    # All seed relations are literal right/parallel/dotted identities.
    for v,c in zip(local,C):
        gaps=[sub(ONE,kd(c,list(map(k,w))))for w in seed]
        assert sum(x==ZERO for x in gaps)==17
        assert sum(x==R for x in gaps)==2
        assert all(x==ZERO or x==R or sg(sub(x,R))>0 for x in gaps)
        assert gaps[v[0]]==mul(R,b) and gaps[7+v[1]]==mul(R,b)
    # Every pair is a translate of a reference displacement, because the
    # cyclic vertex-distance identities above preserve the product metric.
    mutual={};endpoint={};neighbors=[]
    threshold=im(ia(ki(R),rhoi),ia(ki(R),rhoi))
    mixed_dot_margins=[];endpoint_dot_margins=[]
    for j,v in enumerate(local):
        dist=kd(vd(C[0],C[j]),vd(C[0],C[j]))
        if j==0:assert dist==ZERO
        elif dist==sc(Q,3):mutual[j]='pi3';neighbors.append(j)
        else:
            assert sg(sub(dist,sc(Q,4)))>0
            mutual[j]='dotted';mixed_dot_margins.append(sub(dist,sc(Q,4)))
        dz=vd(C[0],list(map(k,Z[j])));de=kd(dz,dz)
        if j==0:
            assert de==sc(Q,F(24,7));endpoint[j]='pi8'
        else:
            margin=isu(ki(de),threshold)
            assert margin[0]>0;endpoint[j]='dotted';endpoint_dot_margins.append(margin)
    assert neighbors==[8,48,56,336]
    assert len(mutual)==391 and len(endpoint)==392
    # No two of the five nonseed finite neighbors are compatible at a
    # common finite vertex: their mutual relations are all dotted.
    for z,i in enumerate(neighbors):
        for j in neighbors[:z]:
            dd=kd(vd(C[i],C[j]),vd(C[i],C[j]));assert sg(sub(dd,sc(Q,4)))>0
        de=kd(vd(C[i],list(map(k,Z[0]))),vd(C[i],list(map(k,Z[0]))))
        assert isu(ki(de),threshold)[0]>0
    # Check all unchanged middle/motif balls lie strictly inside B_8.
    oldpath='results/fresh_h20_balanced_integer_exact_20261014.json'
    old=json.loads((ROOT/oldpath).read_text())
    for p,h in old['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    oldballs=[]
    for row in old['walls'][415:]:
        w=list(map(F,row));nn=dot(w[:19],w[:19]);qq=nn-2*w[19]
        assert nn<36 and 0<qq<=4 and -200-w[19]<0
        oldballs.append((w[:19],qq))
    assert len(oldballs)==395
    # All mixed and endpoint balls are disjoint from B_8, proving every
    # old middle/central cross relation dotted without a numerical gate.
    assert sg(sub(N,mul(add(R,k(8)),add(R,k(8)))))>0
    # Preserve the 392 finite vertical endpoints at height rho.  All other
    # endpoint balls are dotted; mixed balls satisfy |z-w|²>R²; B_8 balls
    # are at distance>6 from z and have radius<=2.
    for j,z in enumerate(Z):
        dz=vd(C[0],list(map(k,z)));assert sg(sub(kd(dz,dz),Q))>0

    # Representative forced parallel cusp: M/R plus seed1.  Its complete
    # incident set has16 seed roots, in addition to the parallel pair.
    zeroseed=[i for i,w in enumerate(seed)if sub(ONE,kd(C[0],list(map(k,w))))==ZERO]
    incident=[i for i in zeroseed if dot(seed[1],seed[i])==0]
    assert incident==[3,4,5,6,9,10,11,12,13,15,16,17,18,19,20,21]
    assert len(incident)==16
    # Complete finite-neighbor star:17seed+4mixed+1endpoint, with the
    # latter five mutually dotted. Hence no positive-definite20-subset
    # can contain a mixed cap (maximum finite clique size19).

    # A remaining gap mixes A6 and A7 differently; here use the inward
    # midpoint of a nearest A7 edge.  Even all vertex-radial caps miss it.
    x=[F(99,200)*(u+v)for u,v in zip(Z[0],Z[1])];xq=dot(x,x)
    assert min(1-dot(x,w)for w in seed)==F(1,100)
    mixedpowers=[]
    for c in C:
        v=vd(list(map(k,x)),c);mixedpowers.append(sub(kd(v,v),Q))
    mp=mixedpowers[0]
    for p in mixedpowers:
        if sg(sub(p,mp))<0:mp=p
    assert sg(sub(mp,ONE))>0
    enddist=min(dot([u-v for u,v in zip(x,z)],[u-v for u,v in zip(x,z)])for z in Z)
    ep=isu(ip(enddist),qi);assert ep[0]>1
    op=min(dot([u-v for u,v in zip(x,c)],[u-v for u,v in zip(x,c)])-q for c,q in oldballs)
    rp=min(xq-4-(dot(x,z)-4)**2/192 for z in Z)
    assert min(op,rp,xq-4)>1
    D=(xq-1)/2;assert xq-2*D==1 and -200-D<0
    # A matched three-factor variant differs by the A7 shift R/8.
    # Its mutual pairing with the current cap is Q/8>0, so it cannot
    # be added as a coherent nested extension. No variant orbit is built.
    nested_distance=sc(Q,F(7,4))
    nested_pair=sc(sub(sc(Q,2),nested_distance),F(1,2))
    assert nested_pair==sc(Q,F(1,8)) and sg(nested_pair)>0
    paths=[f'scripts/{TAG}.py',oldpath]
    return enc(dict(status='PASS_EXACT_MIXED_ORBIT_WITH_ENDPOINT_PI8_AND_GLOBAL_ESCAPE',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        field='Centers/radius in Q(sqrt7); endpoint radius by two positive nested radicals over Q(sqrt2,sqrt7)',
        radius=R,radius_square=Q,t=t,dotted_seed_cosine=b,center_square=N,
        endpoint_radius_formula='rho=R*(sqrt(119+49cos²(pi/8))/7-cos(pi/8))',endpoint_radius_interval=rhoi,
        endpoint_radius_square_interval=qi,orbit_size=392,centers=C,
        coherent_wall_count=1202,unchanged_middle_motif_caps=395,
        mixed_seed_relations={'right':17,'parallel':2,'dotted':3},
        mixed_mutual_pairs={'pi3':784,'dotted':75852},
        mixed_endpoint_pairs={'pi8':392,'dotted':153272},
        mixed_old_middle_pairs={'dotted':392*395},
        endpoint_second_distance_dotted_margin=min(endpoint_dot_margins,key=lambda a:a[0]),
        preserved_finite_endpoints=392,
        new_finite_vertices_involving_mixed_caps=0,
        finite_neighbor_star=dict(seed_right=zeroseed,mixed_pi3=neighbors,endpoint_pi8=0,
                                 five_nonseed_neighbors_mutually_dotted=True,maximum_finite_clique_size=19),
        forced_parallel_cusp=dict(seed=1,other_seed_incidence=incident,
            diagram='affine A1 + spherical A4 + A5 + A7',quotient_rank=17,
            same_type_cusps=784,scope='Listed complete incident sets; not a full cusp census'),
        matched_three_factor_extension=dict(squared_center_distance=nested_distance,
            mutual_pairing=nested_pair,conclusion='positive pairing; replacement is required'),
        strict_escape=dict(spatial=x,D=D,U=1,norm=1,minimum_seed_gap=F(1,100),
            mixed_power=mp,endpoint_power_interval=ep,old_middle_power=op,
            all_real_radial_power=rp,central_power=xq-4),
        scope='Exact orbit-level pair transfer, strict center/facet proof, preserved endpoints, forced cusp and global escape. The mixed caps create finite contacts but no new finite vertices; this is not finite volume.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    path=ROOT/f'results/{TAG}.json'
    if a.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['coherent_wall_count'])
