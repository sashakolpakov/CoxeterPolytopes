#!/usr/bin/env python3
"""A complete binary cap orbit, compressed exact pair proof, and local scout.

The radius-one centers move three unit coordinates from half to full.
Two intersecting triples alternate with unit/golden sign parity. This
cuts every half-probe strictly, outside the previous literal-half family.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_weighted_binary_20261014'
A, B = {0, 1, 2}, {0, 3, 4}
PHI = (F(0), F(1))
def q(a=0, b=0): return F(a), F(b)
def add(a,b): return a[0]+b[0], a[1]+b[1]
def neg(a): return -a[0],-a[1]
def sub(a,b): return add(a,neg(b))
def mul(a,b): return a[0]*b[0]+a[1]*b[1],a[0]*b[1]+a[1]*b[0]+a[1]*b[1]
def sc(a,t): return a[0]*t,a[1]*t
def sign(a):
    x,y=a[0]+a[1]/2,a[1]/2
    if not y:return (x>0)-(x<0)
    if not x:return (y>0)-(y<0)
    if x*y>0:return (x>0)-(x<0)
    d=x*x-5*y*y
    assert d
    return ((x>0)-(x<0))*((d>0)-(d<0))
def total(xs):
    z=q()
    for x in xs:z=add(z,x)
    return z
W=[q(1)]*10+[q(2)]+[q(1,1)]*8
def dot(a,b):return total(sc(w,x*y)for w,x,y in zip(W,a,b,strict=True))
def dist(a,b):
    d=[x-y for x,y in zip(a,b,strict=True)]
    return dot(d,d)
def rel(d2):
    c=sc(sub(d2,q(2)),F(1,2))
    if c==q():return 2
    if c==q(F(1,2)):return 3
    if c==sc(PHI,F(1,2)):return 5
    if sign(sub(c,q(1)))>=0:return 0
    return -1
def center(mask):
    parity=(mask & (((1<<19)-1)^(1<<10))).bit_count()%2
    full=A if parity==0 else B
    return [F((-1 if(mask>>i)&1 else 1)*(2 if i in full else 1),2)
            for i in range(19)]
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(list,tuple)):return [enc(y)for y in x]
    if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
    return x

def run():
    assert len(A)==len(B)==3 and len(A&B)==1
    N=q(20,8);CN=q(F(29,4),2)
    assert dot(center(0),center(0))==CN
    # The central relation is dotted; all 38 seed relations are explicit.
    central=sc(sub(CN,q(5)),F(1,4));assert sign(sub(central,q(1)))>0
    seed={'right':3,'pi3':7,'pi4':1,'pi5':8,'dotted':19}
    # Reference displacements can be compressed into coordinate classes.
    # Same type: squared distance =4h_A+h_other+2h_heavy+(1+phi)h_gold.
    # Different type: base1 +4h_intersection+2h_symmetric_difference+
    # h_neither+2h_heavy+(1+phi)h_gold.
    same=Counter();cross=Counter();rows=[]
    for ha,hn,hh,hg in product(range(4),range(8),range(2),range(9)):
        if (ha+hn+hg)%2:continue
        d=add(q(4*ha+hn+2*hh),sc(q(1,1),hg))
        m=comb(3,ha)*comb(7,hn)*comb(8,hg)
        r=rel(d) if d!=q() else 'identity'
        assert r!=-1
        same[r]+=m
        rows.append(['same',ha,hn,hh,hg,d,r,m])
    for hi,hd,hn,hh,hg in product(range(2),range(5),range(6),range(2),range(9)):
        if (hi+hd+hn+hg)%2!=1:continue
        d=add(q(1+4*hi+2*hd+hn+2*hh),sc(q(1,1),hg))
        m=comb(4,hd)*comb(5,hn)*comb(8,hg)
        r=rel(d);assert r!=-1
        cross[r]+=m
        rows.append(['different',hi,hd,hn,hh,hg,d,r,m])
    assert sum(same.values())==sum(cross.values())==2**18
    assert same['identity']==1
    # Each probe is strictly cut by its matched cap, irrespective of signs.
    matched_power=q(F(-1,4))
    assert sub(dist(center(0),[F(1,2)]*19),q(1))==matched_power
    # Full vertex-centered sqrt2 orbit: H-1 is either 0,1 or >1.
    endpoint_counts=Counter()
    for h1,hh,hg in product(range(11),range(2),range(9)):
        H=add(q(h1+2*hh),sc(q(1,1),hg))
        if H==q():continue
        c=sub(H,q(1)); assert c==q() or sign(sub(c,q(1)))>=0
        endpoint_counts['right' if c==q() else 'parallel_or_dotted']+=comb(10,h1)*comb(8,hg)
    # Every endpoint/binary pair is dotted, proved by a uniform distance.
    min_cross=q(F(17,4),2)
    assert dist(center(0),[F(1)]*19)==min_cross
    # min_cross > (1+sqrt2)^2=3+2sqrt2, using sqrt2<3/2.
    assert sign(sub(min_cross,q(6)))>0
    assert sign(sub(N,q(16)))>0 # central/endpoint certainly dotted
    # A finite neighbor list for the representative binary wall.
    neighbor_masks=[]
    # Finite contacts require <=1 golden, <=3 unit flips; exhaustive
    # formulas above certify that no other displacement can be finite.
    for ku in range(4):
        for uu in combinations(range(10),ku):
            for hh in (0,1):
                for gg in [()]+[(j,)for j in range(11,19)]:
                    mask=sum(1<<j for j in uu+gg)+(hh<<10)
                    if not mask:continue
                    r=rel(dist(center(0),center(mask)))
                    if r in (2,3,5):neighbor_masks.append((mask,r))
    finite_expected=sum(same[k]+cross[k]for k in (2,3,5))
    assert len(neighbor_masks)==finite_expected
    # Exact common escape, despite strict coverage of every half-probe.
    x=[F(0)]*10+[F(99,100)]*9
    xq=dot(x,x)
    binary_power=sub(add(q(F(19,4)),sc(q(10,8),F(49,100)**2)),q(1))
    endpoint_power=sub(add(q(10),sc(q(10,8),F(1,100)**2)),q(2))
    central_power=sub(xq,q(4))
    assert binary_power==sub(dist(x,center(0)),q(1))
    assert all(sign(sub(p,q(1)))>0 for p in (binary_power,endpoint_power,central_power))
    # All binary centers have three full unit coordinates. At x the unit
    # contribution is19/4; matching the other signs minimizes distance.
    D=sc(sub(xq,q(1)),F(1,2))
    return enc(dict(status='PASS_EXACT_PARITY_TRIPLE_CAP_ORBIT_WITH_GLOBAL_ESCAPE',
        sources={f'scripts/{TAG}.py':sha256(Path(__file__).read_bytes()).hexdigest()},
        full_unit_triples=[sorted(A),sorted(B)],type_rule='unit-plus-golden sign parity; heavy sign omitted',
        orbit_size=2**19,radius=1,center_square=CN,seed_relations=seed,
        central_cosine=central,reference_same_type_counts=dict(same),reference_other_type_counts=dict(cross),
        compressed_pair_cases=rows,finite_neighbor_masks=neighbor_masks,
        all_half_probes_strictly_cut=2**19,matched_probe_power=matched_power,
        endpoint_radius_square=2,endpoint_orbit_size=2**19,
        endpoint_reference_pair_counts=dict(endpoint_counts),binary_endpoint_min_distance_square=min_cross,
        coherent_actual_facets=39+2**20,preserved_finite_endpoints=2**19,
        strict_escape=dict(spatial=x,D=D,U=1,norm=1,minimum_seed_gap_normalized=F(1,100),
            binary_power=binary_power,endpoint_power=endpoint_power,central_power=central_power),
        scope='Exact compressed pair proof for the entire parity-triple radius-one orbit and all sqrt2-radius vertex caps, with retained39-wall core. All half-probes are strictly cut and all old finite endpoints remain, but the common escape proves infinite volume. Local finite-vertex scout, if run, is separate.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(p.read_text())==out
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['coherent_actual_facets'],len(out['finite_neighbor_masks']))
