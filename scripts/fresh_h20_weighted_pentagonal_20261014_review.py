#!/usr/bin/env python3
"""Separate Fraction Q(sqrt(5)) replay of new pentagonal pair classes.

No producer/project imports. The old R2/R2 census and vertex-radius bound
are inherited. New pairs are reduced by coordinate permutations before
their arithmetic is checked; distance-two pairs are checked directly.
"""
from fractions import Fraction as F
from itertools import combinations
from collections import Counter
from pathlib import Path
from hashlib import sha256
import argparse
import json

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_weighted_pentagonal_20261014'

def q(a=0, b=0): return (F(a), F(b))
def add(x,y): return (x[0]+y[0],x[1]+y[1])
def neg(x): return (-x[0],-x[1])
def sub(x,y): return add(x,neg(y))
def mul(x,y): return (x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def scale(x,a): return (x[0]*a,x[1]*a)
def sq(x): return mul(x,x)
def sign(x):
    a,b=x
    if not b:return (a>0)-(a<0)
    if not a:return (b>0)-(b<0)
    if a*b>0:return (a>0)-(a<0)
    norm=a*a-5*b*b
    assert norm
    return ((a>0)-(a<0))*((norm>0)-(norm<0))
def summ(xs):
    z=q()
    for x in xs:z=add(z,x)
    return z

PHI=q(F(1,2),F(1,2))
R=scale(sub(PHI,q(1)),4)
W=[q(1)]*10+[q(2)]+[sq(PHI)]*8

def center(kind, profile=(), flips=()):
    if kind==5:
        v=[q(1)]*11+[q(1) if j in profile else sub(scale(PHI,2),q(3)) for j in range(8)]
    elif kind==2:
        full,zero=profile
        v=[q(int(i!=zero)) for i in range(10)]+[q(1)]+[q(1) if j==full else sub(q(2),PHI) for j in range(8)]
    else:
        assert kind==1
        v=[q(F(1,2))]*19
    for j in flips:v[11+j]=neg(v[11+j])
    return v

def dot(v,w):return summ(mul(weight,mul(x,y))for weight,x,y in zip(W,v,w,strict=True))
def dist(v,w):
    d=[sub(x,y)for x,y in zip(v,w,strict=True)]
    return dot(d,d)
def relation(d,r,s):
    base=add(sq(r),sq(s));twors=scale(mul(r,s),2)
    for label,c in [('2',q()),('3',q(F(1,2))),('5',scale(PHI,F(1,2)))]:
        if d==add(base,mul(twors,c)):return label
    v=sign(sub(d,add(base,twors)))
    return 'dotted'if v>0 else'parallel'if v==0 else'forbidden'

def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):encode(v)for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [encode(v)for v in x]
    return x

def run():
    primary=json.loads((ROOT/f'results/{TAG}.json').read_text())
    assert primary['status']=='PASS_PENTAGONAL_MOTIF_AND_FOUR_FAMILY_COVERING_OBSTRUCTION'
    for path,h in primary['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h
    old=json.loads((ROOT/'results/fresh_h20_seed_deformation_20261014.json').read_text())
    assert old['R2_nonvertex_family']['cutters_per_half_probe']==88
    assert old['R2_nonvertex_family']['adjacent_golden_probe_pair_counts']=={'forbidden':6484,'right':1260}
    assert old['R4_vertex_family']['maximum_compatible_caps']==512

    # New all-label coordinate gate, independently isolated radicals.
    plo,phi=F(1618033,10**6),F(1618034,10**6)
    slo,shi=F(1414213,10**6),F(1414214,10**6)
    assert plo*plo-plo-1<0<phi*phi-phi-1
    assert slo*slo<2<shi*shi
    assert sign(sub(sq(R),q(4)))>0
    assert sign(sub(scale(sq(R),F(1,2)),sq(PHI)))>0
    lower=(2*slo)/(4*phi-4)-F(1,2)
    upper=(2*shi)/(4*plo-4)-F(1,2)
    assert F(1,2)<lower<upper<slo/2
    assert sub(scale(PHI,2),scale(R,F(1,2)))==q(2)

    sets=[tuple(c)for k in (2,3)for c in combinations(range(8),k)]
    assert len(sets)==84
    profiles2=[(j,z)for j in range(8)for z in range(-1,10)]
    half=center(1)
    for A in sets:
        power=sub(dist(center(5,A),half),sq(R))
        assert sign(power)<0
        assert sign(add(power,scale(sub(PHI,q(1)),2)))>0

    # Quotient by permutations fixing the changed golden bit. Unit zero
    # positions are symmetric; only their presence matters against R5.
    classes55={};classes52={}
    for A in sets:
        for B in sets:
            key=(len(A),len(B),0 in A,0 in B,len(set(A)&set(B)))
            if key not in classes55:classes55[key]=[0,A,B]
            classes55[key][0]+=1
        for p in profiles2:
            key=(len(A),0 in A,p[0]==0,p[0] in A,p[1]>=0)
            if key not in classes52:classes52[key]=[0,A,p]
            classes52[key][0]+=1
    census={}
    for label,classes,kind in [('5/5',classes55,5),('5/2',classes52,2)]:
        counts=Counter()
        for multiplicity,A,B in classes.values():
            rel=relation(dist(center(5,A),center(kind,B,(0,))),R,R if kind==5 else q(2))
            assert rel=='forbidden'
            counts[rel]+=multiplicity
        census[label]=dict(counts)
    assert census=={'5/5':{'forbidden':7056},'5/2':{'forbidden':7392}}

    c15=Counter();c12=Counter()
    for A in sets:
        rel=relation(dist(half,center(5,A,(0,))),q(1),R)
        assert rel==('forbidden'if len(A)==2 else'5'if 0 in A else'2')
        c15[rel]+=1
    for p in profiles2:
        rel=relation(dist(half,center(2,p,(0,))),q(1),q(2))
        assert (rel=='dotted')==(p[0]==0)
        c12[rel]+=1
    assert relation(dist(half,center(1,(),(0,))),q(1),q(1))=='forbidden'
    census['1/5']=dict(c15);census['1/2']=dict(c12);census['1/1']={'forbidden':1}
    for k,v in census.items():assert primary['adjacent_golden_pair_counts'][k]==v

    # Direct coordinate arithmetic, independently of producer's compact
    # distance-two integer formula.
    triples=[A for A in sets if len(A)==3]
    two=Counter()
    for A in triples:
        for B in triples:
            rel=relation(dist(center(5,A),center(5,B,(0,1))),R,R)
            compatible=not(set(A)&set(B)) and len((set(A)|set(B))&{0,1})==1
            assert (rel=='2')==compatible
            assert rel in ('2','forbidden')
            two[rel]+=1
    assert dict(two)==primary['distance_two_k3']['pair_counts']=={'forbidden':2896,'2':240}

    # The proposed four-cap motif and strict spacelike witness, all in
    # normalized spatial coordinates with diagonal weights W.
    motif=[(1,(),()),(5,(0,2,3),(0,)),(5,(4,5,6),(1,)),(1,(),(0,1))]
    centers=[center(k,p,f)for k,p,f in motif];radii=[q(1),R,R,q(1)]
    pairs=Counter()
    for i in range(4):
        assert relation(dot(centers[i],centers[i]),q(2),radii[i])!='forbidden'
        for j in range(i):pairs[relation(dist(centers[i],centers[j]),radii[i],radii[j])]+=1
    assert dict(pairs)=={'2':4,'5':1,'dotted':1}
    x=[q()]*11+[q(F(99,100))]*8
    assert all(sign(sub(sub(dist(x,c),sq(r)),q(1)))>0 for c,r in zip(centers,radii,strict=True))
    assert sign(sub(dot(x,x),q(5)))>0

    # In each Q8: e(A)<=2|A| and e(A)>=8|A|-1024.
    maximum=1024//6
    assert maximum==170 and 8*171-1024>2*171
    missed=2**19-2**11*maximum-512*275
    assert missed==35328==primary['mixed_covering']['minimum_not_strictly_cut']
    # Tangent R1 spheres differ by one unit bit. At a half-probe p,
    # gradient dotted with its ten signed half-unit entries is 1; the
    # quadratic coefficient is the squared displacement norm 10/4.
    assert 2*F(1,1)*F(1,2)==1 and 10*F(1,2)**2==F(5,2)
    paths=[f'scripts/{TAG}.py',f'results/{TAG}.json',f'scripts/{TAG}_review.py',
           'results/fresh_h20_seed_deformation_20261014.json',
           'results/fresh_h20_weighted_mixed_probe_20261014.json']
    return encode(dict(status='PASS_SEPARATE_NEW_PAIR_CLASSES_AND_PROOF_REVIEW',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        independent_basis='1,sqrt(5); stdlib Fraction only; no producer imports',
        adjacent_R5_R5_class_count=len(classes55),adjacent_R5_R2_class_count=len(classes52),
        adjacent_pair_counts=census,direct_distance_two_pair_counts=dict(two),
        motif_pair_counts=dict(pairs),max_covered_per_Q8=maximum,minimum_not_strictly_cut=missed,
        heavy_complement_cosine_bounds=[lower,upper],
        boundary_power_increment='epsilon+(5/2)*epsilon^2',
        scope='Separate arithmetic for new pair classes, distance-two restrictions, motif and counting. '
              'Source/proof review of complete R5 alphabet, graph charging, boundary cases and finite-volume implication. '
              'Inherits prior R2 classification/7744-pair census and radius4 512-cap/275-probe bounds; '
              'not a second implementation of those prior enumerations.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    result=run();out=ROOT/f'results/{TAG}_review.json'
    if args.verify:assert json.loads(out.read_text())==result
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status']);print(result['adjacent_R5_R5_class_count'],result['adjacent_R5_R2_class_count'])
