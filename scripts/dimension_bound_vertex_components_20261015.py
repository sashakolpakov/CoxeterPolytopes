#!/usr/bin/env python3
"""Exact scalar/partition checks for compact Coxeter dimension necessities.

The universal geometric proof is in dimension_bound_lane_20261015.md.
This finite checker certifies its rational identities, exact partition maxima,
class thresholds and quantitative vertex-distribution consequences.
"""
import argparse,json
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import comb
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='dimension_bound_vertex_components_20261015'

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def U(d):return F(4*(d-1),d-2)if d%2==0 else F(4*d,d-1)
def T(d):return F(d*(d-6),8)if d%2==0 else F((d-1)*(d-5),8)
def maximum(d,r):
    q,s=divmod(d,r)
    return q*comb(r,2)+comb(s,2)
@lru_cache(None)
def dp(d,r):
    if not d:return 0
    return max(comb(i,2)+dp(d-i,r)for i in range(1,min(d,r)+1))
@lru_cache(None)
def partition_count(d,r):
    if d==0:return 1
    if r==0:return 0
    return partition_count(d,r-1)+(partition_count(d-r,r)if d>=r else 0)
def enc(v):
    if isinstance(v,F):return str(v)
    if isinstance(v,list):return[enc(x)for x in v]
    if isinstance(v,dict):return{k:enc(x)for k,x in v.items()}
    return v

def run():
    moments=[];thresholds=[];checked=0
    for d in range(4,101):
        b=comb(d,2);assert b*(F(5)/U(d)-1)==T(d)
        if d<=29:moments.append(dict(dimension=d,nikulin_strict_upper=U(d),
            strict_lower_average_connected_pairs=T(d),
            three_free_strict_upper_average_components=d-T(d),rank_partitions=partition_count(d,d)))
        for r in range(1,30):
            exact=maximum(d,r);assert dp(d,r)==exact;checked+=1
            if d>=4*r+1:assert exact<=T(d)
            if d<=29 and r<=7:
                components=(d+r-1)//r
                thresholds.append(dict(dimension=d,maximum_component_rank=r,maximum_connected_pairs=exact,
                    general_compact_excluded=exact<=T(d),
                    minimum_vertex_components=components,maximum_diagram_edges=d-components,
                    three_free_excluded=(d-components)<=T(d)))
    # The exceptional first dimension needs the exact partition maximum.
    boundaries=[]
    for r in range(1,101):
        d=4*r+1
        assert maximum(d,r)==2*r*(r-1)==T(d)
        boundaries.append(dict(r=r,first_excluded_dimension=d,exact_partition_maximum=maximum(d,r)))
        # Thereafter use the elementary pointwise bound C<=d(r-1)/2.
        for d in (4*r+2,4*r+3,4*r+4,4*r+5):
            if d%2==0:
                assert T(d)-F(d*(r-1),2)==F(d*(d-4*r-2),8)>=0
            else:
                t=d-4*r-3
                assert T(d)-F(d*(r-1),2)==F(d*d-(4*r+2)*d+5,8)==F((4*r+3+t)*(1+t)+5,8)>0
    classes=[]
    for r in range(1,8):
        possible=[d for d in range(4,30)if comb(d,2)>T(d)and maximum(d,r)>T(d)]
        flag=[d for d in possible if d-(d+r-1)//r>T(d)]
        classes.append(dict(maximum_component_rank=r,general_bound=4*r,
            bound_combined_with_Vinberg=min(4*r,29),three_free_bound=max(flag)))
    assert [x['three_free_bound']for x in classes]==[4,8,10,11,11,12,12]
    distributions=[]
    for d,r in [(9,2),(12,2),(20,2),(20,3),(20,4),(29,3),(29,5),(29,7)]:
        c=maximum(d,r);b=comb(d,2)
        if T(d)>=c:
            p=(T(d)-c)/(b-c)
            distributions.append(dict(dimension=d,rank_cutoff=r,low_vertex_maximum=c,
                strict_lower_fraction_with_larger_component=p))
    assert next(x['strict_lower_fraction_with_larger_component']for x in distributions if x['dimension']==20 and x['rank_cutoff']==2)==F(5,36)
    assert next(x['strict_lower_fraction_with_larger_component']for x in distributions if x['dimension']==29 and x['rank_cutoff']==3)==F(4,27)
    assert 12-T(12)==3 and 11-T(11)==F(7,2)
    return enc(dict(status='PASS_EXACT_COMPACT_VERTEX_COMPONENT_DIMENSION_NECESSITIES',
        sources={f'scripts/{TAG}.py':sha(f'scripts/{TAG}.py')},
        arithmetic='Standalone Fraction/integer dynamic programming; no geometry or existing producer imports.',
        theorem_scope='Compact finite-sided hyperbolic Coxeter polytopes, all facets and all labels. The three-free corollary additionally requires every missing face to have size2. No arithmetic, fixed-lattice, facet-count or source assumptions.',
        rigorous_new_derivation='Average connected-pair count>T_d; maximum vertex-component rank r implies d<=4r. In the three-free class average number of components<d-T_d, giving the saved refined class bounds.',
        novelty='Derived use of classical Nikulin/Gauss-Bonnet and Alexandrov quadrilateral charging. No publication-priority claim or improvement of unrestricted29/995 bounds.',
        exact_partition_checks=checked,moments=moments,class_bounds=classes,rank_distribution=distributions,
        dimension_rank_thresholds=thresholds,odd_boundary_checks=boundaries,
        universal_identities=dict(first_boundary='Cmax(4r+1,r)=2r(r-1)=T_(4r+1)',
            even_tail='T_d-d(r-1)/2=d(d-4r-2)/8>=0 for even d>=4r+2',
            odd_tail='T_d-d(r-1)/2=((4r+3+t)(1+t)+5)/8>0 for odd d=4r+3+t,t>=0'),
        limitation='The finite tables do not themselves prove the geometric premise; that proof is in the report. The universal tail follows the displayed algebra, not a finite numerical cutoff. No noncompact extension and no unjustified gadget-occurrence/LP theorem is claimed.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    d=run();out=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text());print('PASS saved exact dimension inequalities')
    else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['exact_partition_checks'],d['class_bounds'])
