#!/usr/bin/env python3
"""Exact opposition moments and a surviving aggregate assignment; no optimizer."""
import argparse
from fractions import Fraction as Q
from itertools import combinations
from math import comb
import json
from pathlib import Path
from check_right_angled_mixed_weighted_triangles_research_20261018 import ds_rows
from check_right_angled_mixed_incidence_feasibility_research_20261018 import macaulay

ROOT=Path(__file__).resolve().parents[1]
CERT=ROOT/'results/right_angled_mixed_opposition_moments_research_20261018.json'


def partitions(total,least=1):
    if not total:yield ()
    for first in range(least,total+1):
        for rest in partitions(total-first,first):yield (first,)+rest


def clique_vector(parts):
    result=[1]+[0]*12
    for size in parts:
        for k in range(12,0,-1):result[k]+=size*result[k-1]
    return result


def generate():
    mu=Q(5,4);p1,p2,p3=Q(1,4),Q(1,22),Q(1,220)
    assert p1-3*p2+2*p3==Q(27,220)
    assert p1-7*p2+12*p3==Q(-3,220)
    assert 3*p2-6*p3==Q(6,55)
    bound=1/mu+Q(9,176)/mu**3+Q(27,1760)/mu**4+2*Q(3,440)/mu**4
    assert bound==Q(28808,34375) and 220*bound==Q(1267552,6875)
    # One rational extremal-range configuration also replays all moments.
    x=[Q(1,4)]*3+[Q(-1,12)]*9
    d=[-sum(x[i] for i in triple) for triple in combinations(range(12),3)]
    average=lambda k:sum(value**k for value in d)/220
    assert sum(x)==0 and sum(t*t for t in x)==Q(1,4)
    assert average(1)==0 and average(2)==Q(9,176)
    assert average(3)==-Q(27,220)*sum(t**3 for t in x)
    assert average(4)==Q(3,440)-Q(3,220)*sum(t**4 for t in x)
    assert min(mu+t for t in d)==Q(1,2)
    assert sum(1/(mu+t) for t in d)/220<bound
    parts=list(partitions(12));assert len(parts)==77
    vectors={','.join(map(str,p)):clique_vector(p) for p in parts}
    assert clique_vector((1,)*12)==[comb(12,k) for k in range(13)]
    assert clique_vector((4,4,4))[2:6]==[48,64,0,0]
    assert clique_vector((6,6))[2:6]==[36,0,0,0]
    # Repeat the full earlier aggregate checks with the two stated changes.
    c=Q(1);free=[Q(4956,5),Q(100),Q(5),Q(1,20),Q(1,100000),Q(1,10**11)]
    known=[Q(0)]*6+free+[c,Q(0)]
    a=[sum(x*y for x,y in zip(row,known)) for row in ds_rows(12)]+free
    f=a[0];T,B,C,D=map(Q,(5590,9000,8000,4000));K=Q(64000)
    scale=396*10**22;checked={};beta=lambda k:comb(11,k-1)*2**(12-k)
    def require(family,slack,strict=False):
        slack=Q(slack);assert slack>0 if strict else slack>=0,(family,slack)
        checked.setdefault(family,[]).append(slack)
    for value in a:require('positive_faces',value,True)
    assert a[1]==6*f+1024 and a[3]==5*a[2]-55*f-14080
    E=66*f+22528-5*a[2];quad=11264-2*T-3*E;pent=a[2]-T-E-quad
    polygons=[T,E,quad,pent]
    assert sum(polygons)==a[2] and T+3*quad+5*pent==66*f
    assert 2*T+3*E+quad==11264
    for value in polygons:require('positive_polygons',value,True)
    minima={3:6,4:10,5:16,6:27,7:46,8:82,9:218,10:1923,11:183968}
    nus={3:9,4:15,5:26,6:45,7:81,8:217,9:1922,10:183967,11:1674688396}
    for k,v in minima.items():require('facet_minima',(13-k)*a[k-1]-v*a[k])
    for k,v in nus.items():require('facet_cusp_minima',(13-k)*a[k-1]+beta(k)-v*a[k])
    for k in range(2,7):
        require('strict_nk',Q(comb(12,k),comb(6,k))*a[k]-comb(12,k)*f-beta(k),True)
        for l in range(1,k):require('strict_nk',Q(comb(6,l),comb(6,k))*a[k]-a[l],True)
    for slack in (Q(22,5)*a[2]-66*f-Q(2123,1280)*11264,
                  f-Q(1474,21),Q(1298,9)-f,T-Q(331,437)*a[2],
                  T-Q(1986,59)*f,B-Q(2515,59)*f,C-Q(627,1475)*f):
        require('strict_weighted_densities',slack,True)
    for slack in (B-(7*a[3]-10*a[2]),B-(10*a[3]-10*a[2]-beta(3)),
                  C-(11*a[4]-9*a[3]),C-(16*a[4]-9*a[3]-beta(4)),
                  D-(17*a[5]-8*a[4])):require('canonical_gaps',slack)
    for k,count,v,w in ((3,B,2,3),(4,C,5,5),(5,D,16,10)):
        for slack in (a[k]-count,comb(12,k)*f-v*count,beta(k)-w*count):
            require('canonical_caps',slack)
    for slack in (10*T-6*B,9*B-10*C,8*C-16*D,45*T-30*C,120*T-160*D,36*B-80*D):
        require('canonical_cofaces',slack)
    weights=[Q(8000,99),Q(125,36),Q(5,2)];weights.append(f-sum(weights))
    for value in weights:require('positive_corner_weights',value,True)
    assert sum(weights)==f
    assert 66*weights[0]+48*weights[1]+36*weights[2]==T
    assert 220*weights[0]+64*weights[1]==2*B
    assert 495*weights[0]==5*C and 792*weights[0]==K==16*D==8*C
    require('corner_fraction',weights[0]-Q(209,1534)*f,True)
    require('incoming_k5',(792*f+160*D)/11-K)
    require('opposition_moment_triangle_cap',220*bound*f-2*B)
    # J=0 is compatible with the refined facet gap because all K5 are canonical.
    require('refined_twenty_facet_gap_J_zero',D-(17*a[5]-8*a[4]))
    for value in a+polygons+weights:assert (value*scale).denominator==1
    b=[scale*(f+2048)]+[scale*(a[k]+2**(11-k)*comb(11,k)) for k in range(1,11)]+[scale*(a[11]+1),Q(1)]
    h=[sum(comb(k,i)*(-1)**(k-i)*b[k] for k in range(i,13)) for i in range(13)]
    assert h==h[::-1] and h[0]==1 and all(value.denominator==1 for value in h)
    g=[h[0]]+[h[i]-h[i-1] for i in range(1,7)]
    assert all(value>0 for value in g)
    growth=[macaulay(int(g[i]),i) for i in range(1,6)]
    assert all(g[i+1]<=growth[i-1] for i in range(1,6))
    assert a[11]*scale>=1674688397 and scale>=127466960740760088
    return dict(status='EXACT_OPPOSITION_MOMENTS_AND_FEASIBLE_AGGREGATES',dimension=12,
        reciprocal_mean_upper=str(bound),triangle_to_finite_vertex_upper=str(220*bound),
        centered_moment_coefficients=dict(second='9/176',third_power_sum='-27/220',
            fourth_constant='3/440',fourth_power_sum='-3/220'),
        multipartite_partitions=77,multipartite_clique_vectors=vectors,
        normalized_counts={name:str(value) for name,value in zip(['f']+[f'a{k}' for k in range(1,12)],a)},
        canonical_and_triangle_counts=dict(T=str(T),B3=str(B),C4=str(C),D5=str(D),K5=str(K),J='0'),
        polygon_counts=[str(value) for value in polygons],corner_graph_weights=[str(value) for value in weights],
        integer_scale=scale,scaled_h=[str(value) for value in h],scaled_g=[str(value) for value in g],
        macaulay_upper_bounds=[str(value) for value in growth],
        checked_families={name:dict(rows=len(values),minimum_slack=str(min(values))) for name,values in checked.items()},
        scope=dict(universal_geometric_arguments_written=True,linear_optimizer_used=False,
            closed_canonical_K12_propagation_realizable=False,actual_face_gluing_verified=False,
            new_example_claim=False,dimension_exclusion=False,lean_compilation_performed=False))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write',action='store_true')
    args=parser.parse_args();record=generate()
    if args.write:CERT.write_text(json.dumps(record,indent=2)+'\n')
    else:
        with CERT.open() as stream:saved=json.load(stream)
        assert record==saved,'Saved opposition-moment certificate differs'
    print('PASS exact moments, 77 clique vectors, all aggregate rows and Macaulay checks; no optimizer, Lean or dimension claim')


if __name__=='__main__':main()
