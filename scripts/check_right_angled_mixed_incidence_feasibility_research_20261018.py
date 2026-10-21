#!/usr/bin/env python3
"""Exact feasible aggregate counts in dimension12; never a geometric example.

Pure rational/integer arithmetic; no optimizer, Lean or large computation.
"""
import argparse
from fractions import Fraction as Q
from math import comb
import json
from pathlib import Path
from check_right_angled_mixed_weighted_triangles_research_20261018 import ds_rows

ROOT=Path(__file__).resolve().parents[1]
CERT=ROOT/'results/right_angled_mixed_incidence_feasibility_research_20261018.json'


def macaulay(value,degree):
    remaining,total,previous=value,0,None
    for j in range(degree,0,-1):
        if not remaining:break
        low,high=j-1,max(j,2)
        while comb(high,j)<=remaining:high*=2
        if previous is not None:high=min(high,previous)
        while low+1<high:
            middle=(low+high)//2
            if comb(middle,j)<=remaining:low=middle
            else:high=middle
        assert low<previous if previous is not None else True
        remaining-=comb(low,j);total+=comb(low+1,j+1);previous=low
    assert remaining==0
    return total


def generate():
    n=12;c=Q(1)
    free=list(map(Q,(990,100,5)))+[Q(1,20),Q(1,100000),Q(1,10**11)]
    known=[Q(0)]*6+free+[c,Q(0)]
    a=[sum(x*y for x,y in zip(row,known)) for row in ds_rows(n)]+free
    f=a[0];T,B,C,D=map(Q,(5590,9000,8000,3000))
    scale=396*10**22
    checked={}
    def require(family,slack,strict=False):
        slack=Q(slack)
        assert slack>0 if strict else slack>=0,(family,slack)
        checked.setdefault(family,[]).append(slack)
    beta=lambda k:comb(11,k-1)*2**(12-k)
    for x in a:require('positive_face_counts',x,True)
    assert a[1]==6*f+1024*c
    assert a[3]==5*a[2]-55*f-14080*c
    # Explicit polygon classes: finite/ideal corners (1,2),(0,3),(3,1),(5,0).
    ideal_tri=66*f+2*11264*c-5*a[2]
    quad=11264*c-2*T-3*ideal_tri
    pent=a[2]-T-ideal_tri-quad
    polygons=[T,ideal_tri,quad,pent]
    for x in polygons:require('positive_polygon_types',x,True)
    assert sum(polygons)==a[2]
    assert T+3*quad+5*pent==66*f
    assert 2*T+3*ideal_tri+quad==11264*c
    # All published facet minima used by the packet, plus sharp 3/4/5 minima.
    minima={3:6,4:10,5:16,6:27,7:46,8:82,9:218,10:1923,11:183968}
    nus={3:9,4:15,5:26,6:45,7:81,8:217,9:1922,10:183967,11:1674688396}
    for k,value in minima.items():require('facet_minima',(13-k)*a[k-1]-value*a[k])
    for k,value in nus.items():require('facet_plus_cusp_minima',(13-k)*a[k-1]+beta(k)*c-value*a[k])
    # Nikulin–Khovanskii in its used range 2<=k<=6, 0<=l<k.
    for k in range(2,7):
        require('strict_nk_vertex',Q(comb(12,k),comb(6,k))*a[k]-comb(12,k)*f-beta(k)*c,True)
        for l in range(1,k):require('strict_nk_positive_faces',Q(comb(6,l),comb(6,k))*a[k]-a[l],True)
    weight=Q(2123,1280)
    require('strict_weighted_two_face',Q(22,5)*a[2]-66*f-weight*11264*c,True)
    require('strict_vertex_ratio_lower',f-Q(1474,21)*c,True)
    require('strict_vertex_ratio_upper',Q(1298,9)*c-f,True)
    require('strict_triangle_face_density',T-Q(331,437)*a[2],True)
    require('strict_triangle_finite_density',T-Q(1986,59)*f,True)
    # Both facet and minimum-sum gaps, using established equality rigidity.
    require('canonical_b3_facet_gap',B-(7*a[3]-10*a[2]))
    require('canonical_b3_sum_gap',B-(10*a[3]-10*a[2]-beta(3)*c))
    require('canonical_c4_facet_gap',C-(11*a[4]-9*a[3]))
    require('canonical_c4_sum_gap',C-(16*a[4]-9*a[3]-beta(4)*c))
    require('canonical_d5_facet_gap',D-(17*a[5]-8*a[4]))
    require('strict_b3_density',B-Q(2515,59)*f,True)
    require('strict_c4_density',C-Q(627,1475)*f,True)
    for k,count,finite,ideal in ((3,B,2,3),(4,C,5,5),(5,D,16,10)):
        require('canonical_face_upper',a[k]-count)
        require('canonical_finite_incidence',comb(12,k)*f-finite*count)
        require('canonical_ideal_incidence',beta(k)*c-ideal*count)
    for slack in (10*T-6*B,9*B-10*C,8*C-16*D,45*T-30*C,120*T-160*D,36*B-80*D):
        require('canonical_coface_incidence',slack)
    # An explicit multiset of finite-corner graphs realizes T,2B,5C.
    graph_weights=[Q(8000,99),Q(125,36),Q(5,2)]
    graph_weights.append(f-sum(graph_weights))
    for x in graph_weights:require('positive_corner_graph_weights',x,True)
    assert sum(graph_weights)==f
    assert 66*graph_weights[0]+48*graph_weights[1]+36*graph_weights[2]==T
    assert 220*graph_weights[0]+64*graph_weights[1]==2*B
    assert 495*graph_weights[0]==5*C
    require('p5_corner_clique_capacity',792*graph_weights[0]-16*D)
    require('strict_k4_corner_fraction',graph_weights[0]-Q(209,1534)*f,True)
    # Scale to actual integers, and check the entire truncated simple g-vector.
    scaled=[x*scale for x in a]
    for x in scaled+[p*scale for p in polygons]+[w*scale for w in graph_weights]:
        assert x.denominator==1
    b=[scale*(f+2048*c)]+[scale*(a[k]+2**(11-k)*comb(11,k)*c) for k in range(1,11)]+[scale*(a[11]+c),Q(1)]
    h=[sum(comb(k,i)*(-1)**(k-i)*b[k] for k in range(i,13)) for i in range(13)]
    assert h==h[::-1] and h[0]==1 and all(x.denominator==1 for x in h)
    g=[h[0]]+[h[i]-h[i-1] for i in range(1,7)]
    assert all(x>0 for x in g)
    macaulay_bounds=[macaulay(int(g[i]),i) for i in range(1,6)]
    assert all(g[i+1]<=macaulay_bounds[i-1] for i in range(1,6))
    assert scaled[11]>=1674688397 and scale>=127466960740760088
    return dict(status='EXACT_FEASIBLE_AGGREGATE_ASSIGNMENT_NO_REALIZATION_CLAIM',dimension=12,
        normalization='c=1; a[0] denotes finite vertices f',integer_scale=scale,
        normalized_counts={name:str(x) for name,x in zip(['f']+[f'a{k}' for k in range(1,12)],a)},
        canonical_and_triangle_counts={'T':str(T),'B3':str(B),'C4':str(C),'D5':str(D)},
        polygon_classes=[dict(finite=p,ideal=q,count=str(x)) for (p,q),x in zip(((1,2),(0,3),(3,1),(5,0)),polygons)],
        finite_corner_graphs=[dict(graph=name,count=str(x)) for name,x in zip(('K12','K4,4,4','K6,6','empty12'),graph_weights)],
        checked_families={name:dict(rows=len(values),minimum_slack=str(min(values))) for name,values in checked.items()},
        facet_minima={str(k):v for k,v in minima.items()},
        facet_plus_cusp_minima={str(k):v for k,v in nus.items()},
        scaled_truncated_h=[str(int(x)) for x in h],scaled_truncated_g=[str(int(x)) for x in g],
        macaulay_upper_bounds=[str(x) for x in macaulay_bounds],
        scope=dict(actual_face_gluing_verified=False,lorentz_wall_realization_verified=False,cusp_box_width_compatibility_verified=False,
            k5_implies_p5_claim=False,new_example_claim=False,dimension_exclusion=False,
            linear_optimizer_used=False,lean_compilation_performed=False))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write',action='store_true')
    args=parser.parse_args();record=generate()
    if args.write:CERT.write_text(json.dumps(record,indent=2)+'\n')
    else:
        with CERT.open() as stream:saved=json.load(stream)
        assert record==saved,'Saved aggregate witness differs from exact replay'
    print('PASS exact aggregate counts, polygon/corner mixtures, stated NK/gaps/incidences, g/Macaulay checks; no realization claim; no LP or Lean')


if __name__=='__main__':main()
