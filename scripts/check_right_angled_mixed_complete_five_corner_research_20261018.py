#!/usr/bin/env python3
"""Exact local models for a complete K5 corner; universal geometry is written.

Replay uses small rational arithmetic, no optimizer or Lean.
"""
import argparse
from fractions import Fraction as Q
from itertools import combinations
import json
from math import comb,isqrt
from pathlib import Path
from check_right_angled_mixed_weighted_triangles_research_20261018 import solve_rows
from check_right_angled_mixed_low_facet_rigidity_research_20261018 import dot,model
from check_right_angled_mixed_incidence_feasibility_research_20261018 import generate as aggregate_model

ROOT=Path(__file__).resolve().parents[1]
CERT=ROOT/'results/right_angled_mixed_complete_five_corner_research_20261018.json'


def generate():
    # A qualifying triangle has a missing complementary edge.  These are
    # precisely the incoming complete-corner positions in the written proof.
    pairs=list(combinations(range(5),2));triangles=list(combinations(range(5),3))
    maximum=0;maximizers=0;geometric_maximum=0;geometric_maximizers=0
    for mask in range(1<<len(pairs)):
        edges={pair for bit,pair in enumerate(pairs) if mask&(1<<bit)}
        qualifying=sum(all(pair in edges for pair in combinations(triple,2))
                       and tuple(i for i in range(5) if i not in triple) not in edges
                       for triple in triangles)
        if qualifying>maximum:maximum,maximizers=qualifying,1
        elif qualifying==maximum:maximizers+=1
        incoming=sum(tuple(i for i in range(5) if i not in triple) not in edges
                     and all(pair in edges for pair in pairs
                             if pair!=tuple(i for i in range(5) if i not in triple))
                     for triple in triangles)
        if incoming>geometric_maximum:geometric_maximum,geometric_maximizers=incoming,1
        elif incoming==geometric_maximum:geometric_maximizers+=1
    assert maximum==4 and geometric_maximum==1 and geometric_maximizers==10
    p5=model();base=[tuple(n) for n in p5['outward_unit_normals']]
    rho=base[-1];old=base[:15]
    reflected=[tuple(a+2*b for a,b in zip(n,rho)) for n in base[:5]]
    walls=old+reflected
    assert len(set(walls))==20
    assert all(dot(n,n)==1 for n in walls)
    assert all(dot(n,m) in (0,-1,-2) for n,m in combinations(walls,2))
    constraints=[([-Q(x) for x in n[1:]],Q(n[0])) for n in walls]
    vertices=set()
    for chosen in combinations(constraints,5):
        try:rows=solve_rows([a+[-b] for a,b in chosen],5)
        except StopIteration:continue
        y=tuple(row[-1] for row in rows)
        if all(sum(a*b for a,b in zip(left,y))<=right for left,right in constraints):vertices.add(y)
    source={tuple(map(Q,y)) for y in p5['coordinates']}
    def reflection(y):
        product=-2+sum(y);time=1-4*product
        return tuple((x-2*product)/time for x in y)
    expected=source|{reflection(y) for y in source}
    expected={y for y in expected if not(sum(y)==2 and sum(x*x for x in y)<1)}
    assert vertices==expected and len(vertices)==37
    finite=[y for y in vertices if sum(x*x for x in y)<1]
    ideal=[y for y in vertices if sum(x*x for x in y)==1]
    assert len(finite)==22 and len(ideal)==15
    assert len(finite)+len(ideal)==len(vertices)
    corner_edges=[];complete=[]
    wallset=set(walls)
    for y in finite:
        square=1-sum(x*x for x in y)
        num,den=isqrt(square.numerator),isqrt(square.denominator)
        assert num*num==square.numerator and den*den==square.denominator
        root=Q(num,den);q=(1/root,)+tuple(-x/root for x in y)
        active=[n for n in walls if dot(q,n)==0]
        assert len(active)==5 and all(dot(n,m)==0 for n,m in combinations(active,2))
        edges=sum(tuple(a-b-c for a,b,c in zip(q,n,m)) in wallset for n,m in combinations(active,2))
        corner_edges.append(edges)
        if edges==10:complete.append(y)
    assert len(complete)==2
    assert all(any(sum(a*b for a,b in zip(left,y))==right for y in vertices)
               for left,right in constraints)
    q=(1,0,0,0,0,0)
    cusp_vectors=[tuple(2*a-sum(n[k] for j,n in enumerate(base[:5]) if j!=i)
                       for k,a in enumerate(q)) for i in range(5)]
    assert all(dot(v,v)==0 for v in cusp_vectors)
    assert all(dot(cusp_vectors[i],base[j])==(-1 if i!=j else 0) for i in range(5) for j in range(5))
    assert all(dot(v,w)==-1 for v,w in combinations(cusp_vectors,2))
    assert all(reflected[i]==tuple(2*a-b for a,b in zip(cusp_vectors[i],base[i])) for i in range(5))
    # Preserve the frozen earlier packet.  Only D changes, so replay it and
    # check every inequality in which D occurs, plus the new incoming bound.
    aggregate=aggregate_model()
    a=[Q(aggregate['normalized_counts'][name]) for name in ['f']+[f'a{k}' for k in range(1,12)]]
    f=a[0];T,B,C,D=map(Q,(5590,9000,8000,4000));k5=Q(64000)
    d_slacks=[D-(17*a[5]-8*a[4]),a[5]-D,comb(12,5)*f-16*D,
              Q(42240)-10*D,8*C-16*D,120*T-160*D,36*B-80*D,k5-16*D]
    assert all(x>=0 for x in d_slacks)
    incoming_slack=(comb(12,5)*f+160*D)/11-k5
    assert incoming_slack>0
    assert k5==16*D==8*C  # no residual noncanonical K5 incidences; J=0
    return dict(status='EXACT_COMPLETE_K5_CORNER_LOCAL_CERTIFICATE',
        canonical_seed_vertices=26,canonical_seed_finite_vertices=16,canonical_seed_cusps=10,
        four_subset_cusp_vectors=[list(v) for v in cusp_vectors],
        width_opposite_product_coefficients={'lambda_mu':-1,'lambda':1,'mu':1},
        double_model=dict(walls=20,wall_choices_examined=comb(20,5),vertices=37,
            finite_vertices=22,ideal_vertices=15,complete_k5_corners=2,
            finite_corner_bad_edge_counts=sorted(corner_edges),
            complete_corner_projective_coordinates=[[str(x) for x in y] for y in sorted(complete)],
            outward_normals=[list(n) for n in walls],
            coordinates=[[str(x) for x in y] for y in sorted(vertices)]),
        conditional_aggregate_gap=dict(noncanonical_complete_corner_face_minimum=20,
            generic_noncanonical_face_minimum=17,canonical_face_minimum=16,
            extra_penalty_per_noncanonical_complete_corner_face=3),
        incoming_complete_corner_bound=dict(graphs_examined=1<<len(pairs),
            maximum_qualifying_triangles=maximum,number_of_maximizers=maximizers,
            maximum_geometrically_possible_incoming=geometric_maximum,
            number_of_geometric_maximizers=geometric_maximizers,
            noncanonical_five_face_finite_fraction='1/11',
            global_canonical_five_face_coefficient='160/11'),
        repaired_aggregate_assignment=dict(base_certificate=aggregate['status'],
            old_D5='3000',new_D5=str(D),total_K5=str(k5),
            changed_inequality_slacks=[str(x) for x in d_slacks],
            new_global_incoming_bound_slack=str(incoming_slack),
            all_complete_corners_assigned_to_canonical_five_faces=True,
            noncanonical_complete_corner_faces_J='0'),
        scope=dict(universal_proof_is_written=True,k5_alone_forces_p5=False,
            universal_upper_bound_two_complete_corners=False,
            dimension_ceiling_improved=False,lean_compilation_performed=False))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write',action='store_true')
    args=parser.parse_args();record=generate()
    if args.write:CERT.write_text(json.dumps(record,indent=2)+'\n')
    else:
        with CERT.open() as stream:saved=json.load(stream)
        assert record==saved,'Saved K5-corner certificate differs from exact replay'
    print('PASS K5 geometry, 15504 doubled-model wall choices, 1024 incoming graphs, repaired aggregate counts; no dimension claim or Lean')


if __name__=='__main__':main()
