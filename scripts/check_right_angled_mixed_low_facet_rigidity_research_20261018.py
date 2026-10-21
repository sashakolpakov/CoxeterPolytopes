#!/usr/bin/env python3
"""Replay the compact arithmetic/local-model certificate for RA4/RA5 rigidity.

The universal geometric proof is in the accompanying report, not this checker.
No Lean, optimizer, large files or numerical approximation are used.
"""
import argparse
from fractions import Fraction as Q
from itertools import combinations
import json
from math import comb
from pathlib import Path
from check_right_angled_mixed_weighted_triangles_research_20261018 import solve_rows

ROOT=Path(__file__).resolve().parents[1]
CERT=ROOT/'results/right_angled_mixed_low_facet_rigidity_research_20261018.json'


def dot(a,b):
    return -a[0]*b[0]+sum(x*y for x,y in zip(a[1:],b[1:]))


def model():
    d=5
    constraints=[]
    for i in range(d):constraints.append(([Q(int(i==j)) for j in range(d)],Q(0)))
    for i,j in combinations(range(d),2):
        constraints.append(([Q(int(k in (i,j))) for k in range(d)],Q(1)))
    constraints.append(([Q(1)]*d,Q(2)))
    vertices=set()
    for selected in combinations(constraints,d):
        try:rows=solve_rows([a+[-b] for a,b in selected],d)
        except StopIteration:continue
        point=tuple(row[-1] for row in rows)
        if (min(point)>=0 and all(point[i]+point[j]<=1 for i,j in combinations(range(d),2))
                and sum(point)<=2):vertices.add(point)
    expected={(Q(0),)*d}
    expected.update(tuple(Q(int(i==j)) for j in range(d)) for i in range(d))
    expected.update(tuple(Q(int(i in subset),2) for i in range(d))
                    for k in (3,4) for subset in combinations(range(d),k))
    expected.update(tuple(Q(2 if i==j else 1,3) for j in range(d)) for i in range(d))
    assert vertices==expected and len(vertices)==26
    norms={point:sum(x*x for x in point) for point in vertices}
    assert all(x<=1 for x in norms.values())
    assert sum(x<1 for x in norms.values())==16
    assert sum(x==1 for x in norms.values())==10
    facets=[]
    for a,b in constraints:
        face=[point for point in vertices if sum(x*y for x,y in zip(a,point))==b]
        finite=sum(norms[point]<1 for point in face)
        ideal=sum(norms[point]==1 for point in face)
        assert len(face)==10 and finite==ideal==5
        facets.append(dict(vertices=len(face),finite=finite,ideal=ideal))
    normals=[(0,)+tuple(int(i==j) for j in range(d)) for i in range(d)]
    normals.extend((1,)+tuple(-int(k in (i,j)) for k in range(d))
                   for i,j in combinations(range(d),2))
    normals.append((2,)+(-1,)*d)
    assert all(dot(n,n)==1 for n in normals)
    assert all(dot(n,m) in (0,-1) for n,m in combinations(normals,2))
    degrees=[sum(dot(n,m)==0 for m in normals) for n in normals]
    assert degrees==[10]*16
    return dict(walls=16,active_constraint_choices_examined=comb(16,5),
        vertices=26,finite_vertices=16,ideal_vertices=10,
        coordinates=[[str(x) for x in point] for point in sorted(vertices)],
        facet_vertex_counts=facets,orthogonal_neighbors_per_facet=degrees,
        outward_unit_normals=[list(n) for n in normals])


def generate():
    cases={str(c):dict(facets=15-c,degree_sum_lower=12*(15-c)-12*c,
                    degree_sum_upper=str(Q(3*(15-c)**2,4))) for c in range(4)}
    assert all(Q(row['degree_sum_lower'])>Q(row['degree_sum_upper'])
               for c,row in cases.items() if c!='3')
    assert cases['3']['degree_sum_lower']==108
    assert Q(cases['3']['degree_sum_upper'])==108
    assert Q(108,12)==9 and 9%2==1
    # Eleven facets with two/three cusps and total incidence24 have exactly two triples.
    assert next(k for k in range(12) if 3*k+2*(11-k)==24)==2
    assert 4*7-(3*(2*7-4-2*3)+4*3)==4
    return dict(status='EXACT_ARITHMETIC_AND_LOCAL_MODEL_CERTIFICATE',
        minsum15_low_cusp_degree_cases=cases,
        four_cusp_case=dict(facets=11,cusp_incidences=24,two_cusp_facets=9,
            three_cusp_facets=2,minimum_triangle_ridges_in_three_cusp_facet=4,
            possible_other_triangle_capable_neighbors=1),
        canonical_five_polytope=model(),
        scope=dict(published_minimum_sum15_is_external=True,
            equality_classification_in_written_proof=True,
            independent_geometric_review_complete=False,lean_compilation_performed=False,
            dimension_bound_improved=False,new_example_claim=False))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args();record=generate()
    if args.write:CERT.write_text(json.dumps(record,indent=2)+'\n')
    else:
        with CERT.open() as stream:saved=json.load(stream)
        assert record==saved,'Saved rigidity certificate differs from exact replay'
    print('PASS degree cases, all4368 P5 wall choices, 16finite/10ideal vertices, all16 P4 facet counts; geometry remains written; no Lean invoked')


if __name__=='__main__':main()
