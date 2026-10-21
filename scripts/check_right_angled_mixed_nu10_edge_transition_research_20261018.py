#!/usr/bin/env python3
"""Tiny exact coefficient checks for the written nu10/edge-transition proof."""
import argparse
from fractions import Fraction as Q
from math import comb
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CERT=ROOT/'results/right_angled_mixed_nu10_edge_transition_research_20261018.json'


def generate():
    # In R7: p4, finite/ideal edges6, no k3 corner can have three bad triangles.
    no_zero=[]
    for a in range(1,4):
        for b in range(a,4):
            for c in range(b,4):
                for d in range(c,4):
                    degrees=(a,b,c,d)
                    if sum(degrees)==6:
                        maximum_bad=sum((0,0,1,2)[k] for k in degrees)
                        assert maximum_bad<3;no_zero.append(list(degrees))
    assert no_zero==[[1,1,1,3],[1,1,2,2]]
    potentials={}
    for n in (11,12):
        values=[]
        for k in range(n+1):
            phi=3*comb(n-k,3)-comb(k,2)*(n-k)
            assert phi==Q(n-2,2)*(n-k)*(n-1-2*k)
            values.append(phi)
        potentials[str(n)]=values
    assert potentials['12'][4]==120 and potentials['12'][7]==-75
    # Uniform face bound: -(6e)<=3A0-A2<=30e for e=nu-10>=1.
    assert Q(2*11-4,3)==6
    assert 6*11-36==30       # one cusp
    assert 6*11-48<=30      # two cusps
    assert 6*11-54<=30      # at least three cusps
    assert 6*12-12==30*2    # compact nu12, sharp upper coefficient
    return dict(status='EXACT_NU10_EDGE_TRANSITION_COEFFICIENTS',
        R7=dict(facets=7,cusps=3,finite_vertices=4,bad_triangles=3,
            all_ideal_triangles=1,quadrilaterals_finite3_ideal1=3,
            finite_vertex_ideal_degrees=[0,2,2,2],no_zero_degree_candidates=no_zero),
        ambient_degree_potential=potentials,
        universal_potential_budget=dict(lower_coefficient=-6,upper_coefficient=30),
        scope=dict(universal_geometric_proof_written=True,
            full_classification_of_Q8_isometry_not_claimed=True,
            new_linear_optimizer_used=False,dimension_exclusion=False,lean_compilation_performed=False))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write',action='store_true')
    args=parser.parse_args();record=generate()
    if args.write:CERT.write_text(json.dumps(record,indent=2)+'\n')
    else:
        with CERT.open() as stream:saved=json.load(stream)
        assert saved==record
    print('PASS nu10 degree cases and potential coefficients; no optimizer or dimension claim')


if __name__=='__main__':main()
