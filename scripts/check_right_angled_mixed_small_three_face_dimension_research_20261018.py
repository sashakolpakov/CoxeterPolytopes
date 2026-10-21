#!/usr/bin/env python3
"""Exact nine base cases for the written dimension-eight carrier theorem."""
import argparse
from fractions import Fraction as Q
from math import comb
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CERT=ROOT/'results/right_angled_mixed_small_three_face_dimension_research_20261018.json'


def sign_sqrt2(a,b):
    a,b=Q(a),Q(b)
    sign=lambda value:(value>0)-(value<0)
    if not b:return sign(a)
    if not a:return sign(b)
    if a*b>0:return sign(a)
    return sign(a)*sign(a*a-2*b*b)


def coefficients(k,m):
    C=comb(k,3);d=comb(k-1,2) if k>=3 else 0
    # Delta=A^2-|spatial|^2 with A=m(sqrt2-1)+C/20.
    return (3*m*m-Q(m*C,10)+Q(C*C-k*d*d,400)-m,
            -2*m*m+Q(m*C,10))


def generate():
    cases=[]
    for k in range(9):
        m=9-k;a,b=coefficients(k,m);assert sign_sqrt2(a,b)>0
        aa,bb=coefficients(k,m+1)
        assert sign_sqrt2(aa-a,bb-b)>0
        cases.append(dict(k=k,m=m,delta=[str(a),str(b)],
                          first_m_increment=[str(aa-a),str(bb-b)]))
    assert coefficients(9,0)==(0,0)
    # Every later m-increment increases by 2(sqrt2-1)^2>0.
    assert sign_sqrt2(6,-4)>0
    # For k>=9, C>=84; (C/10)(sqrt2-1)-1>0.
    assert sign_sqrt2(Q(-47,5),Q(42,5))>0
    assert Q(84,20)>0
    # r_i.r_j=-sinh(t_i)sinh(t_j)=-1 gives all three sinh(t)=1.
    return dict(status='EXACT_SMALL_THREE_FACE_DIMENSION_CARRIER_CERTIFICATE',
        opposition_edge_weight='1/20',finite_edge_weight='1',
        n9_base_cases=cases,full_k9_case=dict(delta='0',positive_time_coefficient='21/5'),
        higher_k_formula='Delta=m^2*t^2+m*((C/10)*t-1)+d^2*k*(k-9)/3600',
        conclusion=dict(if_every_three_face_has_facet_plus_cusp_sum_at_most=10,
                        ambient_dimension_at_most=8),
        scope=dict(universal_lorentz_and_carrier_proof_written=True,
            arbitrary_polytope_dimension_exclusion=False,
            positive_excess_fraction_proved=False,optimizer_used=False,
            lean_compilation_performed=False))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write',action='store_true')
    args=parser.parse_args();record=generate()
    if args.write:CERT.write_text(json.dumps(record,indent=2)+'\n')
    else:
        with CERT.open() as stream:saved=json.load(stream)
        assert record==saved
    print('PASS nine exact causal-vector base cases and monotonicity; conditional dimension<=8, no universal dimension improvement')


if __name__=='__main__':main()
