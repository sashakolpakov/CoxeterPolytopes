#!/usr/bin/env python3
"""Exact finite inequalities for the written polygon/three-face coupling."""
import argparse
from fractions import Fraction as Q
from math import comb
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CERT=ROOT/'results/right_angled_mixed_polygon_three_face_coupling_research_20261018.json'


def generate():
    alpha,beta,weight=Q(331,3840),Q(143,384),Q(2123,1280)
    # These are the minimal regions u=0;u=1;v=0;v=1;u,v>=2.
    minimal=[(0,3,0),(1,3,1),(5,0,0),(3,1,0),(2,2,2)]
    gaps=[Q(u)+weight*v-Q(22,5)-alpha*u-beta*j for u,v,j in minimal]
    assert all(value>=0 for value in gaps) and gaps[-1]==gaps[-2]==0
    assert 1-alpha-beta>0 and weight-beta>0
    assert (Q(22,5)-1-2*weight)*3840==318
    assert alpha*3840==331 and beta*3840==1430
    gammas={11:Q(11,23),12:Q(13,27)};local={}
    for n,gamma in gammas.items():
        N=comb(n,3)
        slacks=[(1-gamma)*(N-comb(k,3))-k*comb(n-k,2) for k in range(n+1)]
        assert all(value>=0 for value in slacks) and slacks[4]==0
        local[str(n)]=dict(gamma=str(gamma),slacks=[str(value) for value in slacks])
    # Exceptional-face budget: p=4F-2nu-4 and cost=22(nu-10)+4(8-F).
    assert 22*(10-10)+4*(8-7)==4  # nu10,F7 has p4
    assert 24*11-184-8*10==0  # nu11 has F<=10
    assert 16*12-184>0       # nu>=12,F<=nu; surplus grows thereafter
    assert Q(14*27+26,4)==101
    assert Q(-760*27,4)==-5130
    assert Q(10340*27-2860,4)==69080
    assert Q(3266560*27,4)==22049280
    assert Q(14*23+22,2)==172
    assert Q(-684*23,2)==-7866
    assert Q(7755*23-1815,2)==88275
    assert Q(1336320*23,2)==15367680
    assert 5130*15-69080==7870
    assert 5130*4246-22049280==-267300
    assert Q(267300*50+880*7870,267300+880*101)==Q(92230,1619)
    assert 7866*Q(25,2)-88275==10050
    assert 7866*Q(21230,11)-15367680==-186300
    assert (186300*Q(75,2)+360*10050)/(186300+360*172)==Q(117825,2758)
    bounds={12:Q(184460,1619),11:Q(117825,1379)}
    C={12:Q(40422,40475),11:Q(50254,62055)}
    fractions={12:Q(6737,21047),11:Q(25127,79982)}
    assert (bounds[12]-64)/50==C[12] and (bounds[12]-64)/156==fractions[12]
    assert (bounds[11]-49)/45==C[11] and (bounds[11]-49)/116==fractions[11]
    candidate_slack=Q(649*41706-1430*1986,1121)-21846
    assert candidate_slack==Q(-262152,1121)
    return dict(status='EXACT_POLYGON_THREE_FACE_COUPLING_COEFFICIENTS',
        polygon=dict(finite_coefficient=str(alpha),ideal_flanked_corner_coefficient=str(beta),
            minimal_cases=[list(case) for case in minimal],slacks=[str(value) for value in gaps],
            old_optimum_candidate_slack=str(candidate_slack)),
        local_three_face_bounds=local,exceptional_face_budget_coefficients=dict(nu_excess=22,facet_defect=4),
        universal_two_B_over_f_lower={str(n):str(value) for n,value in bounds.items()},
        universal_C_over_f_lower={str(n):str(value) for n,value in C.items()},
        universal_K4_corner_fraction_lower={str(n):str(value) for n,value in fractions.items()},
        scope=dict(universal_geometric_proofs_written=True,new_linear_optimizer_used=False,
            new_LP_optimum_claimed=False,dimension_exclusion=False,lean_compilation_performed=False))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write',action='store_true')
    args=parser.parse_args();record=generate()
    if args.write:CERT.write_text(json.dumps(record,indent=2)+'\n')
    else:
        with CERT.open() as stream:saved=json.load(stream)
        assert saved==record
    print('PASS polygon minima, 25 local cases, exception coefficients and stronger B/C densities; no optimizer or dimension claim')


if __name__=='__main__':main()
