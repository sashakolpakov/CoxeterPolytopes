#!/usr/bin/env python3
"""Exact reciprocal-square identity and its opposition triangle cap."""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CERT=ROOT/'results/right_angled_mixed_opposition_power_two_research_20261018.json'


def generate():
    mu=Q(5,4)
    add=lambda a,b:[(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))]
    def mul(a,b):
        result=[Q(0)]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):result[i+j]+=x*y
        return result
    r=[mu,Q(1)];r2=mul(r,r)
    terms=[[mu**3], [0,-2*mu**2], [0,0,3*mu], [0,0,0,-4]]
    polynomial=[Q(0)]
    for term in terms:polynomial=add(polynomial,mul(term,r2))
    polynomial=add(polynomial,[0,0,0,0,5*mu,4])
    assert polynomial==[mu**5,0,0,0,0,0]
    bound=1/mu**2+3*Q(9,176)/mu**4+4*Q(27,1760)/mu**5+13*Q(3,440)/mu**5
    assert bound==Q(129248,171875) and 220*bound==Q(516992,3125)
    return dict(status='EXACT_RECIPROCAL_SQUARE_COEFFICIENT_CERTIFICATE',
        mu=str(mu),identity_after_multiplying_mu5_r2=[str(x) for x in polynomial],
        remainder_multiplier_maximum='13',mean_reciprocal_square_upper=str(bound),
        universal_two_B_over_f_upper=str(220*bound),
        scope=dict(written_opposition_moment_proof=True,dimension_exclusion=False,
            optimizer_used=False,lean_compilation_performed=False))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write',action='store_true')
    args=parser.parse_args();record=generate()
    if args.write:CERT.write_text(json.dumps(record,indent=2)+'\n')
    else:
        with CERT.open() as stream:saved=json.load(stream)
        assert record==saved
    print('PASS exact reciprocal-square identity: universal 2B/f <=516992/3125; no dimension claim')


if __name__=='__main__':main()
