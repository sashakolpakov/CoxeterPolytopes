#!/usr/bin/env python3
"""Exact parity identities for normalized product-link curvature."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_cusp_link_rank_monotonicity_20261017'
INPUTS=['scripts/fresh_cusp_link_shape_curvature_20261017.py',
 'results/fresh_cusp_link_shape_curvature_20261017.json',
 'results/fresh_cusp_link_shape_curvature_20261017.md']

def fr(r,y):
    return sum((abs(y+Q(j)-Q(r,2))for j in range(r+1)),Q())/Q(r+1)

def run():
    k,j,t=s.symbols('k j t',integer=True,nonnegative=True)
    even_den=(2*k+1)*(2*k+2);odd_den=(2*k+2)*(2*k+3)
    even_int=s.cancel(((k+1)**2+j*j)/(2*k+2)-(k*(k+1)+j*j)/(2*k+1))
    even_half=s.cancel(((k+1)**2+j*(j+1))/(2*k+2)-(k*(k+1)+j*(j+1)+s.Rational(1,2))/(2*k+1))
    odd_int=s.cancel(((k+1)*(k+2)+j*j)/(2*k+3)-((k+1)**2+j*j)/(2*k+2))
    odd_half=s.cancel(((k+1)*(k+2)+j*(j+1)+s.Rational(1,2))/(2*k+3)-((k+1)**2+j*(j+1))/(2*k+2))
    expected=[((k+1)**2-j*j)/even_den,
              (k*(k+1)-j*(j+1))/even_den,
              ((k+1)**2-j*j)/odd_den,
              ((k+1)*(k+2)-j*(j+1))/odd_den]
    for a,b in zip((even_int,even_half,odd_int,odd_half),expected):assert s.cancel(a-b)==0
    checks=[]
    for r in range(0,61):
        cap=Q(r+2,4*(r+1))if r%2==0 else Q(r+3,4*(r+2))
        vals=[fr(r+1,Q(z,2))-fr(r,Q(z,2))for z in range(r+4)]
        assert max(vals)==cap and min(vals)==0
        checks.append(dict(rank=r,maximum=str(cap),breakpoint_checks=len(vals)))
    assert Q(3,4)-Q(5,2)*Q(5,19)==Q(7,76)
    assert s.expand(4563*(t+2)-8075)==1051+4563*t
    even=s.expand((30+2*t)*(30+2*t-6)-24*(30+2*t))
    odd=s.expand((31+2*t-1)*(31+2*t-5)-24*(31+2*t))
    assert even==4*t*t+60*t and odd==4*t*t+64*t+36
    paths=[Path(__file__)]+[ROOT/p for p in INPUTS]
    return dict(status='PASS_ALL_RANK_CUSP_CURVATURE_MONOTONICITY',
      sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
      parity_breakpoint_formulas=list(map(str,expected)),literal_rank_checks=checks,
      pointwise_increment_cap_rank_at_least17='5/19',
      normalized_curvature_increment_lower_bound='7/76',
      multiple_factor_CS_squared_margin='4563*s-8075 >= 1051 for s>=2',
      single_factor_slack_numerator_even=str(even),single_factor_slack_numerator_odd=str(odd),
      scope='All-rank normalized product-link curvature dominates3s when all factor ranks>=17 and ambient dimension>=30. A dimension29 implication is conditional on complete face payments; the finite-quadrilateral transfer is not supplied here.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['literal_rank_checks']),'corroborating ranks',flush=True)
