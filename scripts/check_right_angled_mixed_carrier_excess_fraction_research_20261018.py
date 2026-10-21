#!/usr/bin/env python3
"""Exact rational interval proof of a dimension-twelve carrier excess fraction."""
import argparse
import fcntl
from fractions import Fraction as Q
from math import comb,isqrt
import json
from pathlib import Path
import run_lean_batches_20261018 as R

ROOT=Path(__file__).resolve().parents[1]
CERT=ROOT/'results/right_angled_mixed_carrier_excess_fraction_research_20261018.json'
SQRT_DEN=10**7
B=Q(1414213562373095,10**15)  # strictly below sqrt(2)


def sqrt_bounds(value):
    value=Q(value);integer=isqrt(value.numerator*SQRT_DEN**2//value.denominator)
    lo=Q(integer,SQRT_DEN);hi=lo if lo*lo==value else Q(integer+1,SQRT_DEN)
    assert lo*lo<=value<=hi*hi and hi-lo<=Q(1,SQRT_DEN)
    return lo,hi


def row_upper(k,xlo,xhi):
    m=12-k;finite=Q(0);ideal=Q(0)
    if m:
        amin=B-sqrt_bounds(xhi)[1];q=xhi/m
        low=sqrt_bounds(xlo)[0]/m;high=sqrt_bounds(q)[1]
        critical=(amin*B+2*q)/(amin+2*B)
        y=max(low,min(high,critical));mu=B-y
        assert amin>0 and mu>0
        finite=m*(1/mu+(q-y*y)/(amin*mu*mu))
    if k>=3:
        ilo,ihi=1-xhi,1-xlo;c=Q(3*(k-3),k*(k-1))
        amin=2-sqrt_bounds(3*ihi)[1]
        low=sqrt_bounds(ilo)[0]/k;high=sqrt_bounds(ihi/k)[1]
        critical=(6*amin+6*c*ihi)/(9*amin+4*k*c)
        y=max(low,min(high,critical));mu=2-3*y
        assert amin>0 and mu>0
        ideal=Q(comb(k,3),20)*(1/mu+c*(ihi-k*y*y)/(amin*mu*mu))
    return finite+ideal


def generate():
    assert B>0 and B*B<2
    records=[];all_gaps=[];boxes=0
    for k in range(13):
        if k==0:intervals=[(Q(1),Q(1))]
        elif k==12:intervals=[(Q(0),Q(0))]
        else:intervals=[(Q(j,16),Q(j+1,16)) for j in range(16)]
        degree=Q(12-k)+Q(comb(k,3),20)
        gaps=[degree-row_upper(k,lo,hi) for lo,hi in intervals]
        assert all(gap>=Q(1,8) for gap in gaps)
        worst=min(range(len(gaps)),key=lambda i:gaps[i])
        records.append(dict(k=k,boxes=len(intervals),degree=str(degree),
            minimum_gap=str(gaps[worst]),worst_interval=[str(x) for x in intervals[worst]]))
        all_gaps+=gaps;boxes+=len(gaps)
    assert boxes==178
    assert min(all_gaps)>Q(125535,10**6)
    # Crude exceptional row bound: 12sqrt2+11(1+sqrt3)<48.
    assert Q(10,7)**2>2 and Q(7,4)**2>3
    assert 12*Q(10,7)+11*(1+Q(7,4))<48
    epsilon=Q(1,8);bad_fraction=epsilon/(48+epsilon)
    assert bad_fraction==Q(1,385) and bad_fraction/14==Q(1,5390)
    # p<=14(nu-10): nu11 has F<=10; nu>=12 has p<=2nu-4.
    assert 4*10-2*11-4==14 and 14*(12-10)-(2*12-4)>0
    return dict(status='EXACT_RATIONAL_INTERVAL_CARRIER_EXCESS_FRACTION',dimension=12,
        square_root_denominator=SQRT_DEN,sqrt_two_lower=str(B),intervals=boxes,
        guaranteed_small_carrier_row_gap=str(epsilon),global_minimum_certified_gap=str(min(all_gaps)),
        per_k=records,exceptional_row_excess_upper='48',
        exceptional_finite_vertex_fraction_lower=str(bad_fraction),
        universal_nu_three_excess_over_f_lower='1/5390',
        scope=dict(written_geometric_and_calculus_proof=True,all_intervals_rationally_enclosed=True,
            sample_based_proof=False,dimension_eleven_claimed=False,dimension_ceiling_improved=False,
            optimizer_used=False,lean_compilation_performed=False))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    if args.write:
        lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a')
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);R.resources()
    record=generate()
    if args.write:CERT.write_text(json.dumps(record,indent=2)+'\n')
    else:
        with CERT.open() as stream:saved=json.load(stream)
        assert saved==record
    print('PASS 178 exact rational interval boxes: n12 carrier excess E>=f/5390; no optimizer or dimension exclusion')


if __name__=='__main__':main()
