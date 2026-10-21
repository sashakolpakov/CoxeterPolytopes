#!/usr/bin/env python3
"""Exact inner gap for every central-compatible radius-two integer-gap cap."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from fractions import Fraction as F
from fresh_h20_balanced_integer_obstruction_20261014 import factor

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_integer_inner_gap_20261014'


def run():
    ff=[factor(n)for n in(6,6,7)]
    # At x=3z/20, distance squared has denominator11200.
    # Enumerate every cap with power <=1/100, including the boundary.
    threshold=44912;parts=[]
    for rows,_,V in ff:
        live=[]
        for r in rows:
            distance=400*r['norm28']-120*r['dot_vertex28']+9*V
            if distance<=threshold:live.append((r,distance))
        parts.append(live)
    records=[];allowed=[]
    for choice in product(*parts):
        dist=sum(a[1]for a in choice)
        if dist>threshold:continue
        rows=[a[0]for a in choice];N=sum(r['norm28']for r in rows)
        central=N in(224,336)or N>=448
        row=dict(gaps=[list(r['gaps'])for r in rows],norm28=N,
                 dot_vertex28=sum(r['dot_vertex28']for r in rows),
                 power_numerator=dist-44800,central_compatible=central)
        records.append(row)
        if central:allowed.append(row)
    assert allowed and all(r['power_numerator']==112 for r in allowed)
    # Every omitted recipe has power>1/100 by construction.
    x2=F(441,100);square=F(1,200);D=(x2-square)/2
    assert D==F(881,400)
    central_power=x2-4
    endpoint_power=196*(1-F(3,20))**2-11**2
    assert central_power>square and endpoint_power>square
    files=[f'scripts/{TAG}.py','scripts/fresh_h20_balanced_integer_obstruction_20261014.py']
    return dict(status='PASS_ALL_RADIUS_TWO_INTEGER_GAP_CAPS_MISS_INNER_PROBE',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in files},
                original_factor_counts=[len(f[0])for f in ff],
                bounded_factor_counts=list(map(len,parts)),
                small_power_recipes=len(records),central_compatible_small_power_recipes=len(allowed),
                exact_power_lower_bound='1/100',power_denominator=11200,records=records,
                escape=dict(spatial_multiple='3/20',vertex=0,D=str(D),U=1,norm=str(square),
                            central_power=str(central_power),minimum_seed_gap='17/20',
                            relaxed_endpoint_radius_bound=11,
                            relaxed_endpoint_power_lower_bound=str(endpoint_power)),
                scope='All nonnegative integerseedgaps summing7,7,8, radius2, compatible withcentralradius2. The enumeration is complete for power<=1/100 by nonnegative factor distances. Every such cap leaves the same norm1/200 strictescape, irrespective of endpoint walls or mutualcompatibility. Even all vertex-centered caps ofradius<=11 leave it. Other capradii, nonintegerseedgaps, andnonvertexendpoints are outside this result.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run()
    p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['small_power_recipes'],d['central_compatible_small_power_recipes'])
