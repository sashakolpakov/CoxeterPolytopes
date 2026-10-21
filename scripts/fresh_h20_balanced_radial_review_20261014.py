#!/usr/bin/env python3
"""Independent Arb replay of all16,072 radial power bounds at2z/5.

Uses the source's sharper radius bounds; it does not independently prove
their all-label coverage or the secondary half-probe incompatibility.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from flint import fmpq as Q, fmpq_mat as M, arb, ctx

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_balanced_radial_review_20261014'
SRC = 'results/fresh_h20_balanced_half_20261014.json'
CENSUS = 'results/fresh_h20_balanced_endpoint_20261014.json'
ctx.prec = 192


def value(row):
    radicals = (1, 2, 3, 6, 5, 10, 15, 30)
    return sum((arb(Q(c))*arb(d).sqrt() for c, d in zip(row, radicals)), arb(0))


def run():
    source, census = [json.loads((ROOT/p).read_text()) for p in (SRC, CENSUS)]
    for doc in (source, census):
        for p, h in doc['sources'].items():
            assert sha256((ROOT/p).read_bytes()).hexdigest() == h, p
    H = M(19, 19)
    for lo, n in ((0, 6), (6, 6), (12, 7)):
        for i in range(n):
            H[lo+i, lo+i] = 1
            if i+1 < n:
                H[lo+i, lo+i+1] = H[lo+i+1, lo+i] = Q(-1, 2)
    inverse = H.inv()
    profiles = census['types'][0]['profiles']
    assert len(profiles) == len(source['radius_bounds']) == 41
    radii = [value(r['radius_upper_bound']) for r in source['radius_bounds']]
    costs = [value(r['normalized_spatial_cost']) for r in profiles]
    cosines = [[value(c) for c in r['cosine_vector']] for r in profiles]
    trials = 0
    for ti, tr in enumerate(census['literal_edge_transports']):
        gaps = [Q(3, 5)+Q(2, 5)*(7 if w in (0, 7) else 8 if w == 14 else 0)
                for w in tr['edge_support_in_profile_order']]
        dual = [sum(gaps[i]*inverse[i, j] for i in range(19)) for j in range(19)]
        dist = sum(gaps[i]*dual[i] for i in range(19))
        assert all(a > 0 for a in dual)
        for pi, (p, radius, cost, cs) in enumerate(zip(profiles, radii, costs, cosines)):
            assert p['id'] == pi
            linear = sum((arb(d)*c for d, c in zip(dual, cs)), arb(0))
            power = arb(dist)-2*radius*linear+(cost-1)*radius*radius
            if ti == pi == 0:
                assert p['contacts'] == []
                assert source['radius_bounds'][0]['radius_upper_bound'] == ['8']+['0']*7
                assert dist-64 == Q(164, 25)
            else:
                assert power > arb(Q(164, 25)), (ti, pi, str(power))
            trials += 1
    assert trials == 16072
    # Every interior-center central-contact bridge ball is contained in |x|<=5.
    # At |x|=28/5, its power is at least distance-to-ball squared >=9/25.
    assert (Q(28, 5)-5)**2 == Q(9, 25) > Q(1, 4)
    D = (Q(784, 25)-Q(1, 4))/2
    assert Q(784, 25)-2*D == Q(1, 4)
    paths = [SRC, CENSUS, str(Path(__file__).resolve().relative_to(ROOT))]
    return dict(status='PASS_INDEPENDENT_16072_ENDPOINT_RADIAL_POWER_BOUNDS',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                trials=trials, minimum_power_bound='164/25', unique_bound_attainer=[0, 0],
                bridge_power_lower_bound='9/25', common_escape_norm='1/4',
                scope='Separate unit-Gram inverse and192-bit Arb evaluations for all41x392 powers. Assumes sharper all-label radius bounds and radius-five bridge theorem established in the source reports. Does not independently check their classification proofs or the secondary two-half-probe theorem.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    out = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert json.loads(path.read_text()) == out
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], flush=True)
