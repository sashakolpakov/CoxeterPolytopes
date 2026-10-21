#!/usr/bin/env python3
"""Independent endpoint census and radius-bound replay.

Reconstructs the cycle Gram matrices, enumerates every reduced contact vector,
and recomputes bounds in Q(sqrt(2),sqrt(5)) with FLINT rational matrices.
No arithmetic or geometry functions are imported from the producers.
"""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
from flint import fmpq as Q, fmpq_mat, arb, ctx

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_balanced_independent_20261014'
BASE = 'fresh_h20_balanced_endpoint_20261014'
ctx.prec = 192


def q(x=0):
    return (Q(x), Q(0), Q(0), Q(0))


ZERO, ONE = q(), q(1)
SQRT2, SQRT5 = (Q(0), Q(1), Q(0), Q(0)), (Q(0), Q(0), Q(1), Q(0))


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def scale(a, b):
    return tuple(x*b for x in a)


def sub(a, b):
    return add(a, scale(b, -1))


def mul(a, b):
    out = [Q(0)]*4
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i ^ j] += x*y*(2 if i & j & 1 else 1)*(5 if i & j & 2 else 1)
    return tuple(out)


def div(a, b):
    cols = [mul(b, tuple(Q(i == j) for i in range(4))) for j in range(4)]
    matrix = fmpq_mat([[cols[j][i] for j in range(4)] for i in range(4)])
    sol = matrix.solve(fmpq_mat([[x] for x in a]))
    out = tuple(sol[i, 0] for i in range(4))
    assert mul(out, b) == a
    return out


def sign(a):
    if a == ZERO:
        return 0
    value = arb(a[0])+arb(a[1])*arb(2).sqrt()+arb(a[2])*arb(5).sqrt()+arb(a[3])*arb(10).sqrt()
    assert not value.contains(0), str(value)
    return 1 if value > 0 else -1


def read_field(a):
    assert len(a) == 8 and all(Q(a[i]) == 0 for i in (2, 3, 6, 7))
    return tuple(Q(a[i]) for i in (0, 1, 4, 5))


def serial(a):
    return list(map(str, a))


def run():
    paths = [f'results/{BASE}{s}.json' for s in ('', '_bounds', '_central_gap')]
    census, bounds, gap = [json.loads((ROOT/p).read_text()) for p in paths]
    for report in (census, bounds, gap):
        for p, h in report['sources'].items():
            assert sha256((ROOT/p).read_bytes()).hexdigest() == h, p
    factors = [list(range(7)), list(range(7, 14)), list(range(14, 22))]
    G = [[0]*22 for _ in range(22)]
    for factor in factors:
        for k, i in enumerate(factor):
            G[i][i] = 2
            G[i][factor[(k+1) % len(factor)]] = -1
            G[i][factor[(k-1) % len(factor)]] = -1
    assert census['normalized_seed_Gram'] == G
    assert census['factors'] == factors
    assert census['positive_affine_marks'] == [1]*22
    typ = census['types'][0]
    ids = [i for i in range(22) if i not in (0, 7, 14)]
    assert typ['edge_support'] == ids
    H = fmpq_mat([[G[i][j] for j in ids] for i in ids])
    inv = H.inv()
    assert H == fmpq_mat(typ['PD19_Gram'])
    assert inv == fmpq_mat(typ['inverse_Gram'])
    assert all(inv[i, j] >= 0 for i in range(19) for j in range(19))
    assert min(inv[i, i] for i in range(19)) == Q(6, 7)
    # Normalized cost uses the inverse unit Gram, which equals2*inv.
    # Three nonright contacts cost at least3*(1/2)^2*(12/7)=9/7>1.
    assert 3*Q(1, 4)*Q(12, 7) > 1
    # Any label>=5 costs at least(12/7)*cos(pi/5)^2>1.
    assert sign(sub(scale(add(q(3), SQRT5), Q(3, 14)), ONE)) > 0
    c3, c4 = q(Q(1, 2)), scale(SQRT2, Q(1, 2))
    catalog = {}
    trials = 0
    for size in range(3):
        for places in combinations(range(19), size):
            for labels in product((3, 4), repeat=size):
                trials += 1
                c = [ZERO]*19
                for i, label in zip(places, labels):
                    c[i] = c3 if label == 3 else c4
                cost = ZERO
                for i in places:
                    for j in places:
                        cost = add(cost, scale(mul(c[i], c[j]), 2*inv[i, j]))
                if sign(sub(cost, ONE)) > 0:
                    continue
                key = tuple(zip((ids[i] for i in places), labels))
                catalog[key] = (c, cost)
    assert trials == 723 and len(catalog) == 41
    saved = {tuple(map(tuple, p['contacts'])): p for p in typ['profiles']}
    assert set(saved) == set(catalog)
    counts = Counter()
    for key, (c, cost) in catalog.items():
        p = saved[key]
        assert [read_field(x) for x in p['cosine_vector']] == c
        assert read_field(p['normalized_spatial_cost']) == cost
        endpoint = 'ideal_PSD' if cost == ONE else 'finite'
        assert p['endpoint'] == endpoint
        counts[endpoint] += 1
    seen = set()
    for tr in census['literal_edge_transports']:
        perm = tr['wall_permutation']
        assert sorted(perm) == list(range(22))
        assert all(G[i][j] == G[perm[i]][perm[j]] for i in range(22) for j in range(22))
        assert tr['omitted'] == [perm[i] for i in (0, 7, 14)]
        assert tr['edge_support_in_profile_order'] == [perm[i] for i in ids]
        seen.add(tuple(tr['omitted']))
    assert seen == set(product(*factors))
    phi = scale(add(ONE, SQRT5), Q(1, 2))
    minima, exceptional_count = [], 0
    for mode, third in zip(bounds['scales'], (ONE, phi)):
        scales = [ONE, ONE, third]
        A = [q(7), q(7), scale(third, 8)]
        N = add(q(112), scale(mul(third, third), 84))
        assert read_field(mode['all_vertex_square']) == N
        low = None
        for p, row in zip(typ['profiles'], mode['profiles']):
            assert p['id'] == row['profile_id']
            c, cost = catalog[tuple(map(tuple, p['contacts']))]
            sums = [ZERO]*3
            for k, factor in enumerate(factors):
                for i, wall in enumerate(ids):
                    if wall in factor:
                        sums[k] = add(sums[k], c[i])
            assert [read_field(r['subtract']) for r in p['omitted_contact_ratios']] == sums
            radii = [div(a, add(b, c4)) for a, b in zip(A, sums)]
            base = radii[0]
            for r in radii[1:]:
                if sign(sub(r, base)) < 0:
                    base = r
            assert read_field(row['continuous_radius_bound']) == base
            candidates = [base]
            exceptions = []
            for f, (a, b) in enumerate(zip(A, sums)):
                for target in (ZERO, c3):
                    den = add(b, target)
                    if den == ZERO:
                        continue
                    r = div(a, den)
                    ratios = [sub(div(aa, r), bb) for aa, bb in zip(A, sums)]
                    good = all(x in (ZERO, c3) or sign(sub(x, c4)) >= 0 for x in ratios)
                    exceptions.append((f, target, r, ratios, good))
                    if good:
                        candidates.append(r)
            assert len(exceptions) == len(row['exceptional_radii'])
            exceptional_count += len(exceptions)
            for (f, target, r, ratios, good), saved_ex in zip(exceptions, row['exceptional_radii']):
                assert saved_ex['factor'] == f and saved_ex['retained'] == good
                assert read_field(saved_ex['target_cosine']) == target
                assert read_field(saved_ex['radius']) == r
                assert list(map(read_field, saved_ex['ratios'])) == ratios
            maximum = candidates[0]
            for r in candidates[1:]:
                if sign(sub(r, maximum)) > 0:
                    maximum = r
            assert read_field(row['all_radius_upper_bound']) == maximum
            zy = ZERO
            for i, wall in enumerate(ids):
                s = scales[0 if wall < 7 else 1 if wall < 14 else 2]
                for j in range(19):
                    zy = add(zy, scale(mul(s, c[j]), 2*inv[i, j]))
            assert read_field(row['vertex_dot_direction']) == zy
            P = add(sub(N, scale(mul(maximum, zy), 2)), mul(sub(cost, ONE), mul(maximum, maximum)))
            assert read_field(row['incenter_ball_power_lower_bound']) == P
            assert sign(sub(q(14), maximum)) >= 0
            assert sign(sub(P, q(68 if third == ONE else 75))) >= (0 if third == ONE else 1)
            if low is None or sign(sub(P, low)) < 0:
                low = P
        minima.append(serial(low))
    # Independently replay the strict Cauchy margins and all22 seed gaps.
    assert Q(65, 4)**2 > 264 and Q(33, 2)**2 > 271
    assert 68+Q(441, 100)-Q(21, 5)*Q(65, 4) == Q(104, 25)
    assert 75+5-2*Q(23, 10)*Q(33, 2) == Q(41, 10)
    golden_x2 = scale(add(q(49), scale(phi, 21)), Q(1, 16))
    assert sign(sub(golden_x2, q(5))) > 0
    assert sign(sub(q(Q(23, 10)**2), golden_x2)) > 0
    for t, scales in ((Q(3, 20), [ONE]*3), (Q(1, 8), [ONE, ONE, phi])):
        for factor, s in zip(factors, scales):
            assert sign(scale(s, 1-t)) > 0
            assert sign(scale(s, 1+(len(factor)-1)*t)) > 0
    paths += [str(Path(__file__).resolve().relative_to(ROOT))]
    return dict(status='PASS_INDEPENDENT_BALANCED_ENDPOINT_CENSUS_BOUNDS_AND_GAP',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                trials=trials, profiles=len(catalog), endpoint_counts=dict(counts),
                transports=len(seen), radius_bounds=82, exceptional_radii=exceptional_count,
                field_basis=['1', 'sqrt2', 'sqrt5', 'sqrt10'], actual_minimum_power_bounds=minima,
                unit_escape_norm='1/4', golden_escape_norm='1/4',
                scope='Separate rational-matrix and four-dimensional field implementation; complete723 reduced trials and392 transports; all82 radius bounds and every exceptional omitted right/pi3 case. Interval sign decisions use192-bit Arb with exact zero handled algebraically. No finite-volume construction is asserted.')


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
    print(out['status'], out['endpoint_counts'], flush=True)
