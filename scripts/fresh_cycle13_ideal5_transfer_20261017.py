#!/usr/bin/env python3
"""Universal cycle13 distance-six transfer via the strict energy spectrum."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sympy as s
import fresh_finite8_cycle_replacement_20261017 as R

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_cycle13_ideal5_transfer_20261017'
A, B, C = R.A, R.B, R.C


def run():
    stems = ('fresh_cycle14_ideal5_transfer_20261017',
        'fresh_spectator_energy_below_one_20261017',
        'fresh_finite8_cycle_replacement_20261017',
        'fresh_finite8_ideal5_local_replacement_20261017')
    for stem in stems:
        d = json.loads((ROOT / f'results/{stem}.json').read_text())
        assert d['status'].startswith('PASS')
        for path, h in d['sources'].items():
            assert sha256((ROOT / path).read_bytes()).hexdigest() == h, (stem, path)
    P, Q, m = 6, 7, 13
    profiles = B.cycle_profiles()
    single, rectangles, minima = [], [], []
    for p in profiles:
        if p['name'] in ('isolated', 'both_ends3'):
            continue
        side = next(i for i, v in enumerate(p['position']) if v)
        pos = p['position'][side]
        contact = pos if side == 0 else m-pos
        w2 = 2 if p['name'].startswith('end4') else 1
        ids = list(range(1, m+1))
        edges = [(i, i+1, 1) for i in range(1, m-1)]+[(m, contact, w2)]
        adj = {i: set() for i in ids}
        for i, j, _ in edges:
            adj[i].add(j)
            adj[j].add(i)
        leaf = m-1 if side == 0 else 1
        near = next(iter(adj[leaf]))
        assert len(adj[leaf]) == 1 and len(adj[near]) == 2
        distance = C.distances(adj, ids, m)[P]
        assert distance == (P, Q)[side]-pos+1 and distance <= 7
        single.append(dict(profile=p['name'], marked_distance=distance,
            fixed_leaf=leaf, fixed_degree_two_neighbor=near,
            positive_pivots=C.tree_pd(ids, edges)))
    keep = [p for p in profiles if p['name'] in ('isolated', 'both_ends3')]
    for c, d in product(keep, repeat=2):
        values = []
        for ub, vb in product((0, 1), repeat=2):
            value = A.value(A.Q(F(1, P)), A.Q(F(1, Q)), c, d,
                c['umax'] if ub else A.Q(), d['umax'] if vb else A.Q())
            values.append(value)
            rectangles.append(dict(c=c['name'], d=d['name'],
                upper_reciprocals=[ub, vb], value=value.wire()))
        low = R.minimum(values)
        if c['name'] == d['name'] == 'isolated':
            assert low == A.Q(-F(5, 21))
        else:
            assert low.sign() > 0
        minima.append(dict(c=c['name'], d=d['name'], minimum=low.wire()))
    alpha = F(1, P)+F(1, Q)
    gate = 1-3*alpha
    energy = 4-1/alpha
    assert gate == F(1, 14) > 0 and energy == F(10, 13)
    forced_rank = energy/(1-energy)
    assert forced_rank == F(10, 3) and forced_rank.denominator != 1
    assert s.simplify(s.Rational(10, 13)-(2-s.sqrt(2))) != 0
    assert s.simplify(s.Rational(10, 13)-(1-s.sqrt(5)/5)) != 0
    # The isolated 4x4 identity is checked formally, including arbitrary Z.
    h, k, p, q, r = s.symbols('h k p q r', nonzero=True)
    aa = s.Rational(alpha.numerator, alpha.denominator)
    T = s.Matrix([[aa, -aa, -h, 0], [-aa, aa, 0, -k],
                  [-h, 0, 2-p, -q], [0, -k, -q, 2-r]])
    expected = 1-2*aa/h**2-2*aa/k**2+aa*(p/h**2-2*q/(h*k)+r/k**2)
    assert s.factor(T.det()/(h*h*k*k)-expected) == 0
    retained = [[i, j] for i in range(m) for j in range(i+1, m) if min(j-i, m-j+i) <= 5]
    removed = [[i, j] for i in range(m) for j in range(i+1, m) if min(j-i, m-j+i) == 6]
    assert len(single) == 6 and len(rectangles) == 16 and len(retained) == 65 and len(removed) == 13
    paths = [Path(__file__), Path(B.__file__), Path(A.__file__), Path(C.__file__)]
    paths += [ROOT / f'{folder}/{stem}.{ext}' for stem in stems
        for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    paths += [ROOT / 'results/fresh_cusp_factor_payment_framework_20261017.md']
    return dict(status='PASS_UNIVERSAL_CYCLE13_IDEAL5_TRANSFER',
        sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
        single_profiles=single, rectangles=rectangles, rectangle_minima=minima,
        isolated_case=dict(arcs=[P, Q], alpha=str(alpha), one_high_height_gap=str(gate),
            forced_heights=['1', '1'], required_spectator_energy=str(energy),
            forced_Ar_rank=str(forced_rank), outcome='excluded: forced rank is not an integer'),
        cycle_pair_rule=dict(retained=retained, transferred=removed, cost=65),
        counts=dict(single_profiles=6, rectangle_values=16, transferred_pairs=13),
        theorem='Every actual ideal cycle13 marked pair at distance6 transfers to finite distance<=7 or an ideal short path. No new finite coupon is required; arbitrary spectators and all native labels are covered.',
        scope='Universal local cycle13 replacement. No global face-coverage composition or ambient dimension bound is asserted.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    d = json.loads(json.dumps(run()))
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert d == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(d, indent=2)+'\n')
    print(d['status'], d['counts'], flush=True)
