#!/usr/bin/env python3
"""Compose the adaptive depth-six finite-eight corridor with cycle payments.

Exact scalar checks and frozen proof bindings. The report gives the universal
supporting-white and current-rule composition arguments.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import fresh_depth6_finite8_adaptive_pilot_20261018 as H
import fresh_depth8_finite8_zero_transfer_20261018 as D
import fresh_finite8_cycle_replacement_20261017 as C

ROOT = H.ROOT
TAG = 'fresh_depth6_finite8_zero_transfer_20261018'
STEMS = (H.TAG, D.TAG, C.TAG, 'fresh_distance8_ideal_corner_replacement_20261017',
         'fresh_cycle16_half_coupon_transfer_20261017',
         'fresh_depth6_finite7_adaptive_transfer_20261018',
         'fresh_depth11_finite7_zero_transfer_20261018',
         'fresh_depth6_finite6_zero_transfer_20261018',
         'fresh_a15_d15_dimension58_20261018',
         'fresh_a15_d15_ambient59_forest_closure_20261018')


def cycle_boundary(p, q):
    assert (p, q) in ((7, 8), (8, 8))
    classes = []
    for name in D.PROFILES:
        if name in ('isolated', 'both_ends3'):
            kind, distance = 'unrestricted', None
        else:
            side = int(name[-1])
            distance = (p, q)[side]-(2 if name.startswith('next') else 1)+1
            if distance == 8 and name.startswith('end3'):
                kind = 'changed_deep8_requires_neighbor'
            else:
                assert distance <= 8
                kind = 'retained_short_or_feature_unit'
        classes.append(dict(name=name, classification=kind, distance=distance))
    keep = [r for r in classes if r['classification'] != 'retained_short_or_feature_unit']
    expected = (['isolated', 'end3_1', 'both_ends3'] if p == 7 else
                ['isolated', 'end3_0', 'end3_1', 'both_ends3'])
    assert [r['name'] for r in keep] == expected
    a, b = F(1, p), F(1, q)
    rows, pairs, strict = [], [], []
    for cr, dr in product(keep, repeat=2):
        c, d = cr['name'], dr['name']
        side = ('c' if cr['classification'].endswith('requires_neighbor') else
                'd' if dr['classification'].endswith('requires_neighbor') else None)
        constants = [D.PROFILES[c][2], D.PROFILES[d][2]]
        if side:
            constants[0 if side == 'c' else 1] -= F(1, 2)
        assert min(constants) >= 0
        values = []
        for u, v in product((0, 1), repeat=2):
            l = u*D.PROFILES[c][0]+v*D.PROFILES[d][0]
            m = u*D.PROFILES[c][1]+v*D.PROFILES[d][1]
            norm = u*D.PROFILES[c][2]+v*D.PROFILES[d][2]
            raw = 1-2*(a*l+b*m)-(a+b)*norm-a*b*(l-m)**2
            correction = (a+b)*(F(u, 2) if side == 'c' else F(v, 2) if side == 'd' else F())
            corrected = 1-2*(a*l+b*m)-(a+b)*(u*u*constants[0]+v*v*constants[1])-a*b*(l-m)**2
            assert corrected == raw+correction
            if (c, d) != ('isolated', 'isolated'):
                assert corrected > 0
                strict.append(corrected)
            values.append(corrected)
            rows.append(dict(c=c, d=d, upper_reciprocals=[u, v], fixed_cauchy_side=side,
                             raw=str(raw), correction=str(correction), lower=str(corrected)))
        pairs.append(dict(c=c, d=d, fixed_cauchy_side=side,
                          corrected_constants=list(map(str, constants)), minimum=str(min(values))))
    assert min(strict) == {(7, 8): F(3, 112), (8, 8): F(1, 16)}[p, q]
    return dict(P=p, Q=q, classifications=classes, retained_profiles=expected,
                ordered_pairs=pairs, rectangle_values=rows,
                strict_pair_count=len(pairs)-1, strict_minimum=str(min(strict)),
                isolated_case_inherited=True)


def isolated():
    rows = []
    inherited = {(r['P'], r['Q']): r for r in C.isolated_cases()}
    for p, q in ((7, 8), (8, 8)):
        alpha = F(1, p)+F(1, q)
        gate, energy = 1-3*alpha, 4-1/alpha
        assert gate > 0 and 0 <= energy < F(1, 2)
        assert inherited[p, q]['required_unscaled_spectator_difference'] == str(energy)
        rows.append(dict(P=p, Q=q, one_height_at_least_sqrt2_gap=str(gate),
                         required_spectator_energy=str(energy),
                         outcome='two_retained_whole_A1_A15_middle_halves' if energy == 0 else
                                 'impossible_by_disjoint_support_energy_gap'))
    assert [r['required_spectator_energy'] for r in rows] == ['4/15', '0']
    return rows


def run():
    docs = {s: json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem, doc in docs.items():
        assert doc['status'].startswith('PASS'), stem
        for path, digest in doc['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == digest, (stem, path)
    pilot = docs[H.TAG]
    assert pilot['universal_payment_theorem'] is False
    assert pilot['profiles']['candidate_count'] == 99 and pilot['profiles']['accepted_count'] == 13
    schur = pilot['schur']
    assert schur['branch_count'] == schur['positive_branches'] == 576 and not schur['bad_branches']
    assert schur['minimum_positive_constant'] == dict(rational='3/44', sqrt2='0', sign=1)
    assert all(r['q_minus2_lower']['sign'] > 0 and r['r_minus2_lower']['sign'] > 0 and
               r['constant']['sign'] > 0 and all(c['value']['sign'] >= 0 for c in r['coefficients'])
               for r in schur['branches'])
    count = H.counts()
    assert H.enc(count) == pilot['counts']
    normalized = D.profile_check(docs)
    boundaries = [cycle_boundary(*pq) for pq in ((7, 8), (8, 8))]
    tails = [D.cycle_region(q) for q in (9, 10, 11)]
    for new, old in zip(tails, docs[D.TAG]['cycle_regions']):
        assert new == old
    assert sum(len(r['ordered_pairs']) for r in boundaries+tails) == 102
    assert sum(len(r['rectangle_values']) for r in boundaries+tails) == 408
    current = docs['fresh_a15_d15_dimension58_20261018']
    assert current['dimension_bound'] == 58
    paths = {Path(__file__), Path(H.__file__), Path(D.__file__), Path(C.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS
                 for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md')))
    return dict(status='PASS_DEPTH6_FINITE8_ZERO_LOCAL_TRANSFER',
                sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
                current_baseline='fresh_a15_d15_dimension58_20261018',
                feature_depth=6, finite_distance=8, previous_feature_depth=8, new_weight='0',
                inherited_profiles=normalized, cycle_boundaries=boundaries, cycle_tails=tails,
                cycle_ordered_pairs=102, cycle_rectangle_values=408,
                strict_cycle_pairs=100, isolated_cases=isolated(),
                retained_arc6_8_fixed_leaf_distance=5, counts=count,
                triangle_minimum=2, quadrilateral_minimum=1, new_coupons=False,
                all_ideal_weights_unchanged=True, other_finite_weights_unchanged=True,
                new_forest_count=False, new_dimension_bound=None,
                theorem='In the completed dimension58 rule, delete every remaining whole finite A/B/D distance8 unit whose marks both have feature depth at least6. Triangles retain total2 and quadrilaterals total1; no new coupon is needed.',
                scope='Written universal supporting-white/cycle/current-rule composition with exact scalar certificates. Inherited profile classification, actual-star geometry and historical local theorems are not regenerated. No new forest maximum or dimension theorem.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        assert not path.exists()
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], 'new cycle minima', [r['strict_minimum'] for r in result['cycle_boundaries']], flush=True)
