#!/usr/bin/env python3
"""Exact scalar barriers to reaching the requested near-example bounds.

These concern specified proof estimates, not existence or impossibility
of stronger geometric theorems.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_dimension_gap_method_barriers_20261016'


def prokhorov_lower(n):
    if n % 2:
        a32, bad = F(6 * (n - 1), n - 3), F(3 * n + 9, 4 * (n - 2))
        expected = F(n * n - 68 * n + 147, 4 * (n - 3) * (n - 2))
    else:
        a32, bad = F(6 * (n - 2), n - 4), F(3 * n + 6, 4 * (n - 3))
        expected = F(n * n - 70 * n + 216, 4 * (n - 4) * (n - 3))
    lower = 7 - a32 - bad
    assert lower == expected
    return dict(dimension=n, mean_facets_bound=str(a32), bad_face_fraction_bound=str(bad),
                total_weight_lower_per_three_face=str(lower))


def narrow(n, loss):
    N = n - 1
    target = F(2 * N - 2, 9)
    lower = (target - loss) / 2 - 2
    assert lower == F(N - 19, 9) - F(loss, 2)
    return dict(dimension=n, affine_rank=N, target_sum=str(target), loss=str(loss),
                balanced_pairing_lower=str(lower), reaches_strict14=lower >= 14)


def run():
    current = narrow(164, 4)
    assert current['balanced_pairing_lower'] == '14'
    ideal = [narrow(n, 0) for n in range(23, 147)]
    assert all(not r['reaches_strict14'] for r in ideal[:-1])
    assert ideal[-1]['dimension'] == 146 and ideal[-1]['reaches_strict14']
    counts = [prokhorov_lower(n) for n in range(5, 69)]
    assert all(F(r['total_weight_lower_per_three_face']) < 0 for r in counts if r['dimension'] <= 66)
    assert all(F(r['total_weight_lower_per_three_face']) > 0 for r in counts if r['dimension'] >= 67)
    paths = [f'scripts/{TAG}.py', 'results/fresh_nikulin_two_budget_frontier_20261015_dimension163.md',
             'results/fresh_nikulin_two_budget_frontier_20261015.md',
             'results/fresh_nikulin_rank162_local_obstruction_20261016.md',
             'results/dimension_bound_prokhorov_20261015.md',
             'results/fresh_sqrt5_discount_structural_floor_20261016.json']
    return dict(status='PASS_EXACT_SCALAR_BARRIERS_FOR_SPECIFIED_BOUND_METHODS',
                sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
                desired_unrestricted_targets=list(range(23, 27)),
                current_narrow_threshold=current, hypothetical_zero_loss=ideal,
                unchanged_prokhorov_lower=counts,
                scope='Exact scalar implications of the stated frozen proof inputs only. Even hypothetical zero loss cannot make the fixed-target balanced narrow-facet scalar certificate exclude dimensions below146. The unchanged Prokhorov bad-face/mean-facet lower estimate is negative through66, so nonnegative edge weights alone cannot give an exclusion there through those two scalar bounds. These do not assert geometric witnesses, optimal universal bounds, or limitations on changed narrowness constants, targets, incidence information, signed weights, or new arguments.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    a = ap.parse_args()
    data = run()
    path = ROOT / f'results/{TAG}.json'
    if a.verify:
        assert data == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(data, indent=2) + '\n')
    print(data['status'])
