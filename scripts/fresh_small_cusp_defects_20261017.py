#!/usr/bin/env python3
"""Exact small-cusp balance and a nonincident affine-carrier example.

No optimizer or imported project arithmetic. The face-localization lemma is
proved in the companion report, using explicitly bound frozen proofs.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from math import comb, prod
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_small_cusp_defects_20261017'
INPUTS = [
    'results/fresh_distance14_unrestricted_dimension156_20261017.md',
    'results/fresh_distance14_unrestricted_dimension156_20261017.json',
    'results/fresh_distance14_unrestricted_dimension156_20261017_review.md',
    'results/fresh_cusp_minimum_factor_rank9_dimension29_20261017.md',
    'results/fresh_cusp_minimum_factor_rank9_dimension29_20261017.json',
    'results/fresh_affine_order10_quad_payment_20261017.md',
    'results/fresh_affine_carriers_order10_20261017.md',
    'results/fresh_vinberg_finite_quad_affine_cutoff9_20261017.md',
    'results/fresh_cusp_link_shape_curvature_20261017.md',
    'results/fresh_cusp_factor_payment_framework_20261017.md',
    'results/fresh_cusp_spectator_role_payments_20261017.md',
    'results/fresh_compact_dimension_bound_20261015.md',
]


def digest(p):
    h = sha256()
    with p.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


@lru_cache(None)
def polynomial(ranks):
    out = [1]
    for r in ranks:
        new = [0] * (len(out) + r)
        for i, a in enumerate(out):
            for j in range(r + 1):
                new[i+j] += a
        out = new
    return tuple(out)


def moment_independently(ranks):
    # Inclusion-exclusion for bounded compositions, independent of convolution.
    d, k = sum(ranks), len(ranks)
    cut, numerator = (d-1)//2, 0
    for j in range(k+1):
        for subset in combinations(ranks, j):
            shift = sum(r+1 for r in subset)
            m = cut-shift
            if m >= 0:
                numerator += (-1)**j * ((d-2*shift)*comb(m+k, k)
                                        - 2*k*comb(m+k, k+1))
    return Q(numerator, prod(r+1 for r in ranks))


def concentrated(d):
    q, a = divmod(d, 28)
    out = [(28,)*q + ((a,) if a else ())]
    for q in range(d//28+1):
        for a in range(28):
            b = d-28*q-a
            if b >= 29:
                out.append((28,)*q + ((a,) if a else ()) + (b,))
    assert len(out) == len({tuple(sorted(x)) for x in out})
    return out


def profile(ranks):
    ranks = tuple(sorted(ranks))
    d = sum(ranks)
    h = polynomial(ranks)
    v = prod(r+1 for r in ranks)
    assert h == h[::-1] and sum(h) == v
    g = sum(a*min(j, d-j) for j, a in enumerate(h))
    kappa = Q(5*g-d*v, 2*v)
    moment = moment_independently(ranks)
    assert kappa == Q(3*d, 4)-Q(5, 2)*moment
    cost = sum((min(Q(14), Q(r, 2)) for r in ranks), Q())
    return dict(ranks=list(ranks), V=v, G=g, curvature=str(kappa),
                cost=str(cost), margin=str(kappa-cost))


def cusp_balance():
    general = [profile(x) for x in concentrated(155)]
    marked = []
    for a in range(1, 9):
        for tail in concentrated(155-a):
            marked.append(dict(distinguished_small_rank=a, **profile((a,)+tail)))
    lower = min(Q(x['margin']) for x in general)
    small = min(Q(x['margin']) for x in marked)
    assert len(general) == 128
    assert lower == -Q(4172075, 82044596)
    assert small == Q(1395158, 20511149) > 0
    assert all(sum(x['ranks']) == 155 for x in general+marked)
    beta_minus_cost = Q(156*150, 8)-14*156
    assert beta_minus_cost == 741
    denominator = 82044596
    coefficients = dict(large_cusp_vertices=4172075,
                        finite_vertices=int(741*denominator),
                        small_cusp_vertices=int(small*denominator),
                        constant=775*denominator)
    assert coefficients == dict(large_cusp_vertices=4172075,
        finite_vertices=60795045636, small_cusp_vertices=5580632,
        constant=63584561900)
    return dict(general_profiles=general, marked_small_profiles=marked,
        general_minimum=str(lower), small_minimum=str(small),
        general_minimizers=[x['ranks'] for x in general if Q(x['margin']) == lower],
        small_minimizers=[x['ranks'] for x in marked if Q(x['margin']) == small],
        balance_coefficients=coefficients,
        restricted_dimension_bound=155,
        hypothesis='Every cusp has at least one irreducible factor of rank at most8; compact chambers are included vacuously.',
        one_cusp_dimension_bound=155,
        inherited_unrestricted_bound=156)


def dot(x, y):
    return -x[0]*y[0]+sum(a*b for a, b in zip(x[1:], y[1:]))


def solve(a, b):
    rows = [[Q(t) for t in row]+[Q(v)] for row, v in zip(a, b)]
    n = len(b)
    for j in range(n):
        pivot = next((i for i in range(j, n) if rows[i][j]), None)
        if pivot is None:
            return None
        rows[j], rows[pivot] = rows[pivot], rows[j]
        p = rows[j][j]
        rows[j] = [t/p for t in rows[j]]
        for i in range(n):
            if i != j:
                p = rows[i][j]
                rows[i] = [t-p*u for t, u in zip(rows[i], rows[j])]
    return tuple(row[-1] for row in rows)


def encoded(x):
    return [str(t) for t in x]


def nonincident_example():
    # Upper half-space: |X|,|Y|<=1, X^2+Y^2+Z^2>=4.
    names = ['R', 'L', 'U', 'D', 'B']
    roots = [tuple(map(Q, x)) for x in [
        (1, 1, 0, 1), (1, -1, 0, 1), (1, 0, 1, 1),
        (1, 0, -1, 1), (Q(-3, 4), 0, 0, Q(-5, 4))]]
    gram = [[dot(a, b) for b in roots] for a in roots]
    assert all(gram[i][i] == 1 for i in range(5))
    assert gram[0][1] == gram[2][3] == -1
    assert all(gram[i][j] == 0 for i in (0, 1) for j in (2, 3))
    assert all(gram[i][4] == -Q(1, 2) for i in range(4))
    center = tuple(map(Q, (5, 0, 0, 4)))
    assert dot(center, center) == -9
    assert all(dot(center, a) < 0 for a in roots)
    dual = tuple(sum(a[j] for a in roots) for j in range(4))
    assert dot(dual, dual) == -3
    section = [-dual[0], *dual[1:]]
    rays, trials = {}, []
    for support in combinations(range(5), 3):
        a = [[-roots[i][0], *roots[i][1:]] for i in support] + [section]
        v = solve(a, [0, 0, 0, -1])
        if v is None:
            trials.append(dict(support=list(support), status='singular'))
            continue
        if any(dot(v, r) > 0 for r in roots):
            trials.append(dict(support=list(support), status='infeasible'))
            continue
        assert v[0] > 0
        v = tuple(t/v[0] for t in v)
        active = [i for i, r in enumerate(roots) if dot(v, r) == 0]
        rays[v] = active
        trials.append(dict(support=list(support), status='feasible', active=active))
    expected = {(Q(1), Q(x, 5), Q(y, 5), Q(3, 5))
                for x in (-2, 2) for y in (-2, 2)}
    apex = (Q(1), Q(0), Q(0), Q(1))
    assert set(rays) == expected | {apex}
    assert all(dot(v, v) < 0 for v in expected) and dot(apex, apex) == 0
    assert rays[apex] == [0, 1, 2, 3] and dot(apex, roots[4]) == -Q(1, 2)
    faces = {names[i]: [v for v, active in rays.items() if i in active] for i in range(5)}
    assert len(faces['B']) == 4 and all(len(faces[x]) == 3 for x in names[:4])
    for i, name in enumerate(names):
        v = tuple(sum(w[j] for w in faces[name]) for j in range(4))
        assert dot(v, v) < 0 and dot(v, roots[i]) == 0
        assert all(dot(v, r) < 0 for j, r in enumerate(roots) if j != i)
    for i, j in [(0, 1), (2, 3)]:
        assert tuple(roots[i][k]+roots[j][k] for k in range(4)) == tuple(2*x for x in apex)
    # Each base vertex has A3 Gram, with the two side roots at its ends.
    # The published rank3 rule gives weight TWO, not zero, to this pair.
    for active in (rays[v] for v in expected):
        i, j = [x for x in active if x != 4]
        assert gram[i][j] == 0 and gram[i][4] == gram[j][4] == -Q(1, 2)
    return dict(metric=[-1, 1, 1, 1], names=names,
        roots=[encoded(x) for x in roots], gram=[encoded(x) for x in gram],
        strict_center=encoded(center), strict_dual=encoded(dual), trials=trials,
        rays=[dict(ray=encoded(v), active=[names[i] for i in a], norm=str(dot(v, v)))
              for v, a in sorted(rays.items())],
        counts=dict(facets=5, finite_vertices=4, ideal_vertices=1, exterior_vertices=0),
        cusp_factors=[['R', 'L'], ['U', 'D']],
        compact_base_star=names, base_corner_weight=2, base_total_weight=8,
        scope='Actual finite-volume rational H3 square pyramid. Its compact base star contains both affine cusp factors although the apex is not on that face. It is NOT an unpaid-face example.')


def run():
    balance = cusp_balance()
    small_costs = []
    for r in range(1, 9):
        assert Q(2*comb(r+1, 2), r+1) == r
        small_costs.append(dict(rank=r, normalized_weight2_cost=r))
    paths = [f'scripts/{TAG}.py'] + INPUTS
    return dict(status='PASS_SMALL_CUSP_LOCALIZATION_AND_DIMENSION156_BALANCE',
        sources={p: digest(ROOT/p) for p in paths},
        cusp_balance=balance, small_factor_weight2_costs=small_costs,
        nonincident_affine_carrier=nonincident_example(),
        geometry_scope='Written face-localization proof inherits the bound finite-weight/ideal/carrier cases; the checker certifies exact counts and the literal counterexample only.',
        unrestricted_dimension_improvement=False,
        no_claim='No bound on nonincident carrier multiplicity and no reduction of the unrestricted156 ceiling.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = run()
    out = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result, indent=2)+'\n')
    b = result['cusp_balance']
    print(result['status'], len(b['general_profiles']), len(b['marked_small_profiles']), flush=True)
