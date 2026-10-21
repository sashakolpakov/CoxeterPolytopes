#!/usr/bin/env python3
"""Exact data for an all-real-height obstruction to a bare H18 lift.

No height sampling or Coxeter-label cutoff. The report proves the universal
geometry from the literal Gram identities certified here. Arithmetic below
is an independent two-coefficient Fraction implementation of Q(sqrt(2)).
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h18_general_height_lift_20261014'
SOURCE = 'data/enumerations/vinberg_h18_mixture_facets_2026_10_13.json'
PREVIOUS = 'results/fresh_dimension_construction_routes_20261014.json'
Z, O = (F(0), F(0)), (F(1), F(0))

def q(a, b=0): return F(a), F(b)
def add(a, b): return a[0]+b[0], a[1]+b[1]
def neg(a): return -a[0], -a[1]
def sub(a, b): return add(a, neg(b))
def mul(a, b): return a[0]*b[0]+2*a[1]*b[1], a[0]*b[1]+a[1]*b[0]
def inv(a):
    n = a[0]*a[0]-2*a[1]*a[1]
    assert n
    return a[0]/n, -a[1]/n
def div(a, b): return mul(a, inv(b))
def total(xs):
    ans = Z
    for x in xs: ans = add(ans, x)
    return ans
def sign(a):
    x, y = a
    if not y: return (x > 0)-(x < 0)
    if not x: return (y > 0)-(y < 0)
    if x*y > 0: return (x > 0)-(x < 0)
    d = x*x-2*y*y
    assert d
    return ((x > 0)-(x < 0))*((d > 0)-(d < 0))

def inertia(matrix):
    a = [row[:] for row in matrix]
    positive = negative = 0
    pivots = []
    while a:
        n = len(a)
        i = next((i for i in range(n) if a[i][i] != Z), None)
        if i is not None:
            rest = [j for j in range(n) if j != i]
            p = a[i][i]
            positive += sign(p) > 0
            negative += sign(p) < 0
            pivots.append(dict(size=1, value=p, sign=sign(p)))
            a = [[sub(a[j][k], div(mul(a[j][i], a[i][k]), p))
                  for k in rest] for j in rest]
        else:
            pair = next(((i, j) for i in range(n) for j in range(i)
                         if a[i][j] != Z), None)
            if pair is None:
                return (positive, negative, n), pivots
            i, j = pair
            b = a[i][j]
            rest = [k for k in range(n) if k not in pair]
            pivots.append(dict(size=2, off_diagonal=b, inertia=[1, 1]))
            a = [[sub(a[k][l], div(add(mul(a[k][i], a[j][l]),
                                       mul(a[k][j], a[i][l])), b))
                  for l in rest] for k in rest]
            positive += 1
            negative += 1
    return (positive, negative, 0), pivots

def digest(path): return sha256((ROOT/path).read_bytes()).hexdigest()
def encode(x):
    if isinstance(x, F): return str(x)
    if isinstance(x, dict): return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [encode(v) for v in x]
    return x

def run():
    data = json.loads((ROOT/SOURCE).read_text())
    prior = json.loads((ROOT/PREVIOUS).read_text())
    for p, expected in prior['sources'].items(): assert digest(p) == expected, p
    names = data['labels']
    G = [[q(x['a'], x['b']) for x in row] for row in data['gram']]
    assert len(G) == len(names) == 40
    assert all(G[i][j] == G[j][i] for i in range(40) for j in range(40))
    assert all(sign(G[i][j]) <= 0 for i in range(40) for j in range(i))
    signature, pivots = inertia(G)
    assert signature == (18, 1, 21)
    short = [i for i in range(40) if G[i][i] == O]
    long = [i for i in range(40) if G[i][i] == q(2)]
    assert len(short) == 14 and len(long) == 26
    h_pairings = [total(row) for row in G]
    h_norm = total(h_pairings)
    assert sign(h_norm) < 0 and all(sign(x) < 0 for x in h_pairings)
    core = [names.index(name) for name in ('L18', 'L25', 'L26')]
    assert set(core) <= set(short)
    assert all(G[i][j] == q(-1) for i, j in combinations(core, 2))
    anchors = []
    for i in core:
        js = [j for j in long if G[i][j] == q(-1)]
        assert js
        anchors.append(dict(short=i, name=names[i], norm=G[i][i],
                            fixed_long_anchors=js, anchor_names=[names[j] for j in js]))
    null_rows = []
    for i, j in combinations(core, 2):
        pairings = [add(G[i][a], G[j][a]) for a in range(40)]
        delta_norm = add(add(G[i][i], G[j][j]), mul(q(2), G[i][j]))
        delta_h = add(h_pairings[i], h_pairings[j])
        assert delta_norm == Z and sign(delta_h) < 0
        assert pairings[i] == pairings[j] == Z
        assert all(sign(pairings[a]) < 0 for a in short if a not in (i, j))
        assert all(sign(pairings[a]) <= 0 for a in long)
        tau = div(q(-1, 0), mul(q(4), delta_h))
        assert sign(tau) > 0
        leading = add(O, mul(q(2), mul(tau, delta_h)))
        quartic = mul(mul(tau, tau), h_norm)
        assert leading == q(F(1, 2)) and sign(quartic) < 0
        null_rows.append(dict(pair=[i, j], names=[names[i], names[j]],
                              norm=delta_norm, old_wall_pairings=pairings,
                              center_pairing=delta_h, tau=tau,
                              escape_square_quadratic_coefficient=leading,
                              escape_square_quartic_coefficient=quartic,
                              sufficient_squared_epsilon_bound=div(neg(leading), quartic)))
    # For a core lift r+v, its fixed norm-two anchor gives cosine squared
    # 1/[2(1+|v|^2)] in (0,1/2]. Only label3 (1/4) and4 (1/2)
    # are possible, yielding |v|^2=1 and0 respectively.
    core_height_squares = {}
    for label, cosine_square in ((3, q(F(1, 4))), (4, q(F(1, 2)))):
        value = sub(inv(mul(q(2), cosine_square)), O)
        core_height_squares[str(label)] = value
    assert core_height_squares == {'3': O, '4': Z}
    paths = [SOURCE, PREVIOUS, f'scripts/{TAG}.py']
    return encode(dict(
        status='PASS_ALL_REAL_SOURCE_DIRECTION_HEIGHT_LIFTS_DIMENSION_AT_LEAST_20_INFINITE',
        sources={p: digest(p) for p in paths}, source_gram_inertia=signature,
        exact_inertia_pivots=pivots, source_dimension=18,
        added_dimension_range='Every integer k>=2; target H^(18+k).',
        short_indices=short, fixed_long_indices=long,
        core_indices=core, core_anchor_data=anchors,
        strict_center_old_root_coefficients=[1]*40,
        strict_center_norm=h_norm, strict_center_pairings=h_pairings,
        core_height_square_by_anchor_label=core_height_squares,
        feasible_core_height_sets='A singleton {0}, or one antipodal unit pair; uses 0 in their convex hull and all-pair acuteness.',
        core_unit_line_cross_contacts='If both endpoints have antipodal unit pairs, their lines are equal or orthogonal by the gap (0,cos(pi/3)).',
        null_pairs=null_rows,
        universal_escape='Choose a core pair with height span < k, a unit perpendicular u, and X=delta+epsilon*u+epsilon^2*tau*h. Then X^2=epsilon^2/2+epsilon^4*tau^2*h^2>0 and every wall pairing is strict for sufficiently small positive epsilon.',
        hypotheses=[
            'Retain all 26 norm-two source walls at zero added height.',
            'Every other wall has old component one of the 14 norm-one source roots, after positive normalization; each height set is finite.',
            'For each of L18,L25,L26, the height set is nonempty and its convex hull contains zero.',
            'All distinct wall pairs satisfy the acute Coxeter condition.',
        ],
        scope='A universal nonexistence result for the displayed bare fixed-source-direction lift class, allowing arbitrary real, unequal and asymmetric finite height sets. No claim about added walls with new old projections, deleting fixed anchors, or arbitrary high-dimensional Coxeter polytopes. No new finite-volume example.',
    ))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    out = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['source_gram_inertia'], result['core_indices'])
