#!/usr/bin/env python3
"""All-real E5 cutter obstruction from an exact rational outer polytope.

Uses only frozen independent geometry. Numerical discovery is not an input.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from fresh_h20_h4_edge_escape_gate_20261016_review import Interval as I
from fresh_h20_residual_cutter_domain_20261016 import sqrt_interval

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_second_cusp_cutter_domain_20261016'
GEOMETRY = 'results/fresh_h20_residual_h4_corner_20261016_review.json'
PARTNER = 'results/fresh_h20_h4_parallel_partner_20261016_review.json'
CUSP = 'results/fresh_h20_second_cusp_completion_20261016_review.json'


def digest(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()


def solve(rows, rhs):
    a = [list(row) + [b] for row, b in zip(rows, rhs)]
    n = len(a)
    for k in range(n):
        pivot = next((i for i in range(k, n) if a[i][k]), None)
        if pivot is None:
            return None
        a[k], a[pivot] = a[pivot], a[k]
        p = a[k][k]
        a[k] = [v / p for v in a[k]]
        for i in range(n):
            if i != k:
                factor = a[i][k]
                a[i] = [v - factor * w for v, w in zip(a[i], a[k])]
    return tuple(row[-1] for row in a)


def run():
    g, p, c = [json.loads((ROOT / name).read_text()) for name in (GEOMETRY, PARTNER, CUSP)]
    assert g['status'] == 'PASS_INDEPENDENT_145_WALL_RESIDUAL_FINITE_CORNER_AND_ESCAPE_RECTANGLE'
    assert p['status'] == 'PASS_INDEPENDENT_146_WALL_FINITE_PARALLEL_CUTTER_AND_IDEAL_VERTEX'
    assert c['status'] == 'PASS_INDEPENDENT_151_WALL_SECOND_FULL_RANK_CUSP_AND_ESCAPE'
    for data in (g, p, c):
        for name, value in data['sources'].items():
            assert digest(name) == value, name
    q = list(map(F, g['four_norms']))
    C = list(map(F, g['corner']['chart']))
    metric = list(map(F, [-1, 1, 1])) + [F(1, 2)] * 2
    dot = lambda a, b: sum(x * y * t for x, y, t in zip(a, b, metric))
    d = -dot(C, C)
    assert d > 0 and str(-d) == g['corner']['norm']
    assert g['four_gram'] == [[str(q[i] if i == j else F(0)) for j in range(4)] for i in range(4)]
    scales = [sqrt_interval(I(d / v)) for v in q]
    pairs = {tuple(row['walls']): I(*row['bounds']) for row in g['new_pairs']}
    # Old133 is used only to force the coefficient of C to be positive:
    # all four cap products and C.133 are negative.
    assert I(*g['corner']['all145_pairs'][133]).hi < 0
    assert all(pairs[(133, 141 + i)].hi < 0 for i in range(4))
    coefficients = []
    walls = [31, 95, 140, 145, 150]
    for old in walls[:3]:
        cp = I(*g['corner']['all145_pairs'][old])
        assert cp.hi < 0
        coefficients.append([pairs[(old, 141 + i)] / cp * scales[i] for i in range(4)])
    z = I(*p['positive_parameter'])
    assert z.lo > 0
    # From the exact recipes, after dividing by C.n_j:
    # n145=-sqrt(2/q141)n141+sqrt(2)z(C/sqrt(d)-n144/sqrt(q144))
    # n150=-n143+6(C/sqrt(d)-n144/sqrt(q144)).
    assert c['heights'][-1] == 6 and F(c['new_norms'][-1]) == q[2]
    coefficients += [[I(1) / z, I(0), I(0), I(1)],
                     [I(0), I(0), sqrt_interval(I(q[2])) / 6, I(1)]]
    e = [F(1), F(-1, 100), F(-49, 50), F(-3, 10), F(-1, 25)]
    ce = dot(C, e)
    assert ce < 0 and c['strict_escape']['norm'] == '31500061/5000000000'
    escape = [I(*c['strict_escape']['all151_pairs'][141 + i]) / ce * scales[i] for i in range(4)]
    assert all(v.lo > 0 for v in escape)
    denominator = 10**6
    lower = lambda v: F((v.lo * denominator).__floor__(), denominator)
    upper = lambda v: F((v.hi * denominator).__ceil__(), denominator)
    # Since z_i>=0, lowering coefficients of A.z<=1 and raising those
    # of B.z>=1 both enlarge the permitted region.
    A = [[lower(v) for v in row] for row in coefficients]
    rhs = [F(1)] * 5
    B = [upper(v) for v in escape]
    A.append([-v for v in B]); rhs.append(F(-1))
    A += [[F(-int(i == j)) for i in range(4)] for j in range(4)]
    rhs += [F(0)] * 4
    assert all(A[2][i] > 0 for i in range(4))
    # Nonnegativity and positive row140 bound each coordinate, so the
    # outer polyhedron is bounded. Exhaust every possible vertex basis.
    vertices = {}
    nonsingular = feasible = 0
    for ids in combinations(range(len(A)), 4):
        v = solve([A[i] for i in ids], [rhs[i] for i in ids])
        if v is None:
            continue
        nonsingular += 1
        assert all(sum(A[i][j] * v[j] for j in range(4)) == rhs[i] for i in ids)
        slacks = [b - sum(x * y for x, y in zip(row, v)) for row, b in zip(A, rhs)]
        if min(slacks) < 0:
            continue
        feasible += 1
        vertices[v] = [i for i, s in enumerate(slacks) if s == 0]
    assert vertices
    maxima = {omit: max(sum(v[j]**2 for j in range(4) if j != omit) for v in vertices) for omit in (0, 3)}
    bounds = {0: F(99, 100), 3: F(9, 20)}
    assert all(maxima[i] < bounds[i] < 1 for i in bounds)
    rows = []
    for v, active in sorted(vertices.items()):
        rows.append(dict(coordinates=list(map(str, v)), active_constraints=active,
                         projection_squares={str(141 + i): str(sum(v[j]**2 for j in range(4) if j != i)) for i in bounds}))
    paths = [f'scripts/{TAG}.py', GEOMETRY, PARTNER, CUSP,
             'scripts/fresh_h20_h4_edge_escape_gate_20261016_review.py',
             'scripts/fresh_h20_residual_cutter_domain_20261016.py']
    return dict(status='PASS_ALL_REAL_RESIDUAL_E5_CUTTERS_DOTTED_TO_141_AND_144',
                sources={name: digest(name) for name in paths},
                old_walls=walls, coordinate_definition='r=C-sum(u_i*n_i), z_i=u_i*sqrt(q_i/d)>=0, d=-C^2; normalized positive C coefficient',
                exact_constraint_intervals=[[v.saved() for v in row] for row in coefficients],
                exact_E5_constraint_intervals=[v.saved() for v in escape],
                rational_outer_matrix=[[str(v) for v in row] for row in A],
                rational_outer_rhs=list(map(str, rhs)), rounding_denominator=denominator,
                boundedness_coordinate_bounds=[str(1 / v) for v in A[2]],
                bases_checked=210, nonsingular_bases=nonsingular, feasible_bases=feasible,
                vertex_count=len(rows), vertices=rows,
                squared_projection_maxima={str(141 + i): str(maxima[i]) for i in bounds},
                squared_projection_upper_bounds={str(141 + i): str(bounds[i]) for i in bounds},
                physical_projection_norm_upper_bounds={str(141 + i): str(-d * (1 - bounds[i])) for i in bounds},
                proof='Acuteness to caps and old133 forces a positive C coefficient. The five old-wall constraints and E5cut constraint are outer-relaxed with exact directed rational rounding. Nonnegative coordinates and positive row140 make the relaxation bounded. All210four-row bases are exhausted; every vertex has both indicated projected squared lengths strictlybelow1. Convexity extends the bounds to the whole relaxation. Projection of each spacelike cutter to n141-perp andn144-perp is timelike, so acuteness forces both contacts strictlydotted.',
                scope='Every real spacelike normal righttoS16, acute to caps141..144 andold31,95,133,140,145,150, with r.E5>=0. No field,angle,norm,height orretained-center restriction. Doesnot exclude cutters dottedtoboth141and144 orwalls leavingS16. E5cut hypothesis isessential; nofinitevolume claim.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--verify', action='store_true'); args = ap.parse_args()
    result = run(); out = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], flush=True)
