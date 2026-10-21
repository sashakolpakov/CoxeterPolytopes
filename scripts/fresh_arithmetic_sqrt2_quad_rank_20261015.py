#!/usr/bin/env python3
"""A scoped complete-star rank obstruction over Q(sqrt2).

Two disjoint one-4 cycles of orders3/4/5, at least one of order>=4,
with one white bridge (or a unique longer white path), pendant white trees,
and no other contacts. This is not a classification of quadrilateral stars.
"""
import argparse
import json
from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_arithmetic_sqrt2_quad_rank_20261015'


def cycle(k, variables=False, offset=0):
    z = s.symbols('z')
    diag = s.symbols(f'd{offset}:{offset+k}') if variables else [1]*k
    G = s.diag(*diag)
    for i in range(k):
        j = (i+1) % k
        G[i, j] = G[j, i] = -z/2 if i == 0 else -s.Rational(1, 2)
    return G, tuple(diag)


def reduced(expr):
    z = s.symbols('z')
    return s.rem(s.expand(expr), z*z-2, z)


def poly(expr, variables):
    return dict(variables=list(map(str, variables)),
                terms=[[list(m), str(c)] for m, c in s.Poly(expr, *variables).terms()])


def cycle_identities():
    z = s.symbols('z')
    records = []
    for k in (3, 4, 5):
        G, ds = cycle(k, True)
        determinant = reduced(G.det(method='domain-ge'))
        a, b = determinant.coeff(z, 0), -determinant.coeff(z, 1)
        assert b == s.Rational(1, 2**(k-1))
        base_a = a.subs({d: 1 for d in ds})
        assert base_a == (0 if k == 3 else -s.Rational(1, 16))
        cofactors = []
        for i in range(k):
            keep = [j for j in range(k) if i != j]
            minor = reduced(G.extract(keep, keep).det(method='domain-ge'))
            assert not minor.has(z)
            assert s.expand(s.diff(a, ds[i])-minor) == 0
            cofactors.append(poly(minor, ds))
        base = G.subs({d: 1 for d in ds}).subs(z, s.sqrt(2))
        proper = []
        for r in range(1, k):
            for support in combinations(range(k), r):
                value = s.expand(base.extract(support, support).det())
                assert value.is_Rational and value > 0
                proper.append(dict(support=list(support), determinant=str(value)))
        records.append(dict(order=k, rational_part=poly(a, ds), radical_coefficient=str(-b),
                            base_rational_part=str(base_a), cofactors=cofactors,
                            positive_proper_minors=proper))
    return records


def bridge_identities():
    z, x = s.symbols('z x')
    rows = []
    for m, n in [(3, 4), (3, 5), (4, 4), (4, 5), (5, 5)]:
        L, dl = cycle(m, True)
        M, dm = cycle(n, True, m)
        G = s.diag(L, M)
        # The universal bridge/cofactor identity is checked with full
        # symbolic diagonal variables, not a numerical specialization.
        u, v = 2, 2
        G[u, m+v] = G[m+v, u] = -x
        dlm, dmm = reduced(L.det(method='domain-ge')), reduced(M.det(method='domain-ge'))
        li, mi = [i for i in range(m) if i != u], [i for i in range(n) if i != v]
        cl = reduced(L.extract(li, li).det(method='domain-ge'))
        cm = reduced(M.extract(mi, mi).det(method='domain-ge'))
        direct = reduced(G.det(method='domain-ge'))
        assert reduced(direct-dlm*dmm+x*x*cl*cm) == 0
        a, b = dlm.coeff(z, 0), -dlm.coeff(z, 1)
        A, B = dmm.coeff(z, 0), -dmm.coeff(z, 1)
        assert s.expand(direct.coeff(z, 1)+a*B+A*b) == 0
        lower = (s.Rational(1, 16)*B if m >= 4 else 0)+(s.Rational(1, 16)*b if n >= 4 else 0)
        assert lower > 0
        rows.append(dict(orders=[m, n], symbolic_radical_coefficient=poly(direct.coeff(z, 1), dl+dm),
                         exact_lower_bound=str(lower), bridge_square_is_rational=True))
    return rows


def tree_positive(nodes, edges):
    values = {v: F(1) for v in nodes}
    adj = {v: {} for v in nodes}
    for a, b, square in edges:
        if a in adj and b in adj:
            adj[a][b] = adj[b][a] = square
    pivots = []
    while values:
        leaves = [v for v in values if len(adj[v]) <= 1]
        if not leaves:
            return False, pivots, 'cycle'
        v = min(leaves)
        pivots.append([v, str(values[v])])
        if values[v] <= 0:
            return False, pivots, 'nonpositive pivot'
        if adj[v]:
            u, square = next(iter(adj[v].items()))
            values[u] -= square/values[v]
            del adj[u][v]
        del values[v], adj[v]
    return True, pivots, None


def marked_distance(nodes, edges, a, b):
    adj = {v: [] for v in nodes}
    for u, v, _ in edges:
        if u in adj and v in adj:
            adj[u].append(v)
            adj[v].append(u)
    todo, seen = deque([(a, 0)]), {a}
    while todo:
        u, d = todo.popleft()
        if u == b:
            return d
        for v in adj[u]:
            if v not in seen:
                seen.add(v)
                todo.append((v, d+1))
    raise AssertionError('disconnected corner')


def placements():
    counts, survivors = Counter(), []
    failures = Counter()
    for m, n in product((3, 4, 5), repeat=2):
        if m == n == 3:
            continue
        base_edges = [(i, (i+1) % m, F(1, 2) if i == 0 else F(1, 4)) for i in range(m)]
        base_edges += [(m+i, m+(i+1) % n, F(1, 2) if i == 0 else F(1, 4)) for i in range(n)]
        for left in combinations(range(m), 2):
            for right in combinations(range(m, m+n), 2):
                W = set(range(m+n))-set(left)-set(right)
                for u in sorted(set(range(m))-set(left)):
                    for v in sorted(set(range(m, m+n))-set(right)):
                        for square in (F(1, 4), F(1, 2)):
                            counts['total'] += 1
                            key = f'{m},{n}'
                            counts['total '+key] += 1
                            edges = base_edges+[(u, v, square)]
                            corners = []
                            for a, b in product(left, right):
                                nodes = W | {a, b}
                                positive, pivots, reason = tree_positive(nodes, edges)
                                if not positive:
                                    failures[key] += 1
                                    break
                                corners.append(dict(marked=[a, b], pivots=pivots,
                                                    distance=marked_distance(nodes, edges, a, b)))
                            else:
                                counts['four_PD'] += 1
                                counts['four_PD '+key] += 1
                                L, _ = cycle(m)
                                M, _ = cycle(n)
                                G = s.diag(L, M).subs(s.symbols('z'), s.sqrt(2))
                                G[u, v] = G[v, u] = -s.sqrt(s.Rational(square.numerator, square.denominator))
                                determinant = s.expand(G.det(method='domain-ge'))
                                coeff = determinant.coeff(s.sqrt(2))
                                assert coeff > 0
                                survivors.append(dict(orders=[m, n], left_black=list(left), right_black=list(right),
                                                      bridge=[u, v], bridge_square=str(square), corners=corners,
                                                      determinant=str(determinant), radical_coefficient=str(coeff)))
    assert counts['total'] == 4032 and counts['four_PD'] == 28
    assert sum(counts['total '+f'{m},{n}'] for m, n in product((4, 5), repeat=2)) == 3528
    assert all(not counts['four_PD '+f'{m},{n}'] for m, n in product((4, 5), repeat=2))
    return dict(counts=dict(sorted(counts.items())), finite_corner_failures=dict(sorted(failures.items())),
                survivors=survivors,
                minimum_marked_distance_counts=dict(Counter(str(min(x['distance'] for x in r['corners'])) for r in survivors)))


def open_examples():
    # These lie outside the proved cycle-pair scope and show why it is not
    # a complete quadrilateral classification or a distance-two theorem.
    z = s.sqrt(2)
    a = s.sqrt(2+z)/2
    path = s.eye(6)
    for i, w in enumerate([a, s.Rational(1, 2), z/2, s.Rational(1, 2), a]):
        path[i, i+1] = path[i+1, i] = -w
    assert s.simplify(path.det()) == 0
    C, _ = cycle(3)
    pair = s.diag(C, C).subs(s.symbols('z'), z)
    pair[2, 5] = pair[5, 2] = -z/2
    assert s.expand(pair.det()) == 0
    # Each retained-black corner is the finite path3--4--3 (F4).
    f4 = []
    for u, v in product((0, 1), (3, 4)):
        support = [u, 2, 5, v]
        G = pair.extract(support, support)
        assert all(G[:j, :j].det() > 0 for j in range(1, 5))
        f4.append(dict(marked=[u, v], support=support, determinant=str(G.det())))
    return dict(path_labels=[8, 3, 4, 3, 8], path_determinant='0',
                two_C3_bridge_label=4, two_C3_determinant='0', four_F4_corners=f4,
                next_weight_lead='Distance<=2 plus one unit on the two end nodes of an entire F4 component keeps component bound2r-2. No exhaustive face proof is asserted.')


def run():
    paths = [f'scripts/{TAG}.py', 'scripts/fresh_arithmetic_sqrt2_local_20261015.py',
             'results/fresh_arithmetic_sqrt2_local_20261015.json',
             'scripts/fresh_arithmetic_sqrt2_local_20261015_review.py',
             'results/fresh_arithmetic_sqrt2_local_20261015_review.json']
    return dict(status='PASS_SCOPED_SQRT2_TWO_CYCLE_QUADRILATERAL_RANK_OBSTRUCTION',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                scope=dict(arithmetic_ground_field='Q(sqrt2)', compact_two_face=True,
                           disjoint_cycles=[3, 4, 5], at_least_one_order_ge4=True,
                           each_cycle='one4 and all other bonds3',
                           connection='one white bridge or a unique longer white path',
                           allowed_extras='pendant white spherical trees, each attached at one white node; orthogonal white components',
                           excluded_from_scope=['shared white nodes', 'extra black contacts', 'extra cross edges', 'other Lanner triangles', 'two cycles both of order3'],
                           global_dimension_bound=False),
                cycle_identities=cycle_identities(), bridge_identities=bridge_identities(),
                complete_direct_placement_census=placements(), open_examples=open_examples(),
                longer_path_gate='Original disjoint Lanner cycles are orthogonal if their only connector has at least2 bonds; two negative directions.',
                proof_limits='All pendant-tree and diagonal-monotonicity arguments are proved in the report, not inferred from bounded tree samples.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = json.loads(json.dumps(run(), default=str))
    out = ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert json.loads(out.read_text()) == result, 'saved certificate mismatch'
        print('PASS saved replay')
    else:
        out.write_text(json.dumps(result, indent=2)+'\n')
        print(result['status'])
    print(json.dumps(result['complete_direct_placement_census']['counts']))
