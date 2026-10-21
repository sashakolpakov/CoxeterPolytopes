#!/usr/bin/env python3
"""Exact quadrilateral payments and nonincident-carrier multiplicities.

No optimizer; no new dimension bound. Individual triangle payments at E8
are deliberately not reconstructed from the inherited count of four exceptions.
"""
import argparse
from collections import Counter, deque
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
import json
from math import gcd, prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonincident_carrier_incidence_20261017'
OLD = 'fresh_selective_distance13_dimension146_20261017'
NODES = frozenset(range(19))
LEFT, RIGHT = frozenset(range(9)), frozenset(range(10, 19))
RADICAL = [-6, -3, -4, -2, -5, -4, -3, -2, -1, 0,
           6, 3, 4, 2, 5, 4, 3, 2, 1]
EDGES = [(0, 1), (0, 2), (0, 4), (2, 3), (4, 5),
         (5, 6), (6, 7), (7, 8), (8, 9)]
EDGES += [(a + 10 if a < 9 else 9, b + 10 if b < 9 else 9)
          for a, b in EDGES.copy()]
GRAPH = {i: set() for i in NODES}
for a, b in EDGES:
    GRAPH[a].add(b)
    GRAPH[b].add(a)
GRAM = [[F(2 if i == j else -int(j in GRAPH[i]))
         for j in range(19)] for i in range(19)]


def sha(path):
    h = sha256()
    with (ROOT / path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def inverse(a):
    n = len(a)
    m = [[F(x) for x in row] + [F(i == j) for j in range(n)]
         for i, row in enumerate(a)]
    for j in range(n):
        p = next(i for i in range(j, n) if m[i][j])
        m[j], m[p] = m[p], m[j]
        x = m[j][j]
        m[j] = [v / x for v in m[j]]
        for i in range(n):
            if i != j:
                x = m[i][j]
                m[i] = [u - x * v for u, v in zip(m[i], m[j])]
    return [r[n:] for r in m]


def subgram(nodes):
    return [[GRAM[i][j] for j in nodes] for i in nodes]


def bilinear(a, g, b):
    return sum((a[i] * g[i][j] * b[j] for i in range(len(a))
                for j in range(len(b))), F())


def positive_definite(a):
    a = [r.copy() for r in a]
    for k in range(len(a)):
        if a[k][k] <= 0:
            return False
        for i in range(k + 1, len(a)):
            for j in range(k + 1, len(a)):
                a[i][j] -= a[i][k] * a[k][j] / a[k][k]
    return True


def component(nodes, a):
    nodes = set(nodes)
    seen, queue = {a}, [a]
    for i in queue:
        for j in GRAPH[i] & nodes - seen:
            seen.add(j)
            queue.append(j)
    return frozenset(seen)


def components(nodes):
    rest, ans = set(nodes), []
    while rest:
        c = component(rest, min(rest))
        ans.append(c)
        rest -= c
    return ans


def distances(nodes, a):
    nodes = set(nodes)
    out, queue = {a: 0}, deque([a])
    while queue:
        i = queue.popleft()
        for j in GRAPH[i] & nodes - out.keys():
            out[j] = out[i] + 1
            queue.append(j)
    return out


@lru_cache(None)
def finite_type(nodes):
    nodes = set(nodes)
    adj = {i: GRAPH[i] & nodes for i in nodes}
    forks = [i for i in nodes if len(adj[i]) == 3]
    assert sum(map(len, adj.values())) == 2 * len(nodes) - 2
    if not forks:
        name = 'A' + str(len(nodes))
    else:
        assert len(forks) == 1
        f = forks[0]
        lengths = []
        for j in adj[f]:
            previous, current, size = f, j, 1
            while len(adj[current]) == 2:
                previous, current = current, next(x for x in adj[current]
                                                   if x != previous)
                size += 1
            lengths.append(size)
        lengths.sort()
        if lengths[:2] == [1, 1]:
            name = 'D' + str(len(nodes))
        else:
            assert lengths in ([1, 2, 2], [1, 2, 3], [1, 2, 4])
            name = 'E' + str(len(nodes))
    assert positive_definite(subgram(sorted(nodes)))
    return name


@lru_cache(None)
def classical_weights(nodes, ideal=False):
    nodes = frozenset(nodes)
    adj = {i: GRAPH[i] & nodes for i in nodes}
    forks = [i for i in nodes if len(adj[i]) == 3]
    if ideal:
        assert len(nodes) == 17 and len(forks) == 2
    else:
        name = finite_type(tuple(sorted(nodes)))
        assert name[0] in 'AD', name
    if len(nodes) <= 7:
        return {p: 2 if len(nodes) <= 4 else 1
                for p in combinations(sorted(nodes), 2)}
    # Literal definition of an invariant classical end segment, not a
    # graph-distance substitute for Vinberg's finite rule.
    end_pairs = set()
    leaves = {i for i in nodes if len(adj[i]) == 1}
    fork_leaves = [{x for x in adj[f] if len(adj[x]) == 1} for f in forks]
    assert all(len(s) == 2 for s in fork_leaves)
    for k in range(1, 7):
        for row in combinations(sorted(nodes), k):
            s = frozenset(row)
            if not s & leaves or component(s, min(s)) != s:
                continue
            if any(bool(t & s) and not t <= s for t in fork_leaves):
                continue
            end_pairs.update(combinations(row, 2))
    out = {}
    for a, b in combinations(sorted(nodes), 2):
        d = distances(nodes, a)[b]
        out[a, b] = int((d <= 3 if ideal else d == 1) or (a, b) in end_pairs)
    assert sum(out.values()) <= 3 * len(nodes)
    if ideal:
        assert sum(out.values()) == 51
    return out


def corner(omitted, marked):
    nodes = NODES - set(omitted)
    a, b = sorted(marked)
    c = component(nodes, a)
    ideal = set(omitted) == {3, 13}
    if b not in c:
        return dict(omitted=sorted(omitted), marked=[a, b], weight=0,
                    connected=False, ideal=ideal)
    if not ideal:
        name = finite_type(tuple(sorted(c)))
        assert name[0] in 'AD'  # No unspecified E8 marked-pair input.
    else:
        name = 'affine_D16'
    return dict(omitted=sorted(omitted), marked=[a, b], connected=True,
                ideal=ideal, component=sorted(c), type=name,
                distance=distances(c, a)[b],
                weight=classical_weights(tuple(sorted(c)), ideal)[a, b])


def primitive(v):
    g = 0
    for x in v:
        g = gcd(g, abs(x))
    return [x // g for x in v]


def pyramid():
    assert all(sum(GRAM[i][j] * RADICAL[j] for j in NODES) == 0 for i in NODES)
    ids = list(range(1, 19))
    j = inverse(subgram(ids))
    assert positive_definite(subgram(sorted(NODES - {1, 11})))
    ps = []
    p = [0] * 19
    p[9] = -1
    ps.append(p)
    for a, b in product(sorted(LEFT), sorted(RIGHT)):
        p = [0] * 19
        p[a], p[b] = -RADICAL[b], RADICAL[a]
        ps.append(primitive(p))
    rows = []
    for p in ps:
        assert sum(a * b for a, b in zip(p, RADICAL)) == 0
        q = bilinear(p[1:], j, p[1:])
        assert q <= 0
        rows.append(dict(covector=p, square=str(q), ideal=q == 0,
                         support=[i for i, x in enumerate(p) if x]))
    center = [sum(p[i] for p in ps) for i in NODES]
    assert max(center) < 0 and bilinear(center[1:], j, center[1:]) < 0
    assert all(bilinear(p[1:], j, center[1:]) < 0 for p in ps)
    old = json.loads((ROOT / 'results' / (OLD + '.json')).read_text())
    previous = old['known_sharp_distance14_pyramid']['extreme_rays']
    assert len(previous) == len(rows) == 82
    for a, b in zip(rows, previous):
        p = list(map(F, b['covector']))
        i = a['support'][0]
        scale = p[i] / a['covector'][i]
        assert p == [scale * x for x in a['covector']]
        assert F(b['square']) == scale * scale * F(a['square'])
    assert Counter(r['ideal'] for r in rows) == {False: 80, True: 2}
    assert [r['support'] for r in rows if r['ideal']] == [[9], [3, 13]]
    # Full factors of the non-simple ideal apex, with actual positive null
    # combinations giving precisely its covector.
    factors = []
    for c in (LEFT, RIGHT):
        k = [abs(RADICAL[i]) if i in c else 0 for i in NODES]
        pairings = [sum(GRAM[i][h] * k[h] for h in NODES) for i in NODES]
        assert pairings == [F(-int(i == 9)) for i in NODES]
        factors.append(dict(nodes=sorted(c), positive_coefficients=[k[i] for i in sorted(c)],
                            rank=8, type='affine_E8'))
    return rows, factors


def face_ledger(vertices, factors):
    records = []
    for allowed in combinations(range(19), 4):
        allowed = frozenset(allowed)
        if not (allowed & LEFT and allowed & RIGHT):
            continue
        face_vertices = [i for i, r in enumerate(vertices) if set(r['support']) <= allowed]
        assert len(face_vertices) in (3, 4)
        sides = [i for i in allowed
                 if sum(vertices[v]['covector'][i] == 0 for v in face_vertices) >= 2]
        assert len(sides) == len(face_vertices)
        whites = NODES - allowed
        forest = [finite_type(tuple(sorted(c))) for c in components(whites)]
        star = whites | set(sides)
        carriers = [i for i, c in enumerate(factors) if set(c['nodes']) <= star]
        on_apex = 0 in face_vertices
        row = dict(white=sorted(whites), side=sorted(sides), vertices=face_vertices,
                   compact=not any(vertices[v]['ideal'] for v in face_vertices),
                   polygon=len(sides), support_forest=sorted(forest),
                   small_carrier_factors=carriers, carrier_cusp_on_face=on_apex)
        if len(sides) == 4:
            assert 9 not in allowed and len(allowed & LEFT) == len(allowed & RIGHT) == 2
            assert carriers == [0, 1] and not on_apex and star == NODES
            l, r = sorted(allowed & LEFT), sorted(allowed & RIGHT)
            cc = [corner((a, b), sorted(allowed - {a, b})) for a, b in product(l, r)]
            row.update(left=l, right=r, corners=cc, weight=sum(c['weight'] for c in cc))
        else:
            assert len(carriers) == 1
            row['triangle_weight_scope'] = 'Not reconstructed: literal four E8 exclusions are not imported.'
        records.append(row)
    assert len(records) == 3456
    counts = Counter(r['polygon'] for r in records)
    assert counts == {3: 2160, 4: 1296}
    quads = [r for r in records if r['polygon'] == 4]
    weights = Counter(r['weight'] for r in quads)
    assert weights == {0: 1048, 1: 243, 2: 5}
    compact = Counter(r['weight'] for r in quads if r['compact'])
    assert compact == {0: 987, 1: 240, 2: 5}
    specified = next(r for r in quads if r['left'] == [1, 2] and r['right'] == [11, 12])
    assert specified['compact'] and specified['weight'] == 0
    assert [c['type'] for c in specified['corners']] == ['A17', 'A16', 'A16', 'A15']
    carrier_counts = []
    for i in range(2):
        subset = [r for r in records if i in r['small_carrier_factors']]
        carrier_counts.append(dict(factor=i, all_two_faces=len(subset),
            incident_triangles=sum(r['polygon'] == 3 and r['carrier_cusp_on_face'] for r in subset),
            nonincident_triangles=sum(r['polygon'] == 3 and not r['carrier_cusp_on_face'] for r in subset),
            nonincident_quadrilaterals=sum(r['polygon'] == 4 for r in subset),
            unpaid_quadrilaterals=sum(r['polygon'] == 4 and r['weight'] == 0 for r in subset),
            compact_unpaid_quadrilaterals=sum(r['polygon'] == 4 and r['weight'] == 0 and r['compact'] for r in subset)))
    assert all(r == dict(factor=i, all_two_faces=2376, incident_triangles=324,
                        nonincident_triangles=756, nonincident_quadrilaterals=1296,
                        unpaid_quadrilaterals=1048, compact_unpaid_quadrilaterals=987)
               for i, r in enumerate(carrier_counts))
    link_vertices, link_edges, factor_normalization = 81, 648, 162
    assert 1048 > link_edges and 987 > 6 * factor_normalization
    return dict(two_faces=records, counts=dict(triangles=2160, quadrilaterals=1296),
        quad_weight_counts=dict(sorted(weights.items())), compact_quad_weight_counts=dict(sorted(compact.items())),
        specified_quad=specified, small_carrier_counts=carrier_counts,
        apex_link=dict(type='Delta8 x Delta8', vertices=link_vertices, edges=link_edges,
                       triangles=1512, quadrilaterals=1296),
        false_multiplicity_bounds=[
            dict(proposed='Q_bad <= sum edges(link at small carrier cusp)', left=1048, right=648),
            dict(proposed='Q_bad <= 6 sum small_factor_count * link_vertices', left=1048, right=972),
            dict(proposed='compact Q_bad <= 6 sum small_factor_count * link_vertices', left=987, right=972)],
        necessary_factor_budget_coefficient=str(F(1048, 162)),
        necessary_compact_factor_budget_coefficient=str(F(987, 162)))


# A+b sqrt2 field, used only for a separate unbounded all-carrier family.
def add(x, y): return (x[0] + y[0], x[1] + y[1])
def neg(x): return (-x[0], -x[1])
def mul(x, y): return (x[0] * y[0] + 2 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])
def rational(x): return (F(x), F())
ZERO = rational(0)


def lorentz(x, y):
    out = neg(mul(x[0], y[0]))
    for i in range(1, 4): out = add(out, mul(x[i], y[i]))
    return out


def strip_checks():
    rows = []
    for n in (1, 2, 3, 7, 16, 64):
        a = rational(F(1, 2))
        low = (F(), F(1, 2))
        high = (F(), F(2 * n - 1, 2))
        normals = [[a, rational(1), ZERO, a], [a, rational(-1), ZERO, a],
                   [low, ZERO, rational(-1), low], [high, ZERO, rational(1), high]]
        normals += [[rational(k * k), ZERO, (F(), F(k)), rational(k * k - 1)]
                    for k in range(n)]
        assert all(lorentz(x, x) == rational(1) for x in normals)
        assert lorentz(normals[0], normals[1]) == rational(-1)
        assert lorentz(normals[2], normals[3]) == rational(-1)
        for x in range(2):
            for y in (2, 3): assert lorentz(normals[x], normals[y]) == ZERO
        counts = Counter()
        for k in range(n):
            assert lorentz(normals[4 + k], normals[0]) == rational(F(-1, 2))
            assert lorentz(normals[4 + k], normals[1]) == rational(F(-1, 2))
            assert lorentz(normals[4 + k], normals[2]) == (F(), F(-(2 * k + 1), 2))
            assert lorentz(normals[4 + k], normals[3]) == (F(), F(-(2 * (n - 1 - k) + 1), 2))
            for j in range(k):
                p = lorentz(normals[4 + k], normals[4 + j])
                assert p == rational(1 - (k - j) ** 2)
        # Explicit upper-halfspace finite vertices (x,y,z)=(+-1/2,
        # (2j-1)/sqrt2,1/2), normalized hyperboloid lifts.
        for j, sign in product(range(n + 1), (-1, 1)):
            s = F(1, 2) + F((2 * j - 1) ** 2, 2)
            v = [rational(s + 1), rational(sign), (F(), F(2 * j - 1)), rational(s - 1)]
            assert lorentz(v, v) == rational(-1)
            tight = []
            for i, normal in enumerate(normals):
                p = lorentz(v, normal)
                # This family has products either rational or pure sqrt2.
                assert p[0] <= 0 and p[1] <= 0
                if p == ZERO: tight.append(i)
            assert len(tight) == 3
            counts['B3' if j in (0, n) else 'A3'] += 1
        assert counts == ({'B3': 4} if n == 1 else {'B3': 4, 'A3': 2 * n - 2})
        rows.append(dict(caps=n, facets=n + 4, finite_vertices=2 * n + 2,
                         ideal_vertices=1, compact_quadrilateral_caps=n,
                         finite_corner_types=dict(counts), each_cap_weight=8,
                         cusp_factor_ranks=[1, 1], cusp_link_vertices=4))
    return dict(exact_samples=rows,
        all_rank_formula='All integer N>=1: x in[-1/2,1/2], y in[-1/sqrt2,(2N-1)/sqrt2], outside N unit hemispheres centered(0,sqrt2*k).',
        bottom_height_squared_lower_bound='1/4',
        cusp_link='A1 affine x A1 affine, rectangular link, independent of N',
        conclusion='N compact quadrilateral facets share the nonincident x-factor; their weights are8, so this is not an unpaid-face counterexample.')


def run():
    vertices, factors = pyramid()
    faces = face_ledger(vertices, factors)
    strip = strip_checks()
    paths = [f'scripts/{TAG}.py', f'scripts/{OLD}.py', f'results/{OLD}.json',
             f'results/{OLD}.md', f'results/{OLD}_review.json',
             'scripts/fresh_compact_dimension_bound_20261015.py',
             'results/fresh_compact_dimension_bound_20261015.json',
             'results/fresh_small_cusp_defects_20261017.md',
             'results/fresh_cusp_factor_payment_framework_20261017.md',
             'results/fresh_cusp_spectator_role_payments_20261017.md']
    return dict(status='PASS_EXACT_NONINCIDENT_CARRIER_INCIDENCE_AND_SCOPED_COUNTEREXAMPLES',
        sources={p: sha(p) for p in paths}, ambient_dimension=17,
        vertices=vertices, small_apex_factors=factors, face_ledger=faces,
        unbounded_all_carrier_family=strip,
        theorem_inputs=['Known finite-volume H17 pyramid and its complete cone, freshly recomputed here.',
                        'Vinberg finite weights, adjacency/component<=7/invariant end segments<=6; component<=4 weight2.',
                        'Fine affine D16 distance<=3 and invariant end-segment rule.'],
        scope='Complete exact two-face incidence and all1296 quadrilateral weights in one known H17 chamber. No individual E8 triangle weights, new dimension theorem, or global unpaid-face upper bound. The separate all-carrier family is paid and does not refute a defect-sensitive multiplicity theorem.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    out = json.loads(json.dumps(run()))
    target = ROOT / 'results' / f'{TAG}.json'
    if args.verify:
        assert out == json.loads(target.read_text())
    else:
        target.write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
    print(out['status'])
    print('82 vertices;3456 two-faces;1296 quad weights;1048 unpaid,987 compact;6 strip-family checks')
