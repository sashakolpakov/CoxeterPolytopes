#!/usr/bin/env python3
"""Exact local identities for the compact cubic-seven dimension-20 candidate.

The proof of exhaustive geometric reduction is in the matching report.  This
checker verifies the finite graph reductions, two determinant identities,
their cubic-field coefficient consequences, the exceptional marked corners,
and the all-rank weight formulas.  Published arithmetic/face theorems remain
explicit dependencies.  No producer arithmetic is imported.
"""
import argparse
import json
from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from math import comb
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_compact_global_bound_20261015'
VINBERG = 'https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt'


def edge(a, b):
    return tuple(sorted((a, b)))


def distance(vertices, edges, a, b):
    vertices = set(vertices)
    adj = {v: [] for v in vertices}
    for x, y in edges:
        if x in vertices and y in vertices:
            adj[x].append(y)
            adj[y].append(x)
    todo, seen = deque([(a, 0)]), {a}
    while todo:
        x, d = todo.popleft()
        if x == b:
            return d
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                todo.append((y, d+1))
    return None


def forest(vertices, edges):
    parent = {x: x for x in vertices}
    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x
    for a, b in edges:
        if a not in parent or b not in parent:
            continue
        x, y = find(a), find(b)
        if x == y:
            return False
        parent[x] = y
    return True


def connected_graphs(vertices):
    pairs = list(combinations(vertices, 2))
    for mask in range(1 << len(pairs)):
        es = {edge(*p) for i, p in enumerate(pairs) if mask >> i & 1}
        if all(distance(vertices, es, vertices[0], v) is not None for v in vertices):
            yield es


def initial_reduction():
    """All binary supports, including arbitrary cross bonds, not a sample.

    Each order-three minimal indefinite diagram has a distinguished high
    finite bond; any additional high bonds only strengthen the finite-corner
    isolation test. Order-two diagrams are opposite dotted pairs. The test
    deliberately relaxes the other finite-diagram conditions to forests.
    """
    examined, valid = 0, set()
    for wl, wr, shared in [(0, 0, False), (0, 1, False),
                          (1, 0, False), (1, 1, False), (1, 1, True)]:
        left = [0, 2] + ([4] if wl else [])
        right = [1, 3] + ([4 if shared else 5] if wr else [])
        vertices = sorted(set(left+right))
        white = sorted(set(vertices)-{0, 1, 2, 3})
        for el in connected_graphs(left):
            for er in connected_graphs(right):
                internal = {edge(*p) for p in combinations(left, 2)}
                internal |= {edge(*p) for p in combinations(right, 2)}
                cross = sorted(set(combinations(vertices, 2))-internal)
                for hl in (sorted(el) if wl else [None]):
                    for hr in (sorted(er) if wr else [None]):
                        high = {x for x in (hl, hr) if x is not None}
                        for mask in range(1 << len(cross)):
                            bridges = {p for i, p in enumerate(cross) if mask >> i & 1}
                            if not shared and not bridges:
                                continue  # orthogonal indefinite blocks
                            examined += 1
                            es = el | er | bridges
                            good = True
                            for a, b in product((0, 2), (1, 3)):
                                corner = set(white+[a, b])
                                ce = {p for p in es if set(p) <= corner}
                                if not forest(corner, ce):
                                    good = False
                                    break
                                deg = Counter(x for p in ce for x in p)
                                if any(deg[x] != 1 or deg[y] != 1
                                       for x, y in high if {x, y} <= corner):
                                    good = False
                                    break
                                d = distance(corner, ce, a, b)
                                if d is not None and d <= 2:
                                    good = False
                                    break
                            if good:
                                assert wl == wr == 1 and not shared
                                assert hl == (0, 2) and hr == (1, 3)
                                assert bridges == {(4, 5)}
                                valid.add(tuple(sorted(es)))
    assert len(valid) == 9  # three nonempty attachment patterns on each side
    return dict(examined_relaxed_cases=examined,
                remaining_supports=[list(x) for x in sorted(valid)],
                conclusion='two disjoint high opposite bonds, joined only by the white-white bond')


def polynomial_record(expr, variables):
    return dict(variables=[str(x) for x in variables],
                terms=[[list(m), str(c)] for m, c in s.Poly(s.expand(expr), *variables).terms()])


def determinant_checks():
    c, p, q, r = s.symbols('c p q r')
    alpha, beta, gamma, delta, epsilon = s.symbols('alpha beta gamma delta epsilon')
    A, B, C, d, e, L = s.symbols('A B C d e L')
    t, P = 1-c*c, 8*c**3-4*c*c-4*c+1
    # Nodes a,c,b,d,u,v,w; only displayed off-diagonal products are nonzero.
    G = s.eye(7)
    for i, j, z in [(0, 1, c), (2, 3, c), (1, 4, alpha),
                    (2, 5, beta), (3, 6, gamma), (4, 5, delta), (4, 6, epsilon)]:
        G[i, j] = G[j, i] = -z
    G[4, 4], G[5, 5], G[6, 6] = p, q, r
    N = ((p*t-A)*(q*r*t-q*C-r*B+B*C)
         - d*(r*t-C)*t-e*(q*t-B)*t-2*L*c*t)
    substitution = {A: alpha**2, B: beta**2, C: gamma**2,
                    d: delta**2, e: epsilon**2, L: delta*epsilon*beta*gamma}
    raw = G.det(method='domain-ge')
    assert s.expand(raw-N.subs(substitution)) == 0
    rem = s.rem(N, P, c)
    detH = p*q*r-d*r-e*q
    assert s.expand(rem.coeff(c, 1)-(detH-8*L)/8) == 0
    eq = {p: 1, q: 1, r: 1, B: s.Rational(1, 4), C: s.Rational(1, 4),
          d: s.Rational(1, 4), e: s.Rational(1, 4), L: s.Rational(1, 16)}
    reduced = s.factor(rem.subs(eq))
    assert s.expand(reduced-(4*A-1)*(16*c*c-9)/64) == 0
    exceptional = s.factor(N.subs(eq | {A: s.Rational(1, 4)}))
    assert exceptional == (4*c+3)*P/64
    # One white neighbor on either high pair: a genuine degree-two
    # polynomial in z=c^2, whose leading coefficient is detH2>0.
    G6 = s.eye(6)
    for i, j, z in [(0, 1, c), (2, 3, c), (1, 4, alpha),
                    (3, 5, beta), (4, 5, delta)]:
        G6[i, j] = G6[j, i] = -z
    G6[4, 4], G6[5, 5] = p, q
    N6 = (p*q-d)*t*t-(p*B+q*A)*t+A*B
    assert s.expand(G6.det(method='domain-ge')-N6.subs(substitution)) == 0
    z = s.symbols('z')
    squared_minpoly = 64*z**3-80*z*z+24*z-1
    assert s.rem(squared_minpoly.subs(z, c*c), P, c) == 0
    assert s.Poly(P, c).is_irreducible and s.Poly(squared_minpoly, z).is_irreducible
    # Positive-definite H <= Hmax in Loewner order; direct partial
    # derivatives are its positive principal 2x2 minors.
    assert s.diff(detH, p) == q*r
    assert s.diff(detH, q) == p*r-e
    assert s.diff(detH, r) == p*q-d
    return dict(
        cubic=polynomial_record(P, (c,)), squared_cubic=polynomial_record(squared_minpoly, (z,)),
        seven_determinant=polynomial_record(N, (c, p, q, r, A, B, C, d, e, L)),
        seven_remainder=polynomial_record(rem, (c, p, q, r, A, B, C, d, e, L)),
        linear_coefficient=str(s.factor(rem.coeff(c, 1))),
        determinant_H=str(detH), equality_remainder=str(reduced),
        exceptional_factorization=str(exceptional),
        six_determinant=polynomial_record(N6, (c, p, q, A, B, d)),
        direct_matrix=G.tolist(), direct_matrix_six=G6.tolist())


def scalar_checks():
    x = s.symbols('x')
    P = 8*x**3-4*x*x-4*x+1
    assert P.subs(x, s.Rational(9, 10)) == -s.Rational(1, 125)
    assert s.diff(P, x).subs(x, s.Rational(4, 5)) > 0
    assert s.diff(P, x, 2).subs(x, s.Rational(4, 5)) > 0
    assert F(3, 4) > F(4, 5)**2
    # c>.9 gives single7 cost>25/19 and dense7/single14 cost>5.
    assert 1/(4*(1-F(9, 10)**2)) == F(25, 19)
    assert 1/(2*(1-F(9, 10))) == 5
    assert (5-1)*(F(25, 19)-1)-1 == F(5, 19) > 0
    # Schur equality sandwich detH <=1-d-e<=1/2<=8L.
    assert 1-2*F(1, 4) == 8*F(1, 2)**4 == F(1, 2)
    # Odd n=21+2u and even n=22+2u: exact h-vector lower minus upper.
    u = s.symbols('u')
    odd = s.expand(((20+2*u)*(16+2*u))/8-(40+4*u))
    even = s.expand(((22+2*u)*(16+2*u))/8-(42+4*u))
    assert odd == u*u/2+5*u and even == u*u/2+11*u/2+2
    h21 = []
    for i in range(11):
        v = comb(21-i, 2)+comb(i, 2)-100
        assert v == (10-i)*(11-i)
        h21.append(v)
    assert 5*50-comb(21, 2) == 40 == 2*21-2
    return dict(c_lower='9/10', P_at_lower='-1/125',
                single7_cost_lower='25/19', dense7_or_single14_cost_lower='5',
                negative_plane_determinant_lower='5/19', H_sandwich_equality='1/2',
                odd_gap=str(odd), even_gap=str(even), h21_coefficients=h21,
                strict_dimension21_lower=40, dimension21_vertex_upper=40)


def diagram(kind, rank):
    if kind in ('A', 'B', 'F', 'I'):
        return [(i, i+1) for i in range(rank-1)]
    if kind == 'D':
        return [(i, i+1) for i in range(rank-3)] + [(rank-3, rank-2), (rank-3, rank-1)]
    assert kind == 'E' and rank in (6, 7, 8)
    es, nxt = [], 1
    for length in (rank-4, 2, 1):
        prev = 0
        for _ in range(length):
            es.append((prev, nxt))
            prev, nxt = nxt, nxt+1
    return es


def graph_weights():
    rows = []
    types = [(k, n) for k in ('A', 'B') for n in range(1, 34)]
    types += [('D', n) for n in range(4, 34)] + [('E', n) for n in (6, 7, 8)]
    types += [('F', 4), ('I', 2)]
    for kind, rank in types:
        es = diagram(kind, rank)
        ds = {p: distance(range(rank), es, *p) for p in combinations(range(rank), 2)}
        base = sum(d <= 2 for d in ds.values())
        if kind in ('A', 'B', 'F', 'I'):
            assert base == max(rank-1, 0)+max(rank-2, 0)
        else:
            assert base == 2*rank-2
        half_pairs = []
        if (kind, rank) == ('D', 5):
            deg = Counter(x for e in es for x in e)
            half_pairs = [p for p, d in ds.items() if d == 3 and all(deg[x] == 1 for x in p)]
            assert len(half_pairs) == 2
        bonus = 1 if kind == 'I' else 0
        w = F(base+bonus)+F(len(half_pairs), 2)
        deficit = 2*rank-w
        assert deficit >= (1 if (kind, rank) == ('D', 5) else 2)
        rows.append(dict(type=kind, rank=rank, edges=es, distance_two_weight=base,
                         half_pairs=half_pairs, high_I2_bonus=bonus, weight=str(w), deficit=str(deficit)))
    # Structural formulas, not extrapolation: lengthening a long classical
    # arm adds one distance1 pair and one distance2 pair.
    return rows


def exceptional_corners():
    # a,c,b,d,u,v,w =0..6. Opposite high bonds do not occur at corners.
    es = {(0, 1), (2, 3), (1, 4), (2, 5), (3, 6), (4, 5), (4, 6)}
    high = {(0, 1), (2, 3)}
    rows = []
    for a, b in [(0, 2), (2, 1), (1, 3), (3, 0)]:
        nodes = [4, 5, 6, a, b]
        ce = {p for p in es if set(p) <= set(nodes)}
        assert not high & ce and forest(nodes, ce)
        G = s.eye(5)
        for x, y in ce:
            i, j = nodes.index(x), nodes.index(y)
            G[i, j] = G[j, i] = -s.Rational(1, 2)
        minors = [G[:i, :i].det() for i in range(1, 6)]
        assert all(x > 0 for x in minors)
        deg = Counter(x for e in ce for x in e)
        d = distance(nodes, ce, a, b)
        typ = 'D5' if sorted(deg[x] for x in nodes) == [1, 1, 1, 2, 3] else 'A4+A1'
        weight = F(1, 2) if typ == 'D5' and d == 3 and deg[a] == deg[b] == 1 else F(0)
        rows.append(dict(marked=[a, b], nodes=nodes, edges=sorted(ce), type=typ,
                         marked_distance=d, weight=str(weight), leading_minors=list(map(str, minors))))
    assert Counter(x['type'] for x in rows) == {'D5': 2, 'A4+A1': 2}
    assert sum(F(x['weight']) for x in rows) == 1
    return rows


def run():
    arithmetic = json.loads((ROOT/'results/fresh_arithmetic_highdim_20261015.json').read_text())
    catalogue = arithmetic['cubic_triangle_catalogue']
    assert catalogue['allowed_finite_labels'] == [2, 3, 4, 6, 7, 14]
    expected_triangles = [(2, 3, 7), (2, 3, 14), (2, 4, 7),
                          (2, 7, 7), (2, 7, 14), (3, 3, 7), (7, 7, 7)]
    assert [tuple(x['labels']) for x in catalogue['accepted']] == expected_triangles
    assert all(7 in x or 14 in x for x in expected_triangles)
    assert arithmetic['field_obstruction']['excluded_cubic_degree'] == 3
    files = [f'scripts/{TAG}.py', 'scripts/fresh_arithmetic_highdim_20261015.py',
             'results/fresh_arithmetic_highdim_20261015.json',
             'scripts/fresh_compact_dimension_bound_20261015.py',
             'results/fresh_compact_dimension_bound_20261015.json']
    return dict(status='PASS_COMPACT_CUBIC_SEVEN_DIMENSION20_PROOF_CANDIDATE',
                scope=dict(compact=True, arithmetic=True, ground_field='Q(cos(2pi/7))',
                           proposed_upper_dimension=20, unconditional_compact_bound_changed=False,
                           novelty_or_priority_claim=False),
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
                theorem_dependencies=[VINBERG,
                    'Cubic-seven finite labels and Lanner triangle catalogue: frozen arithmetic-highdim package.',
                    'Compact Coxeter face/elliptic correspondence, absence of affine subdiagrams.',
                    'Finite spherical diagram classification; Vinberg arithmetic cyclic-product criterion.',
                    'Nonnegative symmetric h-vector of the dual simplicial polytope.'],
                inherited_field_guards=dict(allowed_finite_labels=catalogue['allowed_finite_labels'],
                                            Lanner_triangles=expected_triangles,
                                            maximum_Lanner_order=3),
                initial_graph_reduction=initial_reduction(),
                determinant_identities=determinant_checks(), scalar=scalar_checks(),
                graph_weights=graph_weights(), exceptional_corners=exceptional_corners(),
                verification_scope='Exact finite checks and identities; the full geometric/weighting proof is in the report.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = json.loads(json.dumps(run(), default=str))
    path = ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert json.loads(path.read_text()) == result, 'saved certificate mismatch'
        print('PASS saved replay')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
        print(result['status'])
    print(json.dumps(dict(local_graph_cases=result['initial_graph_reduction']['examined_relaxed_cases'],
                          graph_records=len(result['graph_weights']), bound=20)))
