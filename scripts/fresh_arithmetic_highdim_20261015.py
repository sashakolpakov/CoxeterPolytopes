#!/usr/bin/env python3
"""Exact checks for a compact Lanner-order<=3 dimension-28 proof candidate.

The geometric charging proof and published theorem inputs are in the matching
report. This program checks its finite graph cases, all-rank formula anchors,
scalar consequence, and a separate complete cubic-seven triangle catalogue.
No unbounded Coxeter-label cutoff is used. No producer arithmetic is imported.
"""
import argparse
import json
from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
from math import comb, gcd, isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = "fresh_arithmetic_highdim_20261015"
VINBERG = "https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt"
NIKULIN = "https://arxiv.org/pdf/0708.3991"


class K:
    """Q[t]/(t^3+t^2-2t-1); signs by rational root isolation."""
    def __init__(self, x=0):
        self.v = x.v if isinstance(x, K) else (
            tuple(map(F, x)) if isinstance(x, (tuple, list)) else (F(x), F(0), F(0)))
    def __add__(self, y):
        return K(tuple(a+b for a, b in zip(self.v, K(y).v)))
    __radd__ = __add__
    def __neg__(self):
        return K(tuple(-a for a in self.v))
    def __sub__(self, y):
        return self + -K(y)
    def __rsub__(self, y):
        return K(y) + -self
    def __mul__(self, y):
        q = [F(0)]*5
        for i, a in enumerate(self.v):
            for j, b in enumerate(K(y).v):
                q[i+j] += a*b
        for i in (4, 3):
            q[i-3] += q[i]
            q[i-2] += 2*q[i]
            q[i-1] -= q[i]
        return K(q[:3])
    __rmul__ = __mul__
    def __bool__(self):
        return any(self.v)
    def __eq__(self, y):
        return self.v == K(y).v
    def encoded(self):
        return [str(x) for x in self.v]
    def sign(self, embedding=0):
        if not self:
            return 0
        a, b = ROOT_INTERVALS[embedding]
        assert cubic(a)*cubic(b) < 0
        for _ in range(300):
            lo = hi = self.v[2]
            for x in reversed(self.v[:2]):
                z = [lo*a, lo*b, hi*a, hi*b]
                lo, hi = min(z)+x, max(z)+x
            if lo > 0:
                return 1
            if hi < 0:
                return -1
            mid = (a+b)/2
            if cubic(a)*cubic(mid) < 0:
                b = mid
            else:
                a = mid
        raise ArithmeticError(self.encoded())


ROOT_INTERVALS = [(F(124, 100), F(125, 100)),
                  (F(-45, 100), F(-44, 100)),
                  (F(-181, 100), F(-180, 100))]


def cubic(x):
    return x**3+x*x-2*x-1


def det3(a):
    return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
            - a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
            + a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))


def norm(x):
    basis = [K((1, 0, 0)), K((0, 1, 0)), K((0, 0, 1))]
    cols = [(x*y).v for y in basis]
    return det3([[cols[j][i] for j in range(3)] for i in range(3)])


def distance(vertices, edges, start, target):
    allowed = set(vertices)
    adj = {v: [] for v in vertices}
    for a, b in edges:
        if a in allowed and b in allowed:
            adj[a].append(b)
            adj[b].append(a)
    todo = deque([(start, 0)])
    seen = {start}
    while todo:
        a, d = todo.popleft()
        if a == target:
            return d
        for b in adj[a]:
            if b not in seen:
                seen.add(b)
                todo.append((b, d+1))
    return None


def connected_graphs(vertices):
    pairs = list(combinations(vertices, 2))
    for mask in range(1 << len(pairs)):
        edges = [p for i, p in enumerate(pairs) if mask >> i & 1]
        if all(distance(vertices, edges, vertices[0], v) is not None for v in vertices):
            yield edges


def local_graph_checks():
    triangles = list(connected_graphs([0, 1, 2]))
    assert len(triangles) == 4 and min(map(len, triangles)) == 2
    rows = []
    # Two minimal indefinite diagrams have black sets {0,2} and {1,3}.
    # Each can use at most one white. If disjoint, one cross edge suffices;
    # if they overlap, their only possible common node is white.
    for wl, wr, overlap in [(0, 0, False), (0, 1, False),
                            (1, 0, False), (1, 1, False), (1, 1, True)]:
        left = [0, 2] + ([4] if wl else [])
        right = [1, 3] + ([4 if overlap else 5] if wr else [])
        white = sorted((set(left) | set(right)) - {0, 1, 2, 3})
        for el in connected_graphs(left):
            for er in connected_graphs(right):
                bridges = [None] if overlap else list(product(left, right))
                for bridge in bridges:
                    edges = el+er+([bridge] if bridge is not None else [])
                    ds = []
                    for a, b in product([0, 2], [1, 3]):
                        d = distance(white+[a, b], edges, a, b)
                        if d is not None:
                            ds.append(d)
                    assert ds and min(ds) <= 3
                    rows.append(dict(left=left, right=right, left_edges=el,
                                     right_edges=er, bridge=bridge,
                                     shortest_actual_corner_path=min(ds)))
    return dict(triangle_shapes=triangles, quadrilateral_cases=rows,
                max_quad_path=max(x['shortest_actual_corner_path'] for x in rows))


def graph(kind, rank):
    if kind in ('A', 'B', 'F', 'H', 'I'):
        return [(i, i+1) for i in range(rank-1)]
    if kind == 'D':
        return [(i, i+1) for i in range(rank-3)] + [(rank-3, rank-2), (rank-3, rank-1)]
    assert kind == 'E' and rank in (6, 7, 8)
    edges, vertex = [], 1
    for length in (rank-4, 2, 1):
        prev = 0
        for _ in range(length):
            edges.append((prev, vertex))
            prev, vertex = vertex, vertex+1
    assert vertex == rank
    return edges


def graph_record(kind, rank):
    edges = graph(kind, rank)
    counts = Counter(distance(list(range(rank)), edges, a, b)
                     for a, b in combinations(range(rank), 2))
    q = sum(counts[d] for d in (1, 2, 3))
    if kind in ('A', 'B', 'F', 'H', 'I'):
        expected = sum(max(rank-d, 0) for d in (1, 2, 3))
    elif kind == 'D':
        expected = 6 if rank == 4 else 3*rank-5
    else:
        expected = 3*rank-4
    assert q == expected
    assert 3*rank-q >= (3 if rank == 1 else 4)
    return dict(type=kind, rank=rank, edges=edges,
                distance_counts=[counts[d] for d in (1, 2, 3)],
                weight=q, deficit=3*rank-q)


def phi(n):
    return sum(gcd(k, n) == 1 for k in range(1, n+1))


def discriminant_cubic(a, b, c):
    return a*a*b*b-4*b**3-4*a**3*c-27*c*c+18*a*b*c


def triangles():
    t = K((0, 1, 0))
    s = t*t+t-1
    assert s*s == t+2 and s*s*s-s*s-2*s+1 == 0
    assert s.sign() == 1
    squares = {2: K(), 3: K(1), 4: K(2), 6: K(3), 7: t+2, 14: t*t+t+1}
    assert all(q.sign(e) >= 0 and (4-q).sign(e) > 0
               for q in squares.values() for e in range(3))
    # If phi(m)|6, every prime p|m has p-1|6: p in {2,3,7}.
    # Exponents are bounded by p^(e-1)|6. Thus this is a complete
    # analytic finite reduction, not a numerical label cutoff.
    candidates = sorted({2**a*3**b*7**c for a in range(3) for b in range(3)
                         for c in range(2)
                         if 2**a*3**b*7**c >= 3 and phi(2**a*3**b*7**c) in (2, 6)})
    assert candidates == [3, 4, 6, 7, 9, 14, 18]
    assert discriminant_cubic(1, -2, -1) == 49
    assert discriminant_cubic(0, -3, 1) == 81
    assert gcd(49, 81) == 1
    # Both cubics are irreducible by the rational-root test. If their
    # fields coincided, its discriminant would divide both 49 and81,
    # impossible for a nontrivial number field (Minkowski).
    assert all(cubic(x) != 0 and x**3-3*x+1 != 0 for x in (-1, 1))
    radical_norms = []
    for mask in range(1, 8):
        q = K(1)
        for j, m in enumerate((4, 6, 14)):
            if mask >> j & 1:
                q *= squares[m]
        nq = norm(q)
        assert nq.denominator == 1 and nq > 0
        assert isqrt(nq.numerator)**2 != nq.numerator
        radical_norms.append(dict(mask=mask, norm=str(nq)))
    rows, accepted = [], []
    for labels in combinations_with_replacement(squares, 3):
        row = dict(labels=labels)
        if sum(F(1, m) for m in labels) >= 1:
            row['reason'] = 'not_hyperbolic'
        elif not any(m in (7, 14) for m in labels):
            row['reason'] = 'no_cubic_field_generator'
        elif 2 not in labels and any(labels.count(m) % 2 for m in (4, 6, 14)):
            row['reason'] = 'three_cycle_outside_K'
        else:
            cycle = K(0 if 2 in labels else 1)
            if 2 not in labels:
                for m in (4, 6, 14):
                    for _ in range(labels.count(m)//2):
                        cycle *= squares[m]
                for _ in range(labels.count(7)):
                    cycle *= s
            determinant = 8-2*sum((squares[m] for m in labels), K())-2*cycle
            signs = [determinant.sign(e) for e in range(3)]
            row.update(determinant=determinant.encoded(), conjugate_signs=signs)
            row['reason'] = 'accepted' if signs == [-1, 1, 1] else 'conjugate_signature'
            if row['reason'] == 'accepted':
                assert determinant.v[1:] != (F(0), F(0))
                accepted.append(row)
        rows.append(row)
    assert len(rows) == 56
    assert [tuple(r['labels']) for r in accepted] == [
        (2, 3, 7), (2, 3, 14), (2, 4, 7), (2, 7, 7),
        (2, 7, 14), (3, 3, 7), (7, 7, 7)]
    return dict(allowed_finite_labels=list(squares), phi_candidates=candidates,
                excluded_other_cubic_labels=[9, 18], cubic_discriminants=[49, 81],
                basis='1,t,t^2; t=2cos(2pi/7)', polynomial=[-1, -2, 1, 1],
                root_intervals=[[str(a), str(b)] for a, b in ROOT_INTERVALS],
                squareclass_norms=radical_norms, all_cases=rows,
                counts=dict(Counter(r['reason'] for r in rows)), accepted=accepted)


def run():
    records = [graph_record(k, n) for k in ('A', 'B') for n in range(1, 33)]
    records += [graph_record('D', n) for n in range(4, 33)]
    records += [graph_record('E', n) for n in (6, 7, 8)]
    records += [graph_record('F', 4), graph_record('H', 3),
                graph_record('H', 4), graph_record('I', 2)]
    # Increasing a classical long arm by one beyond these base ranks adds
    # precisely one pair at each distance1,2,3: the displayed all-rank
    # formulas follow by induction, not extrapolation from this range.
    for k, start in [('A', 3), ('B', 3), ('D', 5)]:
        for n in range(start, 32):
            a, b = graph_record(k, n), graph_record(k, n+1)
            assert [y-x for x, y in zip(a['distance_counts'], b['distance_counts'])] == [1, 1, 1]
    # A connected rank>=2 component has deficit>=4. If every component
    # has rank1 and total rank>=2, its total deficit is at least6.
    # Thus the disconnected all-rank bound is q3<=3n-4.
    assert 2*3 >= 4
    h29 = []
    for i in range(15):
        c = comb(29-i, 2)+comb(i, 2)-2*98
        assert c == (14-i)*(15-i) and c >= 0
        h29.append(c)
    assert 5*98-comb(29, 2) == 84 > 3*29-4 == 83
    # n<8(3-4/n)+6 gives n^2-30n+32<0. At n=29+u:
    assert 29**2-30*29+32 == 3
    assert 2*29-30 == 28
    # A proper three-node path with labels>=6,>=3 is non-PD.
    assert 1-F(3, 4)-F(1, 4) == 0
    bound = [f'scripts/{TAG}.py',
             'scripts/fresh_compact_dimension_bound_20261015.py',
             'results/fresh_compact_dimension_bound_20261015.json',
             'scripts/fresh_compact_dimension_bound_20261015_review.py',
             'results/fresh_compact_dimension_bound_20261015_review.json']
    return dict(
        status='PASS_COMPACT_LANNER_ORDER_LE3_AND_CUBIC_SEVEN_DIMENSION28_CANDIDATE',
        sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in bound},
        scope=dict(compact_only=True, arithmetic_required_for_geometric_lemma=False,
                   geometric_hypothesis='Every Lanner subdiagram has order at most3.',
                   arithmetic_application='Ground field not a subfield of Q(sqrt2,sqrt5), in particular odd degree>1.',
                   high_dimension_cubic_field='Q(cos(2pi/7))', bound=28,
                   surviving_compact_dimension29_fields=['Q(sqrt2)', 'Q(sqrt5)'],
                   novelty_or_current_best_claim=False),
        published_inputs=[
            dict(url=VINBERG, location='printed pp.70--72, Proposition2',
                 statement='Angle weights: vertex total<=cn and k-face total>=5-k imply n<8c+6.'),
            dict(url=VINBERG, location='printed pp.74--77, finite diagrams and Propositions5,7--9',
                 statement='Finite-diagram classification and compact Coxeter face/Lanner correspondence.'),
            dict(url=VINBERG, location='printed p.91, Proposition17',
                 statement='Every Lanner determinant generates the arithmetic ground field.'),
            dict(url=VINBERG, location='printed p.69, Theorem2; section7',
                 statement='Compact arithmetic n>=22 ground fields are Qsqrt2,Qsqrt5,Qcos(2pi/7).'),
            dict(url=NIKULIN, location='PDF p.2, introduction, corrected arXiv version',
                 statement='Explicit n>=22 list includingQ and the three nonrational fields; also real cyclotomic degrees in section5.')],
        graph_checks=records, graph_check_count=len(records),
        all_rank_formulas=dict(path='sum_{d=1}^3 max(r-d,0)',
                               D4=6, D_ge5='3r-5', E6_E7_E8='3r-4',
                               forest_rank_ge2='<=3n-4'),
        local_face_graphs=local_graph_checks(),
        scalar=dict(Vinberg_c='3-4/n', necessary_polynomial='n^2-30n+32<0',
                    at_n_29_plus_u=[3, 28, 1],
                    h29_coefficients=h29, dimension29_lower_strict=84,
                    dimension29_upper=83),
        field_obstruction=dict(proper_three_node_path_boundary=0,
                               allowed_Lanner_order_ge4_labels=[2, 3, 4, 5],
                               determinant_field='Q(sqrt2,sqrt5)', field_degree=4,
                               excluded_cubic_degree=3),
        cubic_triangle_catalogue=triangles(),
        verification_limits=[
            'Geometric charging and all-rank induction are proved in the report; finite checks do not replace them.',
            'No automated external proof verification or priority search completeness is claimed.',
            'The cubic triangle catalogue is auxiliary and not used to derive dimension28.'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    path = ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert json.loads(path.read_text()) == json.loads(json.dumps(result)), 'saved certificate mismatch'
        print('PASS saved replay')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
        print(result['status'])
    print(json.dumps(dict(graphs=result['graph_check_count'],
                          quad_shapes=len(result['local_face_graphs']['quadrilateral_cases']),
                          cubic_triangles=len(result['cubic_triangle_catalogue']['accepted']),
                          bound=28)))
