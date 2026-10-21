#!/usr/bin/env python3
"""Independent Fraction/graph reconstruction of the cubic-seven <=20 candidate.

No producer imports, SymPy, or shared determinant engine. See the report for
the hostile audit of geometric completeness and the entire-D5 assertion.
"""
import argparse
from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_compact_global_bound_20261015_review'
PRIMARY = 'fresh_compact_global_bound_20261015'
NV = 9
NAMES = ['c', 'p', 'q', 'r', 'alpha', 'beta', 'gamma', 'delta', 'epsilon']
ZERO = (0,)*NV


class Poly:
    def __init__(self, value=0):
        if isinstance(value, Poly):
            self.d = value.d.copy()
        elif isinstance(value, dict):
            self.d = {m: F(c) for m, c in value.items() if c}
        else:
            self.d = {} if not value else {ZERO: F(value)}
    def __add__(self, other):
        out = self.d.copy()
        for m, c in Poly(other).d.items():
            out[m] = out.get(m, 0)+c
        return Poly(out)
    __radd__ = __add__
    def __neg__(self):
        return Poly({m: -c for m, c in self.d.items()})
    def __sub__(self, other):
        return self+-Poly(other)
    def __rsub__(self, other):
        return Poly(other)+-self
    def __mul__(self, other):
        out = {}
        for m, c in self.d.items():
            for n, e in Poly(other).d.items():
                k = tuple(a+b for a, b in zip(m, n))
                out[k] = out.get(k, 0)+c*e
        return Poly(out)
    __rmul__ = __mul__
    def __pow__(self, n):
        out = Poly(1)
        for _ in range(n):
            out *= self
        return out
    def __bool__(self):
        return bool(self.d)
    def __eq__(self, other):
        return self.d == Poly(other).d
    def substitute(self, assignments):
        out = {}
        for m, c in self.d.items():
            m = list(m)
            for i, v in assignments.items():
                c *= F(v)**m[i]
                m[i] = 0
            m = tuple(m)
            out[m] = out.get(m, 0)+c
        return Poly(out)
    def coefficient_c(self, n):
        out = {}
        for m, c in self.d.items():
            if m[0] == n:
                out[(0,)+m[1:]] = c
        return Poly(out)
    def reduce_c(self):
        out = self.d.copy()
        while any(m[0] >= 3 for m in out):
            m = max((m for m in out if m[0] >= 3), key=lambda x: x[0])
            a = out.pop(m)
            # c^3 = c^2/2 + c/2 -1/8.
            for power, factor in [(2, F(1, 2)), (1, F(1, 2)), (0, F(-1, 8))]:
                k = (m[0]-3+power,)+m[1:]
                out[k] = out.get(k, 0)+a*factor
                if not out[k]:
                    del out[k]
        return Poly(out)
    def encoded(self):
        return [[list(m), str(c)] for m, c in sorted(self.d.items())]


def variable(i):
    m = [0]*NV
    m[i] = 1
    return Poly({tuple(m): 1})


def determinant(matrix):
    """Sparse Laplace recursion; all operations are in independent Poly."""
    size = len(matrix)
    def recurse(row, columns):
        if row == size:
            return Poly(1)
        total = Poly()
        for j, column in enumerate(columns):
            if matrix[row][column]:
                total += (-1)**j*matrix[row][column]*recurse(row+1, columns[:j]+columns[j+1:])
        return total
    return recurse(0, tuple(range(size)))


def symbolic_checks():
    c, p, q, r, a, b, g, d, e = [variable(i) for i in range(NV)]
    P = 8*c**3-4*c**2-4*c+1
    t = 1-c*c
    Hdet = p*q*r-d*d*r-e*e*q
    G = [[Poly(int(i == j)) for j in range(7)] for i in range(7)]
    for i, j, z in [(0, 1, c), (2, 3, c), (1, 4, a), (2, 5, b),
                    (3, 6, g), (4, 5, d), (4, 6, e)]:
        G[i][j] = G[j][i] = -z
    G[4][4], G[5][5], G[6][6] = p, q, r
    raw = determinant(G)
    formula = ((p*t-a*a)*(q*r*t-q*g*g-r*b*b+b*b*g*g)
               -d*d*(r*t-g*g)*t-e*e*(q*t-b*b)*t-2*d*e*b*g*c*t)
    assert raw == formula
    remainder = raw.reduce_c()
    assert remainder.coefficient_c(1) == F(1, 8)*(Hdet-8*d*e*b*g)
    equality = {1: 1, 2: 1, 3: 1, 5: F(1, 2), 6: F(1, 2),
                7: F(1, 2), 8: F(1, 2)}
    assert remainder.substitute(equality) == F(1, 64)*(4*a*a-1)*(16*c*c-9)
    assert raw.substitute(equality | {4: F(1, 2)}) == F(1, 64)*(4*c+3)*P
    G6 = [[Poly(int(i == j)) for j in range(6)] for i in range(6)]
    for i, j, z in [(0, 1, c), (2, 3, c), (1, 4, a), (3, 5, b), (4, 5, d)]:
        G6[i][j] = G6[j][i] = -z
    G6[4][4], G6[5][5] = p, q
    raw6 = determinant(G6)
    assert raw6 == (p*q-d*d)*t*t-(p*b*b+q*a*a)*t+a*a*b*b
    assert raw6.coefficient_c(4) == p*q-d*d
    assert set(m[0] for m in raw6.d) <= {0, 2, 4}
    squared = 64*c**6-80*c**4+24*c*c-1
    assert squared.reduce_c() == 0
    # Rational-root test for the primitive cubic in c.
    assert all(P.substitute({0: F(sign, den)}) != 0
               for den in (1, 2, 4, 8) for sign in (-1, 1))
    assert P.substitute({0: F(9, 10)}) == F(-1, 125)
    assert 24*F(4, 5)**2-8*F(4, 5)-4 > 0
    assert 48*F(4, 5)-8 > 0
    return dict(variables=NAMES, raw_seven=raw.encoded(), remainder=remainder.encoded(),
                coefficient_c=remainder.coefficient_c(1).encoded(),
                raw_six=raw6.encoded(), six_leading_coefficient=raw6.coefficient_c(4).encoded(),
                equality_remainder=remainder.substitute(equality).encoded(),
                exceptional_factorization=raw.substitute(equality | {4: F(1, 2)}).encoded())


def dist(adj, subset, a, b):
    allowed = sum(1 << x for x in subset)
    seen = 1 << a
    todo = deque([(a, 0)])
    while todo:
        v, d = todo.popleft()
        if v == b:
            return d
        bits = adj[v] & allowed & ~seen
        while bits:
            bit = bits & -bits
            bits -= bit
            w = bit.bit_length()-1
            seen |= bit
            todo.append((w, d+1))
    return None


def connected(adj, vertices):
    return all(dist(adj, vertices, vertices[0], x) is not None for x in vertices)


def is_forest(adj, vertices):
    # A component with |E| >= |V| is cyclic, computed by component search.
    left = set(vertices)
    while left:
        start = next(iter(left))
        component = {x for x in left if dist(adj, left, start, x) is not None}
        count_edges = sum(bool(adj[a] & 1 << b) for a, b in combinations(component, 2))
        if count_edges != len(component)-1:
            return False
        left -= component
    return True


def graph_reduction():
    """Enumerate whole binary graphs, independently of producer shape loops."""
    passing, graphs_examined = set(), 0
    for left, right in [([0, 2], [1, 3]), ([0, 2], [1, 3, 4]),
                        ([0, 2, 4], [1, 3]), ([0, 2, 4], [1, 3, 5]),
                        ([0, 2, 4], [1, 3, 4])]:
        vertices = sorted(set(left+right))
        pairs = list(combinations(vertices, 2))
        white = sorted(set(vertices)-{0, 1, 2, 3})
        for mask in range(1 << len(pairs)):
            graphs_examined += 1
            adj = [0]*6
            edges = []
            for k, (a, b) in enumerate(pairs):
                if mask >> k & 1:
                    adj[a] |= 1 << b
                    adj[b] |= 1 << a
                    edges.append((a, b))
            if not connected(adj, left) or not connected(adj, right):
                continue
            if not set(left) & set(right):
                if not any(adj[a] & 1 << b for a in left for b in right):
                    continue
            corners = [white+[a, b] for a, b in product([0, 2], [1, 3])]
            if any(not is_forest(adj, corner) or
                   (dist(adj, corner, corner[-2], corner[-1]) or 99) <= 2
                   for corner in corners):
                continue
            choices = []
            for nodes in (left, right):
                if len(nodes) == 2:
                    choices.append([None])
                else:
                    possible = []
                    for a, b in combinations(nodes, 2):
                        if not adj[a] & 1 << b:
                            continue
                        if all(not {a, b} <= set(corner) or
                               all(sum(bool(adj[x] & 1 << y) for y in corner) == 1 for x in (a, b))
                               for corner in corners):
                            possible.append((a, b))
                    choices.append(possible)
            if not all(choices):
                continue
            assert left == [0, 2, 4] and right == [1, 3, 5]
            assert choices == [[(0, 2)], [(1, 3)]]
            passing.add(tuple(edges))
    assert graphs_examined == 35904
    assert len(passing) == 9
    return dict(all_binary_graphs=graphs_examined, remaining_supports=sorted(passing))


def graph_family(kind, n):
    if kind in ('A', 'B', 'F', 'I'):
        return [(j, j+1) for j in range(n-1)]
    if kind == 'D':
        return [(j, j+1) for j in range(n-3)]+[(n-3, n-2), (n-3, n-1)]
    edges, k = [], 1
    for length in (n-4, 2, 1):
        prev = 0
        for _ in range(length):
            edges.append((prev, k))
            prev, k = k, k+1
    return edges


def graph_weights():
    rows = []
    kinds = [(k, n) for k in ('A', 'B') for n in range(1, 41)]
    kinds += [('D', n) for n in range(4, 41)]+[('E', n) for n in (6, 7, 8)]
    kinds += [('F', 4), ('I', 2)]
    for kind, n in kinds:
        edges = graph_family(kind, n)
        adj = [0]*n
        for a, b in edges:
            adj[a] |= 1 << b
            adj[b] |= 1 << a
        counts = Counter(dist(adj, range(n), a, b) for a, b in combinations(range(n), 2))
        base = counts[1]+counts[2]
        if kind in ('A', 'B', 'F', 'I'):
            assert base == max(n-1, 0)+max(n-2, 0)
        else:
            assert base == 2*n-2
        halves = []
        if kind == 'D' and n == 5:
            halves = [(a, b) for a, b in combinations(range(n), 2)
                      if adj[a].bit_count() == adj[b].bit_count() == 1
                      and dist(adj, range(n), a, b) == 3]
            assert len(halves) == 2
        value = F(base+int(kind == 'I'))+F(len(halves), 2)
        assert 2*n-value >= (1 if kind == 'D' and n == 5 else 2)
        rows.append(dict(type=kind, rank=n, weight=str(value), deficit=str(2*n-value),
                         half_pairs=halves))
    # Rank>=6 cannot be a lone exceptional D5. Two or more components
    # contribute total deficit>=2; a single other component also does.
    return rows


def rational_ldl(matrix):
    a = [list(map(F, row)) for row in matrix]
    pivots = []
    for k in range(len(a)):
        pivot = a[k][k]
        assert pivot > 0
        pivots.append(pivot)
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a)):
                a[i][j] -= a[i][k]*a[k][j]/pivot
    return pivots


def corner_checks():
    edges = [(0, 1), (2, 3), (1, 4), (2, 5), (3, 6), (4, 5), (4, 6)]
    rows = []
    for x, y in [(0, 2), (2, 1), (1, 3), (3, 0)]:
        nodes = [4, 5, 6, x, y]
        selected = [(a, b) for a, b in edges if a in nodes and b in nodes]
        assert (0, 1) not in selected and (2, 3) not in selected
        matrix = [[F(i == j) for j in range(5)] for i in range(5)]
        adj = [0]*7
        for a, b in selected:
            matrix[nodes.index(a)][nodes.index(b)] = -F(1, 2)
            matrix[nodes.index(b)][nodes.index(a)] = -F(1, 2)
            adj[a] |= 1 << b
            adj[b] |= 1 << a
        pivots = rational_ldl(matrix)
        degree = sorted(adj[a].bit_count() for a in nodes)
        typ = 'D5' if degree == [1, 1, 1, 2, 3] else 'A4+A1'
        value = F(1, 2) if typ == 'D5' else F(0)
        if value:
            assert adj[x].bit_count() == adj[y].bit_count() == 1
            assert dist(adj, nodes, x, y) == 3
        rows.append(dict(marked=[x, y], type=typ, pivots=list(map(str, pivots)), weight=str(value)))
    assert sum(F(x['weight']) for x in rows) == 1
    assert Counter(x['type'] for x in rows) == {'D5': 2, 'A4+A1': 2}
    return rows


def run():
    primary = json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    assert primary['status'] == 'PASS_COMPACT_CUBIC_SEVEN_DIMENSION20_PROOF_CANDIDATE'
    reduced = graph_reduction()
    assert json.loads(json.dumps(reduced['remaining_supports'])) == primary['initial_graph_reduction']['remaining_supports']
    polys = symbolic_checks()
    weights = graph_weights()
    corners = corner_checks()
    assert F(1, 4)/(1-F(9, 10)**2) == F(25, 19)
    assert F(1, 2)/(1-F(9, 10)) == 5
    assert 4*F(6, 19)-1 == F(5, 19) > 0
    assert 1-2*F(1, 4) == 8*F(1, 2)**4 == F(1, 2)
    h = [comb(21-i, 2)+comb(i, 2)-100 for i in range(11)]
    assert h == [(10-i)*(11-i) for i in range(11)]
    assert 5*50-comb(21, 2) == 40 == 2*21-2
    # Exact polynomial coefficient identities for both unbounded parities.
    odd = [F(20*16, 8)-40, F(2*20+2*16, 8)-4, F(4, 8)]
    even = [F(22*16, 8)-42, F(2*22+2*16, 8)-4, F(4, 8)]
    assert odd == [0, 5, F(1, 2)] and even == [2, F(11, 2), F(1, 2)]
    paths = [f'scripts/{TAG}.py']
    for stem in [PRIMARY, 'fresh_arithmetic_highdim_20261015',
                 'fresh_arithmetic_highdim_20261015_review']:
        paths += [f'scripts/{stem}.py', f'results/{stem}.json']
        document = json.loads((ROOT/f'results/{stem}.json').read_text())
        for p, expected in document.get('sources', {}).items():
            assert sha256((ROOT/p).read_bytes()).hexdigest() == expected, p
    return dict(status='PASS_INDEPENDENT_FRACTION_AND_EXHAUSTIVE_QUAD_SCOPE_REVIEW_CUBIC_SEVEN_BOUND20',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                scope=dict(compact_only=True, arithmetic_ground_field='Q(cos(2pi/7))', bound=20,
                           independent_polynomial_arithmetic=True, producer_imports=False,
                           inherited_field_classification=True, novelty_claim=False),
                binary_graph_reduction=reduced, independent_polynomials=polys,
                spherical_graphs=weights, four_exceptional_corners=corners,
                scalar=dict(negative_plane_margin='5/19', determinant_sandwich='1/2',
                            h21_coefficients=h, odd_gap=list(map(str, odd)), even_gap=list(map(str, even))),
                source_proof_audit=[
                    'Every opposite minimal indefinite block is order2/3 and contains both opposite black nodes.',
                    'The graph reduction includes shared whites, every cross-bond subset and every possible distinguished high bond.',
                    'All further white contacts form a complete bipartite forest; each black endpoint has at most one white neighbor.',
                    'All off-core white components attach only once and have rational positive Schur corrections.',
                    'The complete two-face star has at most rank|W|+3, so its reduced6/7node matrix is singular.',
                    'Arithmetic cyclic products make L rational; equality of the determinant sandwich removes all external white branches.',
                    'Exhaustion of black-white contacts then makes each weighted D5 an entire vertex component.',
                    'Strict h-vector inequality and the finite-component deficit prove every dimension>=21 impossible.'],
                limitations=['This is an independent mathematical/arithmetic review, not a formal proof-assistant verification.',
                             'The earlier cubic field classification and published compact geometry are explicit inherited dependencies.'])


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = json.loads(json.dumps(run()))
    path = ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert json.loads(path.read_text()) == result
        print('PASS saved replay')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
        print(result['status'])
    print(json.dumps(dict(binary_graphs=result['binary_graph_reduction']['all_binary_graphs'],
                          remaining_graphs=len(result['binary_graph_reduction']['remaining_supports']),
                          spherical_graphs=len(result['spherical_graphs']), bound=20)))
