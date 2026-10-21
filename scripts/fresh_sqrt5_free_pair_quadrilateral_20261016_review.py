#!/usr/bin/env python3
"""Independent Fraction / Q(sqrt5) review; no producer imports."""
import argparse
from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_free_pair_quadrilateral_20261016_review'
PRIMARY = 'results/fresh_sqrt5_free_pair_quadrilateral_20261016.json'
PRIOR = 'results/fresh_sqrt5_universal_self_edges_20261016.json'


@dataclass(frozen=True)
class K:
    a: F = F(0)
    b: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, 'a', F(self.a))
        object.__setattr__(self, 'b', F(self.b))

    @staticmethod
    def of(v):
        return v if isinstance(v, K) else K(v)

    def __add__(self, v):
        v = K.of(v)
        return K(self.a+v.a, self.b+v.b)

    __radd__ = __add__

    def __neg__(self):
        return K(-self.a, -self.b)

    def __sub__(self, v):
        return self+-K.of(v)

    def __rsub__(self, v):
        return K.of(v)+-self

    def __mul__(self, v):
        v = K.of(v)
        return K(self.a*v.a+5*self.b*v.b, self.a*v.b+self.b*v.a)

    __rmul__ = __mul__

    def __truediv__(self, v):
        v = K.of(v)
        norm = v.a*v.a-5*v.b*v.b
        assert norm
        return self*K(v.a/norm, -v.b/norm)

    def __rtruediv__(self, v):
        return K.of(v)/self

    def __bool__(self):
        return bool(self.a or self.b)

    def conj(self):
        return K(self.a, -self.b)

    def sign(self):
        sign = lambda q: (q > 0)-(q < 0)
        if not self.a:
            return sign(self.b)
        if not self.b or sign(self.a) == sign(self.b):
            return sign(self.a)
        return sign(self.a)*sign(self.a*self.a-5*self.b*self.b)

    def saved(self):
        return [str(self.a), str(self.b)]


def read_phi(v):
    a, b = map(F, v)
    return K(a+b/2, b/2)


def rational_inverse(A):
    """Independent symmetric LDL and triangular solves, no row reduction."""
    n = len(A)
    L = [[F(i == j) for j in range(n)] for i in range(n)]
    D = []
    for j in range(n):
        d = A[j][j]-sum(L[j][k]**2*D[k] for k in range(j))
        assert d > 0
        D.append(d)
        for i in range(j+1, n):
            L[i][j] = (A[i][j]-sum(L[i][k]*L[j][k]*D[k] for k in range(j)))/d
    inv = [[F(0) for _ in range(n)] for _ in range(n)]
    for col in range(n):
        y = []
        for i in range(n):
            y.append(F(i == col)-sum(L[i][k]*y[k] for k in range(i)))
        z = [y[i]/D[i] for i in range(n)]
        x = [F(0)]*n
        for i in reversed(range(n)):
            x[i] = z[i]-sum(L[k][i]*x[k] for k in range(i+1, n))
        for i in range(n):
            inv[i][col] = x[i]
    assert all(sum(A[i][k]*inv[k][j] for k in range(n)) == (i == j)
               for i in range(n) for j in range(n))
    return inv, D


def sha(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def run():
    data = json.loads((ROOT/PRIMARY).read_text())
    prior = json.loads((ROOT/PRIOR).read_text())
    assert data['status'] == 'PASS_EXACT_NO_QUADRILATERAL_FOR_210_FIXED_UNIFORM_EDGE_PAIRS'
    assert prior['exceptional_vertex_type'] == ['D5', 'E8', 'E8']
    A = [[F(2*(i == j)) for j in range(21)] for i in range(21)]
    # D5 has arms1,1,2; E8 has arms1,2,4. Literal producer numbering.
    edges = [(0, 1), (1, 2), (2, 3), (1, 4)]
    E = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)]
    edges += [(a+off, b+off) for off in (5, 13) for a, b in E]
    for i, j in edges:
        A[i][j] = A[j][i] = F(-1)
    assert [[read_phi(v) for v in row] for row in prior['exceptional_full_finite_gram']] == [[K(v) for v in row] for row in A]
    C, D = rational_inverse(A)
    assert all(v >= 0 for row in C for v in row)
    phi = K(F(1, 2), F(1, 2))
    lam = K(1)
    for _ in range(8):
        lam = lam*phi
    assert lam == read_phi(prior['uniform_rule']['lambda_value'])
    T = 2*lam
    assert T == read_phi(data['fixed_tau'])
    p = [2*T-T*T*C[i][i] for i in range(21)]
    v = [[K(i == k)-T*C[k][i] for k in range(21)] for i in range(21)]
    contact = [[sum(A[k][l]*v[i][l] for l in range(21)) for k in range(21)] for i in range(21)]
    for i in range(21):
        assert p[i].sign() < 0 < p[i].conj().sign()
        assert contact[i] == [K(A[k][i])-(T if k == i else 0) for k in range(21)]
        assert sum(v[i][k]*contact[i][k] for k in range(21))+p[i] == K(2)
    assert len(data['pairs']) == 210
    minimum = None
    rows = []
    sectors = 0
    intercomponent = 0
    for (i, j), rec in zip(combinations(range(21), 2), data['pairs']):
        assert rec['ports'] == [i, j]
        b = sum(v[i][k]*contact[j][k] for k in range(21))
        assert b == K(A[i][j])+T*T*C[i][j] and b.sign() >= 0
        assert b == read_phi(rec['b'])
        assert [p[i], p[j]] == list(map(read_phi, rec['new_root_residual_diagonals']))
        product = p[i]*p[j]
        endpoint = []
        for x, key in ((-2, 'determinant_at_minus2'), (0, 'determinant_at_zero')):
            d = product-(K(x)-b)*(K(x)-b)
            assert d.sign() > 0 and d == read_phi(rec[key])
            endpoint.append(d.saved())
            if minimum is None or (d-minimum[0]).sign() < 0:
                minimum = d, i, j, x
        coefficients = [product-b*b, 2*b, K(-1)]
        assert coefficients == list(map(read_phi, rec['determinant_coefficients_constant_linear_quadratic']))
        assert product == read_phi(rec['rank_closure_roots']['radicand'])
        assert rec['rank_on_entire_finite_corner_interval'] == 23
        assert rec['signature_on_entire_finite_corner_interval'] == dict(positive=21, negative=2, zero=0)
        assert (product-(b+2)*(b+2)).sign() > 0
        assert (product-b*b).sign() > 0
        orthogonal = (i < 5) != (j < 5) or (5 <= i < 13) != (5 <= j < 13)
        assert orthogonal == (C[i][j] == 0) == rec['orthogonal_components']
        intercomponent += orthogonal
        local = []
        assert len(rec['real_partial_star_endpoint_sectors']) == 2
        for (a, c), sec in zip(((i, j), (j, i)), rec['real_partial_star_endpoint_sectors']):
            den = T*C[a][a]-1
            aa, ac = -p[a], -p[c]
            margin = ac*den*den-aa*T*T*C[a][c]*C[a][c]
            sq = aa*C[a][a]/(den*den)-1
            assert den.sign() > 0 and aa.sign() > 0 and ac.sign() > 0
            assert margin.sign() > 0 and sq == -1/(den*den)
            assert sq.sign() < 0
            assert sec['released_root'] == a and sec['other_new_root'] == c
            assert den == read_phi(sec['denominator'])
            assert sq == read_phi(sec['endpoint_square'])
            assert margin == read_phi(sec['strict_other_wall_squared_margin'])
            # Dot products of vi with the old roots are zero except port a,
            # where -sqrt(aa)/den<0. New-a product is exactly zero because
            # -sqrt(aa)*(1+(1-T*Caa)/den)=0. New-c is strictly negative by
            # the saved positive-square comparison; Cac>=0 fixes its sign.
            assert 1+(1-T*C[a][a])/den == K()
            local.append(dict(released_root=a, other_new_root=c,
                              square=sq.saved(), margin=margin.saved()))
            sectors += 1
        rows.append(dict(ports=[i, j], endpoint_determinants=endpoint, sectors=local))
    assert sectors == 420 and intercomponent == 144
    assert minimum[0] == read_phi(data['smallest_endpoint_determinant']['value'])
    assert list(minimum[1:3]) == data['smallest_endpoint_determinant']['ports']
    assert minimum[3] == data['smallest_endpoint_determinant']['x']
    assert minimum[0] == read_phi(['20169647', '32635176'])
    seen, queue, references = set(), [PRIMARY, PRIOR], 0
    while queue:
        fn = queue.pop()
        if fn in seen:
            continue
        seen.add(fn)
        for pth, h in json.loads((ROOT/fn).read_text()).get('sources', {}).items():
            assert sha(pth) == h, (fn, pth)
            references += 1
            if pth.endswith('.json'):
                queue.append(pth)
    paths = [f'scripts/{TAG}.py', PRIMARY, PRIOR,
             'scripts/fresh_sqrt5_free_pair_quadrilateral_20261016.py']
    return dict(status='PASS_INDEPENDENT_FRACTION_FREE_PAIR_QUADRILATERAL_AND_REAL_SECTORS',
                sources={p: sha(p) for p in paths}, field_basis=['1', 'sqrt5'],
                old_graph_edges=[list(e) for e in edges], old_positive_pivots=list(map(str, D)),
                cases=210, orthogonal_component_pairs=144, same_component_pairs=66,
                interval_endpoint_checks=420, strict_real_endpoint_sector_checks=420,
                smallest_endpoint_determinant=dict(value_sqrt5_basis=minimum[0].saved(),
                                                  ports=list(minimum[1:3]), x=minimum[3]),
                recursive_manifests_checked=len(seen), recursive_references_checked=references,
                pairs=rows,
                scope='Independent literal D5+2E8 Gram, rational LDL inverse, sqrt5 field signs, all420 determinant endpoints and420 real endpoint-sector checks. Continuous concavity/rank and real sector construction source-reviewed. Only the fixed phi^8 edges are excluded from quadrilateral closure; no arbitrary edge-length, arithmetic partial-star, finite-volume or dimension21 exclusion.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    out = run()
    path = ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'])
