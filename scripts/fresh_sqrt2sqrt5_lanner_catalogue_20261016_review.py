#!/usr/bin/env python3
"""Independent F25 triangle and local-weight review, standard library only.

Uses the opposite quadratic tower, linear-span ground-field degrees, and
connected minimal-support generation instead of the primary algorithms.
No imports from the producer, project arithmetic, or other review scripts.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
from math import comb, gcd, isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt2sqrt5_lanner_catalogue_20261016_review'
PRIMARY = 'results/fresh_sqrt2sqrt5_lanner_catalogue_20261016.json'


def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def sign(x):
    return (x > 0) - (x < 0)


def ratroot(x):
    if x < 0:
        return None
    a, b = isqrt(x.numerator), isqrt(x.denominator)
    return Q(a, b) if a*a == x.numerator and b*b == x.denominator else None


def add5(x, y):
    return tuple(a+b for a, b in zip(x, y))


def mul5(x, y):
    return (x[0]*y[0] + 5*x[1]*y[1], x[0]*y[1] + x[1]*y[0])


def scale5(x, c):
    return tuple(c*a for a in x)


def inv5(x):
    norm = x[0]*x[0] - 5*x[1]*x[1]
    assert norm
    return (x[0]/norm, -x[1]/norm)


def sign5(x):
    a, b = x
    if not a:
        return sign(b)
    if not b or sign(a) == sign(b):
        return sign(a)
    return sign(a)*sign(a*a-5*b*b)


def root5(x):
    """Solve (a+b sqrt5)^2=x by rational trace and norm, exhaustively."""
    a, b = x
    if not b:
        p = ratroot(a)
        if p is not None:
            return (p, Q(0))
        q = ratroot(a/5)
        return None if q is None else (Q(0), q)
    normroot = ratroot(a*a-5*b*b)
    if normroot is None:
        return None
    for signed_norm in (normroot, -normroot):
        p = ratroot((a+signed_norm)/2)
        if p:
            z = (p, b/(2*p))
            assert mul5(z, z) == x
            return z if sign5(z) > 0 else scale5(z, -1)
    return None


class K:
    # Internal order: 1,sqrt5,sqrt2,sqrt10; unlike the primary tower.
    def __init__(self, a=0, b=0, c=0, d=0):
        self.v = tuple(map(Q, (a, b, c, d)))

    @staticmethod
    def of(x):
        return x if isinstance(x, K) else K(x)

    def __add__(self, other):
        return K(*(a+b for a, b in zip(self.v, self.of(other).v)))

    __radd__ = __add__

    def __neg__(self):
        return K(*(-a for a in self.v))

    def __sub__(self, other):
        return self + -self.of(other)

    def __rsub__(self, other):
        return self.of(other) + -self

    def __mul__(self, other):
        a, b = self.v[:2], self.v[2:]
        c, d = self.of(other).v[:2], self.of(other).v[2:]
        return K(*add5(mul5(a, c), scale5(mul5(b, d), 2)),
                 *add5(mul5(a, d), mul5(b, c)))

    __rmul__ = __mul__

    def inverse(self):
        # Product of the other three conjugates divided by rational norm.
        conjugates = [self.emb(s, t) for s, t in ((-1, 1), (1, -1), (-1, -1))]
        numerator = conjugates[0]*conjugates[1]*conjugates[2]
        norm = self*numerator
        assert norm.v[1:] == (0, 0, 0) and norm.v[0]
        return K(*(a/norm.v[0] for a in numerator.v))

    def __truediv__(self, other):
        return self*self.of(other).inverse()

    def __eq__(self, other):
        return self.v == self.of(other).v

    def emb(self, s, t):
        # s changes sqrt2, t changes sqrt5.
        a, b, c, d = self.v
        return K(a, t*b, s*c, s*t*d)

    def sign(self):
        a, b = self.v[:2], self.v[2:]
        sa, sb = sign5(a), sign5(b)
        if not sa:
            return sb
        if not sb or sa == sb:
            return sa
        return sa*sign5(add5(mul5(a, a), scale5(mul5(b, b), -2)))

    def root(self):
        a, b = self.v[:2], self.v[2:]
        if b == (0, 0):
            p = root5(a)
            if p is not None:
                return K(*p)
            q = root5(scale5(a, Q(1, 2)))
            return None if q is None else K(0, 0, *q)
        h = root5(add5(mul5(a, a), scale5(mul5(b, b), -2)))
        if h is None:
            return None
        for hh in (h, scale5(h, -1)):
            p = root5(scale5(add5(a, hh), Q(1, 2)))
            if p is None or p == (0, 0):
                continue
            q = mul5(b, scale5(inv5(p), Q(1, 2)))
            z = K(*p, *q)
            assert z*z == self
            return z if z.sign() > 0 else -z
        return None

    def data(self):
        a, b, c, d = self.v
        return list(map(str, (a, c, b, d)))

    def integral(self):
        a, b, c, d = self.v
        return all(x.denominator == 1 for x in (a-b, 2*b, c-d, 2*d))


E = ((1, 1), (-1, 1), (1, -1), (-1, -1))
S2, S5 = K(0, 0, 1), K(0, 1)
GOLD = (1+S5)/2
D = {2: K(0), 3: K(1), 4: K(2), 5: GOLD+1,
     6: K(3), 8: 2+S2, 10: 2+GOLD}


def rank(rows):
    rows = [list(row) for row in rows]
    r = 0
    for j in range(4):
        pivot = next((i for i in range(r, len(rows)) if rows[i][j]), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        p = rows[r][j]
        rows[r] = [a/p for a in rows[r]]
        for i in range(r+1, len(rows)):
            p = rows[i][j]
            rows[i] = [a-p*b for a, b in zip(rows[i], rows[r])]
        r += 1
    return r


def field_degree(generators):
    # A finite-dimensional rational algebra in a field is already a field.
    basis = [K(1)]
    while True:
        old = len(basis)
        for x in list(basis):
            for y in generators:
                z = x*y
                if rank([a.v for a in basis]+[z.v]) > len(basis):
                    basis.append(z)
        if len(basis) == old:
            return old


def totient(n):
    return sum(gcd(n, j) == 1 for j in range(1, n+1))


def alphabet():
    # If phi(m)|8, each prime divisor has p-1 dividing8; 9 is not prime.
    divisors = [d for d in range(1, 9) if 8 % d == 0]
    primes = [d+1 for d in divisors if all((d+1) % p for p in range(2, isqrt(d+1)+1))]
    assert primes == [2, 3, 5]
    powers = []
    for p in primes:
        choices = [1]
        q = p
        while 8 % (q-q//p) == 0:
            choices.append(q)
            q *= p
        powers.append(choices)
    candidates = sorted(a*b*c for a, b, c in product(*powers)
                        if a*b*c >= 2 and 8 % totient(a*b*c) == 0)
    assert candidates == [2, 3, 4, 5, 6, 8, 10, 12, 15, 16, 20, 24, 30]
    cyclic_witnesses = []
    for m in (15, 16, 20, 30):
        units = [a for a in range(1, m) if gcd(a, m) == 1]
        witness = next(a for a in units if pow(a, 2, m) not in (1, m-1))
        assert pow(witness, 4, m) in (1, m-1)
        cyclic_witnesses.append(dict(modulus=m, order4_witness=witness))
    assert K(3).root() is None  # excludes12 and24, both containing sqrt3.
    for value in D.values():
        assert value.integral()
        assert all(value.emb(*e).sign() >= 0 and (4-value.emb(*e)).sign() > 0 for e in E)
    return dict(candidates=candidates, complete_labels=list(D), cyclic_witnesses=cyclic_witnesses)


def catalogue(primary):
    rows = []
    accepted = []
    for labels in combinations_with_replacement(D, 3):
        ds = [D[m] for m in labels]
        rec = dict(labels=list(labels))
        if sum((Q(1, m) for m in labels), Q()) >= 1:
            rec['status'] = 'NOT_COMPACT_HYPERBOLIC_TRIANGLE'
        else:
            p = (ds[0]*ds[1]*ds[2]).root()
            if p is None:
                rec['status'] = 'CYCLIC_PRODUCT_NOT_IN_F'
            else:
                assert p.sign() >= 0 and p.integral()
                degree = field_degree(ds+[p])
                determinant = 8-2*sum(ds, K())-2*p
                signs = [determinant.emb(*e).sign() for e in E]
                assert signs[0] == -1
                rec.update(cycle=p.data(), determinant=determinant.data(), signs=signs, ground_degree=degree)
                if degree != 4:
                    rec['status'] = 'PROPER_GROUND_SUBFIELD'
                elif signs[1:] != [1, 1, 1]:
                    rec['status'] = 'NONPOSITIVE_CONJUGATE_DETERMINANT'
                else:
                    rec['status'] = 'ARITHMETIC_LANNER_TRIANGLE_FULL_F'
                    accepted.append(list(labels))
        rows.append(rec)
    assert len(rows) == 84 and accepted == [[2, 5, 8], [4, 5, 5]]
    for independent, saved in zip(rows, primary['triangles']['records']):
        assert independent['labels'] == saved['labels']
        assert independent['status'] == saved['status']
        if 'cycle' in independent:
            assert independent['cycle'] == saved['cyclic_product']
            assert independent['determinant'] == saved['determinant']
            assert independent['signs'] == saved['determinant_signs']
    counts = dict(Counter(row['status'] for row in rows))
    assert counts == primary['triangles']['counts']
    return dict(records=rows, accepted=accepted, counts=counts)


def projection_costs(primary):
    out = []
    for labels in ((2, 5, 8), (4, 5, 5)):
        for opposite in sorted(set(labels)):
            contacts = list(labels)
            contacts.remove(opposite)
            cycle = (D[opposite]*D[contacts[0]]*D[contacts[1]]).root()
            assert cycle is not None
            # For unit normals, inverse of [[1,-c],[-c,1]].
            numerator = (D[contacts[0]]+D[contacts[1]]+cycle)/4
            denominator = 1-D[opposite]/4
            cost = numerator/denominator
            margin = cost-2
            assert margin.sign() == (-1 if opposite == 2 else 1)
            out.append(dict(triangle=list(labels), opposite_label=opposite,
                            cost=cost.data(), margin=margin.data(), sign=margin.sign()))
    assert len(out) == len(primary['quadrilateral_costs']) == 5
    for row, saved in zip(out, primary['quadrilateral_costs']):
        assert row['cost'] == saved['cost'] and row['sign'] == saved['sign']
    # At the exceptional white-centered258, its8 bond plus the compulsory
    # white bridge has determinant <=1-D8/4-1/4<0, already a proper corner.
    isolation_margin = 1-D[8]/4-K(Q(1, 4))
    assert isolation_margin.sign() < 0
    return dict(cases=out, white_centered258_corner_determinant_upper=isolation_margin.data())


def edge(a, b):
    return (min(a, b), max(a, b))


def components(vertices, edges):
    parent = {v: v for v in vertices}
    def find(v):
        while parent[v] != v:
            v = parent[v]
        return v
    cycles = 0
    for a, b in sorted(edges):
        x, y = find(a), find(b)
        if x == y:
            cycles += 1
        else:
            parent[x] = y
    return len({find(v) for v in vertices}), cycles


def all_subsets(items):
    for k in range(len(items)+1):
        for subset in combinations(items, k):
            yield set(subset)


def supports(primary):
    shapes = [((0, 2), (1, 3)), ((0, 2), (1, 3, 4)),
              ((0, 2, 4), (1, 3)), ((0, 2, 4), (1, 3, 5)),
              ((0, 2, 4), (1, 3, 4))]
    survivors, examined = set(), 0
    for left, right in shapes:
        vertices = sorted(set(left+right))
        whites = sorted(set(vertices)-{0, 1, 2, 3})
        le, re = set(combinations(left, 2)), set(combinations(right, 2))
        cross = sorted(set(combinations(vertices, 2))-le-re)
        lchoices = [e for e in all_subsets(sorted(le)) if components(left, e)[0] == 1]
        rchoices = [e for e in all_subsets(sorted(re)) if components(right, e)[0] == 1]
        for ls, rs, cs in product(lchoices, rchoices, list(all_subsets(cross))):
            examined += 1
            es = ls | rs | cs
            if not set(left) & set(right) and not cs:
                continue  # orthogonal indefinite blocks are impossible.
            valid = True
            for a, b in ((0, 1), (1, 2), (2, 3), (3, 0)):
                corner = set(whites+[a, b])
                ce = {e for e in es if set(e) <= corner}
                if components(corner, ce)[1]:
                    valid = False
                if edge(a, b) in ce or any(edge(a, w) in ce and edge(b, w) in ce for w in whites):
                    valid = False
            if valid:
                assert left == (0, 2, 4) and right == (1, 3, 5)
                assert cs == {(4, 5)}
                survivors.add(tuple(sorted(es)))
    assert examined == 8976 and len(survivors) == 16
    converted = [[list(e) for e in es] for es in sorted(survivors)]
    assert converted == primary['unmarked_quadrilateral_supports']['surviving_supports']
    return dict(connected_support_candidates=examined, survivors=converted,
                all_white_scope='The obstructing principal six-node matrix suffices regardless of other whites.')


def counting():
    checks = []
    for r in range(1, 41):
        for kind in ('path', 'fork') if r >= 4 else ('path',):
            es = {(i, i+1) for i in range(r-1)}
            if kind == 'fork':
                es.remove((r-2, r-1))
                es.add((r-3, r-1))
            degrees = [sum(v in e for e in es) for v in range(r)]
            total = len(es)+sum(comb(d, 2) for d in degrees)
            assert total == (0 if r == 1 else 2*r-3 if kind == 'path' else 2*r-2)
            assert total <= 2*r-2
            checks.append(dict(rank=r, kind=kind, total=total))
    # Symmetric h-vector: f2/f0 is the weighted mean of these coefficients.
    paired = [Q(comb(21-i, 2)+comb(i, 2), 2) for i in range(22)]
    assert min(paired) == 50 and paired[0] == 105
    assert [p-50 for p in paired[:11]] == [Q((10-i)*(11-i), 2) for i in range(11)]
    assert 5*min(paired)-comb(21, 2) == 40 == 2*21-2
    # All-rank parity formula, checked coefficientwise as a quadratic.
    gaps = {}
    for parity, base in (('odd', 21), ('even', 22)):
        values = []
        for u in range(3):
            n = base+2*u
            lower = Q((n-1)*(n-5), 8) if n % 2 else Q(n*(n-6), 8)
            values.append(lower-(2*n-2))
        c = values[0]
        a = (values[2]-2*values[1]+values[0])/2
        b = values[1]-a-c
        assert (c, b, a) == ((Q(0), Q(5), Q(1, 2)) if parity == 'odd'
                            else (Q(2), Q(11, 2), Q(1, 2)))
        gaps[parity] = list(map(str, (c, b, a)))
    return dict(graph_checks=checks, dimension21_paired_coefficients=list(map(str, paired)),
                strict_lower=40, vertex_upper=40, parity_gap_coefficients=gaps,
                all_rank_proof='Finite connected diagrams are paths or single-trivalent trees; degree identity holds at every rank.')


def run():
    primary = json.loads((ROOT / PRIMARY).read_text())
    for path, digest in primary['sources'].items():
        assert sha(path) == digest, path
    assert S2*S2 == 2 and S5*S5 == 5 and GOLD*GOLD == GOLD+1
    # Exact radical/inverse checks spanning both towers and physical signs.
    for v in product((-1, 0, 1), repeat=4):
        z = K(*v)
        if z == 0:
            continue
        assert z*z.inverse() == 1
        root = (z*z).root()
        assert root is not None and root*root == z*z and root.sign() > 0
    labels = alphabet()
    assert labels['complete_labels'] == primary['alphabet']['finite_labels']
    out = dict(status='PASS_INDEPENDENT_F25_COMPACT_BOUND20_REVIEW',
               sources={p: sha(p) for p in [f'scripts/{TAG}.py', PRIMARY,
                       'results/fresh_sqrt2sqrt5_lanner_catalogue_20261016.md']},
               alphabet=labels, triangles=catalogue(primary),
               costs=projection_costs(primary), supports=supports(primary),
               counting=counting(), compact_dimension_upper_bound_candidate=20,
               external_inputs=['Nikulin0708.3991 section2: hyperbolic subgraphs retain the ground field',
                                'Nikulin0708.3991 section3.1 formula4: arithmetic Lanner order>=4 fields',
                                'Compact acute face correspondence and absence of affine diagrams',
                                'Finite Coxeter diagram classification and simple-polytope h-vector symmetry/nonnegativity'],
               scope='Independent exact finite-angle/triangle catalogue, all position costs, minimal-pair support reduction, and weight threshold; all-rank and all-white arguments are documented in the review report. External validation and priority remain open.')
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    out = run()
    path = ROOT / 'results' / f'{TAG}.json'
    if args.verify:
        assert out == json.loads(path.read_text())
        print('PASS fresh saved independent F25 replay')
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], out['triangles']['accepted'], out['supports']['connected_support_candidates'])
