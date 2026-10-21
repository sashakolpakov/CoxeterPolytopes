#!/usr/bin/env python3
"""Independent exact Gram-inverse audit of the fixed 192/228 cap banks.

Does not import either producer. Uses Golden instead of their Fraction-pair
arithmetic, a literal ten-dimensional Gram inverse instead of their feature
formula, a different exhaustive assignment generator, and set-based clique
enumeration. No floating point arithmetic or numerical geometry is used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import combinations, product
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden as G

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data/enumerations/fresh_h9_split_leaf_caps_20261013_exact.json'
DOTTED = ROOT / 'data/enumerations/fresh_h9_split_leaf_caps_20261013_dotted_exact.json'
OUT = ROOT / 'results/fresh_h9_split_leaf_caps_independent_20261013.json'
Z, O, P = G(0), G(1), G(0, 1)


def parse(v):
    return G.from_json(v) if isinstance(v, dict) else G(*map(Q, v))


def enc(x):
    return [str(x.a), str(x.b)]


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), Z)


def mv(m, x):
    return tuple(dot(r, x) for r in m)


def inverse(m):
    n = len(m)
    a = [list(r) + [G(i == j) for j in range(n)] for i, r in enumerate(m)]
    for k in range(n):
        pivot = next(i for i in range(k, n) if a[i][k])
        a[k], a[pivot] = a[pivot], a[k]
        v = a[k][k]
        a[k] = [x / v for x in a[k]]
        for i in range(n):
            if i != k:
                v = a[i][k]
                a[i] = [x - v * y for x, y in zip(a[i], a[k])]
    return tuple(tuple(r[n:]) for r in a)


@dataclass(frozen=True)
class E:
    """Q(phi,sqrt(phi)); phi' < 0 proves this is a genuine extension."""
    a: G = Z
    b: G = Z

    def __add__(self, other):
        return E(self.a + other.a, self.b + other.b)

    def __neg__(self):
        return E(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if not self.b:
            return E(self.a * other.a, self.a * other.b)
        if not other.b:
            return E(self.a * other.a, self.b * other.a)
        return E(self.a * other.a + P * self.b * other.b,
                 self.a * other.b + self.b * other.a)

    def scale(self, other):
        return E(self.a * other, self.b * other)

    def sign(self):
        sa, sb = self.a.sign(), self.b.sign()
        if not sb:
            return sa
        if not sa or sa == sb:
            return sb
        delta = (self.a * self.a - P * self.b * self.b).sign()
        assert delta != 0  # sqrt(phi) is not in Q(phi).
        return sa if delta > 0 else sb

    def encode(self):
        return [self.a.to_json(), self.b.to_json()]


def edot(a, b):
    return sum((x * y for x, y in zip(a, b)), E())


def emv(m, x):
    return tuple(sum((y.scale(v) for v, y in zip(r, x)), E()) for r in m)


def compatible(x):
    # If b != 0, the embedding phi -> phi' < 0 sends this element to
    # a nonreal number. Every conjugate of 2*cos(pi/m) is real.
    return (not x.b and x.a in (Z, -O, -P)) or (x + E(G(2))).sign() <= 0


def half_assignments():
    """Enumerate pair sums first and all four independent U values second."""
    vals = (Z, O, P)
    sums = {a + b for a, b in product(vals, repeat=2)}
    assert len(sums) == 6
    answer = set()
    for s in sums:
        for x in product(vals, repeat=4):
            y = tuple(s - v for v in x)
            if all(v in vals for v in y):
                answer.add((x, y))
    assert len(answer) == 51
    return answer


def clique_list(neighbors):
    """Independent set-based Bron--Kerbosch, reversing producer branch order."""
    found = []

    def visit(chosen, remaining, earlier):
        if not remaining and not earlier:
            found.append(tuple(sorted(chosen)))
            return
        pivot = max(remaining | earlier,
                    key=lambda i: (len(remaining & neighbors[i]), i))
        for v in sorted(remaining - neighbors[pivot], reverse=True):
            visit(chosen + [v], remaining & neighbors[v], earlier & neighbors[v])
            remaining.remove(v)
            earlier.add(v)

    visit([], set(range(len(neighbors))), set())
    return sorted(found)


def digest(x):
    return hashlib.sha256(json.dumps(x, separators=(',', ':')).encode()).hexdigest()


def audit():
    start = time.monotonic()
    saved = json.loads(BASE.read_text())
    dotted = json.loads(DOTTED.read_text())
    a, b = P + Q(1, 2), -G(5, -2) / 22
    assert b == -a / (12 * a + 1)
    eps = (1,) * 4 + (-1,) * 4
    gram = tuple(tuple(
        G(2) if i == j else Z if i // 8 == j // 8 else
        -(a + b * eps[i % 8] * eps[j % 8] + 2 * (i % 8 == j % 8))
        for j in range(16)) for i in range(16))
    assert gram == tuple(tuple(map(parse, r)) for r in saved['gram'])
    assert enc(a) == saved['a'] and enc(b) == saved['b']

    # Sum/difference invariant spaces give all sixteen eigenvalue signs.
    eigenvalues = [-8*a, -8*b] + [Z]*6 + [G(4)]*6 + [4+8*a, 4+8*b]
    signs = [x.sign() for x in eigenvalues]
    signature = [signs.count(1), signs.count(-1), signs.count(0)]
    assert signature == saved['seed_signature'] == [9, 1, 6]
    assert all((a + b * eps[i] * eps[j] + 2*(i == j) - 2).sign() > 0
               for i in range(8) for j in range(8))
    center = (O,) * 16
    gc = mv(gram, center)
    assert all(x.sign() < 0 for x in gc) and dot(center, gc).sign() < 0
    assert list(map(enc, center)) == saved['strict_center']
    assert enc(dot(center, gc)) == saved['strict_center_square']

    # All U_i plus V_0,V_4 form a ten-dimensional basis. Six exact kernel
    # identities U_i+V_i-U_j-V_j=0 give all constraints on leaf covectors.
    basis = tuple(range(8)) + (8, 12)
    small = tuple(tuple(gram[i][j] for j in basis) for i in basis)
    inv = inverse(small)
    for i in range(10):
        assert mv(small, tuple(inv[j][i] for j in range(10))) == tuple(G(j == i) for j in range(10))
    for j in (0, 4):
        for i in range(j + 1, j + 4):
            v = [Z]*16
            v[i], v[i+8], v[j], v[j+8] = O, O, -O, -O
            assert mv(gram, v) == (Z,)*16

    # Solve the literal Gram system for each compatible leaf covector.
    half = half_assignments()
    expected = {}
    for (x, y), (u, v) in product(half, repeat=2):
        magnitudes = x + u + y + v
        cov = tuple(-magnitudes[i] for i in basis)
        coeff = mv(inv, cov)
        if dot(coeff, cov) != 2:
            continue
        literal = [Z]*16
        for i, z in zip(basis, coeff):
            literal[i] = z
        assert mv(gram, literal) == tuple(-z for z in magnitudes)
        expected[magnitudes] = coeff
    assert len(expected) == 192 == saved['cap_count']
    assert saved['ordered_assignments_exhausted'] == 2601
    coeffs, covs = [], []
    seen = set()
    for row in saved['caps']:
        mags = tuple(map(parse, row['x'] + row['y']))
        assert mags not in seen and mags in expected
        seen.add(mags)
        r = tuple(map(parse, row['normal']))
        gr = mv(gram, r)
        assert gr == tuple(-z for z in mags)
        assert dot(r, gr) == 2 and sum(gr, Z).sign() < 0
        coeffs.append(tuple(E(z) for z in expected[mags]))
        covs.append(tuple(E(gr[i]) for i in basis))
    assert seen == set(expected)

    finite_edges, neighbors = [], [set() for _ in coeffs]
    for i, j in combinations(range(192), 2):
        val = edot(coeffs[i], covs[j])
        assert not val.b
        if compatible(val):
            finite_edges.append([i, j, enc(val.a)])
            neighbors[i].add(j)
            neighbors[j].add(i)
    assert finite_edges == saved['compatibility_edges']
    finite_cliques = clique_list(neighbors)
    assert len(finite_cliques) == saved['maximal_clique_count'] == 8478
    assert max(map(len, finite_cliques)) == saved['maximum_compatible_caps'] == 6
    assert finite_cliques == [tuple(c['caps']) for c in saved['all_maximal_clique_escapes']]

    witness_cov = []
    for row in saved['strict_spacelike_witnesses']:
        w = tuple(map(parse, row['coefficients']))
        gw = mv(gram, w)
        assert enc(dot(w, gw)) == row['square'] and dot(w, gw).sign() > 0
        assert sum(gw, Z).sign() < 0 and all(z.sign() < 0 for z in gw)
        witness_cov.append(tuple(E(gw[i]) for i in basis))
    assert len(witness_cov) == 36
    finite_coverage = [[edot(r, gw).sign() < 0 for r in coeffs] for gw in witness_cov]
    for row in saved['all_maximal_clique_escapes']:
        assert all(finite_coverage[row['strict_escape']][j] for j in row['caps'])
    print({'finite_bank': 'PASS', 'caps': 192, 'cliques': len(finite_cliques)}, flush=True)

    # Add the 36 declared dotted patterns and independently solve the same
    # ten-dimensional Gram system, with the target covectors now in E.
    d = E(P, P)
    assert (d - E(G(2))).sign() > 0
    meta = dotted['metadata']
    assert meta['dotted_magnitude'] == d.encode()
    patterns = [list(i+j) for i, j in product(combinations(range(4), 2), combinations(range(4, 8), 2))]
    assert patterns == meta['new_cap_patterns']
    for pattern, raw in zip(patterns, meta['new_cap_normals']):
        x = tuple(E(O) if i in pattern else d for i in range(8))
        y = tuple(d if i in pattern else E(O) for i in range(8))
        full_cov = tuple(-z for z in x+y)
        cov = tuple(full_cov[i] for i in basis)
        coeff = emv(inv, cov)
        assert edot(coeff, cov) == E(G(2))
        r = tuple(E(parse(z[0]), parse(z[1])) for z in raw)
        assert emv(gram, r) == full_cov
        assert edot(r, full_cov) == E(G(2))
        assert sum(full_cov, E()).sign() < 0
        coeffs.append(coeff)
        covs.append(cov)
    assert len(coeffs) == meta['cap_count'] == 228
    edges = [[i, j, E(parse(v)).encode()] for i, j, v in finite_edges]
    neighbors += [set() for _ in range(36)]
    for i, j in combinations(range(228), 2):
        if j < 192:
            continue
        val = edot(coeffs[i], covs[j])
        if compatible(val):
            edges.append([i, j, val.encode()])
            neighbors[i].add(j)
            neighbors[j].add(i)
    edges.sort(key=lambda e: e[:2])
    assert len(edges) == meta['compatibility_edges'] == 8460
    assert digest(edges) == meta['edge_digest']
    cliques = clique_list(neighbors)
    assert len(cliques) == meta['maximal_clique_count'] == 156990
    assert max(map(len, cliques)) == meta['maximum_caps'] == 14
    assert digest(cliques) == meta['clique_digest']
    selected = meta['selected_fourteen_caps']
    assert len(selected) == 14 and tuple(selected) in cliques
    print({'combined_graph': 'PASS', 'caps': 228, 'cliques': len(cliques)}, flush=True)

    cut_sets = []
    for row in dotted['witnesses']:
        w = tuple(G(Q(z)) for z in row['coefficients'])
        gw = mv(gram, w)
        norm, orientation = dot(w, gw), sum(gw, Z)
        assert norm.sign() > 0 and orientation.sign() < 0
        assert all(z.sign() < 0 for z in gw)
        assert E(norm).encode() == row['square']
        assert E(orientation).encode() == row['center_pairing']
        signs = [edot(c, tuple(E(gw[i]) for i in basis)).sign() for c in coeffs]
        assert 0 not in signs
        cuts = {i for i, s in enumerate(signs) if s > 0}
        assert sorted(cuts) == row['cut_caps']
        cut_sets.append(cuts)
    assert len(cut_sets) == 12
    assert all(any(not cuts.intersection(c) for cuts in cut_sets) for c in cliques)
    assert dotted['remaining_maximal_cliques'] == 0

    # A malformed cap is rejected by the underlying Gram equation, not by
    # comparing a previously produced result with itself.
    r = list(map(parse, saved['caps'][0]['normal']))
    target_cov = mv(gram, r)
    r[0] += 1
    assert mv(gram, r) != target_cov
    result = {
        'status': 'PASS_INDEPENDENT_GRAM_INVERSE_AND_BOTH_FINITE_BANKS',
        'sources': [str(BASE.relative_to(ROOT)), str(DOTTED.relative_to(ROOT))],
        'seed_signature': signature, 'independent_gram_basis': list(basis),
        'finite_assignments': 2601, 'finite_caps': 192,
        'finite_edges': len(finite_edges), 'finite_maximal_cliques': len(finite_cliques),
        'finite_max_caps': 6, 'finite_strict_witnesses': len(witness_cov),
        'additional_dotted_caps': 36, 'combined_caps': 228,
        'combined_edges': len(edges), 'combined_maximal_cliques': len(cliques),
        'combined_max_caps': 14, 'combined_strict_witnesses': len(cut_sets),
        'forged_cap_gram_equation_rejected': True,
        'scope': ('Every Coxeter completion retaining this literal seed and using only '
                  'the complete leaf-{0,-1,-phi} cap bank, optionally supplemented by '
                  'these 36 specific dotted caps, has infinite volume. This is not '
                  'an obstruction to arbitrary other caps or all dotted leaf values.'),
        'elapsed_seconds': round(time.monotonic() - start, 3),
    }
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = audit()
    if args.verify:
        old = json.loads(OUT.read_text())
        assert {k: v for k, v in old.items() if k != 'elapsed_seconds'} == {
            k: v for k, v in result.items() if k != 'elapsed_seconds'}
    else:
        OUT.write_text(json.dumps(result, indent=2) + '\n')
    print(result)


if __name__ == '__main__':
    main()
