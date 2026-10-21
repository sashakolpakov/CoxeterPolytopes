#!/usr/bin/env python3
"""Independent pair arithmetic and directed-capacity checks, without producer imports.

This does not reconstruct the producer's 555 physical projections. It checks
their resulting eleven tau classes with independent symbolic minimal
polynomials, and checks all five directed-incidence capacity classes.
"""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_multispectator_cap_capacity_20261017_review'
PRIMARY = 'fresh_sqrt5_multispectator_cap_capacity_20261017'
KINDS = ('disjoint', 'common_head', 'chain', 'swapped')
PATTERNS = ((False, False, False, False), (True, False, False, False),
            (True, False, False, True), (False, False, True, False),
            (False, True, False, False))


def capacity(mask, k):
    assert k >= 2 and mask in PATTERNS
    return {PATTERNS[0]: 1, PATTERNS[1]: k // 2,
            PATTERNS[2]: 2 * (k // 2), PATTERNS[3]: 1 if k == 2 else 3,
            PATTERNS[4]: k - 1}[mask]


def compatible(a, b, mask):
    # Ordered pairs are arrows: first entry is the tau-contact spectator.
    q, p = a
    r, s = b
    assert a != b
    if q == r:
        return False
    if p == s:
        return mask[1]
    if p == r and q == s:
        return mask[3]
    if p == r or q == s:
        return mask[2]
    assert not set(a) & set(b)
    return mask[0]


def clique_maximum(k, mask):
    arrows = list(permutations(range(k), 2))
    adj = [set(j for j, b in enumerate(arrows)
               if j != i and compatible(a, b, mask)) for i, a in enumerate(arrows)]
    best = 0

    def visit(chosen, remaining):
        nonlocal best
        best = max(best, chosen)
        while chosen + len(remaining) > best:
            i = remaining.pop()
            visit(chosen + 1, remaining & adj[i])

    visit(0, set(range(len(arrows))))
    return best


def run():
    source = ROOT / 'scripts' / (PRIMARY + '.py')
    data = ROOT / 'results' / (PRIMARY + '.json')
    d = json.loads(data.read_text())
    assert d['status'] == 'PASS_EXACT_TWO_A1_MULTISPECTATOR_CAP_CAPACITIES'
    for path, digest in d['sources'].items():
        assert sha256((ROOT / path).read_bytes()).hexdigest() == digest, path
    X = sp.Symbol('X')
    independent = []
    patterns = {}
    for t in d['tau_pair_tests']:
        a, b = map(sp.Rational, t['tau'])
        tau = a + b * sp.sqrt(5)
        assert sp.sign(tau) == 1
        h = (3 - tau) / 2
        assert sp.simplify(h + tau / 2) == sp.Rational(3, 2)
        products = (h, h + sp.Rational(1, 2), h + sp.sqrt(tau) / 2,
                    h + sp.sqrt(tau))
        mask = []
        for kind, product in zip(KINDS, products):
            sign = int(sp.sign(product))
            # A root is an algebraic integer iff its monic minimal polynomial
            # over Q has integer coefficients. No relative-norm implementation.
            polynomial = sp.Poly(sp.minimal_polynomial(product ** 2, X), X).monic()
            coefficients = polynomial.all_coeffs()
            integral = all(coefficient.q == 1 for coefficient in coefficients)
            allowed = sign <= 0 and integral
            saved = t['pairs'][kind]
            assert saved['product_sign'] == sign
            assert saved['two_cycle_integral'] == integral
            assert saved['allowed_by_necessary_tests'] == allowed
            mask.append(allowed)
            independent.append(dict(tau=t['tau'], kind=kind, product_sign=sign,
                square_minimal_polynomial_descending=list(map(str, coefficients)),
                two_cycle_integral=integral, allowed_by_necessary_tests=allowed))
        mask = tuple(mask)
        assert mask in PATTERNS and tuple(t['pattern']) == mask
        patterns[tuple(t['tau'])] = mask
    assert len(patterns) == 11 and len(independent) == 44

    occurrence_patterns = Counter()
    for o in d['occurrences']:
        mask = patterns[tuple(o['tau'])]
        assert tuple(o['pair_pattern']) == mask
        single = capacity(mask, o['A1_count'])
        assert o['single_orientation_capacity'] == single
        assert o['orientation_count'] == len(o['orientation_orbit'])
        assert o['capacity'] == single * o['orientation_count']
        occurrence_patterns[mask] += 1
    assert sum(occurrence_patterns.values()) == 555
    for c in d['columns']:
        single = capacity(tuple(c['pair_pattern']), c['A1_count'])
        assert c['single_orientation_capacity'] == single
        assert c['capacity'] == single * len(c['orientation_orbit'])
        assert Fraction(c['slack']) == c['capacity'] * Fraction(c['initial_simplex_mass']) - Fraction(c['load'])

    small = []
    for mask in PATTERNS:
        for k in range(2, 8):
            exact = clique_maximum(k, mask)
            assert exact == capacity(mask, k)
            small.append(dict(pattern=list(mask), copies=k, maximum=exact))
    return dict(status='PASS_INDEPENDENT_MULTISPECTATOR_PAIR_ARITHMETIC_AND_CAPACITIES',
        sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest()
                 for p in (Path(__file__), source, data)},
        counts=dict(tau_classes=11, independent_minimal_polynomials=44,
                    occurrence_capacity_checks=555, column_capacity_checks=len(d['columns']),
                    small_directed_graph_checks=len(small)),
        pair_checks=independent, small_directed_graph_checks=small,
        occurrence_pattern_counts=[dict(pattern=list(p), count=occurrence_patterns[p])
                                   for p in PATTERNS],
        scope='Independent pair signs and minimal-polynomial integrality, capacity identities, and small graph checks. All-k graph bounds have written proofs. Producer physical projection/maps/orbits received a source read, not a duplicate 555-projection reconstruction. No optimization or dimension claim.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    path = ROOT / 'results' / (TAG + '.json')
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['counts'], flush=True)
