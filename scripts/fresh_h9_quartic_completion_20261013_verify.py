#!/usr/bin/env python3
"""Independent Fraction replay of full quartic-anchor root certificates."""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from math import comb
import json
from pathlib import Path

from scripts.fresh_h9_quartic_completion_20261013_arithmetic import (
    K, Quadratic, Z, O, P, add, sub, mul, div, scale, parse, enc,
    dot, pairing, matvec, finite_label, old_pair_gate, new_pair_gate)

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads((ROOT/path).read_text())


def matrix(rows):
    return tuple(tuple(map(parse, row)) for row in rows)


def node(row, G):
    a = tuple(map(parse, row['coefficients_rational']))
    b = tuple(map(parse, row['coefficients_radical']))
    d = parse(row['discriminant'])
    assert K.sign(d) > 0
    return {'a': a, 'b': b, 'd': d, 'q': scale(O, 2),
            'pa': matvec(K, G, a), 'pb': matvec(K, G, b)}


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def polynomial_division(p, q):
    p, q = trim(p), trim(q)
    out = [Q(0)]*max(1, len(p)-len(q)+1)
    while len(p) >= len(q) and any(p):
        k = len(p)-len(q)
        c = p[-1]/q[-1]
        out[k] += c
        for i, v in enumerate(q):
            p[i+k] -= c*v
        p = trim(p)
    return trim(out), p


@lru_cache(None)
def cyclotomic(n):
    p = [Q(-1)]+[Q(0)]*(n-1)+[Q(1)]
    for d in range(1, n):
        if n % d == 0:
            p, rem = polynomial_division(p, cyclotomic(d))
            assert not any(rem)
    return tuple(p)


def evaluate(p, x):
    out = Q(0)
    for c in reversed(p):
        out = out*x+c
    return out


def variations(chain, x):
    signs = [1 if y > 0 else -1 for p in chain if (y := evaluate(p, x))]
    return sum(a != b for a, b in zip(signs, signs[1:]))


def angle_pool(path):
    def phi(m):
        return sum(__import__('math').gcd(m, k) == 1 for k in range(1, m+1))
    expected = [m for m in range(2, 2049) if m == 2 or phi(m) in (2, 4, 8, 16, 32)]
    rows = read(path)
    assert [r['label'] for r in rows] == expected
    out = []
    for r in rows:
        m = r['label']
        coefficients = tuple(map(Q, r['coefficients_descending']))
        lo, hi = map(Q, r['largest_root_interval'])
        degree = len(coefficients)-1
        p = list(reversed(coefficients))
        if m == 2:
            assert p == [0, 1] and lo == hi == 0
        else:
            assert degree == phi(m)//2
            substituted = [Q(0)]*(2*degree+1)
            for k, c in enumerate(p):
                for j in range(2*k+1):
                    substituted[j+degree-k] += c*comb(2*k, j)
            assert tuple(substituted) == cyclotomic(m)
            if lo == hi:
                assert degree == 1 and evaluate(p, lo) == 0
            else:
                derivative = [i*p[i] for i in range(1, len(p))]
                chain = [p, derivative]
                while len(chain[-1]) > 1:
                    _, remainder = polynomial_division(chain[-2], chain[-1])
                    chain.append([-x for x in remainder])
                assert variations(chain, lo)-variations(chain, hi) == 1
                assert variations(chain, hi)-variations(chain, Q(4)) == 0
        out.append((m, coefficients, lo, hi))
    return out


def load_old(G):
    graph = read('data/enumerations/fresh_h8_alternate_selection_20261013.json')
    base, polar, coupled = [read(graph[key]) for key in ('base_source', 'polar_source', 'coupled_source')]
    old = []
    for record in graph['nodes']:
        kind, index = record['source'], record['index']
        if kind == 'coupled':
            n = node(coupled['caps'][index], G)
        else:
            row = base['caps'][index] if kind == 'base' else polar['polar_caps'][index]
            assert kind in ('base', 'polar')
            a = tuple(map(parse, row['seed_coefficients']))
            n = {'a': a, 'b': (Z,)*10, 'd': Z, 'q': parse(row['norm']),
                 'pa': matvec(K, G, a), 'pb': (Z,)*10}
        old.append(n)
    assert len(old) == 493
    for path in ('data/enumerations/fresh_h8_threefree_combined_caps_20261013.json',
                 'data/enumerations/fresh_h8_cap_completion_20261013_bottleneck_caps.json'):
        old += [node(c, G) for c in read(path)['caps']]
    assert len(old) == 1084
    return old


def check(source):
    data = read(source)
    base = read(data['base_source'])
    G, I = matrix(base['seed_gram']), matrix(base['seed_gram_inverse'])
    center = tuple(map(parse, base['strict_center_y']))
    old = load_old(G)
    anchor = data['anchor']
    assert anchor == 1039 and parse(data['anchor_radicand']) == old[anchor]['d']
    E = Quadratic(K, old[anchor]['d'])
    IE = tuple(tuple(E.lift(x) for x in row) for row in I)
    ar = tuple(zip(old[anchor]['a'], old[anchor]['b']))
    chosen = data['selected_configuration']
    assert len(set(chosen)) == len(chosen) == 22 and anchor in chosen
    target = read(data['target_source'])
    ys = [tuple(map(parse, r['certificate']['y'])) for r in target['targets']]
    for y in ys:
        assert K.sign(pairing(K, y, I, y)) > 0
        assert K.sign(pairing(K, y, I, center)) < 0
        assert all(K.sign(x) > 0 for x in y)
    polynomials = angle_pool(data['angle_source'])
    required = []
    if data.get('required_cap_source'):
        rr = read(data['required_cap_source'])
        assert parse(rr['anchor_radicand']) == E.d
        for index in data['required_cap_indices']:
            c = rr['caps'][index]
            f = Quadratic(E, E.parse(c['discriminant_over_E']))
            required.append((index, f.d, tuple(map(f.parse, c['seed_coefficients']))))
    counts = {'caps': 0, 'seed_pairs': 0, 'old_pairs': 0, 'target_signs': 0}
    for c in data['caps']:
        w0, w1 = (tuple(map(E.parse, c[name])) for name in ('line_constant', 'line_direction'))
        A = pairing(E, w1, IE, w1)
        B = E.mul(E.rational(2), pairing(E, w0, IE, w1))
        C = E.sub(pairing(E, w0, IE, w0), E.rational(2))
        assert [A, B, C] == list(map(E.parse, c['norm_polynomial_descending']))
        D = E.sub(E.mul(B, B), E.mul(E.rational(4), E.mul(A, C)))
        assert E.sign(D) > 0 and D == E.parse(c['discriminant_over_E'])
        field = Quadratic(E, D)
        x = field.parse(c['parameter'])
        value = field.add(field.add(field.mul(field.lift(A), field.mul(x, x)),
                                    field.mul(field.lift(B), x)), field.lift(C))
        assert field.sign(value) == 0
        pa = tuple(map(field.parse, c['seed_pairings']))
        ca = tuple(map(field.parse, c['seed_coefficients']))
        assert len(pa) == len(ca) == len(w0) == len(w1) == 10
        for p, a, b in zip(pa, w0, w1):
            assert field.sign(field.add(p, field.add(field.lift(a),
                                                     field.mul(x, field.lift(b))))) == 0
        GF = tuple(tuple(field.lift(E.lift(t)) for t in row) for row in G)
        assert all(field.sign(field.sub(p, q)) == 0
                   for p, q in zip(pa, matvec(field, GF, ca)))
        assert field.sign(field.sub(dot(field, pa, ca), field.rational(2))) == 0
        center_value = dot(field, ca, tuple(field.lift(E.lift(t)) for t in center))
        assert field.sign(center_value) > 0
        other = c['rational_anchor']
        assert (other in chosen and old[other]['d'] == Z) or (
            data.get('include_seed_anchors') and -10 <= other < 0)
        M, N = map(parse, c['anchor_pairing_magnitudes'])
        br = (tuple(E.lift(t) for t in old[other]['a']) if other >= 0 else
              tuple(E.one if j == -1-other else E.zero for j in range(10)))
        for row, value in ((ar, M), (br, N)):
            assert E.sign(E.sub(dot(E, row, w0), E.lift(value))) == 0
            assert E.sign(dot(E, row, w1)) == 0
            assert field.sign(field.add(dot(field, pa, tuple(field.lift(t) for t in row)),
                                        field.lift(E.lift(value)))) == 0
        seed_rows = []
        for j, p in enumerate(pa):
            assert field.sign(p) <= 0
            square = field.mul(p, p)
            comparison = field.sign(field.sub(square, field.rational(4)))
            label = None if comparison >= 0 else finite_label(field, square, polynomials)
            kind = 'dotted' if comparison > 0 else 'parallel' if comparison == 0 else 'finite' if label is not None else 'noncoxeter'
            assert kind in ('finite', 'dotted')
            seed_rows.append({'seed': j, 'kind': kind, 'label': label})
        assert seed_rows == c['seed_relations']
        rows = [{'old_node': j, **old_pair_gate(field, pa, old[j], polynomials)}
                for j in chosen]
        assert rows == c['pair_checks']
        kept = [r['old_node'] for r in rows if r['compatible']]
        dropped = [r['old_node'] for r in rows if not r['compatible']]
        assert kept == c['retained_nodes'] and dropped == c['dropped_nodes']
        assert anchor in kept and (other < 0 or other in kept) and len(kept) >= data['minimum_retained']
        rchecks = [{'cap_index': j, **new_pair_gate(field, pa, dd, cc, polynomials)}
                   for j, dd, cc in required]
        assert rchecks == c.get('required_cap_pair_checks', [])
        assert all(r['compatible'] for r in rchecks)
        cuts = [j for j, y in enumerate(ys)
                if field.sign(dot(field, ca, tuple(field.lift(E.lift(t)) for t in y))) < 0]
        assert cuts == c['cut_target_indices'] and set(cuts) & set(data.get('focus_target_indices', [5, 22]))
        counts['caps'] += 1
        counts['seed_pairs'] += 10
        counts['old_pairs'] += len(chosen)
        counts['target_signs'] += len(ys)
    return {'status': 'PASS_FULL_QUARTIC_ANCHOR_ROOT_REPLAY',
            'source': source, 'field_degree_bound': 8, 'pair_field_degree_bound': 16,
            'safe_labels': [r[0] for r in polynomials], **counts,
            'scope': 'Accepted roots, exact full anchor conditions, all saved seed and old-cap '
                     'relations, center orientations and target signs; no volume claim.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--source', default='data/enumerations/fresh_h9_quartic_completion_20261013_targeted.json')
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    result = check(a.source)
    output = ROOT / a.source.replace('.json', '_verify.json')
    if a.verify:
        assert json.loads(output.read_text()) == result
    else:
        output.write_text(json.dumps(result, indent=2)+'\n')
    print(result, flush=True)
