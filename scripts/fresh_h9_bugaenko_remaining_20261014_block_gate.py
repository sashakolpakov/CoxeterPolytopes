#!/usr/bin/env python3
"""Exact two-Lanner-block obstruction and literal catalogue application.

No classification of whole products is used. The hypothesis is that every
union of proper subsets of the two blocks is spherical. Complete internal
Lanner tables and all permitted cross-edge categories are enumerated.
"""
import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path
import sympy as s
from fresh_h9_highbase_construction_20261014_cell import codes, finite

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_bugaenko_remaining_20261014_block_gate'
KINDS = [(4, 3), (5, 3), (4, 4), (5, 4), (5, 5)]
X = s.Symbol('x', real=True)


@lru_cache(None)
def canonical_codes(n):
    pairs = list(combinations(range(n), 2))
    index = {p: i for i, p in enumerate(pairs)}
    maps = [[index[tuple(sorted((p[i], p[j])))] for i, j in pairs]
            for p in permutations(range(n))]
    return sorted({min(tuple(c[i] for i in m) for m in maps) for c in codes(n)})


def gram(n, labels):
    G = s.eye(n)
    for (i, j), m in labels.items():
        G[i, j] = G[j, i] = 0 if m == 2 else (-X if m == 7 else -s.cos(s.pi / m))
    return G


def sign_proof_negative(expr):
    p = s.Poly(s.expand(s.radsimp(s.cancel(expr))), X)
    if p.degree() == 0:
        c = s.simplify(p.nth(0))
        assert c.is_negative is True, c
        return dict(polynomial=str(c), degree=0, negative_bernstein_coefficients=[str(c)])
    t = s.Symbol('t', real=True)
    shifted = s.Poly(s.expand(p.as_expr().subs(X, s.Rational(9, 10) + t / 10)), t)
    degree = shifted.degree()
    bernstein = [s.simplify(sum(shifted.nth(j) * s.binomial(i, j) / s.binomial(degree, j)
                                 for j in range(i + 1))) for i in range(degree + 1)]
    assert all(c.is_negative is True for c in bernstein), (p, bernstein)
    return dict(polynomial=str(p.as_expr()), degree=degree, interval=['9/10', '1'],
                negative_bernstein_coefficients=list(map(str, bernstein)))


def check_internal_tables():
    rows = []
    for n in [4, 5]:
        for code in canonical_codes(n):
            G = gram(n, dict(zip(combinations(range(n), 2), code)))
            det = s.simplify(G.det())
            assert det.is_negative is True
            for removed in range(n):
                keep = [i for i in range(n) if i != removed]
                B = G.extract(keep, keep)
                assert all(s.simplify(B[:k, :k].det()).is_positive is True for k in range(1, n))
            rows.append(dict(order=n, code=code, determinant=str(det)))
    # cos(pi/7) is the largest real root of this elementary triple-angle factor.
    p = s.Poly(8 * X**3 - 4 * X**2 - 4 * X + 1, X)
    assert p.count_roots(s.Rational(9, 10), 1) == 1
    assert p.count_roots(1, s.oo) == 0
    assert p.eval(s.Rational(9, 10)) < 0 and p.eval(1) > 0
    return rows


def table(n, m):
    A_pairs = list(combinations(range(n), 2))
    B_pairs = list(combinations(range(n, n + m), 2))
    cross = list(product(range(n), range(n, n + m)))
    proper = [[i for i in range(n + m) if i not in [a, b]]
              for a in range(n) for b in range(n, n + m)]
    proper_pairs = [tuple(combinations(v, 2)) for v in proper]
    def spherical(d):
        return all(finite(n + m - 2, tuple(d.get(p, 2) for p in pp)) for pp in proper_pairs)
    records = []
    internal_count = forced_right_count = 0
    for ca in canonical_codes(n):
        A = gram(n, dict(zip(A_pairs, ca)))
        inverse = A.inv()
        for cb in canonical_codes(m):
            internal_count += 1
            d = dict(zip(A_pairs, ca))
            d.update(zip(B_pairs, cb))
            active = []
            for pair in cross:
                d[pair] = 3
                if spherical(d):
                    active.append(pair)
                del d[pair]
            if not active:
                forced_right_count += 1
                continue
            completions = []
            for values in product(range(2, 8), repeat=len(active)):
                if max(values) == 2:
                    continue  # two unchanged orthogonal negative blocks
                trial = dict(d)
                trial.update(zip(active, values))
                if not spherical(trial):
                    continue
                assert all(a <= 6 for a in values)
                assert sum(a == 7 for a in ca + cb) <= 1
                G = gram(n + m, trial)
                S = G[n:, n:] - G[n:, :n] * inverse * G[:n, n:]
                proof = sign_proof_negative(S.det())
                completions.append(dict(cross_values=values, schur=proof))
            records.append(dict(first_code=ca, second_code=cb,
                                potentially_nonright_cross_pairs=active, completions=completions))
    return dict(orders=[n, m], internal_type_pairs=internal_count,
        all_cross_pairs_forced_right=forced_right_count,
        internal_pairs_allowing_nonright=records,
        nonright_completions=sum(len(r['completions']) for r in records))


def applications(inc):
    records = []
    allowed = {tuple(sorted(k)) for k in KINDS}
    for row in inc['cases']:
        masks = row['canonical_minimal_nonfaces']
        chosen = None
        for i, a in enumerate(masks):
            for j in range(i + 1, len(masks)):
                b = masks[j]
                if a & b or tuple(sorted([a.bit_count(), b.bit_count()])) not in allowed:
                    continue
                union = a | b
                if any(k not in [i, j] and c & union == c for k, c in enumerate(masks)):
                    continue
                A = [k for k in range(14) if a >> k & 1]
                B = [k for k in range(14) if b >> k & 1]
                # Direct actual-face check, not just a clutter heuristic.
                for x in A:
                    for y in B:
                        face = union ^ (1 << x) ^ (1 << y)
                        assert any(v & face == face for v in row['canonical_vertices'])
                chosen = dict(case=row['case'], first_missing_face=A, second_missing_face=B,
                              minimal_nonface_indices=[i, j], proper_union_checks=len(A) * len(B))
                break
            if chosen:
                break
        if chosen:
            records.append(chosen)
    return records


def run():
    tables = [table(n, m) for n, m in KINDS]
    assert [r['nonright_completions'] for r in tables] == [10, 27, 0, 2, 3]
    ip = 'results/fresh_h9_bugaenko_surgery_20261014_incidence.json'
    inc = json.loads((ROOT / ip).read_text())
    for p, h in inc['sources'].items():
        assert sha256((ROOT / p).read_bytes()).hexdigest() == h, p
    applied = applications(inc)
    assert len(applied) == 2365
    assert set([22, 23, 318, 582, 629, 633, 636, 637, 995]) <= {r['case'] for r in applied}
    files = [ip, str(Path(__file__).resolve().relative_to(ROOT)),
        'scripts/fresh_h9_highbase_construction_20261014_cell.py',
        'scripts/independent_h9_14_record30_a1_labels.py']
    return dict(status='PASS_EXACT_TWO_LANNER_BLOCK_OBSTRUCTION',
        sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in files},
        internal_table_counts={str(n): len(canonical_codes(n)) for n in [3, 4, 5]},
        internal_signature_checks=check_internal_tables(), tables=tables,
        excluded_catalogue_types=len(applied), applications=applied,
        scope='Two disjoint Lanner blocks of orders (4,3),(5,3),(4,4),(5,4),(5,5), '
              'all unions of proper block subsets spherical. Least-label tests force every '
              'unlisted cross pair right; every remaining finite-label category is exhausted. '
              'Orthogonal blocks give two negatives. Otherwise a negative exact Schur '
              'determinant gives at least one further negative direction beyond the first '
              'Lanner block. Category7 covers all finite m>=7; Bernstein coefficients prove '
              'the Schur determinant strictly negative on [9/10,1], containing cos(pi/m). '
              'The literal catalogue face checks establish the hypothesis. No whole-polytope '
              'product classification, numerical inertia, or dotted-distance bound is used. '
              'Complete Lanner/finite-diagram classifications are source-bound dependencies.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    a = ap.parse_args()
    data = json.loads(json.dumps(run()))
    target = ROOT / f'results/{TAG}.json'
    if a.verify:
        assert data == json.loads(target.read_text())
    else:
        target.write_text(json.dumps(data, indent=2) + '\n')
    print(data['status'], data['excluded_catalogue_types'])
