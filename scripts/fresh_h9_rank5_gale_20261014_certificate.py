#!/usr/bin/env python3
"""Terminal ledger and an explicit rational primal realization; no search.

This uses stdlib Fraction arithmetic and imports no proposal or geometry engine.
The separate all-case audit remains the source of the full Bareiss replay.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_rank5_gale_20261014'
OUT = ROOT / f'results/{TAG}_certificate.json'


@lru_cache(None)
def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def read(path):
    d = json.loads((ROOT / path).read_text())
    for p, h in d.get('sources', {}).items():
        assert sha(p) == h, p
    return d


def inverse(A):
    n = len(A)
    B = [[F(x) for x in row] + [F(i == j) for j in range(n)]
         for i, row in enumerate(A)]
    for k in range(n):
        j = next(j for j in range(k, n) if B[j][k])
        B[j], B[k] = B[k], B[j]
        a = B[k][k]
        B[k] = [x / a for x in B[k]]
        for i in range(n):
            if i != k:
                a = B[i][k]
                B[i] = [x - a*y for x, y in zip(B[i], B[k])]
    return [row[n:] for row in B]


def rank(A):
    B = [[F(x) for x in row] for row in A]
    row = 0
    for col in range(len(B[0])):
        j = next((j for j in range(row, len(B)) if B[j][col]), None)
        if j is None:
            continue
        B[j], B[row] = B[row], B[j]
        a = B[row][col]
        B[row] = [x/a for x in B[row]]
        for i in range(row+1, len(B)):
            a = B[i][col]
            B[i] = [x-a*y for x, y in zip(B[i], B[row])]
        row += 1
        if row == len(B):
            break
    return row


def enc(x):
    if isinstance(x, (list, tuple)):
        return [enc(y) for y in x]
    x = F(x)
    return [x.numerator, x.denominator]


def explicit_primal(g):
    V, weights = g['gale'], g['positive_kernel']
    # H=diag(weights)*V has column sums zero. Its transpose kernel
    # is the affine evaluation space of the desired rational vertices.
    H = [[F(weights[i]*V[i][j]) for j in range(5)] for i in range(15)]
    A = [list(row) for row in zip(*H)]
    assert all(sum(row) == 0 for row in A)
    inv = inverse([row[:5] for row in A])
    K = [[-sum(inv[i][k]*A[k][j] for k in range(5))
          for j in range(5, 15)] for i in range(5)]
    K += [[F(i == j) for j in range(10)] for i in range(10)]
    assert all(sum(row) == 1 for row in K)
    assert all(sum(A[i][k]*K[k][j] for k in range(15)) == 0
               for i in range(5) for j in range(10))
    points = [row[:9] for row in K]
    assert rank([[F(1)] + row for row in points]) == 10
    assert len(set(map(tuple, points))) == 15
    supports = []
    for face in g['facets']:
        C = [i for i in range(15) if not face >> i & 1]
        # Fix the final positive-circuit coordinate to one, solve the other
        # five using the scaled Gale rows, and check all six are positive.
        B = [[H[i][j] for i in C[:5]] for j in range(5)]
        J = inverse(B)
        vals = [-sum(J[i][j]*H[C[5]][j] for j in range(5)) for i in range(5)] + [F(1)]
        assert all(x > 0 for x in vals)
        y = [F(0)]*15
        for i, x in zip(C, vals):
            y[i] = x
        assert all(sum(H[i][j]*y[i] for i in range(15)) == 0 for j in range(5))
        functional = [y[14]] + [y[i]-y[14] for i in range(5, 14)]
        values = [functional[0] + sum(a*b for a, b in zip(functional[1:], p)) for p in points]
        assert values == y
        assert sum(1 << i for i, x in enumerate(values) if x == 0) == face
        assert rank([[F(1)] + points[i] for i in range(15) if face >> i & 1]) == 9
        supports.append(dict(facet=face, nonnegative_affine_functional=enc(functional)))
    return dict(case=0, points=enc(points), affine_dimension=9,
                simplex_facet_dimension=8, facets=len(supports),
                supporting_functionals=supports,
                convention='functional[0] + dot(functional[1:], point) >= 0')


def run():
    primary_path = f'results/{TAG}.json'
    audit_path = f'results/{TAG}_independent.json'
    state_path = f'results/{TAG}_state.json'
    p, a, s = map(read, (primary_path, audit_path, state_path))
    assert p['status'] == 'FINISHED_BOUNDED_DIRECT_RANK5_GALE_SEARCH'
    assert a['status'] == 'PASS_ALL_SAVED_DIRECT_RANK5_GALE_CASES'
    assert s['status'] == p['status'] and s['phase'] == p['phase'] == 'complete'
    assert s['options'] == p['options'] and s['counts'] == p['counts']
    assert p['options'] == dict(random_seed=2026101406, random_trials=3000,
                              walks=4, steps_per_walk=10000, label_node_limit=3000)
    assert p['first_full_category_survivor'] is None and s['survivor'] is None
    assert s['walk'] == p['walk_index'] == 4
    count = len(p['cases'])
    assert count == p['distinct_exact_nerves'] == a['cases_checked'] == a['saved_primary_cases']
    audits = {r['case']: r for r in a['records']}
    assert len(audits) == count and set(audits) == set(range(count))
    canonical = set()
    literals = set()
    gates, origins = Counter(), Counter()
    vertices = []
    first = None
    for i, ref in enumerate(p['cases']):
        assert ref['case'] == i and sha(ref['path']) == ref['sha256']
        d = read(ref['path'])
        assert d['case'] == i and d['sources'] == p['sources']
        r = audits[i]
        assert r['input'] == ref['path'] and r['input_sha256'] == ref['sha256']
        assert r['status'] == 'PASS_SEPARATE_INTEGER_GALE_INCIDENCE_AND_LITERAL_GATE'
        saved_audit = read(f'results/{TAG}_audit_cases/case_{i:05d}.json')
        assert saved_audit['result'] == r
        g = d['geometry']
        key = tuple(g['canonical_minimal_nonfaces'])
        literal = tuple(g['minimal_nonfaces'])
        assert key not in canonical and literal not in literals
        canonical.add(key)
        literals.add(literal)
        assert g['simplex_boundary_link_factors'] == []
        assert r['vertices'] == len(g['facets'])
        assert r['gate'] == d['gate']['gate'] == ref['gate']
        gates[r['gate']] += 1
        origin = d['origin']
        origins[origin['phase'] + (str(origin['walk']) if 'walk' in origin else '')] += 1
        vertices.append(len(g['facets']))
        if i == 0:
            first = g
    assert gates == {'literal_induced_pattern': 902, 'two_lanner_blocks': 161}
    assert sum(vertices) == a['vertices_checked']
    assert a['five_minors_checked'] == 3003*count
    assert a['six_circuits_checked'] == 5005*count
    assert a['vertex_links_checked'] == 15*count
    c = p['counts']
    assert c['random_proposals']+c['walk_proposals'] == 43000
    assert sum(c[k] for k in ('large_minimal_nonface','unused_vertex',
                              'admissible_novel_family_proposals','simplex_link_factor')) == 43000
    files = [f'scripts/{TAG}_certificate.py', primary_path, audit_path, state_path]
    return dict(status='PASS_TERMINAL_LEDGER_AND_EXPLICIT_RATIONAL_PRIMAL',
                sources={f: sha(f) for f in files},
                saved_nerves=count, distinct_literal_minimal_nonface_sets=len(literals),
                distinct_saved_canonical_keys=len(canonical), origins=dict(origins),
                exact_exclusions=dict(gates), unresolved_saved_nerves=0,
                full_label_searches=0, metric_solves=0,
                facet_count=15, dimension=9,
                actual_vertex_count_range=[min(vertices),max(vertices)],
                independent_audit_counts={k:a[k] for k in ('vertices_checked','five_minors_checked','six_circuits_checked','vertex_links_checked')},
                explicit_rational_primal=explicit_primal(first),
                scope='All saved incidence/gate applications checked; inherited local obstruction proofs are dependencies, not recomputed. Distinct saved canonical keys use the existing colored-clutter pynauty routine; no separate isomorphism census. Bounded proposal search, not all fifteen-facet polytopes.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    out = run()
    if args.verify:
        assert out == json.loads(OUT.read_text())
    else:
        OUT.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], out['saved_nerves'], out['exact_exclusions'])
