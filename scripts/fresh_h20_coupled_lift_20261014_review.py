#!/usr/bin/env python3
"""Independent small Fraction geometry audit, not a second all-pair census."""
import argparse
import hashlib
import json
from fractions import Fraction as F
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_coupled_lift_20261014'


def read(path):
    return json.loads((ROOT / path).read_text())


def mv(A, v):
    return [sum((x*y for x, y in zip(row, v)), F(0)) for row in A]


def scalar(x, y):
    return sum((a*b for a, b in zip(x, y)), F(0))


def sign(a, b, d):
    assert d > 0
    if not b:
        return (a > 0) - (a < 0)
    if not a:
        return (b > 0) - (b < 0)
    if a*b > 0:
        return (a > 0) - (a < 0)
    s = a*a-d*b*b
    return ((a > 0)-(a < 0))*((s > 0)-(s < 0))


def positive_pivots(A):
    A = [row[:] for row in A]
    out = []
    while A:
        p = A[0][0]
        assert p > 0
        out.append(p)
        A = [[A[i][j]-A[i][0]*A[0][j]/p
              for j in range(1, len(A))] for i in range(1, len(A))]
    return out


def run():
    doc = read(f'results/{TAG}.json')
    old = read('results/fresh_dimension_construction_routes_20261014.json')
    for d in (doc, old):
        for name, digest in d['sources'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    Q = [[F(x) for x in row] for row in doc['form']]
    assert len(Q) == 21 and all(len(row) == 21 for row in Q)
    assert all(Q[i][j] == Q[j][i] for i in range(21) for j in range(21))
    assert all(F(x['b']) == 0 for row in old['form'] for x in row)
    assert Q == [[F(x['a']) for x in row] for row in old['form']]

    def inner(v, w):
        return scalar(v, mv(Q, w))

    roots = [(list(map(F, r['a'])), list(map(F, r['b'])), F(r['d']))
             for r in doc['normals']]
    assert len(roots) == 92
    covs = [(mv(Q, a), mv(Q, b), d) for a, b, d in roots]

    def at(i, p):
        a, b, d = covs[i]
        return sign(scalar(p, a), scalar(p, b), d)

    H = list(map(F, doc['strict_center']))
    assert inner(H, H) == F(doc['center_norm']) == -4296
    for i, ((a, b, d), (qa, qb, _)) in enumerate(zip(roots, covs)):
        assert scalar(a, qb) == 0
        assert scalar(a, qa)+d*scalar(b, qb) == F(doc['norms'][i]) > 0
        assert at(i, H) == -1
    for i, r in enumerate(old['normals']):
        assert roots[i] == ([F(x['a']) for x in r], [F(x['b']) for x in r], F(2))

    # Four exact rational support/frame reconstructions establish rank 20.
    for rec in doc['base_recipes']:
        S = rec['support']
        assert len(S) == len(set(S)) == 18
        assert all(not any(roots[i][1]) for i in S)
        Gram = [[scalar(roots[i][0], covs[j][0]) for j in S] for i in S]
        assert positive_pivots(Gram) == list(map(F, rec['pd18_pivots']))
        h, w, m, p = [list(map(F, rec[k])) for k in ('h0', 'w', 'm', 'finite_vertex')]
        v = [F(x['a']) for x in old['strict_escape']['boundary_point']]
        assert all(F(x['b']) == 0 for x in old['strict_escape']['boundary_point'])
        assert m == [x+F(rec['t'])*y for x, y in zip(v, h)]
        d, q = F(rec['radicand']), inner(m, m)
        assert q > 0 and inner(w, w) > 0 and d*inner(w, w) == q
        assert inner(m, w) == 0
        assert inner(h, h) < 0 and inner(p, p) < 0 and inner(p, H) < 0
        assert p == [x-inner(h, m)*y/q for x, y in zip(h, m)]
        for i in S:
            assert scalar(h, covs[i][0]) == scalar(m, covs[i][0]) == scalar(w, covs[i][0]) == 0
        assert inner(p, m) == inner(p, w) == 0
        assert F(rec['cap_norm']) == 2*q

    # The source's height sign changes really are form isometries, not merely
    # set permutations. Recheck all old strict escape images independently.
    ea = [F(x['a']) for x in old['strict_escape']['point']]
    eb = [F(x['b']) for x in old['strict_escape']['point']]
    flips = 0
    old_target_signs = 0
    for block, signs in enumerate(product((1, -1), repeat=2)):
        D = [1]*19+list(signs)
        assert all(D[i]*D[j]*Q[i][j] == Q[i][j] for i in range(21) for j in range(21))
        flip = lambda v: [x*s for x, s in zip(v, D)]
        assert flip(H) == H
        for i, j in enumerate(doc['height_sign_permutations'][block]):
            assert (flip(roots[i][0]), flip(roots[i][1]), roots[i][2]) == roots[j]
        for j, rec in enumerate(doc['base_recipes']):
            m, w = [list(map(F, rec[k])) for k in ('m', 'w')]
            d = F(rec['radicand'])
            assert roots[60+8*block+2*j] == (flip(m), flip(w), d)
            assert roots[61+8*block+2*j] == (flip(m), flip([-x for x in w]), d)
        a, b = flip(ea), flip(eb)
        assert sign(inner(a, a)+2*inner(b, b), 2*inner(a, b), F(2)) > 0
        assert sign(inner(a, H), inner(b, H), F(2)) < 0
        for qa, qb, d in covs[:60]:
            assert d == 2
            assert sign(scalar(a, qa)+2*scalar(b, qb), scalar(a, qb)+scalar(b, qa), d) < 0
            old_target_signs += 1
        flips += 1

    zeros = set()
    vertex_tests = 0
    for row in doc['finite_vertices']:
        p = list(map(F, row['point']))
        assert inner(p, p) == F(row['norm']) < 0 and inner(p, H) < 0
        s = [at(i, p) for i in range(92)]
        assert all(x <= 0 for x in s)
        inc = tuple(i for i, x in enumerate(s) if x == 0)
        assert len(inc) == 20 and inc == tuple(sorted(row['incident']))
        zeros.add(inc)
        vertex_tests += len(s)
    assert len(zeros) == 16

    cusp_tests = 0
    for rec, saved in zip(old['repaired_cusps'], doc['preserved_three_cusp_links']):
        p = [F(x['a']) for x in rec['delta']]
        assert all(F(x['b']) == 0 for x in rec['delta'])
        assert inner(p, p) == 0 and inner(p, H) < 0
        assert all(at(i, p) <= 0 for i in range(60))
        assert [at(i, p) for i in range(60, 92)] == saved['new_signs'] == [-1]*32
        cusp_tests += 32
    assert cusp_tests == 96
    p = list(map(F, doc['strict_escape']['point']))
    assert inner(p, p) == F(doc['strict_escape']['norm']) > 0
    assert inner(p, H) == F(doc['strict_escape']['center_pairing']) < 0
    assert all(at(i, p) < 0 for i in range(92))

    files = [f'scripts/{TAG}.py', f'results/{TAG}.json',
             'results/fresh_dimension_construction_routes_20261014.json',
             'data/enumerations/vinberg_h18_mixture_facets_2026_10_13.json',
             f'scripts/{TAG}_review.py']
    return dict(status='PASS_SCOPED_INDEPENDENT_FRACTION_GEOMETRY_AND_SOURCE_REVIEW',
                sources={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
                root_norms_and_center=92, rational_support_frames=4,
                sign_flip_isometries=flips, old_escape_image_signs=old_target_signs,
                actual_finite_vertices=16, vertex_pairing_tests=vertex_tests,
                new_cusp_strict_pairings=cusp_tests, new_escape_strict_pairings=92,
                new_escape_norm=doc['strict_escape']['norm'],
                scope='Separate Fraction geometry checks from frozen literals; no producer imports. '
                      'All 2416 new pair gates and 320 open germs additionally received source/proof '
                      'review and a fresh primary saved replay, not a second full pair census. '
                      'Old signature/old-old pairs and bounded links remain inherited dependencies.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    out = run()
    path = ROOT / f'results/{TAG}_review.json'
    if args.verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], out['vertex_pairing_tests'], out['new_escape_norm'])
