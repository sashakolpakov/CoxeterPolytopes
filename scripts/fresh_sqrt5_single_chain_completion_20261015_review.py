#!/usr/bin/env python3
"""Independent Fraction Q(sqrt5) replay of all-length single-chain completion.

Imports only the previously frozen independent field/matrix checker, never a
producer. Gram entries come from independent DSU colors and exact square roots.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import fresh_sqrt5_quad_filters_20261015_review as R

ROOT = Path(__file__).resolve().parents[1]
BASE = 'fresh_sqrt5_single_chain_completion_20261015'
TAG = BASE + '_review'


def dot(a, b):
    out = R.Z
    for x, y in zip(a, b):
        out = R.add(out, R.mul(x, y))
    return out


def chain(core, contacts, length):
    n = core['nodes']
    old = dict(zip(combinations(range(n), 2), core['labels']))
    g = {p: old.get(p, 2) for p in combinations(range(n + length), 2)}
    for i, m in enumerate(contacts):
        g[i, n] = m
    for i in range(n, n + length - 1):
        g[i, i + 1] = 3
    return dict(nodes=n + length, white=list(range(4, n + length)),
                labels=[g[p] for p in combinations(range(n + length), 2)])


def check_unpaid(row):
    n = row['nodes']
    g = dict(zip(combinations(range(n), 2), row['labels']))
    edge = lambda i, j: g[tuple(sorted((i, j)))]
    for a, b in product((0, 2), (1, 3)):
        assert edge(a, b) == 2
        assert not any(edge(a, w) != 2 and edge(b, w) != 2 for w in row['white'])
        ids = row['white'] + [a, b]
        adj = {i: [j for j in ids if j != i and edge(i, j) != 2] for i in ids}
        seen = {a}
        todo = [a]
        while todo:
            for j in adj[todo.pop()]:
                if j not in seen:
                    seen.add(j)
                    todo.append(j)
        if b not in seen or len(seen) != 4:
            continue
        degrees = sorted(len(adj[i]) for i in seen)
        labels = sorted(edge(i, j) for i, j in combinations(sorted(seen), 2)
                        if edge(i, j) != 2)
        assert not (degrees == [1, 1, 2, 2] and labels == [3, 3, 5]
                    and len(adj[a]) == len(adj[b]) == 1)


def run():
    primary = json.loads((ROOT / 'results' / (BASE + '.json')).read_text())
    source = json.loads((ROOT / 'results/fresh_sqrt5_quad_contacts_20261015.json').read_text())
    union = json.loads((ROOT / 'results/fresh_sqrt5_quad_exhaustion_20261015.json').read_text())
    expected = {}
    for c in source['records']:
        if c['status'] == 'UNPAID_FINITE_CONTACT_DOMAIN' and not c['singular'] and c['accepted']:
            for j, t in enumerate(c['trials']):
                if t['status'] == 'PASS_EXTRA_WHITE_NECESSARY_TESTS':
                    expected[c['class_index'], j] = (c, t)
    assert len(expected) == 458
    assert {(r['source_class_index'], r['trial_index']) for r in primary['records']} == set(expected)
    assert len(primary['records']) == 458
    cache, records = {}, []
    counts, lengths = Counter(), Counter()
    corners = 0
    for saved in primary['records']:
        key = saved['source_class_index'], saved['trial_index']
        c, trial = expected[key]
        row = union['labelled_union_supports'][c['union_index']]
        assert trial['contacts'] == saved['contacts']
        if c['class_index'] not in cache:
            colors = R.coloring(row)
            assert colors is not None
            H = R.matrix(row, colors)
            Hi = R.inverse(H)
            n = len(H)
            assert all(dot(H[i], [Hi[k][j] for k in range(n)]) == (R.ONE if i == j else R.Z)
                       for i in range(n) for j in range(n))
            assert R.inertia(H) == [n - 1, 1, 0]
            assert R.inertia([[R.conj(x) for x in r] for r in H]) == [n, 0, 0]
            cache[c['class_index']] = H, Hi, colors
        H, Hi, oldcolors = cache[c['class_index']]
        one = chain(row, trial['contacts'], 1)
        colors = R.coloring(one)
        assert colors is not None and colors[:-1] == oldcolors
        G = R.matrix(one, colors)
        v, d = [r[-1] for r in G[:-1]], G[-1][-1]
        cost = dot(v, [dot(r, v) for r in Hi])
        denom = R.sub(R.mul(R.TWO, cost), d)
        assert d == R.unpack(saved['first_norm'])
        assert cost == R.unpack(saved['projection_cost'])
        assert denom == R.unpack(saved['denominator'])
        out = dict(source_class_index=key[0], trial_index=key[1],
                   first_norm=R.enc(d), projection_cost=R.enc(cost), denominator=R.enc(denom))
        if denom == R.Z:
            status = 'NO_LENGTH_ZERO_DENOMINATOR'
        else:
            length = R.div(d, denom)
            assert length == R.unpack(saved['length_in_field'])
            out['length_in_sqrt5_basis'] = R.enc(length)
            if length[1] or length[0].denominator != 1 or length[0] < 1:
                status = 'NO_POSITIVE_INTEGER_LENGTH'
            else:
                ell = int(length[0])
                assert R.mul((F(ell + 1), F(0)), d) == R.mul((F(2 * ell), F(0)), cost)
                whole = chain(row, trial['contacts'], ell)
                assert saved['length'] == ell and saved['row'] == whole
                colors = R.coloring(whole)
                assert colors is not None
                M = R.matrix(whole, colors)
                N = len(M)
                phy = R.inertia(M)
                other = R.inertia([[R.conj(x) for x in r] for r in M])
                assert phy == saved['physical_inertia'] == [N - 2, 1, 1]
                assert other == saved['conjugate_inertia'] == [N - 1, 0, 1]
                for a, b in product((0, 2), (1, 3)):
                    ids = whole['white'] + [a, b]
                    C = [[M[i][j] for j in ids] for i in ids]
                    assert R.inertia(C) == [len(ids), 0, 0]
                    assert R.inertia([[R.conj(x) for x in r] for r in C]) == [len(ids), 0, 0]
                    corners += 1
                # Independent direct chain inverse, not the closed formula used
                # by the producer, verifies the rank-one determinant factor.
                n = row['nodes']
                Q = [[M[i][j] for j in range(n, N)] for i in range(n, N)]
                Qi = R.inverse(Q)
                assert Qi[0][0] == R.div((F(2 * ell), F(0)), R.mul(d, (F(ell + 1), F(0))))
                assert R.sub(R.ONE, R.mul(cost, Qi[0][0])) == R.Z
                check_unpaid(whole)
                status = 'UNIQUE_LENGTH_UNPAID_PRINCIPAL_STAR'
                lengths[ell] += 1
                out.update(length=ell, physical_inertia=phy, conjugate_inertia=other)
        assert status == saved['status']
        counts[status] += 1
        out['status'] = status
        records.append(out)
    assert dict(counts) == primary['counts']
    assert {str(k): v for k, v in lengths.items()} == primary['unique_positive_integer_lengths']
    assert corners == 4 * 163
    files = [f'scripts/{TAG}.py', f'scripts/{BASE}.py', f'results/{BASE}.json',
             'scripts/fresh_sqrt5_quad_filters_20261015_review.py',
             'results/fresh_sqrt5_quad_filters_20261015_review.json',
             'results/fresh_sqrt5_quad_contacts_20261015.json',
             'results/fresh_sqrt5_quad_contacts_20261015_review.json',
             'results/fresh_sqrt5_quad_exhaustion_20261015.json']
    return dict(status='PASS_INDEPENDENT_ALL_LENGTH_SINGLE_CHAIN_REVIEW',
                sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in files},
                profiles=len(records), nonsingular_cores=len(cache), counts=dict(counts),
                valid_lengths={str(k): v for k, v in sorted(lengths.items())},
                finite_corners_in_both_embeddings=corners, records=records,
                arithmetic='Independent Fraction basis1,sqrt5, DSU gauge, exact-square-root Gram entries, inverse identity checks and full congruence inertia; no producer imports.',
                scope='All458 selected finite profiles and all163 resulting principal stars. The all-length proof is for one simply-laced path with only its first node touching the core, plus orthogonal positive-definite spectators. No arbitrary branching, multiple attached components, opposite-end attachments, full polytope construction or global dimension bound.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    out = run()
    p = ROOT / 'results' / (TAG + '.json')
    if args.verify:
        assert out == json.loads(p.read_text())
    else:
        p.write_text(json.dumps(out, indent=2) + '\n')
    print(out['status'], out['profiles'], out['counts'], out['valid_lengths'])
