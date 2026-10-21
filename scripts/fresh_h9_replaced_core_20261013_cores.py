#!/usr/bin/env python3
"""Exact pair and connected affine census for two changed39-wall cores."""
import argparse
import json
from itertools import combinations
from scripts.fresh_h9_replaced_core_20261013_base import *

OUT = ROOT/'data/enumerations/fresh_h9_replaced_core_20261013_cores.json'


def normalized_field(labels=()):
    e = Quadratic(K, scale(O, 2))
    f = Quadratic(e, e.lift(add(scale(O, 2), P)))
    weights = {2: f.zero, 3: f.rational(-1), 4: f.lift((Z, K.neg(O))),
               5: lift_K(f, K.neg(P)), 10: (e.zero, e.neg(e.one))}
    if 6 in labels:
        g = Quadratic(f, f.rational(3))
        weights = {m: g.lift(x) for m, x in weights.items()}
        weights[6] = f.zero, f.neg(f.one)
        f = g
    return f, weights


def inertia(f, matrix):
    a = [list(row) for row in matrix]
    positive = negative = zero = 0
    while a:
        pivot = next((i for i in range(len(a)) if f.sign(a[i][i])), None)
        if pivot is not None:
            order = [pivot]+[i for i in range(len(a)) if i != pivot]
            a = [[a[i][j] for j in order] for i in order]
            d = a[0][0]
            positive += f.sign(d) > 0
            negative += f.sign(d) < 0
            a = [[f.sub(a[i][j], f.div(f.mul(a[i][0], a[0][j]), d))
                  for j in range(1, len(a))] for i in range(1, len(a))]
            continue
        pair = next(((i, j) for i, j in combinations(range(len(a)), 2)
                     if f.sign(a[i][j])), None)
        if pair is None:
            zero += len(a)
            break
        order = list(pair)+[i for i in range(len(a)) if i not in pair]
        a = [[a[i][j] for j in order] for i in order]
        b = a[0][1]
        positive += 1
        negative += 1
        a = [[f.sub(a[i][j], f.div(f.add(f.mul(a[i][0], a[1][j]),
                                         f.mul(a[i][1], a[0][j])), b))
              for j in range(2, len(a))] for i in range(2, len(a))]
    return positive, negative, zero


def affine_census(names, pairs):
    labels = {r['label'] for r in pairs.values() if r['kind'] == 'finite'}
    f, weights = normalized_field(labels)
    def rel(i, j):
        return pairs['|'.join(sorted((i, j)))]
    seen, pending, affine, counts = set(), [(name,) for name in names], [], {}
    rows = []
    while pending:
        ids = pending.pop()
        if ids in seen:
            continue
        seen.add(ids)
        matrix = [[f.rational(2) if a == b else weights[rel(a, b)['label']]
                   for b in ids] for a in ids]
        sig = inertia(f, matrix)
        rows.append({'walls': list(ids), 'inertia': list(sig)})
        counts[str(len(ids))] = counts.get(str(len(ids)), 0)+1
        if sig[1]:
            continue
        if sig[2]:
            assert sig == (len(ids)-1, 0, 1)
            affine.append(list(ids))
            continue
        assert sig == (len(ids), 0, 0)
        assert len(ids) <= 9
        for name in names:
            if name in ids:
                continue
            relations = [rel(name, other) for other in ids]
            if any(r['kind'] != 'finite' for r in relations):
                continue
            if not any(r['label'] != 2 for r in relations):
                continue
            pending.append(tuple(sorted((*ids, name))))
    return {'connected_tested_subsets': len(rows), 'subset_counts_by_size': counts,
            'affine_circuits': sorted(affine), 'exact_inertias': sorted(rows, key=lambda r: r['walls']),
            'scope': 'Complete connected affine principal-subdiagram census. Every proper '
                     'principal subset of an irreducible affine Coxeter Gram matrix is '
                     'positive definite; connected PD subsets are exhaustively extended. '
                     'Indefinite or already-affine subsets need no extension.'}


def run():
    G, I, h, ps = profiles()
    pool = angle_pool(ANGLE_SOURCE)
    base_deleted = ['old31', 'old51', 'old197', 'old496', 'old1028', 'new4']
    union = [name for name in ps if name not in base_deleted]
    pairs = {}
    for a, b in combinations(union, 2):
        row = relation(ps[a], ps[b], pool)
        assert row['kind'] in ('finite', 'dotted')
        pairs['|'.join(sorted((a, b)))] = row
    cases = []
    for extra in ('old310', 'old837'):
        names = [name for name in union if name != extra]
        assert len(names) == 39
        census = affine_census(names, pairs)
        cases.append({'removed': base_deleted+[extra], 'walls': names,
                      'facet_count': len(names), **census})
        print({'removed_extra': extra, 'facets': len(names),
               'affine': census['affine_circuits'],
               'tested': census['connected_tested_subsets']}, flush=True)
    return {'status': 'EXACT_CHANGED_CORE_AUDIT',
            'source': 'data/enumerations/fresh_h9_quartic_completion_20261013_opposite.json',
            'base_source': BASE_SOURCE, 'common_timelike_center': list(map(enc, h)),
            'pair_checks': pairs, 'cases': cases,
            'scope': 'Two genuine39-wall replacements, preserving the ten-wall seed. '
                     'All pair angles are exact, with no parallel pairs. Absence of '
                     'affine circuits alone is not a compactness or volume certificate.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    args = p.parse_args()
    result = run()
    if args.verify:
        assert json.loads(OUT.read_text()) == result
    else:
        OUT.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'])
