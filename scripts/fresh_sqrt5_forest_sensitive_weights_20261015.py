#!/usr/bin/env python3
"""Exact A1-sensitive partial weight rule and all rank29 vertex budgets.

Covers the112 previously sealed quadrilateral cores. The130 finite and24
symbolic dotted unpaid extension frontiers remain outside this result.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import fresh_sqrt5_sealed_weights_20261015 as T
import fresh_sqrt5_uniform_face_weights_20261015 as U

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_forest_sensitive_weights_20261015'


def forest(H):
    remaining = set(H)
    out = []
    while remaining:
        ids = T.component(H, min(remaining))
        remaining.difference_update(ids)
        out.append({i: {j: m for j, m in H[i].items() if j in ids} for i in ids})
    return out


def dynamic(components, totals, bound, allow_a1):
    best = [F(0)] + [None] * bound
    witnesses = [()] + [None] * bound
    for n in range(1, bound+1):
        options = [(best[n-r['rank']] + v, tuple(sorted(witnesses[n-r['rank']] + (r['type'],))))
                   for r, v in zip(components, totals)
                   if r['rank'] <= n and best[n-r['rank']] is not None
                   and (allow_a1 or r['rank'] != 1)]
        if options:
            best[n], witnesses[n] = max(options)
    return best, witnesses


def run():
    sealed = json.loads((ROOT/'results/fresh_sqrt5_sealed_weights_20261015.json').read_text())
    contacts = json.loads((ROOT/'results/fresh_sqrt5_quad_contacts_20261015.json').read_text())
    unions = json.loads((ROOT/'results/fresh_sqrt5_quad_exhaustion_20261015.json').read_text())
    triangle = json.loads((ROOT/'results/fresh_sqrt5_triangle_weights_20261015.json').read_text())
    selected = {r['class_index']: r for r in contacts['records']
                if r['status']=='UNPAID_FINITE_CONTACT_DOMAIN' and r['singular'] and not r['accepted']}
    assert len(selected) == 112
    pure_d7 = set()
    for row in sealed['face_rows']:
        if row['class_index'] in (535, 563):
            assert len(row['coefficients']) == 1
            pure_d7.add(row['coefficients'][0]['orbit'])
    assert len(pure_d7) == 2
    half_d7 = {c['pair_orbit'] for r in sealed['face_rows'] if r['class_index'] in (785, 860)
               for c in r['corners'] if c['connected'] and c['type']=='D7'}
    assert len(half_d7) == 2 and not half_d7.intersection(pure_d7)
    extra = {o['key']: F(o['weight']) for r in sealed['component_rows'] for o in r['pair_orbits']}

    def weight(H, a, b, a1):
        pk = T.keystr(T.tree_key(H, (a, b)))
        bonus = extra.get(pk, F(0))
        if T.typename(H) == 'D7' and not a1 and pk not in pure_d7:
            bonus = F(1, 2) if pk in half_d7 else F(0)
        return F(T.baseline(H, a, b)) + bonus

    face_rows = []
    for ci, rec in sorted(selected.items()):
        row = unions['labelled_union_supports'][rec['union_index']]
        g = T.edges(row)
        corners = []
        for a, b in product((0, 2), (1, 3)):
            H = T.graph(row['white'] + [a, b], g)
            components = forest(H)
            a1 = any(len(C) == 1 for C in components)
            C = next(C for C in components if a in C)
            vals = [weight(C, a, b, a1 or outside_a1) if b in C else F(0)
                    for outside_a1 in (False, True)]
            corners.append(dict(marked=[a, b], core_forest=[T.typename(J) for J in components],
                                core_has_a1=a1, weights=list(map(str, vals))))
        sums = [sum(F(c['weights'][i]) for c in corners) for i in (0, 1)]
        face_rows.append(dict(class_index=ci, union_index=rec['union_index'],
                              sums_without_with_external_a1=list(map(str, sums)), corners=corners))
    failures = [r for r in face_rows if min(map(F, r['sums_without_with_external_a1'])) < 1]
    if failures:
        print('unpaid', json.dumps(failures), flush=True)
    assert not failures

    component_rows = []
    totals = [[], []]
    for name, H in U.models(29):
        row = dict(type=name, rank=len(H), pair_weights=[])
        sums = [F(0), F(0)]
        for a, b in combinations(H, 2):
            vals = [weight(H, a, b, bool(flag)) for flag in (0, 1)]
            sums = [x+y for x, y in zip(sums, vals)]
            row['pair_weights'].append(dict(marked=[a, b], weights=list(map(str, vals))))
        row['totals_without_with_a1'] = list(map(str, sums))
        component_rows.append(row)
        for i in (0, 1):
            totals[i].append(sums[i])
    no_a1, no_a1_forests = dynamic(component_rows, totals[0], 29, False)
    with_a1_tail, with_a1_forests = dynamic(component_rows, totals[1], 28, True)
    maximum = max(no_a1[29], with_a1_tail[28])
    assert maximum == 84
    assert with_a1_tail[28] == 84
    assert with_a1_forests[28] == ('D7', 'D7', 'D7', 'D7')
    d7 = next(r for r in component_rows if r['type']=='D7')
    assert d7['totals_without_with_a1'] == ['16', '21']
    # The rule only increases the previously checked distance-two/H4 triangle
    # rule, so all triangle extensions and all short-paid quads stay paid.
    assert all(r['h4_only_weight'] >= 0 for f in triangle['markings'] for r in f['corners'])
    assert all(sum(r['h4_only_weight'] for r in f['corners']) >= 2 for f in triangle['markings'])
    files = [f'scripts/{TAG}.py', 'scripts/fresh_sqrt5_uniform_face_weights_20261015.py',
             'scripts/fresh_sqrt5_sealed_weights_20261015.py',
             'results/fresh_sqrt5_sealed_weights_20261015.json',
             'results/fresh_sqrt5_quad_exhaustion_20261015.json',
             'results/fresh_sqrt5_quad_contacts_20261015.json',
             'results/fresh_sqrt5_quad_contacts_20261015_review.json',
             'results/fresh_sqrt5_triangle_weights_20261015.json',
             'results/fresh_sqrt5_triangle_weights_20261015_review.json']
    return dict(status='PASS_EXACT_A1_SENSITIVE_PARTIAL_FACE_RULE_AND_RANK29_VERTEX_BUDGET',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
                rule='Keep the saved nonnegative pair bonuses. For entireD7 without an entireA1 in the full vertex forest, retain bonus1 on the535/563 marked-pair orbits, bonus1/2 on the785/860 orbits, and zero other bonuses. Keep all saved D7 bonuses if an entireA1 exists. Base rule is distance<=2 plus the established H4 bonuses.',
                pure_d7_orbits=sorted(pure_d7), half_d7_orbits=sorted(half_d7), sealed_quad_count=112, sealed_quad_cases=224,
                sealed_quad_rows=face_rows, component_rows=component_rows,
                rank29_without_a1=dict(maximum=str(no_a1[29]), forest=list(no_a1_forests[29])),
                rank29_with_a1=dict(maximum=str(with_a1_tail[28]), forest=['A1']+list(with_a1_forests[28])),
                rank29_maximum=str(maximum), rank29_target=84,
                no_a1_rank_maxima=[None if v is None else str(v) for v in no_a1],
                with_a1_tail_rank_maxima=list(map(str, with_a1_tail)),
                sealing_argument='For each of the112 cores, every first external white node touching the core either fails a necessary arithmetic/finite-corner/compactness test or creates an existing short-path or entireH4 payment. Thus, unless already paid by the baseline rule, every external white node is orthogonal to the entire core. The corner forest is the displayed core corner forest plus the same external white forest at every corner. This gives exactly the two tested cases according to whether that external forest contains an entireA1. Nonnegative bonuses preserve all baseline payments.',
                covered='All triangle faces using the previously reviewed arbitrary-extension H4-only triangle rule; quadrilaterals paid by a short marked path or an entireH4; and the112 sealed singular finite principal classes. The44 unpaid-sealed nonsingular cores and2 certified affine dotted cores remain excluded by the inherited contact argument.',
                unresolved='The130 nonsingular finite unpaid classes with458 accepted single-white profiles and24 symbolic dotted unpaid classes with369 profiles. Their arbitrary external-white completions must satisfy this same weight rule before the rank29 face-count contradiction follows. No new dimension bound is established.',
                scope='Partial arithmetic face-weight theorem, pending independent review of this new forest-sensitive step. Exact rank29 vertex budget succeeds; complete quadrilateral coverage remains open. The old component-only barrier is bypassed by dependence on the entire vertex forest.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    out = run()
    path = ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], out['rank29_without_a1'], out['rank29_with_a1'])
