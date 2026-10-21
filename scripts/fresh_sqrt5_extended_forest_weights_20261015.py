#!/usr/bin/env python3
"""Extend the A1-sensitive rule to the46 newly sealed singular stars."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import fresh_sqrt5_sealed_weights_20261015 as T
import fresh_sqrt5_uniform_face_weights_20261015 as U
import fresh_sqrt5_forest_sensitive_weights_20261015 as P

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_extended_forest_weights_20261015'


def run():
    original = json.loads((ROOT/'results/fresh_sqrt5_forest_sensitive_weights_20261015.json').read_text())
    frontier = json.loads((ROOT/'results/fresh_sqrt5_singular_extension_frontier_20261015.json').read_text())
    unions = json.loads((ROOT/'results/fresh_sqrt5_quad_exhaustion_20261015.json').read_text())
    models = U.models(29)
    weights = {}
    old_components = {r['type']: r for r in original['component_rows']}
    for name, H in models:
        for r in old_components[name]['pair_weights']:
            a, b = r['marked']
            key = T.keystr(T.tree_key(H, (a, b)))
            value = tuple(map(F, r['weights']))
            assert key not in weights or weights[key] == value
            weights[key] = value
    sealed = [r for r in frontier['records'] if not r['accepted']]
    unpaid = [r for r in sealed if not r['payment']['paid']]
    assert len(sealed) == 46 and len(unpaid) == 14
    upgraded = {}
    for r in unpaid:
        for c in r['payment']['corners']:
            if c['pair_orbit'] is not None:
                key = c['pair_orbit']
                assert key in weights
                upgraded[key] = dict(type=c['marked_component'], previous=list(map(str, weights[key])))
                weights[key] = tuple(max(v, F(1)) for v in weights[key])
    face_rows = []
    inputs = [('old', r['class_index'], unions['labelled_union_supports'][r['union_index']])
              for r in original['sealed_quad_rows']]
    inputs += [('new', r['index'], r['row']) for r in sealed]
    for family, ident, row in inputs:
        g = T.edges(row)
        corners = []
        for a, b in product((0, 2), (1, 3)):
            H = T.graph(row['white'] + [a, b], g)
            components = P.forest(H)
            has_a1 = any(len(C)==1 for C in components)
            C = next(C for C in components if a in C)
            if b in C:
                key = T.keystr(T.tree_key(C, (a, b)))
                vals = [weights[key][int(has_a1)], weights[key][1]]
            else:
                vals = [F(0), F(0)]
            corners.append(dict(marked=[a, b], core_has_a1=has_a1, weights=list(map(str, vals))))
        sums = [sum(F(c['weights'][i]) for c in corners) for i in (0, 1)]
        assert min(sums) >= 1, (family, ident, sums)
        face_rows.append(dict(family=family, index=ident, corners=corners,
                              sums_without_with_external_a1=list(map(str, sums))))
    totals = [[], []]
    components = []
    for name, H in models:
        sums = [F(0), F(0)]
        for a, b in combinations(H, 2):
            vals = weights[T.keystr(T.tree_key(H, (a, b)))]
            sums = [x+y for x, y in zip(sums, vals)]
        components.append(dict(type=name, rank=len(H), totals=list(map(str, sums))))
        for i in (0, 1):
            totals[i].append(sums[i])
    no_a1, no_forests = P.dynamic(components, totals[0], 29, False)
    with_a1, with_forests = P.dynamic(components, totals[1], 28, True)
    maximum = max(no_a1[29], with_a1[28])
    assert maximum == 84, (no_a1[29], with_a1[28])
    files = [f'scripts/{TAG}.py', 'scripts/fresh_sqrt5_sealed_weights_20261015.py',
             'scripts/fresh_sqrt5_uniform_face_weights_20261015.py',
             'scripts/fresh_sqrt5_forest_sensitive_weights_20261015.py',
             'results/fresh_sqrt5_forest_sensitive_weights_20261015.json',
             'results/fresh_sqrt5_forest_sensitive_weights_20261015_review.json',
             'results/fresh_sqrt5_singular_extension_frontier_20261015.json',
             'results/fresh_sqrt5_quad_exhaustion_20261015.json']
    return dict(status='PASS_EXACT_EXTENDED_FOREST_RULE_WITH_RANK29_BUDGET84',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
                upgraded_orbits=[dict(key=k, **r, new=list(map(str, weights[k]))) for k, r in sorted(upgraded.items())],
                old_sealed_cores=112, new_sealed_cores=46, face_cases=len(face_rows)*2,
                face_rows=face_rows, component_rows=components,
                rank29_without_a1=dict(maximum=str(no_a1[29]), forest=list(no_forests[29])),
                rank29_with_a1=dict(maximum=str(with_a1[28]), forest=['A1']+list(with_forests[28])),
                rank29_maximum=str(maximum), no_a1_rank_maxima=[None if v is None else str(v) for v in no_a1],
                with_a1_tail_rank_maxima=list(map(str, with_a1)),
                rule='Starting with the reviewed A1-sensitive rule, raise each connected marked-pair orbit occurring in one of the14 newly unpaid sealed faces to weight at least1, in both A1 contexts. All other weights are unchanged.',
                scope='Pays all112 old and46 new sealed cores under arbitrary outside-white extension, conditional on their respective sealing certificates. All baseline triangle/short-path/H4 payments are preserved. Remaining nonsingular, nonsealed singular and dotted frontiers still prevent a dimension conclusion. This expanded rule is pending independent review.')


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
    print(out['status'], 'orbits', len(out['upgraded_orbits']), out['rank29_without_a1'], out['rank29_with_a1'])
