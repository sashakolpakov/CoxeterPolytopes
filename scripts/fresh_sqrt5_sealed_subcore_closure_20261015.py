#!/usr/bin/env python3
"""Find reviewed sealed subcores in every accepted finite first extension.

All white subsets are searched; this finite containment reduction is not an
arbitrary-depth attachment closure.
"""
import argparse
import json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_all_quad_unions_20261015_classes as C
import fresh_sqrt5_quad_contacts_20261015 as A

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_sealed_subcore_closure_20261015'


def run():
    contacts = json.loads((ROOT/'results/fresh_sqrt5_quad_contacts_20261015.json').read_text())
    base = json.loads((ROOT/'results/fresh_sqrt5_quad_exhaustion_20261015.json').read_text())
    lookup = {}
    sizes = set()
    for rec in contacts['records']:
        if rec['status']=='UNPAID_FINITE_CONTACT_DOMAIN' and rec['singular'] and not rec['accepted']:
            row = base['labelled_union_supports'][rec['union_index']]
            key, order = C.canonical(row)
            assert key not in lookup
            lookup[key] = (rec['class_index'], rec['union_index'], row, order)
            sizes.add(len(row['white']))
    assert len(lookup) == 112
    records = []
    counts = Counter()
    subsets = 0
    for rec in contacts['records']:
        if rec['status']!='UNPAID_FINITE_CONTACT_DOMAIN' or rec['singular'] or not rec['accepted']:
            continue
        source = base['labelled_union_supports'][rec['union_index']]
        for ti, trial in enumerate(rec['trials']):
            if trial['status']!='PASS_EXTRA_WHITE_NECESSARY_TESTS':
                continue
            row = A.extension(source, trial['contacts'])
            g = A.edges(row)
            matches = []
            singular = trial['physical']['determinant'] == ['0', '0']
            for size in sorted(sizes):
                for white in combinations(row['white'], size):
                    ids = list(range(4)) + list(white)
                    sub = dict(nodes=len(ids), white=list(range(4, len(ids))),
                               labels=[g[tuple(sorted((a, b)))] for a, b in combinations(ids, 2)])
                    key, order = C.canonical(sub)
                    subsets += 1
                    if key not in lookup:
                        continue
                    ci, ui, target, target_order = lookup[key]
                    mp = [None] * len(ids)
                    for a, b in zip(order, target_order):
                        mp[a] = b
                    C.check_map(sub, target, mp)
                    matches.append(dict(class_index=ci, union_index=ui, nodes_in_extension=ids,
                                        to_sealed_representative=mp))
            assert not matches or singular
            status = 'PAID_BY_REVIEWED_SEALED_SUBCORE' if matches else 'OPEN_SINGULAR' if singular else 'OPEN_NONSINGULAR'
            counts[status] += 1
            records.append(dict(source_class_index=rec['class_index'], source_union_index=rec['union_index'],
                                trial_index=ti, contacts=trial['contacts'], singular=singular,
                                status=status, sealed_subcores=matches))
    assert len(records) == 458
    assert sum(r['singular'] for r in records) == 147
    files = [f'scripts/{TAG}.py', 'scripts/fresh_sqrt5_all_quad_unions_20261015_classes.py',
             'scripts/fresh_sqrt5_quad_contacts_20261015.py',
             'results/fresh_sqrt5_quad_contacts_20261015.json',
             'results/fresh_sqrt5_quad_exhaustion_20261015.json',
             'results/fresh_sqrt5_forest_sensitive_weights_20261015.json']
    return dict(status='PASS_EXACT_FIRST_EXTENSION_SEALED_SUBCORE_CONTAINMENT',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
                extensions=len(records), white_subsets_checked=subsets,
                counts=dict(sorted(counts.items())), records=records,
                scope='Exact complete marked induced-subdiagram search among the458 accepted one-white extensions of130 nonsingular finite unpaid cores. Any found sealed subcore inherits its arbitrary-extension payment from the forest-sensitive rule. Unmatched extensions and all24 symbolic dotted families remain open. This is not a search over arbitrary collections of external white vertices and establishes no new dimension bound.')


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
    print(out['status'], out['counts'], out['white_subsets_checked'])
