#!/usr/bin/env python3
"""Compose the dense-third partition and transport it on twelve literal cores."""
import argparse
import json
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_star_third34_20261014_gates'
OLD = 'results/fresh_h9_star_core_transfers_20261014_gates.json'
TRANSFER = 'results/fresh_h9_star_core_transfers_20261014.json'
TRANSFER_REVIEW = 'results/fresh_h9_star_core_transfers_20261014_review.json'
GEOM = 'results/fresh_h9_tree_cube_wedges_20261014.json'
NEW = 'results/fresh_h9_star_third34_20261014.json'
NEW_REVIEW = 'results/fresh_h9_star_third34_20261014_review.json'
OLD_THEOREMS = [
    'results/fresh_h9_star_mixed_singles_20261014.json',
    'results/fresh_h9_star_opposite_core_20261014.json',
]
CORE = {3, 4, 6, 7, 9, 10, 12, 13, 11, 14, 17}


def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def key(family):
    return tuple(tuple(sorted(map(tuple, family[k]))) for k in ('fixed', 'bounds'))


def applies(family, labels):
    return (all(labels[a, b] == m for a, b, m in family['fixed'])
            and all(labels[a, b] >= m for a, b, m in family['bounds']))


def run():
    seen = {}

    def read(path):
        if path in seen:
            return seen[path]
        data = json.loads((ROOT / path).read_text())
        seen[path] = data
        for p, digest in data.get('sources', {}).items():
            assert sha(p) == digest, p
            if p.endswith('.json'):
                read(p)
        return data

    old, transfer, tr, new, nr = [read(p) for p in
                                 (OLD, TRANSFER, TRANSFER_REVIEW, NEW, NEW_REVIEW)]
    assert old['status'] == 'PASS_REVIEWED_LOCAL_CORE_GATE_COMPOSITION'
    assert transfer['status'] == 'PASS_TWELVE_LITERAL_CORE_TRANSFERS'
    assert tr['status'] == 'PASS_INDEPENDENT_TWELVE_INDUCED_CORE_TRANSFERS'
    assert tr['sources'][TRANSFER] == sha(TRANSFER)
    assert new['status'] == 'PASS_ALL_LABEL_DENSE_THIRD_SINGLE_ANCHOR_CORE'
    assert nr['status'].startswith('PASS_') and nr['sources'][NEW] == sha(NEW)
    assert old['case'] == transfer['case'] == new['case'] == 0
    assert len(old['families']) == 253
    assert set(new['rank_core']) == CORE
    assert set(transfer['source_core']) == CORE
    assert new['complete_maximal_core_faces'] == transfer['source_maximal_faces']

    predecessors = []
    paths = [OLD, TRANSFER, TRANSFER_REVIEW, NEW, NEW_REVIEW, GEOM]
    for p in OLD_THEOREMS:
        q = p[:-5] + '_review.json'
        data, review = read(p), read(q)
        assert data['status'].startswith('PASS_')
        assert review['status'].startswith('PASS_') and review['sources'][p] == sha(p)
        assert set(data['rank_core']) == CORE
        assert sorted(map(tuple, data['complete_maximal_core_faces'])) == sorted(map(tuple, transfer['source_maximal_faces']))
        predecessors.extend((p, f) for f in data['families'])
        paths.extend((p, q))
    assert len(predecessors) == 4

    # The partition uses the integer alternatives (3,3), (>=4,3),
    # (3,>=4), (>=4,>=4). No upper bound is introduced.
    source_families, partitions = [], []
    for p, f in predecessors:
        assert len(f['fixed']) == 6 and len(f['bounds']) == 4
        assert [3, 14, 3] in f['fixed'] and [4, 14, 3] in f['fixed']
        fixed = [r for r in f['fixed'] if r[:2] not in ([3, 14], [4, 14])]
        assert len(fixed) == 4
        matching = [r for r in new['families'] if r['fixed'] == fixed]
        assert len(matching) == 2
        bounds = sorted(f['bounds'] + [[3, 14, 3], [4, 14, 3]])
        assert {tuple(map(tuple, r['bounds'])) for r in matching} == {
            tuple(map(tuple, sorted(f['bounds'] + [[3, 14, 3], [4, 14, 4]]))),
            tuple(map(tuple, sorted(f['bounds'] + [[3, 14, 4], [4, 14, 3]]))),
        }
        source = dict(id='dense_' + f['id'], origin='reviewed_dense_third_partition',
                      fixed=fixed, bounds=bounds)
        for a, b in product((3, 4), repeat=2):
            labels = {tuple(r[:2]): r[2] for r in source['fixed'] + source['bounds']}
            labels[3, 14], labels[4, 14] = a, b
            assert applies(source, labels)
            assert applies(f, labels) or any(applies(r, labels) for r in matching)
        source_families.append(source)
        partitions.append(dict(source_family=source['id'], exact_33_theorem=p,
                               exact_33_family=f['id'], larger_own_theorem=NEW,
                               larger_own_families=[r['id'] for r in matching],
                               integer_partition=['3,3', '>=4,3', '3,>=4', '>=4,>=4']))
    assert len({key(f) for f in source_families}) == 4

    # Recheck that every newly needed domain witness is entirely in the
    # certified induced core; inherited theorem witnesses are already bound
    # by the independent transfer review.
    new_supports = ([r['support'] for r in new['positive_kernel_guards']]
                    + [r['support'] for r in new['projection_budgets']])
    assert len(new_supports) == 34
    assert all(set(s) <= CORE for s in new_supports)
    geom = read(GEOM)['representatives'][0]['geometry']
    vertices = [set(v['incident']) for v in geom['vertices']]
    maximal = list(map(set, transfer['source_maximal_faces']))
    assert all(any(set(s) <= v for v in maximal) for s in new_supports)
    families = old['families'][:]
    keys = {key(f) for f in families}
    images, witnessed = [], []
    for placement in transfer['placements']:
        mp = dict(placement['core_map'])
        assert set(mp) == CORE and len(set(mp.values())) == 11
        assert all(any({mp[x] for x in pair} <= v for v in vertices)
                   for pair in combinations(CORE, 2))
        supports = [sorted(mp[x] for x in s) for s in new_supports]
        assert all(any(set(s) <= v for v in vertices) for s in supports)
        witnesses = [dict(support=s, actual_vertex=next(i for i, v in enumerate(vertices)
                                                      if set(s) <= v)) for s in supports]
        ids = []
        for f in source_families:
            image = {k: sorted([sorted((mp[a], mp[b])) + [m] for a, b, m in f[k]])
                     for k in ('fixed', 'bounds')}
            ident = (f'dense_anchor{placement["anchor_block"]}_'
                     f'owner{placement["owner_block"]}_' + f['id'])
            row = dict(id=ident, origin='reviewed_dense_core_transfer',
                       source_family=f['id'], **image)
            assert key(row) not in keys
            keys.add(key(row))
            families.append(row)
            images.append(row)
            ids.append(ident)
        witnessed.append(dict(anchor_block=placement['anchor_block'],
                              owner_block=placement['owner_block'],
                              core_map=placement['core_map'], family_ids=ids,
                              new_domain_witnesses=witnesses))
    assert len(images) == 48 and len(families) == 301 and len(witnessed) == 12
    batch_path = 'results/fresh_h9_star_core_transfers_20261014_search_batch_011.json'
    batch = read(batch_path)
    labels = {tuple(r[:2]): r[2] for r in batch['solution']['pair_labels']}
    hits = [f['id'] for f in images if applies(f, labels)]
    assert hits and not any(applies(f, labels) for f in old['families'])
    paths += [batch_path, f'scripts/{TAG}.py']
    return dict(status='PASS_REVIEWED_DENSE_THIRD_CORE_COMPOSITION', case=0,
                sources={p: sha(p) for p in paths}, families=families,
                inherited_predicates=253, new_predicates=48, total_predicates=301,
                source_families=source_families, partition=partitions,
                placements=witnessed, new_witness_occurrences=408,
                saved_cursor=batch['cumulative_additional_nodes'],
                saved_survivor_new_gate_hits=hits,
                proof='For finite integer own labels >=3, either both equal3 '
                      '(the independently reviewed shared/opposite core theorems), '
                      'or at least one is>=4 (the independent six-polynomial theorem). '
                      'All four single-anchor orientations are covered. The twelve '
                      'independently reviewed full induced-core maps transport this '
                      'union; all408 new domain-witness occurrences are actual faces. '
                      'This is a local transfer, not a full-nerve symmetry claim.',
                scope='301 reviewed predicates for star0; 48 new broad dense-third '
                      'families. No whole-incidence exclusion and no search performed '
                      'by this artifact.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    data = json.loads(json.dumps(run()))
    out = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(data, indent=2) + '\n')
    print(data['status'], data['total_predicates'], data['new_witness_occurrences'])
