#!/usr/bin/env python3
"""Require the independent single-pi4-third proof and transport on twelve exact cores."""
import argparse
import json
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_star_pi4_resume_20261015_gates'
OLD = 'results/fresh_h9_star_single_third_opposite_20261015_gates.json'
NEW = 'results/fresh_h9_star_third_four_20261015.json'
REVIEW = 'results/fresh_h9_star_third_four_20261015_review.json'
MAPS = 'results/fresh_h9_star_core_transfers_20261014.json'
MAP_REVIEW = 'results/fresh_h9_star_core_transfers_20261014_review.json'
GEOM = 'results/fresh_h9_tree_cube_wedges_20261014.json'
BATCH = 'results/fresh_h9_star_single_third_opposite_20261015_search_batch_014.json'


def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def key(f):
    return tuple(tuple(sorted(map(tuple, f[k]))) for k in ('fixed', 'bounds'))


def run():
    seen = {}

    def read(path):
        if path in seen:
            return seen[path]
        d = json.loads((ROOT / path).read_text())
        seen[path] = d
        for p, h in d.get('sources', {}).items():
            assert sha(p) == h, p
            if p.endswith('.json'):
                read(p)
        return d

    old, new, review, maps, mr = [read(p) for p in (OLD, NEW, REVIEW, MAPS, MAP_REVIEW)]
    assert old['status'] == 'PASS_REVIEWED_OPPOSITE_SINGLE_THIRD_CORE_TRANSFERS'
    assert new['status'] == 'PASS_SINGLE_PI4_THIRD_SHARED_OR_OPPOSITE_ANCHORS'
    assert review['status'].startswith('PASS_') and review['sources'][NEW] == sha(NEW)
    assert maps['status'] == 'PASS_TWELVE_LITERAL_CORE_TRANSFERS'
    assert mr['status'] == 'PASS_INDEPENDENT_TWELVE_INDUCED_CORE_TRANSFERS'
    assert mr['sources'][MAPS] == sha(MAPS)
    assert old['case'] == new['case'] == maps['case'] == 0
    assert new['rank_core'] == maps['source_core']
    assert new['complete_maximal_core_faces'] == maps['source_maximal_faces']
    core = set(new['rank_core'])
    supports = [r['support'] for r in new['positive_kernel_guards'] + new['projection_budgets']]
    assert len(supports) == 114 and all(set(s) <= core for s in supports)
    assert all(any(set(s) <= set(f) for f in maps['source_maximal_faces']) for s in supports)
    assert len(new['families']) == 8
    for f in new['families']:
        assert len(f['fixed']) == 6 and len(f['bounds']) == 4
        assert all(set(r[:2]) <= core for k in ('fixed', 'bounds') for r in f[k])
    geom = read(GEOM)['representatives'][0]['geometry']
    vertices = [set(v['incident']) for v in geom['vertices']]
    minimal = list(map(set, geom['minimal_nonfaces']))
    families = old['families'][:]
    keys = {key(f) for f in families}
    assert len(families) == len(keys) == 397
    placements, images = [], []
    for row in maps['placements']:
        mp = dict(row['core_map'])
        assert set(mp) == core and len(set(mp.values())) == 11
        target = set(mp.values())
        induced = {tuple(sorted(v & target)) for v in vertices}
        maximal = sorted(f for f in induced if not any(set(f) < set(g) for g in induced))
        assert maximal == sorted(tuple(sorted(mp[x] for x in f)) for f in maps['source_maximal_faces'])
        assert [list(m) for m in sorted(tuple(sorted(m)) for m in minimal if m <= target)] == row['induced_minimal_nonfaces']
        assert all(any({mp[x] for x in pair} <= v for v in vertices)
                   for pair in combinations(core, 2))
        witness_rows = []
        for support in supports:
            image = sorted(mp[x] for x in support)
            witness_rows.append(dict(support=image, actual_vertex=next(i for i, v in enumerate(vertices)
                                                                      if set(image) <= v)))
        ids = []
        for f in new['families']:
            image = {k: sorted([sorted((mp[a], mp[b])) + [m] for a, b, m in f[k]])
                     for k in ('fixed', 'bounds')}
            ident = f'pi4_third_anchor{row["anchor_block"]}_owner{row["owner_block"]}_' + f['id']
            lifted = dict(id=ident, origin='reviewed_pi4_third_core_transfer',
                          source_family=f['id'], **image)
            assert key(lifted) not in keys
            keys.add(key(lifted))
            families.append(lifted)
            images.append(lifted)
            ids.append(ident)
        placements.append(dict(core_map=row['core_map'], anchor_block=row['anchor_block'],
                               owner_block=row['owner_block'], family_ids=ids,
                               new_domain_witnesses=witness_rows))
    assert len(images) == 96 and len(families) == 493
    assert len(placements) == 12 and sum(len(p['new_domain_witnesses']) for p in placements) == 1368
    batch = read(BATCH)
    labels = {tuple(r[:2]): r[2] for r in batch['solution']['pair_labels']}

    def hit(f):
        return all(labels[a, b] == m for a, b, m in f['fixed']) and \
            all(labels[a, b] >= m for a, b, m in f['bounds'])

    hits = [f['id'] for f in images if hit(f)]
    assert hits and not any(hit(f) for f in old['families'])
    paths = [OLD, NEW, REVIEW, MAPS, MAP_REVIEW, GEOM, BATCH, f'scripts/{TAG}.py']
    return dict(status='PASS_REVIEWED_PI4_THIRD_CORE_TRANSFERS', case=0,
                sources={p: sha(p) for p in paths}, families=families,
                inherited_predicates=397, new_predicates=96, total_predicates=493,
                placements=placements, new_witness_occurrences=1368,
                saved_cursor=batch['cumulative_additional_nodes'], saved_survivor_new_gate_hits=hits,
                proof='The eight independently checked shared/opposite anchor and owner/last-root orientations are transported '
                      'by all twelve previously independently checked full induced-core maps. '
                      'Each of the1368 new path/cycle/budget support occurrences is an actual face; '
                      'all55 core pairs meet and the complete induced face complex is rechecked. '
                      'No full-nerve automorphism is assumed for these twelve placements.',
                scope='493 reviewed star0 predicates.96 new single-pi4-third/shared-or-opposite-anchor families. '
                      'No whole-incidence exclusion and no search performed here.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    data = json.loads(json.dumps(run()))
    out = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(data, indent=2)+'\n')
    print(data['status'], data['total_predicates'], data['new_witness_occurrences'])
