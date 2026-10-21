#!/usr/bin/env python3
"""Bind the independently reviewed four three-large-block predicates."""
import argparse
import json
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_star_large_block_resume_20261015_gates'
OLD = 'results/fresh_h9_star_pi4_resume_20261015_gates.json'
OLD_REVIEW = 'results/fresh_h9_star_pi4_resume_20261015_gates_review.json'
NEW = 'results/fresh_h9_star_three_large_blocks_20261015.json'
REVIEW = 'results/fresh_h9_star_three_large_blocks_20261015_review.json'
BATCH = 'results/fresh_h9_star_pi4_resume_20261015_search_batch_015.json'


def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()


def run():
    checked = {}

    def read(p):
        if p in checked:
            return checked[p]
        d = json.loads((ROOT / p).read_text())
        checked[p] = d
        for q, h in d.get('sources', {}).items():
            assert sha(q) == h, q
            if q.endswith('.json'):
                read(q)
        return d

    old, old_review, new, review, batch = [read(p) for p in (OLD, OLD_REVIEW, NEW, REVIEW, BATCH)]
    assert old['status'] == 'PASS_REVIEWED_PI4_THIRD_CORE_TRANSFERS'
    assert new['status'] == 'PASS_FOUR_THREE_LARGE_BLOCK_LORENTZ_OBSTRUCTIONS'
    assert old_review['status'].startswith('PASS_') and old_review['sources'][OLD] == sha(OLD)
    assert review['status'].startswith('PASS_') and review['sources'][NEW] == sha(NEW)
    assert old['case'] == new['case'] == batch['case'] == 0
    assert len(old['families']) == 493 and len(new['families']) == 4
    assert new['path_guard_occurrences'] == 96 and new['own_Lanner_occurrences'] == 12
    assert new['actual_leaf_triangles'] == 4
    families = old['families'] + new['families']
    keys = {tuple(tuple(sorted(map(tuple, f[k]))) for k in ('fixed', 'bounds')) for f in families}
    assert len(keys) == len(families) == 497
    blocks = [(3, 4), (6, 7), (9, 10), (12, 13)]
    assert {tuple(tuple(r) for r in f['bounds']) for f in new['families']} == {
        tuple(tuple(b) + (6,) for j, b in enumerate(blocks) if j != i) for i in range(4)}
    assert all(f['fixed'] == [] for f in new['families'])
    assert batch['status'] == 'FULL_MANDATORY_BRIDGE_SURVIVOR'
    assert batch['cumulative_additional_nodes'] == 119877
    labels = {tuple(r[:2]): r[2] for r in batch['solution']['pair_labels']}

    def hit(f):
        return all(labels[a, b] == m for a, b, m in f['fixed']) and all(labels[a, b] >= m for a, b, m in f['bounds'])

    assert not any(hit(f) for f in old['families'])
    hits = [f['id'] for f in new['families'] if hit(f)]
    assert len(hits) == 4
    paths = [f'scripts/{TAG}.py', OLD, OLD_REVIEW, NEW, REVIEW, BATCH]
    return dict(status='PASS_REVIEWED_THREE_LARGE_BLOCK_BINDER', case=0,
                sources={p: sha(p) for p in paths}, families=families,
                inherited_predicates=493, new_predicates=4, total_predicates=497,
                saved_cursor=119877, saved_survivor_new_gate_hits=hits,
                proof='The four exact placements were independently reconstructed against the full literal incidence. '
                      'No additional symmetry or contact hypothesis is used. Three internal labels at least6 are forbidden; '
                      'the unbounded category7 remains all finite labels at least7.',
                scope='497 compact-star0 predicates; no search or whole-incidence exclusion in this binder.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    out = run()
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out, indent=2) + '\n')
    print(out['status'], out['total_predicates'])
