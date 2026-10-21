#!/usr/bin/env python3
"""Close the literal high-block orientation of the reviewed sparse theorem."""
import argparse
import json
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_t8_sparse_search_20261014_twin'
PROOF = 'results/fresh_h9_t8_sparse_search_20261014_all_blocks.json'
REVIEW = 'results/fresh_h9_t8_all_blocks_20261014_review.json'
GEOM = 'results/fresh_h9_tree_cube_wedges_20261014.json'
STATE = 'results/fresh_h9_t8_sparse_search_20261014_final_case_08.json'

def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()

def run():
    data = [json.loads((ROOT / p).read_text()) for p in (PROOF, REVIEW, GEOM, STATE)]
    for obj in data:
        for p, h in obj['sources'].items():
            assert sha(p) == h, p
    proof, review, geom, state = data
    assert review['status'] == 'PASS_INDEPENDENT_T8_ALL_INTERNAL_BLOCK_LABELS_EXTENSION'
    g = geom['representatives'][8]['geometry']
    perm = list(range(18))
    perm[1], perm[2] = 2, 1
    def mapped(s):
        return sorted(perm[i] for i in s)
    vertices = {tuple(sorted(v['incident'])) for v in g['vertices']}
    missing = {tuple(m) for m in g['minimal_nonfaces']}
    assert {tuple(mapped(v)) for v in vertices} == vertices
    assert {tuple(mapped(m)) for m in missing} == missing
    fixed = sorted([mapped((a, b)) + [m] for a, b, m in proof['fixed_labels']])
    bounds = sorted([mapped((a, b)) + [m] for a, b, m in proof['finite_label_lower_bounds']])
    labels = {(a, b): m for a, b, m in state['solution']['pair_labels']}
    assert all(labels[a, b] == m for a, b, m in fixed)
    assert all(labels[a, b] >= m for a, b, m in bounds)
    paths = [PROOF, REVIEW, GEOM, STATE, f'scripts/{TAG}.py']
    return dict(status='PASS_LITERAL_TWIN_TRANSFER_AND_SAVED_STATE_EXCLUSION',
        sources={p: sha(p) for p in paths}, permutation=perm,
        complete_vertex_count=len(vertices), complete_missing_face_count=len(missing),
        fixed_labels=fixed, finite_label_lower_bounds=bounds,
        saved_state_excluded=True,
        scope='Exact incidence automorphism swapping roots1 and2 transfers the entire '
              'independently reviewed all-block sparse theorem, whose other labels are free. '
              'All462 vertices and21 minimal nonfaces are checked, not just selected supports. '
              'No new arithmetic or metric search. This closes both high-root orientations '
              'of that theorem, not the whole T8 incidence.')

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    d = run()
    out = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert d == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(d, indent=2) + '\n')
    print(d['status'], d['complete_vertex_count'], d['complete_missing_face_count'])
