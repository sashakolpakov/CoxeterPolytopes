#!/usr/bin/env python3
"""Logical composition of the exhaustive raw fixed-source attachment branches."""
import argparse
import json
from hashlib import sha256
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_h8_raw_attachment_20261014'
NAMES = {
    'gap': 'fresh_h9_h8_facet_cover_20261014_gap_gate',
    'pair': 'fresh_h9_h8_corner_anchors_20261014_pair_gate',
    'pair_review': 'fresh_h9_h8_corner_anchors_20261014_pair_review',
    'radial': 'fresh_h9_h8_global_cap_gap_20261014_radial',
    'radial_review': 'fresh_h9_h8_global_cap_gap_20261014_radial_review',
    'polar': 'fresh_h9_h8_corner_anchors_20261014_polar_gate',
    'polar_review': 'fresh_h9_h8_corner_anchors_20261014_polar_review',
    'high': 'fresh_h9_h8_all_high_20261014',
    'high_review': 'fresh_h9_h8_all_high_20261014_review',
}
def sha(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    paths = {key: f'results/{name}.json' for key, name in NAMES.items()}
    docs = {key: json.loads((ROOT/path).read_text()) for key, path in paths.items()}
    pending, seen = list(paths.values()), set()
    while pending:
        p = pending.pop()
        if p in seen: continue
        seen.add(p); doc = json.loads((ROOT/p).read_text())
        for path, h in doc.get('sources', {}).items():
            assert sha(path) == h, path
            if path.endswith('.json'): pending.append(path)
    assert all(d['status'].startswith('PASS_') for d in docs.values())
    for key in ('pair', 'radial', 'polar', 'high'):
        assert docs[key+'_review']['sources'][paths[key]] == sha(paths[key])
    assert docs['gap']['ordinary_minimum_norm'] == ['-1', '1']
    assert docs['gap']['exceptional_minimum_norm'] == ['-1/4', '1/4']
    assert sorted(r['source_vertex'] for r in docs['gap']['exceptional_rays']) == [0, 32, 33]
    assert docs['pair']['status'] == 'PASS_ALL_LABEL_TWO_EXCEPTIONAL_CENTER_ANCHORS_EXCLUDED'
    assert docs['radial']['eliminated_low_old_walls'] == list(range(11))
    assert docs['polar']['source_anchor_vertices'] == [0, 32, 33]
    assert docs['polar']['status'] == 'PASS_ALL_LOW_POLAR11_SOURCE_ANCHORS_EXCLUDED'
    assert docs['high']['status'] == 'PASS_ALL_HIGH_FIXED_H8_RAW_ATTACHMENT_EXCLUDED'
    assert docs['high']['high_threshold'] == [['-1', '1'], ['0', '0']]
    assert docs['high']['required_dotted_pairs'] == [[4, 11], [5, 6], [8, 11]]
    eliminated_single_indices = docs['radial']['eliminated_low_old_walls']+[11]
    assert sorted(eliminated_single_indices) == list(range(12))
    sources = list(paths.values())+[f'scripts/{TAG}.py']
    return dict(status='PASS_ALL_REAL_RAW_FIXED_H8_ATTACHMENT_BRANCHES_EXCLUDED',
                sources={p: sha(p) for p in sources},
                branch_cover=[dict(low_center_count='at least two',excluded_by='pair'),
                              dict(low_center_count=1,old_indices=list(range(11)),excluded_by='radial'),
                              dict(low_center_count=1,old_indices=[11],excluded_by='polar'),
                              dict(low_center_count=0,excluded_by='high')],
                theorem='For every real all-dotted cap direction, the frozen mixed twelve-wall H9 chamber cannot be completed to finite volume by the fixed Bugaenko eleven-facet compact H8 source attached with perpendicular source-side walls alone.',
                scope='Logical composition, not a new independent replay of the source-gap census. All new pair, radial and projection obstructions have separate arithmetic reviews. Additional nonperpendicular covering walls, changed retained walls and other source polytopes are outside scope; no general H9 nonexistence statement.')
if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--verify', action='store_true'); args = parser.parse_args()
    result = run(); out = ROOT/f'results/{TAG}.json'
    if args.verify: assert result == json.loads(out.read_text())
    else: out.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'])
