#!/usr/bin/env python3
"""Transfer a reviewed proof only after checking every literal face premise.

No global incidence automorphism is asserted or used.  The four maps merely
rename the variables of the same negative-plane proof; all its incidence
premises are reconstructed in the fixed T8 geometry.
"""
import argparse
import json
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_t8_sparse_search_20261014_bridge_variants'
BASE = 'results/fresh_h9_t8_sparse_search_20261014_two_high_blocks.json'
REVIEW = 'results/fresh_h9_t8_sparse_search_20261014_two_high_blocks_review.json'
GEOM = 'results/fresh_h9_tree_cube_wedges_20261014.json'

def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()

def bound(p, seen=None):
    if seen is None:
        seen = set()
    obj = json.loads((ROOT / p).read_text())
    if p in seen:
        return obj
    seen.add(p)
    for q, h in obj.get('sources', {}).items():
        assert sha(q) == h, q
        if q.endswith('.json'):
            bound(q, seen)
    return obj

def run():
    base, review, data = [bound(p) for p in (BASE, REVIEW, GEOM)]
    assert base['status'] == 'PASS_T8_TWO_HIGH_BLOCKS_ALL_OTHER_LABELS_FREE_EXCLUDED'
    assert 'PASS' in review['status']
    geo = data['representatives'][8]['geometry']
    vertices = [set(v['incident']) for v in geo['vertices']]
    missing = {tuple(sorted(m)) for m in geo['minimal_nonfaces']}
    def actual(s):
        return next(i for i, v in enumerate(vertices) if set(s) <= v)
    records = []
    for bridge in (14, 15):
        for anchor in (11, 17):
            mapping = list(range(18))
            if bridge == 15:
                for a, b in ((4, 7), (5, 8), (14, 15)):
                    mapping[a], mapping[b] = mapping[b], mapping[a]
            if anchor == 17:
                mapping[11], mapping[17] = mapping[17], mapping[11]
            assert sorted(mapping) == list(range(18))
            def mapped(s):
                return [mapping[x] for x in s]
            fields = {}
            for key, support_key in (
                ('all45_core_meeting_pairs', 'pair'),
                ('cross_right_guards', 'support'),
                ('anchor_high_right_budgets', 'support'),
                ('third_other_blocks_right_budgets', 'support'),
                ('third_anchor_upper_stars', 'support'),
            ):
                rows = []
                for index, row in enumerate(base[key]):
                    s = mapped(row[support_key])
                    rows.append(dict(base_row=index, support=s,
                                     actual_vertex=actual(s)))
                fields[key] = rows
            triples = [sorted(mapped(t)) for t in base['missing_triples']]
            assert all(tuple(t) in missing for t in triples)
            support = mapped(base['positive_support'])
            bounds = [sorted(mapped((a, b))) + [m]
                      for a, b, m in base['finite_label_lower_bounds']]
            records.append(dict(bridge=bridge, anchor=anchor, mapping=mapping,
                finite_label_lower_bounds=sorted(bounds), fixed_labels=[],
                positive_support=support, positive_support_actual_vertex=actual(support),
                missing_triples=triples, incidence_checks=fields))
    assert len(records) == 4
    distinct_bounds = sorted({tuple(tuple(b) for b in r['finite_label_lower_bounds'])
                             for r in records})
    assert len(distinct_bounds) == 2
    sources = [BASE, REVIEW, GEOM, f'scripts/{TAG}.py']
    return dict(status='PASS_FOUR_LITERAL_TWO_HIGH_BRIDGE_PROOF_TRANSFERS',
        sources={p: sha(p) for p in sources}, variants=records,
        distinct_condition_sets=distinct_bounds,
        checks_per_variant=dict(meeting_pairs=45, cross_right_guards=24,
            anchor_budgets=2, third_budgets=4, upper_stars=4,
            positive_support=1, missing_triples=3),
        scope='Four literal relabelings of one independently reviewed negative-plane proof. '
              'Every selected actual face and missing triple is checked anew. No incidence '
              'automorphism is assumed; arithmetic and all-label inequalities are inherited '
              'from the bound primary and independent certificates. The two distinct sets '
              'of four internal label inequalities each exclude compact H9 realizations '
              'of T8 with all other labels/distances free. This is not a whole-incidence exclusion.')

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = json.loads(json.dumps(run()))
    out = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], len(result['variants']))
