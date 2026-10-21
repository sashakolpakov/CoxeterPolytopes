#!/usr/bin/env python3
"""Bind the independently checked case31 metric transfers before cursor continuation."""
import argparse
import json
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_tnext_20261015_case31_gates'
OLD = 'results/fresh_h9_tnext_20261015_gates_case_31.json'
NEW = 'results/fresh_h9_tnext_20261015_case31_transfers.json'
REVIEW = 'results/fresh_h9_tnext_20261015_case31_transfers_review.json'
ENGINE = 'scripts/fresh_h9_tnext_20261015_search.py'


def run():
    def sha(p):
        return sha256((ROOT / p).read_bytes()).hexdigest()
    seen = set()
    def read(p):
        d = json.loads((ROOT / p).read_text())
        if p not in seen:
            seen.add(p)
            for q, h in d.get('sources', {}).items():
                assert sha(q) == h, q
                if q.endswith('.json'):
                    read(q)
        return d
    old, new, review = map(read, (OLD, NEW, REVIEW))
    assert old['total_predicates'] == 153
    assert new['status'] == 'PASS_LITERAL_T31_TWO_HIGH_METRIC_TRANSFERS'
    assert review['status'].startswith('PASS_') and review['sources'][NEW] == sha(NEW)
    families = old['families'] + new['new_conditions']
    keys = {(tuple(map(tuple, f['fixed'])), tuple(sorted(map(tuple, f['bounds'])))) for f in families}
    assert len(families) == len(keys) == 155
    paths = [OLD, NEW, REVIEW, ENGINE, f'scripts/{TAG}.py']
    return dict(status='PASS_TNEXT_LITERAL_REVIEWED_CONDITIONS', case=31,
                sources={p: sha(p) for p in paths},
                inherited_predicates=153, total_predicates=155, families=families,
                scope='The original153 individually checked T31 predicates plus two reviewed '
                      'all-label negative-plane transfers. The saved first survivor is excluded; '
                      'whole-incidence status belongs to the separate resumed traversal.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    d = run()
    out = ROOT / f'results/{TAG}.json'
    if a.verify:
        assert d == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(d, indent=2) + '\n')
    print(d['status'], d['total_predicates'])
