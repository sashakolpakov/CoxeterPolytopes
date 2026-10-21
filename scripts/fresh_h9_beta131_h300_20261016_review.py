#!/usr/bin/env python3
"""Fresh height300 replay through the independent Golden census."""
import argparse
import json
import fresh_h9_beta131_pilot_20261016_review as I

TAG = 'fresh_h9_beta131_h300_20261016_review'


def run():
    I.TAG = TAG
    I.PRIMARY = 'results/fresh_h9_beta131_h300_20261016.json'
    result = I.run()
    old_path = 'results/fresh_h9_beta131_h100_20261016_review.json'
    old = json.loads((I.ROOT / old_path).read_text())
    assert result['height_bound'] == '300'
    assert result['levels'][:old['controls']] == old['levels']
    paths = [old_path, 'scripts/fresh_h9_beta131_pilot_20261016_review.py']
    result['sources'].update({p: I.sha(p) for p in paths})
    result['preserved_independent_prefix'] = dict(
        height=100, controls=old['controls'], walls=old['walls'])
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    path = I.ROOT / 'results' / f'{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['controls'], result['walls'], result['counts'])
