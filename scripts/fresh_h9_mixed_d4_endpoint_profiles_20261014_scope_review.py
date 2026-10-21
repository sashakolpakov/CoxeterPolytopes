#!/usr/bin/env python3
"""Narrow source/proof audit; not a second 10,753-profile arithmetic census."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = 'fresh_h9_mixed_d4_endpoint_profiles_20261014'
TAG = BASE + '_scope_review'


def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()


def run():
    data = json.loads((ROOT / ('results/' + BASE + '.json')).read_text())
    for p, h in data['sources'].items():
        assert sha(p) == h, p
    rows = list(data['profiles'].values())
    assert len(rows) == len({tuple(r['support']) for r in rows}) == 52
    accepted = sum(len(r['finite_profiles']) for r in rows)
    tails = sum(len(r['dihedral_tail_profiles']) for r in rows)
    trials = sum(r['finite_trial_count'] for r in rows)
    assert (accepted, tails, trials) == (1825, 72, 10753)
    for r in rows:
        assert r['finite_trial_count'] == len(r['finite_profiles']) + sum(r['rejected_trial_counts'].values())
        assert all(len(p['contacts']) <= 3 for p in r['finite_profiles'])
        assert all(len(p['contacts']) != 3 or all(m == 3 for _, m in p['contacts']) for p in r['finite_profiles'])
    # Independent elementary principal-minor bounds, with roots of square2.
    four_star_form = 2 * (2**2 + 4) - 2 * 4 * 2
    heavy_three_star_schur = F(2) - F(2 + 1 + 1, 2)
    high_path_schur = F(2) - F(3 + 1, 2)
    assert four_star_form == heavy_three_star_schur == high_path_schur == 0
    paths = ['scripts/' + BASE + '.py', 'results/' + BASE + '.json', 'scripts/' + TAG + '.py']
    return dict(status='PASS_NARROW_SOURCE_PROOF_AND_METADATA_AUDIT',
                sources={p: sha(p) for p in paths},
                counts=dict(edges=52, finite_trials=trials, accepted_profiles=accepted, isolated_dihedral_tails=tails),
                independent_checks=dict(four_contact_star_form=four_star_form,
                    three_contact_heavy_schur=str(heavy_three_star_schur),
                    label_at_least_six_path_schur=str(high_path_schur)),
                proof_review=['cycle and star completeness', 'isolated arbitrary-label I2 tails',
                    'block diagonal inverse Schur cost', 'oriented Lorentz frame and exponential parameter',
                    'strict positive parameter for a genuine added facet', 'finite-endpoint-only scope'],
                limitations='Source/proof audit plus independent elementary bounds and saved-data checks. No second exact support-inverse/profile enumeration and no new global cap feasibility test.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    out = ROOT / ('results/' + TAG + '.json')
    if args.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['counts'])
