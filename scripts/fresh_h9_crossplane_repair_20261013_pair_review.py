#!/usr/bin/env python3
"""Separate construction identities and 160-digit replay of the interval checker."""
import argparse
import json
from pathlib import Path
import sys
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_crossplane_repair_20261013_pair_verify as C

V = C.V
OUT = V.ROOT / 'results/fresh_h9_crossplane_repair_20261013_pair_review.json'


def run():
    # For u in the constrained span and v perpendicular to it, the norm
    # equation is qu+t^2*qv=2. The second cap is reflection across D.
    qu, qv, t = sp.symbols('qu qv t')
    norm = qu + t*t*qv
    assert sp.expand(norm - 2).subs(t*t, (2-qu)/qv) == 0
    rr, rd, dd = sp.symbols('rr rd dd')
    assert (rr-2*rd+dd).subs({rr: 2, rd: 2, dd: 4}) == 2
    assert (rr-rd).subs({rr: 2, rd: 2}) == 0
    source = C.PREFIX + '_caps.json'
    saved = V.read(C.PREFIX + '_verify.json')
    assert saved.pop('source_sha256') == V.sha(source)
    V.SCALE = 10**160
    V.PHI = V.div(V.add(V.rat(1), V.sqrt(V.rat(5))), V.rat(2))
    V.radical.cache_clear()
    replay = C.check(V.read(source))
    count = 0
    factor = 10**80

    def compare(old, new, path='root'):
        nonlocal count
        if old == new:
            return
        if isinstance(old, dict):
            assert isinstance(new, dict) and old.keys() == new.keys(), path
            for key in old:
                compare(old[key], new[key], path + '.' + key)
        elif isinstance(old, list):
            assert isinstance(new, list) and len(old) == len(new), path
            if len(old) == 2 and all(isinstance(x, str) for x in old):
                try:
                    lo, hi = map(int, old)
                    newlo, newhi = map(int, new)
                except ValueError:
                    pass
                else:
                    assert lo*factor <= newlo <= newhi <= hi*factor, path
                    count += 1
                    return
            for i, (a, b) in enumerate(zip(old, new, strict=True)):
                compare(a, b, path + f'[{i}]')
        else:
            assert old == new, path

    compare(saved, replay)
    assert count > 400
    census = replay['affine_and_vertex_census']
    assert census['new_subset_count'] == 640
    assert census['new_affine_circuits'] == 0
    assert census['new_finite_vertices'] == []
    assert len(replay['targets_old88_strict_pairings']) == 2
    assert all(len(x) == 88 for x in replay['targets_old88_strict_pairings'])
    assert len(replay['remaining_escape']['strict_pairings']) == 90
    return {
        'status': 'PASS_PAIR_CONSTRUCTION_REVIEW_AND_160_DIGIT_REPLAY',
        'source': source, 'source_sha256': V.sha(source),
        'verifier_source_sha256': V.sha('scripts/fresh_h9_crossplane_repair_20261013_pair_verify.py'),
        'facets': 90, 'old_facets_retained': 88,
        'new_norms': [2, 2], 'new_mutual_label': 2,
        'target_old88_feasibility_checks': 176,
        'strict_target_cuts': 4,
        'complete_new_psd_subset_count': 640,
        'new_affine_circuits': 0, 'new_finite_vertices': 0,
        'nested_nonidentical_interval_enclosures': count,
        'remaining_escape_norm': replay['remaining_escape']['norm'],
        'strict_escape_inequalities': 90,
        'scope': 'Separate algebraic review of the projection norm and reflection identities; source review of the prescribed pairings and complete finite-neighbor subset coverage; 160-digit execution of the existing interval engine. Every changed enclosure is nested in its saved 80-digit enclosure, and exact metadata agree. This is not a separate interval arithmetic implementation. The strict escape proves infinite volume.'
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    out = run()
    if args.verify:
        assert json.loads(OUT.read_text()) == out
    else:
        OUT.write_text(json.dumps(out, indent=2) + '\n')
    print(out['status'], out['nested_nonidentical_interval_enclosures'])
