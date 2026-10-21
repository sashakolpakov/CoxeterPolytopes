#!/usr/bin/env python3
"""Separate source audit / precision replay of the remote83 and85 certificates.

This deliberately reuses the constructors' exact interval engine. It is a
different-precision replay, not a second independent field implementation.
Every archived strict enclosure is checked against a tighter enclosure.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import fresh_h9_remote_caps_20261013_verify as V
from scripts import fresh_h9_remote_caps_20261013_coupled_verify as C

OUT = 'results/fresh_h9_remote_caps_independent_20261013.json'


def intervals(d, coupled):
    """Return all numerical enclosures, leaving exact metadata untouched."""
    paths = []
    for i in range(len(d['dotted_checks'])):
        paths.append(('dotted_checks', i, 'pairing'))
    for i in range(len(d['target_cut_pairings'])):
        paths.append(('target_cut_pairings', i))
    for i in range(len(d['remaining_escape']['strict_pairings'])):
        paths.append(('remaining_escape', 'strict_pairings', i))
    if not coupled:
        for i in range(len(d['vertex_nonincident_checks'])):
            paths.append(('vertex_nonincident_checks', i, 'pairing'))
    else:
        for k in d['boost_equation']:
            paths.append(('boost_equation', k))
        for vertex in ('joint_vertex', 'attached_vertex', 'original_vertex'):
            for i in range(len(d[vertex]['nonincident_checks'])):
                paths.append((vertex, 'nonincident_checks', i, 'pairing'))
            for i in range(len(d[vertex].get('point_enclosures', []))):
                paths.append((vertex, 'point_enclosures', i))
    result = {}
    cleaned = copy.deepcopy(d)
    for path in paths:
        a = d
        b = cleaned
        for k in path[:-1]:
            a = a[k]
            b = b[k]
        x = a[path[-1]]
        assert len(x) == 2 and all(isinstance(z, str) for z in x)
        lo, hi = map(int, x)
        assert lo <= hi
        result[path] = (lo, hi)
        b[path[-1]] = 'CHECKED_INTERVAL'
    cleaned.pop('interval_scale', None)
    return result, cleaned


def run():
    digits = 100
    V.SCALE = 10**digits
    V.PHI = V.div(V.add(V.rat(1), V.sqrt(V.rat(5))), V.rat(2))
    V.radical.cache_clear()
    rows = []
    hashes = {}
    for module, coupled in ((V, False), (C, True)):
        prefix = module.PREFIX
        source = prefix + '_caps.json'
        saved = prefix + '_verify.json'
        for p in (source, saved):
            hashes[p] = V.sha(p)
        data = V.read(source)
        old = V.read(saved)
        new = module.check(data)
        new['source_sha256'] = V.sha(source)
        oi, om = intervals(old, coupled)
        ni, nm = intervals(new, coupled)
        assert om == nm, 'Non-interval metadata changed during replay'
        assert set(oi) == set(ni)
        old_scale = 10**80 if coupled else int(old['interval_scale'])
        assert V.SCALE % old_scale == 0
        factor = V.SCALE // old_scale
        for p, (lo, hi) in ni.items():
            a, b = oi[p]
            assert a*factor <= lo <= hi <= b*factor, p
        rows.append({
            'walls': new.get('total_facets'),
            'status': new['status'],
            'strict_enclosures_replayed_and_nested': len(ni),
            'new_old_dotted_pairs': new['new_old_dotted_pairings'],
            'new_old_right_pairs': new['new_old_right_pairings'],
            'remaining_escape_norm': new['remaining_escape']['norm'],
            'remaining_escape_center_pairing': new['remaining_escape']['future_center_pairing'],
            'actual_finite_vertices_checked': 3 if coupled else 1,
            'finite_edges_checked': len(new.get('finite_edges', [])),
        })
    for module in (V, C, V.K):
        path = Path(module.__file__).resolve().relative_to(ROOT).as_posix()
        hashes[path] = V.sha(path)
    return {
        'status': 'PASS_SOURCE_REVIEW_AND_100_DIGIT_OUTWARD_INTERVAL_REPLAY',
        'decimal_precision': digits,
        'source_hashes': hashes,
        'stages': rows,
        'scope': ('Source-level proof review plus a100-digit replay of the same '
                  'exact outward-rounded integer interval implementation. This '
                  'is not a second independent field implementation. Exact '
                  'frame/boost identities prove norms and zero/finite pairings; '
                  'strict interval signs prove every remaining relation, actual '
                  'vertex inequality, and saved spacelike escape. The frozen79 '
                  'norm/angle proof remains an explicitly hash-anchored dependency. '
                  'Both83 and85 stages have infinite volume; no global frontier '
                  'completeness or completion obstruction is asserted.')
    }


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    args = p.parse_args()
    out = run()
    if args.verify:
        assert out == json.loads((ROOT/OUT).read_text())
    else:
        (ROOT/OUT).write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'status': out['status'], 'stages': out['stages']}))
