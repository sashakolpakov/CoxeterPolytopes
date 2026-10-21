#!/usr/bin/env python3
"""Small symbolic geometry audit, not a second full golden-wall census."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_intrinsic_caps_20261014_golden_review'
BASE = 'fresh_h20_intrinsic_caps_20261014_golden_'


def run():
    paths = [f'{folder}/{BASE}{stage}.{ext}' for stage in ('seed', 'extension')
             for folder, ext in (('scripts', 'py'), ('results', 'json'))]
    data = [json.loads((ROOT / f'results/{BASE}{stage}.json').read_text())
            for stage in ('seed', 'extension')]
    for d in data:
        for p, h in d['sources'].items():
            assert hashlib.sha256((ROOT / p).read_bytes()).hexdigest() == h, p
    seed, ext = data
    assert len(seed['walls']) == 275 and len(ext['new_walls']) == 252
    assert len(seed['finite_endpoints']) == len(ext['new_actual_facet_points']) == 252
    assert sum(seed['pair_counts'].values()) == 275 * 274 // 2
    assert sum(ext['new_pair_counts'].values()) == 527 * 526 // 2 - 275 * 274 // 2
    assert sum(ext['all527_pair_counts'].values()) == 527 * 526 // 2

    # D.U=-1: a vector with spatial square n and null coordinates A,B
    # has square n-2AB. Use arbitrary spatial products, not saved numerics.
    n, m, p, R, T, h2 = S.symbols('n m p R T h2')
    sphere_D = (n - R**2) / 2
    point_D = (n + h2) / 2
    sphere_norm = S.expand(n - 2 * sphere_D)
    point_norm = S.expand(n - 2 * point_D)
    sphere_pair = S.expand(p - sphere_D - (m - T**2) / 2)
    distance2 = n + m - 2 * p
    assert sphere_norm == R**2
    assert point_norm == -h2
    assert S.expand(sphere_pair - (R**2 + T**2 - distance2) / 2) == 0
    point_sphere_pair = S.expand(p - point_D - (m - T**2) / 2)
    assert S.expand(point_sphere_pair - (T**2 - h2 - distance2) / 2) == 0
    assert S.expand(point_sphere_pair.subs({h2: 1, T: 1}) + distance2 / 2) == 0
    # Central radius-two sphere at zero and radius-one cap centered at a.
    assert S.expand(sphere_pair.subs({n: 0, p: 0, R: 2, T: 1}) - (5-m)/2) == 0
    assert point_norm.subs(h2, -1) == 1
    paths.append(str(Path(__file__).resolve().relative_to(ROOT)))
    return dict(
        status='PASS_SCOPED_SYMBOLIC_GEOMETRY_AND_SOURCE_REVIEW',
        sources={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
        symbolic_identities=['sphere_norm', 'point_norm', 'sphere_distance_pair',
                             'point_sphere_distance_pair', 'height_one_facet_pair',
                             'central_cap_pair', 'negative_height_escape_norm'],
        frozen_counts=dict(seed_facets=275, extension_facets=527,
                           certified_old_finite_vertices=252,
                           seed_pairs=37675, extension_new_pairs=100926),
        scope='Independent symbolic identities plus source/proof review and manifest/count checks. '
              'Not a second full implementation or replay of all scalar pairings; no complete vertex census. '
              'Both chambers have certified strict spacelike escapes and therefore infinite volume.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--verify', action='store_true')
    args = ap.parse_args(); out = run(); path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert json.loads(path.read_text()) == out
    else:
        path.write_text(json.dumps(out, indent=2) + '\n')
    print(out['status'])
