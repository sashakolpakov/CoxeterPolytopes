#!/usr/bin/env python3
"""Numerical monotone scores in one fixed residual affine chart.

These diagnostics guide search; they are not exact geometry certificates.
"""
import argparse
import json
from hashlib import sha256
from pathlib import Path
import numpy as np
from scipy.spatial import ConvexHull

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_residual_monotone_scores_20261016'
INPUTS = [
    'results/fresh_h20_second_cusp_pi8_cap_20261016_residual_scout.json',
    'results/fresh_h20_pi8_triangle_face_20261016_residual_scout.json',
]


def run():
    a, b = [json.loads((ROOT / p).read_text()) for p in INPUTS]
    seed, samples = 20261016, 65536
    rng = np.random.default_rng(seed)
    directions = np.abs(rng.normal(size=(samples, 4)))
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    records = []
    previous_radii = None
    for name, pack in [('151', a['old151']), ('152_pi8', b['old152']), ('153_triangle', b['new153'])]:
        points = np.array([row['z'] for row in pack['vertices']])
        hull = ConvexHull(points)
        coefficients, rhs = hull.equations[:, :4], -hull.equations[:, 4]
        # The discarded origin facets are the nonnegative-coordinate
        # walls; they impose no radial upper bound inside this orthant.
        mask = rhs > 1e-10
        coefficients, rhs = coefficients[mask], rhs[mask]
        radii = []
        for part in np.array_split(directions, 32):
            products = coefficients @ part.T
            bounds = np.divide(rhs[:, None], products, out=np.full_like(products, np.inf), where=products > 1e-14)
            radii.extend(np.min(bounds, axis=0))
        radii = np.array(radii)
        assert np.isfinite(radii).all() and min(radii) > 0
        if previous_radii is not None:
            assert np.max(radii - previous_radii) < 1e-8
        previous_radii = radii
        outside = radii > 1
        # S3 has area2*pi^2; its positive orthant has1/16of that.
        # Radial4-volume integrates r^3dr, giving prefactor pi^2/32.
        observations = np.pi**2 / 32 * np.maximum(radii**4 - 1, 0)
        fraction = float(np.mean(outside))
        worst = sorted(pack['vertices'], key=lambda row: -row['squared_radius'])[:8]
        records.append(dict(chamber=name, affine_Euclidean_volume=float(hull.volume),
            maximum_vertex_squared_radius=float(np.max(np.sum(points * points, axis=1))),
            sampled_exterior_direction_count=int(np.sum(outside)),
            estimated_exterior_direction_fraction=fraction,
            fraction_sampling_standard_error=float(np.sqrt(fraction * (1 - fraction) / samples)),
            estimated_exterior_Euclidean_volume=float(np.mean(observations)),
            volume_sampling_standard_error=float(np.std(observations, ddof=1) / np.sqrt(samples)),
            largest_exterior_vertices=worst))
    assert all(records[i + 1]['affine_Euclidean_volume'] < records[i]['affine_Euclidean_volume'] for i in (0, 1))
    removed = records[0]['affine_Euclidean_volume'] - records[-1]['affine_Euclidean_volume']
    paths = [f'scripts/{TAG}.py'] + INPUTS
    return dict(status='PASS_NUMERICAL_FIXED_CHART_MONOTONE_RESIDUAL_SCORES',
        sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
        seed=seed, shared_sample_count=samples, records=records,
        combined_affine_volume_removed=removed,
        combined_fraction_of_starting_affine_volume_removed=removed / records[0]['affine_Euclidean_volume'],
        method='Floating Qhull hulls of saved numerical vertices in the same C+four-caps chart; common absolute-Gaussian uniformS3-positive-orthant directions; radial integrals. Standarderrors describe sampling only, not geometry/rounding errors.',
        interpretation='Each cut shrinks the region although exteriorvertexcount can rise. These cuts leave the largest z142-region and maximumradius untouched; prioritize walls meeting its old133/149 boundaries.',
        scope='Numerical search guidance only. AffineEuclideanvolume is not hyperbolicvolume. No exactcone census, probability guarantee, finitevolume proof orclaim that one residualface captures the wholeambientchamber.')


def close(a, b):
    if isinstance(a, float):
        return isinstance(b, (int, float)) and abs(a - b) <= 1e-7 * (1 + abs(a))
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(close(v, b[k]) for k, v in a.items())
    if isinstance(a, list):
        return len(a) == len(b) and all(close(x, y) for x, y in zip(a, b))
    return a == b


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--verify', action='store_true'); args = ap.parse_args()
    result = run(); out = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert close(result, json.loads(out.read_text()))
    else:
        out.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], flush=True)
