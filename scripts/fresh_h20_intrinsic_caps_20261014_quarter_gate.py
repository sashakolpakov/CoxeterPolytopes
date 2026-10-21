#!/usr/bin/env python3
"""Exact two-probe obstruction for the 252 uniform quarter-vertex caps.

Independent Fraction arithmetic; does not import the construction or MILP code.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / 'results/fresh_h20_intrinsic_caps_20261014_golden_seed.json'
OUTPUT = ROOT / 'results/fresh_h20_intrinsic_caps_20261014_quarter_gate.json'
Z = (F(0), F(0))


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def scale(x, r):
    return x[0] * r, x[1] * r


def mul(x, y):
    a, b = x
    c, d = y
    return a*c+b*d, a*d+b*c+b*d


def sign(x):
    # 2*(a+b*phi) = (2*a+b) + b*sqrt(5).
    a, b = 2*x[0]+x[1], x[1]
    if not b:
        return (a > 0) - (a < 0)
    if not a:
        return (b > 0) - (b < 0)
    if (a > 0) == (b > 0):
        return (a > 0) - (a < 0)
    d = a*a-5*b*b
    assert d  # sqrt(5) is irrational.
    return ((a > 0) - (a < 0)) * ((d > 0) - (d < 0))


def pack(x):
    return [str(t) for t in x]


def run():
    p = json.loads(INPUT.read_text())
    gram = [[F(x) for x in row] for row in p['spatial_Gram']]
    sparse = [(i, j, x) for i, row in enumerate(gram)
              for j, x in enumerate(row) if x]
    vertices = [[tuple(F(x) for x in c) for c in row[:19]]
                for row in p['walls'][23:]]
    assert len(vertices) == 252

    def dist(i, j):
        v = [add(a, scale(b, -1)) for a, b in zip(vertices[i], vertices[j])]
        out = Z
        for k, l, b in sparse:
            out = add(out, scale(mul(v[k], v[l]), b))
        return out

    probes = [0, 3]
    covers = [[j for j in range(252) if sign(add(dist(i, j), (F(-36), F(0)))) <= 0]
              for i in probes]
    assert covers == [[0, 7, 14], [3, 10, 17]]
    conflicts = []
    for i in covers[0]:
        for j in covers[1]:
            d = dist(i, j)
            c = add(scale(d, F(1, 72)), (F(-1), F(0)))
            assert d in [(F(48), F(0)), (F(60), F(12))]
            if sign(c) < 0:
                reason = 'obtuse'
            else:
                assert sign(c) > 0 and sign(add(c, (F(-1, 2), F(0)))) < 0
                reason = 'strictly_between_cos_pi_2_and_cos_pi_3'
            conflicts.append({'caps': [i, j], 'old_vertex_distance_square': pack(d),
                              'normalized_cosine': pack(c), 'reason': reason})
    # The old central, vertex and half-vertex balls cannot cover either probe.
    old_margins = []
    for i in probes:
        q = [scale(v, F(1, 4)) for v in vertices[i]]

        def power(center, radius_square):
            v = [add(a, scale(b, -1)) for a, b in zip(q, center)]
            out = (F(-radius_square), F(0))
            for k, l, b in sparse:
                out = add(out, scale(mul(v[k], v[l]), b))
            assert sign(out) > 0
            return out

        powers = [power([Z]*19, 4)]
        for v in vertices:
            powers.extend([power(v, 1), power([scale(x, F(1, 2)) for x in v], 1)])
        minimum = powers[0]
        for v in powers[1:]:
            if sign(add(v, scale(minimum, -1))) < 0:
                minimum = v
        old_margins.append(pack(minimum))
    sources = {}
    for f in [INPUT, Path(__file__)]:
        sources[str(f.relative_to(ROOT))] = hashlib.sha256(f.read_bytes()).hexdigest()
    return {'status': 'PASS_EXACT_TWO_PROBE_OBSTRUCTION', 'sources': sources,
            'field': 'Q(phi), phi^2=phi+1, phi>1',
            'bank': {'centers': 'z_v/4 for all 252 product-alcove vertices',
                     'radius': '3/2', 'coverage_distance_square_threshold': '36'},
            'probe_vertex_indices': probes,
            'probe_omitted_source_walls': [p['finite_endpoints'][i]['omitted_source_walls'] for i in probes],
            'covering_cap_indices': covers, 'cross_conflicts': conflicts,
            'old_527_spherical_minimum_powers': old_margins,
            'conclusion': 'No pairwise Coxeter-compatible subset of this 252-cap bank covers both probes, even with all 527 old walls.',
            'scope': 'Exact exclusion only of the uniform quarter-vertex bank. Face-barycenter caps and other placements are not excluded.'}


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = run()
    if args.verify:
        assert json.loads(OUTPUT.read_text()) == result
    else:
        OUTPUT.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'])
