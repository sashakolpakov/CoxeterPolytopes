#!/usr/bin/env python3
"""Exact mixed weighted-box covering obstruction, with every finite label.

Replays the frozen R2 cutter/R4 compatibility precursor, then checks the
all-label extension, a feasible golden Q2, and the global Q8 covering count.
This is a scoped obstruction for two cap templates, not for arbitrary caps.
"""
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import argparse
import hashlib
import json

import fresh_h20_seed_deformation_20261014 as P

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_weighted_mixed_probe_20261014'
PRIOR = 'fresh_h20_seed_deformation_20261014'
q, add, sub, mul, scale, sign, total = P.q, P.add, P.sub, P.mul, P.scale, P.sign, P.total
O, Z, PHI = P.O, P.Z, P.PHI


def run():
    prior = json.loads((ROOT / f'results/{PRIOR}.json').read_text())
    # This replays the 88 recipes, every adjacent-probe pair, and the R4 gates.
    assert P.run() == prior
    assert prior['R4_vertex_family']['maximum_compatible_caps'] == 512
    assert prior['R2_nonvertex_family']['adjacent_golden_probe_pair_counts'] == {
        'forbidden': 6484, 'right': 1260}

    # For radius two, opposite physical seed gaps sum to 2b. The smaller
    # one is <=b<2, hence right or finite. Its label is at most 3,4,5 for
    # b=1,sqrt2,phi. Only two complementary gaps escaped the prior 2..6 test.
    sqrt2 = (F(1414213, 10**6), F(1414214, 10**6))
    phi = (F(1618033, 10**6), F(1618034, 10**6))
    assert sqrt2[0]**2 < 2 < sqrt2[1]**2
    assert phi[0]**2-phi[0]-1 < 0 < phi[1]**2-phi[1]-1
    assert 1 < sqrt2[0] < sqrt2[1] < phi[0] < phi[1]
    assert phi[1]**2 < 3 and phi[1] < 2
    gaps = [(2*sqrt2[0]-1, 2*sqrt2[1]-1),
            (2*phi[0]-sqrt2[1], 2*phi[1]-sqrt2[0])]
    assert all(F(91, 50) < lo < hi < F(46, 25) for lo, hi in gaps)
    # pi>3, cosine decreases on (0,pi), and cos(x)<1-x²/2+x⁴/24 at x=3/7.
    x = F(3, 7)
    twice_cos7_upper = 2*(1-x*x/2+x**4/24)
    assert twice_cos7_upper < F(91, 50)
    # 2cos(pi/8)=sqrt(2+sqrt2)>46/25.
    assert 2+sqrt2[0] > F(46, 25)**2
    # Remaining smaller-gap choices, with their complementary relations.
    one_coordinate_cases = {
        'unit': [['0', 'parallel'], ['1', 'pi3']],
        'heavy': [['0', 'dotted'], ['1', 'forbidden_between_pi7_pi8'],
                  ['sqrt2', 'pi4']],
        'golden': [['0', 'dotted'], ['1', 'dotted'],
                   ['sqrt2', 'forbidden_between_pi7_pi8'], ['phi', 'pi5']]}

    weights = [O]*10 + [q(2)] + [add(O, PHI)]*8
    N = total(weights)
    assert N == (F(20), F(8))
    def dot(a, b):
        return total(mul(w, mul(x, y)) for w, x, y in zip(weights, a, b, strict=True))
    def dist(a, b):
        d = [sub(x, y) for x, y in zip(a, b, strict=True)]
        return dot(d, d)
    def center(a, zero_unit, pattern):
        return ([q(i != zero_unit) for i in range(10)] + [O] +
                [scale(O if j == a else sub(q(2), PHI),
                       -1 if (pattern >> j) & 1 else 1) for j in range(8)])
    def half(pattern):
        return [q(F(1, 2))]*11 + [q(F(-1 if (pattern >> j) & 1 else 1, 2))
                                      for j in range(8)]

    # Q2 itself IS feasible: the diagonal pairs are pi/3, not obstructions.
    square = [(0, 0, 0), (1, 1, 1), (3, 0, 0), (2, 1, 1)]
    centers = [center(a, i, pattern) for pattern, a, i in square]
    q2_pairs = []
    for i, c in enumerate(centers):
        assert dot(c, c) == (F(26), F(-6))
        assert sign(sub(dot(c, c), q(16))) > 0  # central dotted
        for j, (pattern, _, _) in enumerate(square):
            power = sub(dist(c, half(pattern)), q(4))
            if i == j:
                assert power == (F(8), F(-5)) and sign(power) < 0
            else:
                assert sign(power) > 0
        for j in range(i):
            d2 = dist(c, centers[j])
            assert d2 in (q(8), q(12))
            q2_pairs.append([j, i, 'right' if d2 == q(8) else 'pi3', d2])
    assert [r[2] for r in q2_pairs].count('right') == 4
    assert [r[2] for r in q2_pairs].count('pi3') == 2

    # |vertex-halfprobe|²=N/4+2H. No radius-four boundary equality occurs.
    threshold = sub(q(F(11, 2)), PHI)
    cover_rows, count_by_golden, boundary = [], [0]*9, 0
    for h1, h2, h5 in product(range(11), range(2), range(9)):
        H = add(q(h1+2*h2), scale(add(O, PHI), h5))
        s = sign(sub(H, threshold))
        boundary += s == 0
        if s < 0:
            count = comb(10, h1)*comb(8, h5)
            cover_rows.append([h1, h2, h5, count])
            count_by_golden[h5] += count
    assert boundary == 0
    assert count_by_golden == [187, 88, 0, 0, 0, 0, 0, 0, 0]
    per_R4 = sum(count_by_golden)
    assert per_R4 == 275

    # The induced Q8 edges among R2-covered probes can be oriented toward
    # their unique full-coordinate endpoint. Each vertex accepts <=1 edge.
    # Thus e<=v, whereas degree eight implies e>=8v-1024.
    max_R2_per_fiber = 1024 // 7
    assert max_R2_per_fiber == 146
    assert 8*146-1024 <= 146 and 8*147-1024 > 147
    fibers, probes = 2**11, 2**19
    max_R2, max_R4 = fibers*146, 512*per_R4
    missed = probes-max_R2-max_R4
    assert (max_R2, max_R4, missed) == (299008, 140800, 84480)
    required_R4 = fibers*(256-146)
    assert required_R4 == 225280 and required_R4 > max_R4

    # At a vertex the opposite unit-wall gap is two. R=4 is the only
    # seed-compatible radius with R²>8. Every other radius misses all probes.
    min_vertex_probe_d2 = scale(N, F(1, 4))
    assert min_vertex_probe_d2 == (F(5), F(2))
    assert sign(sub(min_vertex_probe_d2, q(8))) > 0
    assert sign(sub(min_vertex_probe_d2, q(4))) > 0  # central also misses

    paths = [f'scripts/{TAG}.py', f'scripts/{PRIOR}.py', f'results/{PRIOR}.json']
    return P.enc(dict(
        status='PASS_ALL_LABEL_MIXED_R2_AND_ALL_VERTEX_RADII_OBSTRUCTION',
        sources={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        inherited_cutter_replay='PASS', retained_core_facets=39,
        vertex_square=N, total_half_probes=probes,
        all_label_R2_gate=dict(smaller_gap_max_labels=[3, 4, 5],
            coordinate_cases=one_coordinate_cases, complementary_gap_bounds=gaps,
            twice_cos_pi7_upper=twice_cos7_upper, twice_cos_pi8_lower=F(46, 25),
            resulting_center_alphabet_physical=prior['one_coordinate_center_options_physical']),
        feasible_Q2=dict(recipes=[dict(pattern=p, full_golden=a, zero_unit=i)
                                 for p, a, i in square], centers=centers,
                        pairs=q2_pairs, own_probe_power=(F(8), F(-5))),
        R4_half_probe_cover=dict(hamming_trials=198, cover_rows=cover_rows,
            counts_by_golden_difference=count_by_golden, boundary_cases=boundary,
            probes_per_cap=per_R4, maximum_compatible_caps=512,
            maximum_covered_half_probes=max_R4),
        R2_induced_cube_bound=dict(fibers=fibers, probes_per_fiber=256,
            induced_edge_upper='e <= v', induced_edge_lower='e >= 8v-1024',
            maximum_covered_per_fiber=max_R2_per_fiber,
            maximum_covered_half_probes=max_R2),
        mixed_union=dict(required_R4_for_full_coverage=required_R4,
            maximum_total_covered=max_R2+max_R4, minimum_uncovered_half_probes=missed),
        arbitrary_vertex_radius_gate=dict(allowed_large_radius=4,
            other_radius_square_upper=8, minimum_half_probe_distance_square=min_vertex_probe_d2),
        scope='Retain the 38 weighted-box walls and central radius-two cap. '
              'For any finite pairwise Coxeter-compatible collection of arbitrary radius-two '
              'spherical caps and arbitrary vertex-centered seed-compatible spherical caps '
              'of any radius, at least 84480 half-vertex probes miss every closed cap ball. '
              'Every finite seed label, parallel and dotted relations are permitted. '
              'Such a collection cannot give finite volume. Arbitrary nonvertex radii, '
              'changed core, and arbitrary H20 Coxeter polytopes remain outside this theorem.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    data = run()
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert json.loads(path.read_text()) == data
    else:
        path.write_text(json.dumps(data, indent=2)+'\n')
    print(data['status'])
    print(data['mixed_union'])
