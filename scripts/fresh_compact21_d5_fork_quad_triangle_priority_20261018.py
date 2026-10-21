#!/usr/bin/env python3
"""Exact small checks for a D5+2E8 quadrilateral-only priority payment.

All-distance completeness and payment preservation are proved in the
companion report. There is no distance search or optimizer in this file.
"""
import argparse
import fcntl
from fractions import Fraction as Q
from hashlib import file_digest
from itertools import product
import json
from pathlib import Path
import resource
import signal
import sys
import fresh_compact21_d5_fork_quad_single_basic_unit_20261018 as B

P = B.P
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_compact21_d5_fork_quad_triangle_priority_20261018'


def distances(g, start):
    answer, frontier = {start: 0}, {start}
    while frontier:
        following = set()
        for i in frontier:
            for j in range(len(g)):
                if i != j and g[i][j] and j not in answer:
                    answer[j] = answer[i]+1
                    following.add(j)
        frontier = following
    return answer


def marked_fork(g, i, j):
    degrees = [sum(k != t and bool(g[t][k]) for k in range(len(g))) for t in range(len(g))]
    return (len(g) == 5 and len(distances(g, 0)) == 5
            and sorted(degrees) == [1, 1, 1, 2, 3]
            and all(g[u][v] == -1 for u in range(5) for v in range(u) if g[u][v])
            and degrees[i] == degrees[j] == 1 and distances(g, i).get(j) == 2)


def marked_half(g, i, j):
    degrees = [sum(k != t and bool(g[t][k]) for k in range(len(g))) for t in range(len(g))]
    return (len(g) == 5 and sorted(degrees) == [1, 1, 1, 2, 3]
            and all(g[u][v] == -1 for u in range(5) for v in range(u) if g[u][v])
            and degrees[i] == degrees[j] == 1 and distances(g, i).get(j) == 3)


def run():
    s, sqrt2 = P.K((0, 0, 1, 0)), P.K((0, 1, 0, 0))
    phi = (1+s)/2
    endpoint_types = [(None, 0), (0, 3), (0, 4), (1, 3), (2, 3)]
    permitted = []
    label_value = {3: P.K(1), 4: sqrt2}
    for a, b in product(endpoint_types, repeat=2):
        g = [[P.K(2*(i == j)) for j in range(5)] for i in range(5)]
        for i, j in ((0, 1), (1, 2)):
            g[i][j] = g[j][i] = P.K(-1)
        for new, (contact, label) in ((3, a), (4, b)):
            if contact is not None:
                g[contact][new] = g[new][contact] = -label_value[label]
        if P.psd_status(g) == 1 and (not B.basic_payment(g, 3, 4) or marked_fork(g, 3, 4)):
            permitted.append((a, b))
    assert len(permitted) == 19
    unordered = {tuple(sorted(pair, key=lambda z: (-1 if z[0] is None else z[0], z[1])))
                 for pair in permitted}
    assert len(unordered) == 11
    half = [[P.K(2*(i == j)) for j in range(5)] for i in range(5)]
    for i, j in ((0, 1), (1, 2), (1, 3), (2, 4)):
        half[i][j] = half[j][i] = P.K(-1)
    assert P.psd_status(half) == 1 and marked_half(half, 3, 4)
    # The other half corner is the same marked diagram, with the long
    # endpoint named the old side instead of the new side.

    x, y = B.Poly({(1, 0): 1}), B.Poly({(0, 1): 1})
    rn23 = [[Q(5, 4)-x*Q(1, 2), Q(-3, 4), Q(-3, 4)],
            [Q(-3, 4), Q(5, 4), y-Q(3, 4)],
            [Q(-3, 4), y-Q(3, 4), Q(5, 4)]]
    fn23 = (5-2*x)*(2*y+1)-9
    assert 8*B.det_poly(rn23) == (2-y)*fn23
    r0323 = [[Q(5, 4), Q(-3, 4), x-Q(1, 4), Q(-3, 4)],
             [Q(-3, 4), Q(5, 4), Q(-1, 4), y-Q(3, 4)],
             [x-Q(1, 4), Q(-1, 4), Q(5, 4), Q(-1, 4)],
             [Q(-3, 4), y-Q(3, 4), Q(-1, 4), Q(5, 4)]]
    f0323 = (2*y+1)*(3+x-2*x*x)+3*x-7
    assert 4*B.det_poly(r0323) == (2-y)*f0323
    r0423 = [[Q(5, 4), Q(-3, 4), x-Q(1, 2), Q(-3, 4)],
             [Q(-3, 4), Q(5, 4), Q(-1, 2), y-Q(3, 4)],
             [x-Q(1, 2), Q(-1, 2), 1, Q(-1, 2)],
             [Q(-3, 4), y-Q(3, 4), Q(-1, 2), Q(5, 4)]]
    f0423 = (2*y+1)*(1+x-x*x)+3*x-5
    assert 2*B.det_poly(r0423) == (2-y)*f0423
    r2323 = [[Q(5, 4), Q(-3, 4), x-Q(3, 4), Q(-3, 4)],
             [Q(-3, 4), Q(5, 4), Q(-3, 4), y-Q(3, 4)],
             [x-Q(3, 4), Q(-3, 4), Q(5, 4), Q(-3, 4)],
             [Q(-3, 4), y-Q(3, 4), Q(-3, 4), Q(5, 4)]]
    f2323 = (2*x+1)*(2*y+1)-9
    assert 4*B.det_poly(r2323) == (2-x)*(2-y)*f2323
    r0303 = [[Q(5, 4), Q(-3, 4), x-Q(1, 4), Q(-1, 4)],
             [Q(-3, 4), Q(5, 4), Q(-1, 4), y-Q(1, 4)],
             [x-Q(1, 4), Q(-1, 4), Q(5, 4), Q(-3, 4)],
             [Q(-1, 4), y-Q(1, 4), Q(-3, 4), Q(5, 4)]]
    f0303 = (x+1)*(y+1)-1
    assert 4*B.det_poly(r0303) == f0303*((3-2*x)*(3-2*y)-1)

    # Complete unit branches for N-P23, from the earlier proved UV=9 lemma.
    n23_branches = []
    for u in (s-2, 18+9*s):
        v = 9/u
        xx, yy = (u+5)/2, -(v+1)/2
        assert (5-2*xx)*(2*yy+1) == 9
        assert B.integral(xx) and B.integral(yy)
        n23_branches.append(P.enc([xx, yy]))
    assert (s-2+5)/2 == phi*phi
    invalid_y = (1-s)/2
    assert invalid_y.sign() < 0 and (invalid_y+1).sign() > 0

    # P03-P23: finite x=-phi forces a nonintegral y. The two possible
    # finite y's force nonsquare discriminants.
    bad_y = -(3+12*phi)/11
    assert (2*bad_y+1)*(3-phi-2*phi*phi)-3*phi-7 == 0
    assert not B.integral(bad_y)
    ob0323 = [B.certify_nonsquare(v) for v in (P.K(21), 134+50*s)]
    # P04-P23: only finite x=-2 is field-compatible, and y=-8/5 fails.
    bad_y4 = P.K(Q(-8, 5))
    assert (2*bad_y4+1)*(-5)-11 == 0 and not B.integral(bad_y4)
    ob0423 = [B.certify_nonsquare(v) for v in (P.K(28), 34+14*s)]
    # P23-P23: finite x=-1 gives y=-5 and a rational indefinite Gram;
    # finite x=-phi gives a nonintegral y.
    bad_yy = -P.K(Q(1, 2))-Q(9, 10)*s
    assert (1-2*phi)*(2*bad_yy+1) == 9 and not B.integral(bad_yy)
    assert (2*(-1)+1)*(2*(-5)+1) == 9
    rational = [[P.K(2*(i == j)) for j in range(6)] for i in range(6)]
    for i, j in ((0, 1), (1, 2), (2, 3), (2, 4), (2, 5), (3, 5)):
        rational[i][j] = rational[j][i] = P.K(-1)
    assert P.det(rational) == -21

    # The P03-P03 unit branches are U=phi or 1/phi. Both yield the
    # opposite pair {-phi,-phi^2}; the finite side excludes both E8 types.
    h_l = (-1+5*s)/8
    h_p = (1+9*phi)/4
    e8 = [4+2*s, 1+s]
    assert (1+phi+Q(5,4)*phi*phi-2) == h_p
    assert h_p*e8[0] == (19+29*phi)/2
    assert h_p*e8[1] == (9+10*phi)/2
    assert all(not B.integral(h*h8) for h in (h_l, h_p) for h8 in e8)
    assert (-phi+1)*(-phi*phi+1) == 1

    sources = [f'scripts/{TAG}.py',
               'scripts/fresh_compact21_d5_fork_quad_single_basic_unit_20261018.py',
               'results/fresh_compact21_d5_fork_quad_single_basic_unit_20261018.json',
               'results/fresh_compact21_d5_fork_quad_local_witness_20261018.json',
               'scripts/fresh_compact_global_bound_20261015.py',
               'results/fresh_compact21_single_field_20261016.json']
    pins = {}
    for path in sources:
        with (ROOT/path).open('rb') as stream:
            pins[path] = file_digest(stream, 'sha256').hexdigest()
    return {'status': 'PASS_D5_FORK_QUAD_TRIANGLE_PRIORITY_EXACT_CHECKS',
            'sources': pins, 'ordered_finite_patterns': permitted,
            'ordered_count': 19, 'unordered_count': 11,
            'new_corank_polynomials': {'N-P23': fn23.encode(), 'P03-P23': f0323.encode(),
                                      'P04-P23': f0423.encode(), 'P23-P23': f2323.encode(),
                                      'P03-P03': f0303.encode()},
            'N_P23_complete_necessary_branches': n23_branches,
            'P03_P23_bad_finite_x_opposite': P.enc(bad_y),
            'P03_P23_discriminant_obstructions': ob0323,
            'P04_P23_bad_finite_x_opposite': P.enc(bad_y4),
            'P04_P23_discriminant_obstructions': ob0423,
            'P23_P23_bad_finite_x_opposite': P.enc(bad_yy),
            'P03_P03_opposite_products': P.enc([-phi, -phi*phi]),
            'forced_root_heights': P.enc([h_l, h_p]),
            'E8_cyclic_products': P.enc([[h*h8 for h8 in e8] for h in (h_l, h_p)]),
            'retained_half_bonus_in_P13_P23': '1/2 + 1/2 = 1',
            'quad_only_payment_preservation': 'proved in companion report; not inferred from finite sampling',
            'uniform_vertex_allowance_reduction_proved': False, 'new_dimension_bound': None}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    signal.alarm(59)
    rss = lambda: resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024)
    assert rss() < 1024**3
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        result = run()
        assert rss() < 1024**3
        target = ROOT/f'results/{TAG}.json'
        if args.verify:
            assert json.loads(target.read_text()) == json.loads(json.dumps(result))
        else:
            assert not target.exists()
            target.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'ordered_patterns': 19,
                      'unordered_patterns': 11, 'peak_rss_bytes': rss(),
                      'new_dimension_bound': None}))
