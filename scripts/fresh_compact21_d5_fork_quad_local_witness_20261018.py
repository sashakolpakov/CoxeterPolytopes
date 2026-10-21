#!/usr/bin/env python3
"""Exact D5 fork quadrilateral star and a coupled D5 triangle.

This is an actual local acute arithmetic configuration, not a compact
finite-volume polytope or a proof that it extends to dimension21.
The seven roots are field-valued; root6 is scaled to norm2Y to avoid
loading a quadratic extension or using numeric eigensolvers.
"""
import argparse
import fcntl
from fractions import Fraction as Q
from hashlib import file_digest
from itertools import combinations
import json
from pathlib import Path
import resource
import signal
import sys
import fresh_compact29_small_faces_20261015 as P

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_compact21_d5_fork_quad_local_witness_20261018'


def sigma(x):
    a, b, c, d = x.v
    assert not b and not d
    return P.K((a, 0, -c, 0))


def integral(x):
    a, b, c, d = x.v
    assert not b and not d
    return (2*c).denominator == 1 and (a-c).denominator == 1


def run():
    s = P.K((0, 0, 1, 0))
    phi = (1+s)/2
    x = phi*phi
    y = (23+9*s)/2
    assert (2*x-5)*(2*y-5) == 9
    assert y+sigma(y) == 23 and y*sigma(y) == 31
    assert y.sign() > 0 and sigma(y).sign() > 0
    assert (y-4).sign() > 0
    g = [[P.K(2*(i == j)) for j in range(7)] for i in range(7)]
    g[6][6] = 2*y
    for i, j in ((0, 1), (1, 2), (2, 3), (2, 4)):
        g[i][j] = g[j][i] = P.K(-1)
    g[4][5] = g[5][4] = -phi
    g[3][6] = g[6][3] = -y
    assert all(g[i][j].sign() <= 0 for i in range(7) for j in range(i))
    assert P.psd_status(P.sub(g, list(range(5)))) == 1
    assert P.det(P.sub(g, list(range(6)))).sign() < 0
    assert not P.det(g)
    conjugate = [[sigma(z) for z in row] for row in g]
    assert P.psd_status(conjugate) == 0
    assert P.psd_status(P.sub(conjugate, list(range(6)))) == 1
    cartan = [[2*g[i][j]/g[j][j] for j in range(7)] for i in range(7)]
    assert all(integral(z) for row in cartan for z in row)
    # The support graph is a tree. Every closed walk uses every edge an even
    # number of times, so cyclic products belong toZ[phi]. The 5bond givesk.
    nonzero = [(i, j) for i in range(7) for j in range(i) if g[i][j]]
    assert len(nonzero) == 6
    assert x == (3+s)/2 and integral(x) and integral(y)
    white, black = [0, 1, 2], [4, 3, 5, 6]
    adjacent = [(0, 1), (1, 2), (2, 3), (3, 0)]
    corners = []
    for i, j in adjacent:
        ids = white+[black[i], black[j]]
        assert P.psd_status(P.sub(g, ids)) == 1
        corners.append({'roots': ids, 'marked': [black[i], black[j]],
                        'type': ['D5', 'A4+A1', 'A3+2A1', 'A4+A1'][len(corners)],
                        'baseline_distance_two_payment': 1 if not corners else 0})
    r = P.schur(g, white, black)
    basis = P.sub(r, [0, 1, 2])
    assert P.det(basis).sign() < 0 and not P.det(r)
    inverse = P.inverse(basis)
    last = P.mv(inverse, [r[i][3] for i in range(3)])
    assert P.dot(last, last, basis) == r[3][3]
    rays, pairings = [], []
    for i, j in adjacent:
        rows = [basis[k] if k < 3 else [r[t][3] for t in range(3)] for k in (i, j)]
        k = next(k for k in range(3) if k not in (i, j))
        rows.append(basis[k])
        ray = P.mv(P.inverse(rows), [P.K(), P.K(), P.K(-1)])
        dots = P.mv(basis, ray)+[sum((r[t][3]*ray[t] for t in range(3)), P.K())]
        assert all(not dots[k] if k in (i, j) else dots[k].sign() < 0 for k in range(4))
        assert P.dot(ray, ray, basis).sign() < 0
        rays.append(ray)
        pairings.append(dots)
    vertex_gram = [[P.dot(a, b, basis) for b in rays] for a in rays]
    assert all(z.sign() < 0 for row in vertex_gram for z in row)
    tested = skipped = 0
    for size in range(1, 8):
        for ids in combinations(range(7), size):
            if 3 in ids and 6 in ids:
                skipped += 1  # Their two-root minor is strictly negative.
                continue
            assert P.psd_status(P.sub(g, ids)) != 0
            tested += 1
    assert tested+skipped == 127
    # Root5 is precisely theL4 necessary third wall forD5 marks(1,4).
    for omitted in (1, 4):
        assert P.psd_status(P.sub(g, [i for i in range(6) if i != omitted])) == 1
    h = Q(5,4)*x-2
    assert h == P.K((Q(-1,8), 0, Q(5,8), 0))
    # Coupling to the two actualE8 triangle types fails cyclic integrality.
    e8_heights = [4+2*s, 1+s]
    assert all(not integral(h*e8) for e8 in e8_heights)
    sources = [f'scripts/{TAG}.py', 'scripts/fresh_compact29_small_faces_20261015.py',
               'results/fresh_compact21_d5_two_e8_triangle_coupling_20261018.json',
               'results/fresh_compact21_single_field_20261016.json']
    pins = {}
    for path in sources:
        with (ROOT/path).open('rb') as stream:
            pins[path] = file_digest(stream, 'sha256').hexdigest()
    return P.enc({'status': 'PASS_D5_FORK_QUAD_LOCAL_ARITHMETIC_WITNESS',
                  'sources': pins, 'field': 'Q(sqrt5)', 'Gram': g,
                  'root6_norm': 2*y, 'opposite_normalized_squares': [x, y],
                  'physical_signature': [5, 1, 1], 'nonidentity_signature': [6, 0, 1],
                  'integral_Cartan_matrix': cartan, 'corners': corners,
                  'quad_white': white, 'quad_black_cyclic': black, 'residual_Gram': r,
                  'face_rays': rays, 'face_wall_pairings': pairings, 'face_vertex_Gram': vertex_gram,
                  'principal_subsets_tested': tested, 'negative_pair_skips': skipped,
                  'baseline_quad_payment': 1, 'forced_D5_triangle': [1, 4],
                  'E8_triangle_coexistence_excluded': True,
                  'global_compact_extension_proved': False, 'finite_volume_example': False,
                  'new_dimension_bound': None, 'charge_saving_proved': False,
                  'scope': 'Actual compactquadrilateral in a seven-wall local acute chamber, '
                           'with all local arithmetic conditions. OrthogonalE8 spectators '
                           'preserve the star, but no completecompact H21 chamber is produced.'})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    signal.alarm(59)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024)
    assert rss < 1024**3
    lock_path = ROOT/'formal/.proof-batches/exclusive.lock'
    lock_path.parent.mkdir(exist_ok=True)
    with lock_path.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        result = run()
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024)
        assert rss < 1024**3
        target = ROOT/f'results/{TAG}.json'
        if args.verify:
            assert json.loads(target.read_text()) == result
        else:
            assert not target.exists()
            target.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'quad_corners': 4,
                      'principal_subsets': 127, 'baseline_quad_payment': 1,
                      'forced_D5_triangle': [1, 4], 'peak_rss_bytes': rss,
                      'new_dimension_bound': None}))
