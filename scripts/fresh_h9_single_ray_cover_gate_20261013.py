#!/usr/bin/env python3
"""Exact data for a covering-ball obstruction to a single H4 weight ray.

This is an all-center/all-label theorem for each of four separate families.
It does not exclude mixing the families or mixed transverse profiles.
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

from fresh_h9_makarov_face_extension_20261013 import K, build, inv, mv, pair

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_single_ray_cover_gate_20261013'
SOURCE = ROOT/'results/fresh_h9_makarov_face_extension_20261013.json'
Z = K()


def enc(xs):
    return [x.enc() for x in xs]


def run(verify=False):
    s = json.loads(SOURCE.read_text())
    # Replay compactness, the complete old vertex census, and the source gaps.
    assert s == build()
    Q = [[K(*a) for a in row] for row in s['form']]
    G = [row[:6] for row in Q[:6]]
    W = [[K(*a) for a in row[:6]] for row in s['seed_normals'][:7]]
    V = [[K(*a) for a in row] for row in s['compact_H5_vertices']]
    B = [[K(*a) for a in row] for row in s['H4_inverse']]
    A = [mv(G, w) for w in W]
    # Exact rational relaxation: the true polar-wall minimum is sqrt(q6/2).
    low = [K(1)]*6 + [K(F(2, 5))]
    assert (pair(W[6], W[6], G)/2-low[6]*low[6]).sign() > 0
    inverse_by_omission = {}
    singular = []
    for omitted in range(7):
        active = [i for i in range(7) if i != omitted]
        try:
            inverse_by_omission[omitted] = inv([A[i] for i in active])
        except StopIteration:
            singular.append(omitted)

    target = B[3][3]-2
    assert target == K(1, 1)
    supports = []
    minimum_by_zero_count = {}
    admissible_zero_sets = []
    for bits in range(1, 128):
        nonzero = [i for i in range(7) if bits >> i & 1]
        zeros = [i for i in range(7) if i not in nonzero]
        rhs = [-low[i] if i in nonzero else Z for i in range(7)]
        candidates = []
        # All forced equalities must be among the six active constraints.
        for omitted in nonzero:
            if omitted not in inverse_by_omission:
                continue
            active = [i for i in range(7) if i != omitted]
            u = mv(inverse_by_omission[omitted], [rhs[i] for i in active])
            values = mv(A, u)
            if any((a-b).sign() > 0 for a, b in zip(values, rhs)):
                continue
            assert all(not values[i] for i in zeros)
            square = -pair(u, u, G)
            assert square.sign() > 0
            candidates.append({'omitted': omitted, 'normal_old_part': enc(u),
                               'positive_lorentz_norm_square': square.enc()})
            n = len(zeros)
            if n not in minimum_by_zero_count or (square-minimum_by_zero_count[n]).sign() < 0:
                minimum_by_zero_count[n] = square
            if n <= 4:
                assert (square-target).sign() > 0
        if candidates:
            values = [K(*r['positive_lorentz_norm_square']) for r in candidates]
            if any((x-target).sign() <= 0 for x in values):
                assert len(zeros) == 5 and zeros in s['vertex_active_sets']
                admissible_zero_sets.append(zeros)
        supports.append({'zero_indices': zeros, 'nonzero_indices': nonzero,
                         'relaxed_polyhedron_vertices': candidates})

    # A strict old interior point outside every j=3,m=3 vertex-center ball.
    x = [K(a) for a in s['all_label_vertex_center_family_obstruction']['strict_integer_escape'][:6]]
    xx = pair(x, x, G)
    assert xx.sign() < 0 and all(v.sign() < 0 for v in mv(A, x))
    C3 = B[3][0]*B[3][0]/B[0][0]/target
    gaps = []
    for v in V:
        px = pair(v, x, G)
        assert px.sign() < 0
        gap = px*px-C3*pair(v, v, G)*xx
        assert gap.sign() > 0
        gaps.append(gap)

    # Two old vertices farther apart than twice either isolated-ball radius.
    v, w = V[1], V[8]
    vw = pair(v, w, G)
    diameter_lower_cosh_square = vw*vw/(pair(v, v, G)*pair(w, w, G))
    ray_records = []
    expected_D = [K(F(9, 4), F(3, 4)), K(3, 1),
                  K(F(5, 4), F(1, 4)), K(F(3, 4))]
    for j, i in enumerate([3, 3, 3, 0]):
        D = B[j][j]-B[j][i]*B[j][i]/B[i][i]
        assert D == expected_D[j]
        # i is a worst transverse axis for this ray.
        assert all((B[i][i]*B[j][k]*B[j][k]
                    -B[k][k]*B[j][i]*B[j][i]).sign() >= 0 for k in range(4))
        record = {'ray': j, 'worst_axis': i, 'D': D.enc()}
        if j < 2:
            assert (D-2).sign() > 0
            record['reason'] = 'No cap with t>=1 cuts any ideal direction on this axis.'
        else:
            min_isolated_t_square = 1 if j == 2 else 2
            # Strict separation D(t^2+s^2)>2, with every other cap s>=1.
            separation_margin = D*(min_isolated_t_square+1)-2
            assert separation_margin.sign() > 0
            C = B[j][i]*B[j][i]/B[i][i]/(B[j][j]-K(2)/min_isolated_t_square)
            assert (C-1).sign() > 0
            diameter_margin = diameter_lower_cosh_square-(2*C-1)*(2*C-1)
            assert diameter_margin.sign() > 0
            # Coxeter labels m>=4 (j2), or m>=6 (j3), cannot cut the axis.
            first_ineffective_t_square = 2 if j == 2 else 3
            assert (D-K(2)/first_ineffective_t_square).sign() > 0
            record.update({'minimum_isolated_t_square': min_isolated_t_square,
                           'strict_ball_separation_margin': separation_margin.enc(),
                           'largest_isolated_ball_cosh_radius_square': C.enc(),
                           'diameter_margin': diameter_margin.enc(),
                           'first_ineffective_label': 4 if j == 2 else 6})
        ray_records.append(record)
    out = {
        'status': 'PASS_EXACT_FOUR_SEPARATE_SINGLE_RAY_NO_COMPLETION_GATES',
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'source_geometry_fresh_replay': True,
        'support_minimum_relaxation': 'Nonzero old pairings have magnitudes >=1 (first six) or >=2/5 (polar).',
        'singular_six_constraint_omissions': singular,
        'support_count': len(supports),
        'empty_relaxed_polyhedra': sum(not r['relaxed_polyhedron_vertices'] for r in supports),
        'minimum_lorentz_norm_square_by_zero_count': {str(n): x.enc() for n, x in sorted(minimum_by_zero_count.items())},
        'j3_m3_required_lorentz_norm_square': target.enc(),
        'j3_m3_possible_zero_sets_after_relaxation': admissible_zero_sets,
        'support_vertex_census': supports,
        'strict_old_interior_gap_point': enc(x),
        'j3_m3_ball_cosh_radius_square': C3.enc(),
        'ten_strict_vertex_ball_gaps': enc(gaps),
        'diameter_vertex_indices': [1, 8],
        'diameter_lower_cosh_square': diameter_lower_cosh_square.enc(),
        'four_separate_ray_gates': ray_records,
        'scope': 'Fix one j in 0..3. Any finite Coxeter-compatible set of added walls with H4 projection -t*omega_j, t>=1, leaves infinite volume in the 11-wall seed. Old centers and allowed labels are unrestricted. Mixing different rays, using mixed profiles, or deleting seed walls is NOT excluded.',
    }
    path = ROOT/'results'/f'{TAG}.json'
    if verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], 'supports', len(supports), 'empty', out['empty_relaxed_polyhedra'])


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    run(ap.parse_args().verify)
