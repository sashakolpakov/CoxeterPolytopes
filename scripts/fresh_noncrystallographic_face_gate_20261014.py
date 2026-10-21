#!/usr/bin/env python3
"""Exact projected seed geometry for the noncrystallographic ridge gate.

The general compact-face and dimension deductions use the cited theorems;
these finite calculations are their application to the weighted H20 core.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_noncrystallographic_face_gate_20261014'


def q(a=0, b=0): return F(a), F(b)
def add(x, y): return x[0]+y[0], x[1]+y[1]
def neg(x): return -x[0], -x[1]
def sub(x, y): return add(x, neg(y))
def mul(x, y):
    return x[0]*y[0]+x[1]*y[1], x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
def div(x, y):
    yc = (y[0]+y[1], -y[1])
    d = mul(y, yc)
    assert not d[1] and d[0]
    z = mul(x, yc)
    return z[0]/d[0], z[1]/d[0]
def sign(x):
    a, b = x[0]+x[1]/2, x[1]/2
    if not b: return (a>0)-(a<0)
    if not a: return (b>0)-(b<0)
    if a*b>0: return (a>0)-(a<0)
    d = a*a-5*b*b
    assert d
    return ((a>0)-(a<0))*((d>0)-(d<0))
def enc(v):
    if isinstance(v, F): return str(v)
    if isinstance(v, (tuple, list)): return [enc(x) for x in v]
    if isinstance(v, dict): return {str(k): enc(x) for k, x in v.items()}
    return v


def run():
    phi = q(0, 1)
    phi2 = mul(phi, phi)
    k = sub(q(4), phi2)
    assert k == q(3, -1) and sign(k)>0
    # S consists of the norm-four central wall and a UNIT golden seed.
    S = [[q(4), neg(phi)], [neg(phi), q(1)]]
    inverse = [[div(q(1), k), div(phi, k)],
               [div(phi, k), div(q(4), k)]]
    for i in range(2):
        for j in range(2):
            t = add(mul(S[i][0], inverse[0][j]),
                    mul(S[i][1], inverse[1][j]))
            assert t == q(int(i == j))
    unit_square = sub(q(1), inverse[0][0])
    heavy_square = sub(q(1), mul(q(2), inverse[0][0]))
    other_gold_square = sub(q(1), mul(phi2, inverse[0][0]))
    opposite_gold_square = div(mul(q(-4), phi2), k)
    assert sign(unit_square)>0
    assert all(sign(x)<0 for x in
               (heavy_square, other_gold_square, opposite_gold_square))
    distinct_unit_cosine = div(inverse[0][0], unit_square)
    opposite_unit_cosine = div(add(q(1), inverse[0][0]), unit_square)
    assert distinct_unit_cosine == phi2
    assert opposite_unit_cosine == q(3, 2)
    assert sign(sub(distinct_unit_cosine, q(1)))>0
    assert sign(sub(opposite_unit_cosine, q(1)))>0
    # A new wall right to S retains its square on restriction. An old
    # right contact stays right; an old parallel/dotted cosine increases
    # since the positive unit-seed projection has square strictly below1.
    assert sign(sub(q(1), unit_square))>0
    total_seed_weight = q(20, 8)
    prospective_central_vertex_square = sub(total_seed_weight, q(4))
    actual_ridge_point_square = sub(phi2, q(4))
    assert sign(prospective_central_vertex_square)>0
    assert sign(actual_ridge_point_square)<0
    return enc(dict(
        status='PASS_EXACT_WEIGHTED_CORE_I2_5_FACE_GEOMETRY',
        sources={f'scripts/{TAG}.py': sha256(Path(__file__).read_bytes()).hexdigest()},
        spherical_pair_gram=S, determinant=k, inverse=inverse,
        projected_seed_squares=dict(unit=unit_square, heavy=heavy_square,
            other_gold=other_gold_square, opposite_target_gold=opposite_gold_square),
        actual_face_facets=20,
        unit_face_pair_cosines=dict(distinct_coordinates=distinct_unit_cosine,
                                   opposite_sides=opposite_unit_cosine),
        all_old_face_pairs='strictly dotted',
        central_vertex_square=prospective_central_vertex_square,
        explicit_ridge_point_square=actual_ridge_point_square,
        compact_face_dimension_if_finite_volume=18,
        right_angled_compact_dimension_upper_bound=4,
        general_scope='The compact Coxeter (n-2)-face implication for finite labels outside2,3,4,6 uses Borcherds/Allcock plus the Euclidean cusp trace argument. This is not a computational nonexistence theorem for general H20 completions.',
        restricted_completion_obstruction='If all extra walls restricting to this face are right to the chosen I2(5), and all finite contacts among them and with the twenty old unit seed walls are right, the induced18-face is right-angled. A finite-volume completion would make it compact, contradicting the compact right-angled dimension bound. General changed-angle or changed-facet completions are outside this obstruction.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert json.loads(path.read_text()) == result
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'])
