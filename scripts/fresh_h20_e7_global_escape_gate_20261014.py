#!/usr/bin/env python3
"""Every Coxeter facet added to the saved E7 face retains the global escape.

Exact all-magnitude proof. Uses independently implemented recursive field
arithmetic, no numerical optimizer, no finite angle or height cutoffs.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_h20_e8_interior_caps_20261014_exact_review as K

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_e7_global_escape_gate_20261014'
CHAMBER = 'results/fresh_h20_e8_interior_caps_20261014_exact.json'
FACE = 'results/fresh_h20_e7_face_20261014.json'
Z, c, add, sub, mul, scale, sign = K.Z, K.scalar, K.add, K.sub, K.multiply, K.scale, K.sign


def total(xs):
    out = Z
    for x in xs:
        out = add(out, x)
    return out


def encode(x):
    if isinstance(x, F): return str(x)
    if isinstance(x, dict): return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [encode(v) for v in x]
    return x


def run():
    d, f = [json.loads((ROOT/p).read_text()) for p in [CHAMBER, FACE]]
    Q = [[F(x) for x in r] for r in d['form']]
    def dot(v, w):
        qw = [tuple(sum(Q[i][j]*w[j][k] for j in range(21)) for k in range(8)) for i in range(21)]
        return total(mul(a, b) for a, b in zip(v, qw))
    N = list(map(K.vector, d['normals']))
    E = K.vector(d['strict_escape']['point'])
    T = K.vector(f['intrinsic_normalization']['timelike_basis'])
    omega = K.vector(f['unique_nonright_meeting_profile']['fundamental_weight'])
    G = [[F(x) for x in r] for r in f['face_gram']]
    assert K.pd(G)
    J = K.inverse(G)
    A = f['face_walls']
    assert len(A) == 13
    corner = f['E7_support']+A
    corner_gram = [[K.rational(dot(N[i], N[j])) for j in corner] for i in corner]
    assert len(corner) == 20 and K.pd(corner_gram)
    corner_inverse = K.inverse(corner_gram)
    assert all(x >= 0 for row in corner_inverse for x in row)
    assert all(dot(T, N[i]) == Z for i in corner) and dot(T, T) == c(-2)
    assert sign(dot(T, K.vector(d['strict_center']))) < 0
    W = [tuple(total(scale(N[j][k], q) for j, q in zip(A, row)) for k in range(21)) for row in J]
    assert all(dot(w, N[j]) == c(int(i == k)) for i, w in enumerate(W) for k, j in enumerate(A))
    a = [scale(dot(w, E), -1) for w in W]
    B, w = scale(dot(T, E), -1), scale(dot(omega, E), -1)
    assert all(sign(x) > 0 for x in a+[B, w])
    coords = [list(map(K.field, row)) for row in f['intrinsic_normalization']['projected_wall_coordinates']]
    assert all(dot(T, n) == scale(row[13], -2) for n, row in zip(N, coords))
    assert all(dot(lam, n) == row[i] for n, row in zip(N, coords) for i, lam in enumerate(W))
    ow = [dot(omega, n) for n in N]
    # A normalized facet has r=-sum d_i lambda_i+tT-lift*omega,
    # with lift=0 or1; norm2; d_i in {0,1} union [sqrt2,infinity).
    # Positive t follows from existence of any strict future center p:
    # write p=uT-sum b_i Omega_i, u>0,b_i>0, and
    # r=tT-sum c_i Omega_i, c_i>=0. Then
    # r.p=-2tu+c'corner_inverse*b cannot be negative if t<=0.
    # Thus a valid acute addition automatically keeps the old finite apex.
    pair_coeff = [[scale(v, -1) for v in row[:13]]+[scale(row[13], -2)] for row in coords]
    mu = {30: F('0.002038242873333'), 75: F('0.002461273481333'), 76: F('1.639210640401')}
    residual = [sub(v, total(scale(pair_coeff[j][i], q) for j, q in mu.items())) for i, v in enumerate(a+[scale(B, -1)])]
    assert all(sign(v) < 0 for v in residual)
    # Hence every acute right lift pairs negatively with E, even without
    # imposing its norm or the discrete angle restrictions.
    rows, rhs = {}, {}
    for j in [16, 76]:
        assert sign(coords[j][13]) > 0
        rows[j] = [sub(scale(mul(B, coords[j][i]), -1), scale(mul(coords[j][13], a[i]), 2)) for i in range(13)]
        rhs[j] = add(mul(B, ow[j]), scale(mul(coords[j][13], w), 2))
    # If the nonright lift cuts E, t < (a.d+w)/B. Substitution in its
    # old-wall inequalities gives rows[j].d <= rhs[j].
    mixed = [add(scale(rows[16][i], F(194, 1000)), scale(rows[76][i], F(843, 1000))) for i in range(13)]
    mixed_rhs = add(scale(rhs[16], F(194, 1000)), scale(rhs[76], F(843, 1000)))
    assert all(sign(v) > 0 for v in mixed)
    assert sign(mixed_rhs) > 0 and sign(sub(mixed[12], mixed_rhs)) > 0
    # d12<1, so its Coxeter gap forces d12=0.
    assert sign(rows[16][0]) < 0 and sign(sub(rows[16][0], rhs[16])) > 0
    assert all(sign(rows[16][i]) > 0 for i in range(1, 12))
    # With d12=0, the wall16 inequality forces d0>1, hence d0>=sqrt2.
    B2 = mul(B, B)
    aGa = total(scale(mul(a[i], a[j]), G[i][j]) for i in range(13) for j in range(13))
    positive_gap = sub(B2, scale(aGa, 2))
    assert sign(positive_gap) > 0
    # f(d)=d'Jd-2(a.d+w)^2/B^2 is strictly convex because
    # J-2aa'/B^2 is positive definite (rank-one update criterion).
    root2 = (F(0), F(1))+(F(0),)*6
    C = add(mul(a[0], root2), w)
    D = sub(B2, scale(mul(a[11], a[11]), 3))
    assert sign(D) > 0 and sign(C) > 0
    # On d>=0,d12=0,d0>=sqrt2 the minimizer has d0=sqrt2,
    # d11=3*a11*C/D, and all other entries zero. The remaining
    # gradient coordinates are positive, as checked without division.
    gradients = [sub(scale(mul(D, root2), J[i][0]), scale(mul(a[i], C), 2)) for i in range(11)]
    assert all(sign(v) > 0 for v in gradients)
    assert J[0][0] == F(11, 12) and J[11][11] == F(2, 3)
    assert all(J[i][11] == 0 for i in range(11))
    margin = sub(scale(D, F(4, 3)), scale(mul(C, C), 2))
    assert sign(margin) > 0
    # min f=11/6-2C^2/D >1/2. But cut(E) implies y^2>f(d),
    # contradicting the nonright lift's required intrinsic norm1/2.
    paths = [CHAMBER, FACE, f'scripts/{TAG}.py',
             'scripts/fresh_h20_e8_interior_caps_20261014_exact_review.py']
    return encode(dict(status='PASS_ALL_E7_FACE_FACET_LIFTS_RETAIN_GLOBAL87_ESCAPE',
        sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        field_radicands=[1, 2, 79, 158, 11, 22, 869, 1738],
        finite_corner=dict(walls=corner, inverse=corner_inverse,
            orientation_argument='Nonnegative corner inverse plus any strict future center forces every added acute normal to have positive timelike coefficient; the old finite apex is automatically retained.'),
        escape_norm=dot(E, E), target_coefficients=dict(a=a, B=B, w=w),
        right_lift=dict(multipliers={str(i):v for i,v in mu.items()}, residual_coefficients=residual),
        nonright_lift=dict(cut_boundary_rows={str(j):row for j,row in rows.items()},
            cut_boundary_rhs={str(j):v for j,v in rhs.items()},
            mixed_coefficients=mixed, mixed_rhs=mixed_rhs,
            rank_one_positive_gap=positive_gap, denominator=D,
            C=C, positive_gradient_numerators=gradients,
            minimum_norm_minus_half_times_D=margin),
        conclusion='Any valid wall that is an actual additional facet of the retained E7 face and is Coxeter with the87 old walls retains the saved global spacelike escape.',
        scope='No arbitrary-completion exclusion. Walls not producing an actual facet on this face, or constructions replacing old walls, are outside this theorem.'))


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
