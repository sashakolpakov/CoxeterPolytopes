#!/usr/bin/env python3
"""Complete all-label coherent endpoint locus on the H4+H4 [7,30] conic.

Exact physical recipes use the frozen outward-interval engine. The only new
cosine bounds use rational Machin and alternating Taylor certificates.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

import fresh_nonarithmetic_h3_endpoint_domains_20261017 as E

P, U, I = E.P, E.U, E.I
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_h3_edge730_domain_20261017'


def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def arctan_bounds(denominator):
    z = F(1, denominator)
    value = sum(((-1) ** k * z ** (2*k+1) / (2*k+1) for k in range(32)), F(0))
    # The last term is negative; the next positive term is an upper bound.
    return value, value + z ** 65 / 65


def pi_bounds():
    a, b = arctan_bounds(5), arctan_bounds(239)
    lo, hi = 16*a[0]-4*b[1], 16*a[1]-4*b[0]
    scale = 1 << 128
    return F((lo*scale).__floor__(), scale), F((hi*scale).__ceil__(), scale)


PI = pi_bounds()


def cosine_bounds(m):
    assert m >= 3
    def at(x):
        term = F(1)
        value = term
        for k in range(1, 32):
            term *= -x*x/F((2*k-1)*(2*k))
            value += term
        next_term = -term*x*x/F(63*64)
        assert next_term > 0
        return value, value + next_term
    # cos decreases on these positive intervals; alternating terms decrease.
    lo = at(PI[1]/m)[0]
    hi = at(PI[0]/m)[1]
    return I(2*lo, 2*hi)


def theta(m):
    if m in (2, 3, 4, 5, 6, 10):
        return U.theta(m)
    if m == 8:
        return U.sqrt(2 + U.sqrt(I(2)))
    return cosine_bounds(m)


def run():
    endpoint = E.B.read(f'results/{E.TAG}.json')
    domain = next(row for row in endpoint['endpoint_domains'] if row['sides'] == [7, 30])
    assert domain['profiles']['positive_residual_profiles'] == 1
    assert domain['profiles']['finite_label_profiles'][0]['labels'] == [2]*4
    assert domain['profiles']['tails'] == []
    _, qi, normals, h, _, _, _, _ = P.physical_data()
    pair = lambda x, y: U.idot(x, U.imv(qi, y))
    k0 = [I(*x) for x in domain['kernel_future_unit']]
    k1 = [I(*x) for x in domain['kernel_positive_unit']]
    v = [I(*x) for x in endpoint['common_interior_vertex']]
    frame = domain['full_support']
    assert sorted(frame) == [0, 1, 2, 3, 5, 6, 7, 30]
    # k0^2=-1, k1^2=1 and k0.k1=0 are exact construction identities.
    # Thus A,B are null, A.B=1 and r(t)=t*A+B/t has square2.
    sq2 = U.sqrt(I(2))
    a = [U.rnd((x-y)/sq2) for x, y in zip(k0, k1)]
    b = [U.rnd((-x-y)/sq2) for x, y in zip(k0, k1)]
    other_a = [U.rnd((x+y)/sq2) for x, y in zip(k0, k1)]
    other_b = [U.rnd((-x+y)/sq2) for x, y in zip(k0, k1)]
    other_47 = [pair(other_a, normals[47]), pair(other_b, normals[47])]
    assert all(x.l > 0 for x in other_47)
    coefficients = {i: (pair(a, n), pair(b, n)) for i, n in enumerate(normals) if i not in frame}
    assert len(coefficients) == 47
    assert all(c[0].h < 0 for c in coefficients.values())
    alpha, beta = coefficients[51]
    assert alpha.h < 0 < beta.l

    def parameter(c):
        # alpha*t^2+c*t+beta=0 has exactly one positive root.
        return U.rnd((-c-U.sqrt(c*c-4*alpha*beta))/(2*alpha))

    def root(t):
        return [U.rnd(t*x+y/t) for x, y in zip(a, b)]

    def endpoint_point(t):
        return [U.rnd(t*x-y/t) for x, y in zip(a, b)]

    t8, tinf = parameter(theta(8)), parameter(I(2))
    assert t8.l > 5 and t8.h < tinf.l < 6
    certified_tail = []
    for wall, (aa, bb) in coefficients.items():
        # At t8 every wall except51 is already strictly dotted. The numerator
        # of g(t)+2 is decreasing thereafter, as certified by its derivative.
        if wall != 51:
            dotted_value = U.rnd(aa*t8*t8+2*t8+bb)
            derivative = U.rnd(2*aa*t8+2)
            assert dotted_value.h < 0 and derivative.h < 0
        else:
            dotted_value = derivative = None
        # The actual finite endpoint x(t)=t*A-B/t has square -2. Every old
        # nonincident inequality is strict there for every t>=t8.
        vertex_value = U.rnd(aa*t8*t8-bb)
        assert vertex_value.h < 0
        certified_tail.append(dict(wall=wall, coefficients=P.enc([aa, bb]),
                                   dotted_numerator_at_t8=None if dotted_value is None else dotted_value.enc(),
                                   dotted_derivative_at_t8=None if derivative is None else derivative.enc(),
                                   endpoint_numerator_at_t8=vertex_value.enc()))
    va, vb = pair(a, v), pair(b, v)
    interior_at_t8 = U.rnd(va*t8*t8+vb)
    assert va.h < 0 and interior_at_t8.h < 0
    ah, bh = pair(a, h), pair(b, h)
    assert ah.h < 0 < bh.l

    # The remaining bounded part has only the possibilities m=2,...,7 at51.
    # The coefficient signs make g51 strictly decreasing, so these exhaust it.
    gaps = []
    blockers = {2: (49, None), 3: (49, (2, 3)), 4: (49, (4, 5)),
                5: (50, (6, 7)), 6: (50, (10, 11)), 7: (50, (24, 25))}
    for m, (wall, bracket) in blockers.items():
        t = parameter(theta(m))
        assert t.h < t8.l
        aa, bb = coefficients[wall]
        product = U.rnd(aa*t+bb/t)
        if bracket is None:
            assert product.l > 0
            decision = 'OBTUSE'
        else:
            left, right = map(theta, bracket)
            assert (product+left).h < 0 and (product+right).l > 0
            decision = 'CONSECUTIVE_LABEL_GAP'
        gaps.append(dict(label_to51=m, parameter=t.enc(), wall=wall,
                         product=product.enc(), decision=decision,
                         consecutive_labels=bracket))

    # The parallel endpoint is a genuine coherent cap but cannot occur in a
    # finite-volume acute completion: its positive null sum with51 is incident
    # to the retained pi/5 pair0,1.
    rinf = root(tinf)
    nu = [U.rnd(x+y) for x, y in zip(rinf, normals[51])]
    matrix = E.labels()
    assert matrix[0][1] == 5 and matrix[0][51] == matrix[1][51] == 2
    assert [matrix[i][51] for i in P.S] == [2, 2, 2, 2, 6, 4]
    wall51_residual_square = U.rnd(2-I(F(2, 3))*(5+U.sqrt(I(6))))
    assert wall51_residual_square.h < 0
    null_future = pair(nu, v)
    assert null_future.h < 0
    # Null norm and orthogonality are exact: rinf.51=-2, both norms2,
    # and both are perpendicular to old0,1 by the literal recipes.
    paths = [f'scripts/{TAG}.py', f'scripts/{E.TAG}.py', f'results/{E.TAG}.json',
             f'scripts/{P.TAG}.py', f'results/{P.TAG}.json', P.U.CAP, P.C53.CL,
             P.CAP52, P.CAPS55, P.CONE55,
             'results/fresh_h9_fixed_core_cusp_obstruction_20261013.md']
    return dict(status='PASS_COMPLETE_ALL_LABEL_H4_SQUARED_ENDPOINT_LOCUS',
                sources={p: sha(p) for p in paths}, retained_walls=55, old_edge=[7, 30],
                right_support=frame, A=P.enc(a), B=P.enc(b),
                exact_normal_formula='r(t)=t*A+B/t, t>0; A^2=B^2=0,A.B=1.',
                exact_vertex_formula='x(t)=t*A-B/t; x(t)^2=-2, x(t).r(t)=0.',
                other_component_obtuse_wall=47, other_component_positive_coefficients=P.enc(other_47),
                native_coefficient_certificates=certified_tail,
                controlling_wall=51, controlling_alpha_beta=P.enc([alpha, beta]),
                parameter_formula='t(c)=(-c-sqrt(c^2-4*alpha*beta))/(2*alpha), c=2cos(pi/m); t_infinity=t(2).',
                first_allowed_label=8, t8=t8.enc(), t_infinity=tinf.enc(),
                low_label_rejections=gaps, pi_interval=[str(x) for x in PI],
                cosine_isolators={str(m): theta(m).enc() for m in (7, 8, 11, 24, 25)},
                cosine_proof='Machin pi=16atan(1/5)-4atan(1/239),32 alternating terms with next-term brackets, rounded outwards to128bits. Cosine uses32 alternating terms, next-term remainder, and monotonicity on[0,pi/3].',
                common_interior_vertex=P.enc(v), interior_alpha_beta=P.enc([va, vb]),
                interior_numerator_at_t8=interior_at_t8.enc(), endpoint_future_alpha_beta=P.enc([ah, bh]),
                parallel_positive_null=P.enc(nu), parallel_null_future_product=null_future.enc(),
                parallel_forbidden_angle_edge=[0, 1],
                wall51_S6_projected_square=wall51_residual_square.enc(),
                incidence_caution='Old51 has labels6,4 to A2 nodes5,6, and negative S6-projected square. Therefore r.51 finite is not an edge of the actual S6 face. H4+[r,51] is a different positive support for finite m.',
                exact_coherent_locus='Exactly t=t_m for integers m>=8, or t>=t_infinity. At t_m only51 is nonright finite (labelm); at t_infinity only51 is parallel; above it all47 nonframe old products are strictly dotted. Every branch is right to the eight support roots.',
                finite_volume_compatible_locus='Remove t=t_infinity, whose affineA1 positive-null vector is perpendicular to the retained I2(5). This gives candidates, not a promise of finite-volume completion.',
                actual_vertex='For every admitted parameter x(t) is future timelike, right to the eight old roots and r(t), and strict against all other old walls. Thus each root gives an actual finite endpoint of the selected old edge.',
                common_interior='All admitted roots have negative product with the same persistent old finite vertex. Any finite collection therefore has a common strict future interior by sufficiently small perturbation toward the old center. However two distinct roots on this conic have product t/s+s/t>=2, so no two can belong to an acute collection.',
                scope='Complete all-label endpoint locus for the unique H4+H4 old edge family; no finite label cutoff and no restriction to an arithmetic field. No new wall is appended. Coordinated closure needs other supports; compactness of the entire H3 remains unresolved.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    data = run()
    # Preserve ordinary JSON arrays in the reproducible result.
    data = json.loads(json.dumps(data))
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(path.read_text())
        print('PASS saved complete endpoint-locus replay')
    else:
        path.write_text(json.dumps(data, indent=2) + '\n')
    print(data['status'], 'all m>=8 plus strict dotted tail')
