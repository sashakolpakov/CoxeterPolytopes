#!/usr/bin/env python3
"""Universal depth17 finite-distance8 corridor: exact scalar certificate.

All-label minimal-core and actual-face facts are explicit inherited inputs.
The universal inequalities use Fraction only, with no project imports.
No payment budget or ambient dimension conclusion is asserted.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_universal_deep_finite8_corridor_20261018'
DEPTH, OUTWARD, N, LEFT, RIGHT = 17, 15, 39, 15, 23


def j(n, i, k):
    return F((min(i, k)+1)*(n-max(i, k)), n+1)


def inverse_checks():
    rows = []
    for n in (DEPTH, N):
        determinants = [F(1), F(2)]
        for size in range(2, n+1):
            determinants.append(2*determinants[-1]-determinants[-2])
        assert determinants == [F(k+1) for k in range(n+1)]
        for i in range(n):
            for k in range(n):
                value = 2*j(n, i, k)
                if i:
                    value -= j(n, i-1, k)
                if i+1 < n:
                    value -= j(n, i+1, k)
                assert value == int(i == k)
                assert j(n, i, k) > 0
        rows.append(dict(order=n, inverse_product_entries_checked=n*n,
                         positive_inverse_entries=n*n,
                         positive_leading_principal_determinants=list(map(str, determinants[1:]))))
    return rows


def run():
    # The order-ten ceiling is a published/inherited classification input;
    # it is not inferred from the size of a bounded catalogue below.
    core_order_bound = 10
    after_omitted_black_order = core_order_bound-1
    root_distance_bound = after_omitted_black_order-1
    assert root_distance_bound == 8 < OUTWARD
    assert RIGHT-LEFT == 8 and N == 2*OUTWARD+9
    edge_gates = []
    for distance in range(3, OUTWARD+1):
        diagonal = j(DEPTH, distance-1, distance-1)
        assert diagonal > 2
        edge_gates.append(dict(outward_distance=distance,
                               A17_inverse_diagonal=str(diagonal),
                               minimum_native_projection=str(diagonal)))
    assert min(F(r['A17_inverse_diagonal']) for r in edge_gates) == F(5, 2)
    short_routes = [dict(inner_distance=distance, route_length=9-distance)
                    for distance in range(2, 8)]
    assert all(r['route_length'] <= 7 for r in short_routes)

    left = list(range(13, 17))
    right = list(range(22, 26))
    assert left == [LEFT-2, LEFT-1, LEFT, LEFT+1]
    assert right == [RIGHT-1, RIGHT, RIGHT+1, RIGHT+2]
    alpha = F(23, 680)
    block_bounds = []
    for side, positions in [('left', left), ('right', right)]:
        coordinate = (lambda i: i+1) if side == 'left' else (lambda i: N-i)
        for i in positions:
            for k in positions:
                entry = j(N, i, k)
                lower = alpha*coordinate(i)*coordinate(k)
                assert entry >= lower > 0
                block_bounds.append(dict(side=side, i=i, k=k,
                    inverse_entry=str(entry), rank_one_lower_bound=str(lower),
                    difference=str(entry-lower)))
    cross = []
    for i in left:
        for k in right:
            value = F((i+1)*(N-k), N+1)
            assert value == j(N, i, k)
            cross.append(dict(left=i, right=k, inverse_entry=str(value),
                              left_factor=i+1, right_factor=N-k))
    diagonal_rows = [dict(position=i, inverse_diagonal=str(j(N, i, i)))
                     for i in left+right]
    q0 = min(F(r['inverse_diagonal']) for r in diagonal_rows)
    assert q0 == F(91, 10) > 2
    rho2 = 1/(F((N+1)**2)*alpha**2)
    kappa = 1-rho2
    assert rho2 == F(289, 529) and kappa == F(240, 529)
    # At q=q0+u, r=q0+v, the determinant lower bound is this
    # polynomial with strictly positive coefficients for u,v>=0.
    coefficients = dict(constant=kappa*q0*q0-4*q0+4,
                        u=kappa*q0-2, v=kappa*q0-2, uv=kappa)
    assert coefficients == dict(constant=F(13674, 2645),
                                 u=F(1126, 529), v=F(1126, 529), uv=F(240, 529))
    assert all(value > 0 for value in coefficients.values())

    paths = [Path(__file__)]
    for stem in ('fresh_high_rank_quasilanner_census_20261017',
                 'fresh_full_a39_corridor_high_domain_20261017',
                 'fresh_distance9_low_core_20261017',
                 'fresh_selective_distance9_dimension108_20261017',
                 'fresh_selective_distance10_dimension116_20261017',
                 'fresh_distance8_low_core_20261017',
                 'fresh_distance8_high_disjoint_20261017'):
        paths.append(ROOT/f'results/{stem}.md')
    paths.append(ROOT/'data/references/fresh_finite_edge_cusp_tradeoff_20261017_chein.pdf')
    return dict(status='PASS_UNIVERSAL_DEPTH17_FINITE8_SHORT7_ALTERNATIVE',
        sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
        parameters=dict(feature_depth=DEPTH, retained_outward_whites=OUTWARD,
                        path_order=N, marked_positions=[LEFT, RIGHT],
                        left_contact_support=left, right_contact_support=right),
        minimal_core_input=dict(order_bound=core_order_bound,
            after_black_deletion_component_order_at_most=after_omitted_black_order,
            distance_from_retained_mark_at_most=root_distance_bound,
            primary_location='Chein1969 printedpage4, PDFpage3: irreducible hyperbolic order<=10 and the finite/affine-component criterion for every vertex deletion.',
            primary_url='https://www.numdam.org/item/M2AN_1969__3_3_3_0.pdf',
            classification_bridge='Order2 dotted cores already satisfy the ceiling. At order>=3 no dotted edge is permitted by minimality; finite/parallel Coxeter labels remain. A strictPD black deletion and interlacing give signature(order-1,1), while all other deletions have finite/affine connected components, exactly the stated Chein criterion.',
            scope='Published all-label order-ten ceiling inherited from the bound reports; no enumeration in this source establishes that ceiling.'),
        inverse_formula_checks=inverse_checks(), actual_edge_gates=edge_gates,
        inner_contact_shortcuts=short_routes, diagonal_lower_bounds=diagonal_rows,
        positive_block_bounds=block_bounds, cross_block_factorizations=cross,
        inequalities=dict(rank_one_block_factor=str(alpha), q_r_lower_bound=str(q0),
            correlation_squared_upper_bound=str(rho2), determinant_qr_coefficient=str(kappa),
            diagonal_schur_upper_bound=str(2-q0),
            positive_shifted_determinant_coefficients={k:str(v) for k,v in coefficients.items()},
            strict_determinant_lower_bound=str(coefficients['constant'])),
        theorem='For any actual quadrilateral star, a marked distance8 pair in a whole finite A/B/D component with both marks at feature depth>=17 forces a full-star path of length<=7 between the two opposite black pairs, with white internal roots. This covers low/high, disjoint/overlapping cores, all native labels and arbitrary supporting roots.',
        scope='Universal local short-path alternative only; no changed payment system, finite budget or global dimension bound.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    data = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(data, indent=2)+'\n')
    print(data['status'], 'strict determinant bound',
          data['inequalities']['strict_determinant_lower_bound'], flush=True)
