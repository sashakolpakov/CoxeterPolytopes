#!/usr/bin/env python3
"""Resolve the A24 equality and lower deep finite7 feature depth to10."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import fresh_depth11_finite7_zero_transfer_20261018 as B

ROOT = B.ROOT
TAG = 'fresh_depth10_finite7_zero_transfer_20261018'
FRONT = 'fresh_ambient60_payment_frontier_20261018'
STEMS = (B.TAG, FRONT)


def corridor():
    n, a, b = 24, 8, 15
    assert b-a == 7 and a == n-1-b == 8
    left, right = list(range(6, 10)), list(range(14, 18))
    alpha = F(1, 25)
    arm = [dict(outward_distance=j, inverse_diagonal=str(B.inverse(10, j-1, j-1)))
           for j in range(3, 9)]
    assert min(F(r['inverse_diagonal']) for r in arm) == F(24, 11) > 2
    rows, equal, cross = [], {}, []
    for side, ids in (('left', left), ('right', right)):
        equality = []
        for i, j in product(ids, repeat=2):
            u, v = (i+1, j+1) if side == 'left' else (n-i, n-j)
            value = B.inverse(n, i, j)-alpha*u*v
            assert value > 0
            if i == j:
                assert value >= 2
                if value == 2:
                    equality.append(i)
            rows.append(dict(side=side, i=i, j=j, factors=[u, v], residual_entry=str(value)))
        equal[side] = equality
    assert equal == dict(left=[9], right=[14])
    for i, j in product(left, right):
        assert B.inverse(n, i, j) == F((i+1)*(n-j), 25)
        cross.append(dict(i=i, j=j, inverse_entry=str(B.inverse(n, i, j))))
    q, r, s = B.inverse(n, 9, 9), B.inverse(n, 14, 14), B.inverse(n, 9, 14)
    assert (q, r, s) == (F(6), F(6), F(4))
    assert (q-2)*(r-2)-s*s == 0 and 2-q == 2-r == -4
    # The kernel is in the full principal Gram, not merely the Schur block.
    gram = [[F(2 if i == j else -1 if i < n and j < n and abs(i-j) == 1 else 0)
             for j in range(n+2)] for i in range(n+2)]
    for i, j in ((9, 24), (14, 25)):
        gram[i][j] = gram[j][i] = -1
    radical = [B.inverse(n, i, 9)-B.inverse(n, i, 14) for i in range(n)] + [F(1), F(-1)]
    assert all(sum(x*y for x, y in zip(row, radical)) == 0 for row in gram)
    assert radical[0] == F(1, 5) and radical[-2:] == [F(1), F(-1)]
    # The genuine ninth outward white has only the old endpoint contact.
    # Strict actual-edge positivity permits only zero/simple c,d contacts.
    projection = B.inverse(10, 8, 8)
    assert projection == F(18, 11) and 2*projection == F(36, 11) > 2
    contacts = []
    for x, y in product((0, 1), repeat=2):
        value = -radical[0]-x+y
        assert value != 0
        contacts.append(dict(c_contact=x, d_contact=y, radical_pairing=str(value)))
    return dict(feature_depth=10, retained_outward_whites=8,
        minimal_indefinite_order_bound=10, retained_core_contact_radius=8,
        path_order=n, marked_positions=[a, b], left_support=left, right_support=right,
        strict_arm_gates=arm, inverse_checks=[B.path_check(10), B.path_check(n)],
        block_factor=str(alpha), residual_blocks=rows, equality_positions=equal,
        cross_factorizations=cross, determinant_lower_bound='0',
        unique_equality_contacts=dict(c={9:'1'}, d={14:'1'}),
        equality_Schur=[['-4','-4'],['-4','-4']],
        equality_principal_Gram=[list(map(str, row)) for row in gram],
        equality_principal_inertia=[24,1,1], exact_radical=list(map(str, radical)),
        omitted_white_outward_distance=9, omitted_white_old_contact=[0,'-1'],
        omitted_white_edge_inverse_diagonal=str(projection),
        nonsimple_projection_lower=str(2*projection), allowed_omitted_white_contact_magnitudes=[0,1],
        required_contact_difference='1/5', complete_contact_pairings=contacts,
        conclusion='Strict Schur determinant gives two negative directions; the sole equality contradicts the actual ninth white. Thus every newly zero finite7 corner gives a short6 actual corner.')


def run():
    docs = {s: json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem, d in docs.items():
        assert d['status'].startswith('PASS'), stem
        for p, h in d['sources'].items():
            assert sha256((ROOT/p).read_bytes()).hexdigest() == h, (stem, p)
    old = docs[B.TAG]
    assert old['feature_depth'] == 11 and old['finite_distance'] == 7
    assert old['triangle_minimum'] == 2 and old['quadrilateral_minimum'] == 1
    protected = old['triangle_interface']['retained_distance7_elementary_payments']
    assert len(protected) == 14
    assert all(len(r['feature_witness']['path'])-1 <= 5 < 10 for r in protected)
    special = old['triangle_interface']['exceptional_distance7_marked_leaf']
    assert special['fixed_simple_leaf'] in special['alternate_pair']
    front = docs[FRONT]
    assert front['finite']['exact_maximum'] == '413'
    assert ['D17','D17','D26'] in front['finite']['actual_attaining_forests']
    assert all(max(0, r-29) == 0 for r in (17,17,26))
    paths = [Path(__file__)] + [ROOT/f'{folder}/{s}.{ext}' for s in STEMS
        for folder, ext in (('scripts','py'), ('results','json'), ('results','md'))]
    return dict(status='PASS_DEPTH10_FINITE7_ZERO_LOCAL_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        current_baseline=old['current_baseline'], finite_distance=7, feature_depth=10, new_weight='0',
        corridor=corridor(), inherited_triangle_protected_pairs=14,
        triangle_minimum=2, quadrilateral_minimum=1,
        ideal_weights_unchanged=True, all_other_finite_weights_and_coupons_unchanged=True,
        deleted_pair_counts={k:'max(0,n-'+str(c)+')' for k,c in dict(A=27,B=28,D=29).items()},
        rank60_finite_maximum_by_monotonicity_and_attainment='413',
        unchanged_rank60_attainer=['D17','D17','D26'], new_dimension_bound=None,
        new_finite_volume_example=False,
        scope='Universal depth10 corridor equality proof, exact A10/A24 inverse/Schur/kernel/native-contact checks and inherited noncircular short6/triangle interfaces. No triangle census rerun, finite DP, ideal convolution, new ambient bound or physical example.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    data = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(path.read_text())
    else:
        assert not path.exists()
        path.write_text(json.dumps(data, indent=2)+'\n')
    print(data['status'], 'depth', data['feature_depth'], 'finite_maximum',
          data['rank60_finite_maximum_by_monotonicity_and_attainment'], flush=True)
