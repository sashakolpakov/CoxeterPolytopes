#!/usr/bin/env python3
"""All-label obstruction to one-wall triangular closure of the r8/S6 face."""
import argparse
import itertools
import json
from hashlib import sha256
from pathlib import Path

import sympy as sy
import fresh_nonarithmetic_pi8_tetrahedron_completion_20261017 as C

P, U, I, D = C.P, C.U, C.I, C.D
E = D.E
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nonarithmetic_pi8_triangle_facet_gate_20261017'


def run():
    z = C.geometry()
    n, pair = z['normals'], z['pair']
    frame_ids = [0, 1, 2, 3, 5, 6, 55, 7, 30]
    frame = [n[i] for i in frame_ids]
    gram = [[pair(a, b) for b in frame] for a in frame]
    support_pd = P.pd(gram)
    inverse, pivots = P.C53.inverse(gram)
    duals = [P.combine([inverse[j][i] for j in range(9)], frame) for i in range(4, 9)]
    kernel = P.sub(n[51], P.combine(U.imv(inverse, [pair(n[51], x) for x in frame]), frame))
    kernel_square = pair(kernel, kernel)
    assert kernel_square.h < 0
    if pair(kernel, z['h']).l > 0:
        kernel = [-x for x in kernel]
    assert pair(kernel, z['h']).h < 0
    scale = U.sqrt(-kernel_square)
    kernel = [U.rnd(x/scale) for x in kernel]

    # The two endpoint supports are H4+(A3+A1) and H4+(A2+A1+A1).
    # H4 cannot be attached in a finite diagram, so only four contacts remain.
    old_labels = E.labels()
    def label(i, j):
        return 2 if 55 in (i, j) else old_labels[i][j]
    domains = []
    for side in (7, 30):
        ids = [5, 6, 55, side]
        small = sy.Matrix([[2 if i == j else -E.THETA[label(i, j)] for j in ids] for i in ids])
        profiles = E.profiles(small)
        domains.append(dict(side=side, nodes=ids, profiles=profiles))
    finite = sorted({tuple(a['labels']+[b['labels'][3]])
                     for a in domains[0]['profiles']['finite_label_profiles']
                     for b in domains[1]['profiles']['finite_label_profiles']
                     if a['labels'][:3] == b['labels'][:3]})
    assert len(finite) == 46
    # All tails in a rank4 positive support are one-contact isolated-A1
    # tails. Matching the shared entries leaves one at55, and four at30.
    assert domains[0]['profiles']['isolated_A1_slots'] == [2]
    assert domains[1]['profiles']['isolated_A1_slots'] == [2, 3]
    assert sorted(r['labels'][3] for r in domains[0]['profiles']['finite_label_profiles']
                  if r['labels'][:3] == [2, 2, 2]) == [2, 3, 4, 5]
    tails = []
    for side, domain in enumerate(domains):
        for tail in domain['profiles']['tails']:
            slot = tail['slot']
            if slot < 3:
                assert slot == 2
            else:
                assert side == 1
    finite_records = []
    for ls in finite:
        projection = P.combine([-D.theta(m) for m in ls], duals)
        residual = U.rnd(pair(projection, projection)-2)
        if residual.h < 0:
            finite_records.append(dict(labels=ls, residual=residual.enc(), status='NO_REAL_NORM2_ROOT'))
            continue
        assert residual.l > 0
        branches = []
        for sign in (-1, 1):
            normal = [U.rnd(x+sign*U.sqrt(residual)*y) for x, y in zip(projection, kernel)]
            blockers = [(i, U.rnd(pair(normal, v))) for i, v in enumerate(n) if i not in frame_ids]
            wall, value = next((i, v) for i, v in blockers if v.l > 0)
            branches.append(dict(kernel_sign=sign, obtuse_wall=wall, product=value.enc()))
        finite_records.append(dict(labels=ls, residual=residual.enc(), status='BOTH_BRANCHES_OBTUSE', branches=branches))
    assert sum(r['status'] == 'NO_REAL_NORM2_ROOT' for r in finite_records) == 4

    # Tail at55: p=-c*n55/2, p^2=c^2/2<2, so the negative one-dimensional
    # kernel can never raise its square to2. This is exact for every finite m.
    tails.append(dict(nonright_node=55, all_other_contacts=2,
                      parameter='sqrt(3)<=c=2cos(pi/m)<2',
                      projection_square='c^2/2', status='NO_REAL_NORM2_ROOT'))

    # For the four tails at30, p=a*du+c*dv. The future branch is blocked by49.
    # Squaring is used only after certifying the linear term is positive.
    du, dv = [[-x for x in duals[k]] for k in (3, 4)]
    uu, uv, vv = [U.rnd(x) for x in (pair(du, du), pair(du, dv), pair(dv, dv))]
    assert uu.l > 0 and uv.l > 0 and vv.l > 0
    lower, upper = U.sqrt(I(3)), I(2)
    def polynomial(q, x):
        return U.rnd(q[0]+q[1]*x+q[2]*x*x)
    def bernstein(q):
        a = polynomial(q, lower)
        b = U.rnd(a+(upper-lower)*(q[1]+2*q[2]*lower)/2)
        c = polynomial(q, upper)
        assert all(x.l > 0 for x in (a, b, c))
        return P.enc([a, b, c])
    for m in (2, 3, 4, 5):
        a = D.theta(m)
        rho = [U.rnd(a*a*uu-2), U.rnd(2*a*uv), vv]
        assert polynomial(rho, lower).l > 0
        assert rho[1].l >= 0 and rho[2].l > 0
        proofs = []
        for sign, wall in ((1, 49), (-1, 4)):
            linear = [U.rnd(a*pair(du, n[wall])), U.rnd(pair(dv, n[wall]))]
            k = U.rnd(pair(kernel, n[wall]))
            assert k.h < 0
            assert (linear[0]+linear[1]*lower).l > 0
            assert (linear[0]+linear[1]*upper).l > 0
            if sign == 1:
                difference = [U.rnd(linear[0]*linear[0]-k*k*rho[0]),
                              U.rnd(2*linear[0]*linear[1]-k*k*rho[1]),
                              U.rnd(linear[1]*linear[1]-k*k*rho[2])]
                bounds = bernstein(difference)
            else:
                difference, bounds = None, None
            proofs.append(dict(kernel_sign=sign, obtuse_wall=wall,
                               linear_term=P.enc(linear), kernel_product=k.enc(),
                               squared_difference=None if difference is None else P.enc(difference),
                               positive_Bernstein_coefficients=bounds))
        tails.append(dict(nonright_node=30, label_to7=m, right_nodes=[0,1,2,3,5,6,55],
                          parameter='sqrt(3)<=c<=2, a=2cos(pi/m7)',
                          norm_residual_polynomial=P.enc(rho), branches=proofs,
                          status='BOTH_BRANCHES_OBTUSE_FOR_ENTIRE_TAIL'))
    paths = [f'scripts/{TAG}.py', f'scripts/{C.TAG}.py', f'results/{C.TAG}.json',
             f'scripts/{E.TAG}.py', f'results/{E.TAG}.json', f'scripts/{D.TAG}.py',
             f'results/{D.TAG}.json', f'scripts/{P.TAG}.py', f'results/{P.TAG}.json']
    return dict(status='PASS_NO_SINGLE_WALL_TRIANGULAR_CLOSURE_OF_R8_S6_FACET',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                face_support=[0,1,2,3,5,6,55], existing_sides=[7,30], full_frame=frame_ids,
                frame_positive_pivots=support_pd, frame_inverse_pivots=pivots,
                kernel_future_unit=P.enc(kernel), kernel_unnormalized_square=kernel_square.enc(),
                endpoint_domains=domains, finite_contact_profile_order=[5,6,55,7,30],
                finite_profile_count=46, finite_records=finite_records, high_label_tail_count=5,
                high_label_tails=tails,
                completeness='A compact triangle with the two retained sides7,30 and one additional actual facet has two finite new endpoints. Each full incident rank9 diagram is finite, so the H4 factor is orthogonal and the remaining four contacts belong to the complete positive-residual domains for A3+A1 or A2+A1+A1. In each domain a label>=6 must be the sole contact to an isolatedA1 (inverse diagonal1/2); every other diagonal is>=2/3 and any second nonzero contact makes cost>=2. Matching both endpoint domains gives exactly46 finite profiles and five isolated high-label tails. The nine fixed roots are PD, so the orthogonal kernel is negative one-dimensional; the norm equation has no root when p^2<2 and exactly the two displayed roots otherwise.',
                continuous_proof='For the four30-tails, p=a*du+c*dv and rho=p^2-2>0 on[sqrt3,2]. The opposite branch has positive pairing with4 because its linear term is positive and the future kernel product is negative. For the future branch, L49>0 and L49^2-k49^2*rho>0 follow from three strictly positive degree2 Bernstein coefficients. Therefore L49+k49*sqrt(rho)>0 throughout, including every finite label without a cutoff.',
                scope='No one added wall, in any field and with arbitrary finite labels, can make the entire actual S6+r8 two-face a compact triangle with its existing sides7,30. Nonright A2 and r8 contacts are included. A polygon using two or more additional sides, or a different retained-wall construction, remains open. This does not obstruct the already certified compact H4+I2(8) tetrahedron.')


if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data==json.loads(out.read_text())
        print('PASS saved all-label triangle-facet replay')
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'])
