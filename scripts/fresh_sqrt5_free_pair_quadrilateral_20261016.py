#!/usr/bin/env python3
"""Exact continuous free-pairing obstruction for 210 fixed two-edge witnesses.

The two uniform edges are held fixed. The new-to-new product is arbitrary
real. Every value compatible even with a finite two-root corner gives rank23
and negative index2, so no quadrilateral closure in H21 is possible.
"""
import argparse
import json
from hashlib import sha256
from itertools import combinations
from pathlib import Path

from fresh_sqrt5_d4_two_e8_self_edge_20261016 import Q, dot, mv, inverse, pivots, save_vector

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_free_pair_quadrilateral_20261016'
PRIOR = 'results/fresh_sqrt5_universal_self_edges_20261016.json'


def component(i):
    return 0 if i < 5 else 1 if i < 13 else 2


def run():
    prior = json.loads((ROOT / PRIOR).read_text())
    assert prior['status'] == 'PASS_EXACT_UNIVERSAL_SAME_TYPE_EDGES_AND_BALANCED_DUAL_PORT_FLOWS'
    assert prior['exceptional_vertex_type'] == ['D5', 'E8', 'E8']
    A = [[Q(*x) for x in row] for row in prior['exceptional_full_finite_gram']]
    assert len(A) == 21 and all(A[i][i] == Q(2) for i in range(21))
    assert all(A[i][j] in (Q(), Q(-1)) for i in range(21) for j in range(i))
    assert all(x.sign() > 0 for x in pivots(A))
    C = inverse(A)
    lam = Q(*prior['uniform_rule']['lambda_value'])
    assert lam == Q(13, 21)
    T = 2 * lam
    ps = [2 * T - T * T * C[i][i] for i in range(21)]
    assert all(p.sign() < 0 and p.conj().sign() > 0 for p in ps)
    assert all(C[i][j].sign() >= 0 for i in range(21) for j in range(21))
    spatial = [[Q(k == i) - T * C[k][i] for k in range(21)] for i in range(21)]
    contacts = [mv(A, v) for v in spatial]
    assert len(prior['exceptional_all21_edge_and_escape_witnesses']) == 21
    for i, rec in enumerate(prior['exceptional_all21_edge_and_escape_witnesses']):
        assert rec['root'] == i and Q(*rec['norm']) == Q(2)
        assert Q(*rec['tau']) == T and Q(*rec['dotted_product']) == 2 - T
        assert contacts[i] == [A[k][i] - (T if k == i else Q()) for k in range(21)]
        assert 2 - dot(spatial[i], contacts[i]) == ps[i]
    rows = []
    minimum = None
    intercomponent = 0
    for i, j in combinations(range(21), 2):
        p, q = ps[i], ps[j]
        b = A[i][j] + T * T * C[i][j]
        assert dot(spatial[i], contacts[j]) == b
        assert b.sign() >= 0
        product = p * q
        assert product.sign() > 0 and product.conj().sign() > 0
        # After eliminating the 21 old finite roots, the two-root residual
        # is [[p,x-b],[x-b,q]], for free new-to-new product x.
        coefficients = [product - b * b, 2 * b, Q(-1)]
        def determinant(x):
            return coefficients[0] + coefficients[1] * x + coefficients[2] * x * x
        at_minus_two, at_zero = determinant(Q(-2)), determinant(Q())
        assert at_minus_two.sign() > 0 and at_zero.sign() > 0
        # On x=-2+2t, 0<=t<=1, the exact determinant equals its chord
        # plus4t(1-t). Coefficient comparison proves the continuous bound.
        direct = [coefficients[0] - 2 * coefficients[1] + 4 * coefficients[2],
                  2 * coefficients[1] - 8 * coefficients[2], 4 * coefficients[2]]
        chord = [at_minus_two, at_zero - at_minus_two + 4, Q(-4)]
        assert direct == chord
        for x, value in [(-2, at_minus_two), (0, at_zero)]:
            if minimum is None or (value - minimum[0]).sign() < 0:
                minimum = (value, i, j, x)
        # The two real rank22 roots are b+-sqrt(p*q). These inequalities
        # put the negative root strictly below-2 and the positive root above0.
        assert (b + 2).sign() > 0 and (product - (b + 2) * (b + 2)).sign() > 0
        assert (product - b * b).sign() > 0
        orthogonal = component(i) != component(j)
        assert orthogonal == (not C[i][j])
        if orthogonal:
            assert b == Q() and A[i][j] == Q()
            intercomponent += 1
        # Real partial-star check. Let t be a time vector of square-1 and
        # omega_i the old Gram-dual vector. Set a_i=sqrt(-p_i)>0 and
        # new_i=e_i-T*omega_i+a_i*t. The negative rank root is their mutual
        # product. The old vertex is t, and its two finite neighbors are
        # t-a_i/(T*Cii-1)*omega_i. Verify their strict other-new-wall sectors
        # without evaluating any radical, by positive square comparison.
        sectors = []
        for a, c in ((i, j), (j, i)):
            denom = T * C[a][a] - 1
            assert denom.sign() > 0
            aa, ac = -ps[a], -ps[c]
            margin = ac * denom * denom - aa * T * T * C[a][c] * C[a][c]
            assert margin.sign() > 0
            endpoint_square = aa * C[a][a] / (denom * denom) - 1
            assert endpoint_square == -1 / (denom * denom) and endpoint_square.sign() < 0
            sectors.append(dict(released_root=a, other_new_root=c,
                denominator=denom.saved(), endpoint_square=endpoint_square.saved(),
                strict_other_wall_squared_margin=margin.saved()))
        rows.append(dict(
            ports=[i, j], component_ids=[component(i), component(j)],
            orthogonal_components=orthogonal,
            old_product=A[i][j].saved(), inverse_diagonals=save_vector([C[i][i], C[j][j]]),
            inverse_cross=C[i][j].saved(),
            new_root_residual_diagonals=save_vector([p, q]),
            conjugate_new_root_residual_diagonals=save_vector([p.conj(), q.conj()]),
            residual_offdiagonal='x-b', b=b.saved(),
            determinant_coefficients_constant_linear_quadratic=save_vector(coefficients),
            determinant_at_minus2=at_minus_two.saved(), determinant_at_zero=at_zero.saved(),
            rank_on_entire_finite_corner_interval=23,
            signature_on_entire_finite_corner_interval=dict(positive=21, negative=2, zero=0),
            rank_closure_roots=dict(formula='x=b+-sqrt(radicand)', radicand=product.saved(),
                negative_root_strictly_below_minus2=True, positive_root_strictly_above_zero=True),
            real_partial_star_endpoint_sectors=sectors,
        ))
    assert len(rows) == 210 and intercomponent == 144
    # Input source chains are checked without mutating or re-running producers.
    queue = [PRIOR]
    seen = set()
    references = 0
    while queue:
        name = queue.pop()
        if name in seen:
            continue
        seen.add(name)
        for path, h in json.loads((ROOT / name).read_text()).get('sources', {}).items():
            assert sha256((ROOT / path).read_bytes()).hexdigest() == h, (name, path)
            references += 1
            if path.endswith('.json'):
                queue.append(path)
    paths = [f'scripts/{TAG}.py', PRIOR, 'scripts/fresh_sqrt5_d4_two_e8_self_edge_20261016.py']
    return dict(
        status='PASS_EXACT_NO_QUADRILATERAL_FOR_210_FIXED_UNIFORM_EDGE_PAIRS',
        sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
        recursive_manifests_checked=len(seen), recursive_source_references_checked=references,
        vertex_type=['D5', 'E8', 'E8'], ambient_hyperbolic_dimension=21,
        fixed_tau=T.saved(), root_norm='2', cases=210,
        orthogonal_component_pairs=144, same_component_pairs=66,
        interval_endpoint_checks=420, strict_real_endpoint_sector_checks=420,
        smallest_endpoint_determinant=dict(value=minimum[0].saved(), ports=list(minimum[1:3]), x=minimum[3]),
        reduction=dict(
            fixed_data='The old21 Gram A and all contacts from each chosen duplicate to the old roots are retained. Only the mutual product x of the two new roots is free.',
            spherical_corners='The old-old and the two old-new intended corners are copies of A. A fourth spherical corner requires its two norm-two new roots to have -2<x<=0. The checker relaxes this to the closed real interval[-2,0].',
            rank='There are23 normals in a22-dimensional ambient vector space, so rank<=22 is necessary. Since A has rank21, the residual2x2 block must have rank<=1.',
            schur='p_i=2T-T^2*Cii, p_j=2T-T^2*Cjj, b=Aij+T^2*Cij, C=A^-1; residual=[[p_i,x-b],[x-b,p_j]].',
            continuous_certificate='For x=-2+2t and0<=t<=1, D(x)=(1-t)D(-2)+tD(0)+4t(1-t)>0. Both residual diagonals are negative. Thus the residual is negative definite throughout the entire interval.',
            conjugate='Each fixed edge retains positive conjugate residual diagonal. Physical rank23 and negative index2 already exclude every finite-corner value, even before imposing conjugate PSD, algebraicity, integrality or Coxeter angle restrictions.',
            longer_faces='The negative real rank-closure root lies below-2, so the two new walls are ultraparallel. It supplies a real three-vertex partial star with the saved strict endpoint sectors, but cannot be the fourth corner of a quadrilateral. A completion preserving these two edges must use further side walls; longer-face or finite-volume completion is not excluded.',
        ),
        pairs=rows,
        scope='Complete continuous exclusion for all210 pairs of the specific frozen uniform long-edge witnesses. No exclusion of other edge lengths, attachments, vertex types, arithmetic polytopes or dimension21 follows. The real partial-star branch uses declared square roots; no arithmeticity, common integral lattice, complete edge star or finite-volume realization is asserted for it.',
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    out = json.loads(json.dumps(run()))
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out, indent=2) + '\n')
    print(out['status'])
