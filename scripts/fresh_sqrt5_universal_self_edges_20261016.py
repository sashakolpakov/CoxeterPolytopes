#!/usr/bin/env python3
"""Exact obstruction to vertex/edge-port compatibility refinements in rank 21.

All 724 rooted positions in the 68 frozen finite component types admit a
same-type V-arithmetic edge, with one uniform dotted-product rule. This does
not realize simultaneous edge stars or a finite-volume reflection chamber.
The latest restricted dual is extended by explicit balanced self-port flows.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

from fresh_sqrt5_d4_two_e8_self_edge_20261016 import (
    Q, phi, dot, mv, inverse, pivots, save_matrix, save_vector,
)
import fresh_sqrt5_forest_sensitive_weights_20261015_review as C

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_universal_self_edges_20261016'
DUAL = 'results/fresh_sqrt5_dimension21_a4_source_orbits_20261016_cap2.json'


def power(x, n):
    z = Q(1)
    for _ in range(n):
        z *= x
    return z


LAMBDA = power(phi, 8)


def integral(x):
    return x.a.denominator == x.b.denominator == 1


def absq(x):
    return x if x.sign() >= 0 else -x


def finite_gram(g):
    """Integral symmetrization; also covers 4, 6 and 10 bonds over Q(phi)."""
    n = len(g)
    weights = {0: Q(1)}
    queue = [0]
    oriented = []
    scale = {3: Q(1), 4: Q(2), 5: Q(1), 6: Q(3), 10: phi + 2}
    contact = {3: Q(1), 4: Q(2), 5: phi, 6: Q(3), 10: phi + 2}
    cyclic = {3: Q(1), 4: Q(2), 5: phi + 1, 6: Q(3), 10: phi + 2}
    for i in queue:
        for j, m in sorted(g[i].items()):
            if j in weights:
                continue
            weights[j] = weights[i] * scale[m]
            queue.append(j)
            oriented.append((i, j, m))
    assert len(weights) == n
    A = [[2 * weights[i] if i == j else Q() for j in range(n)] for i in range(n)]
    for i, j, m in oriented:
        A[i][j] = A[j][i] = -weights[i] * contact[m]
        assert 4 * A[i][j] * A[i][j] / (A[i][i] * A[j][j]) == cyclic[m]
    assert all(integral(x) for row in A for x in row)
    assert all(integral(2 * A[i][j] / A[i][i]) for i in range(n) for j in range(n))
    assert all(x.sign() > 0 and x.conj().sign() > 0 for x in pivots(A))
    return A


def edge(A, J, i):
    """Duplicate root i, keeping all its contacts except its new dotted pair."""
    n = len(A)
    d, c = A[i][i], J[i][i]
    s, tau = 1 / c, d * LAMBDA
    g = d - tau
    assert (g + d).sign() < 0  # normalized pairing g/d < -1
    assert (2 * s - tau).sign() < 0
    assert tau.sign() > 0 and tau.conj().sign() > 0
    assert (2 * s.conj() - tau.conj()).sign() > 0
    B = [row + [A[j][i] if j != i else g] for j, row in enumerate(A)]
    B.append([A[i][j] if j != i else g for j in range(n)] + [d])
    # Exact block inverse from elimination of the entire original finite Gram.
    q = [Q(j == i) - tau * J[j][i] for j in range(n)]
    delta = tau * (2 - tau * c)
    assert delta.sign() < 0
    K = [[J[j][k] + q[j] * q[k] / delta for k in range(n)] + [-q[j] / delta] for j in range(n)]
    K.append([-x / delta for x in q] + [1 / delta])
    rays = [[-K[j][omitted] for j in range(n + 1)] for omitted in (i, n)]
    for omitted, v in zip((i, n), rays):
        assert mv(B, v) == [-Q(j == omitted) for j in range(n + 1)]
        assert dot(v, mv(B, v)).sign() < 0
    square = dot(rays[0], mv(B, rays[0]))
    mutual = dot(rays[0], mv(B, rays[1]))
    assert square == dot(rays[1], mv(B, rays[1])) == 1 / delta
    assert mutual.sign() < 0 and mutual / square == tau / s - 1
    midpoint = [a + b for a, b in zip(*rays)]
    assert dot(midpoint, mv(B, midpoint)).sign() < 0
    assert all(integral(2 * B[j][k] / B[j][j]) for j in range(n + 1) for k in range(n + 1))
    # Length-two cyclic product of the dotted pair alone generates Q(phi).
    cycle = 4 * g * g / (d * d)
    assert integral(cycle) and cycle.b
    # Any other common-wall deletion leaves a Lorentz principal matrix:
    # its finite Schur diagonal is <= d, while tau > 2d. Its cofactor is
    # negative, as is det(B), so the common-wall inverse diagonal is positive.
    assert (LAMBDA - 2).sign() > 0
    assert all(K[j][j].sign() > 0 for j in range(n) if j != i)
    return B, K, dict(
        root=i, norm=d.saved(), finite_inverse_diagonal=c.saved(),
        finite_residual_diagonal=s.saved(), dotted_product=g.saved(),
        tau=tau.saved(), physical_black_eigenvalues=save_vector([tau, 2 * s - tau]),
        conjugate_black_eigenvalues=save_vector([tau.conj(), (2 * s - tau).conj()]),
        endpoint_rays=[save_vector(v) for v in rays],
        endpoint_square=square.saved(), endpoint_mutual_pairing=mutual.saved(),
        edge_length_cosh=(mutual / square).saved(),
        length_two_cyclic_product=cycle.saved(),
        integral_reflection_coefficients=True,
        common_wall_inverse_diagonals=[dict(root=j, value=K[j][j].saved()) for j in range(n) if j != i],
    )


def block_diagonal(blocks):
    n = sum(len(A) for A in blocks)
    B = [[Q() for _ in range(n)] for _ in range(n)]
    offset = 0
    for A in blocks:
        for i, row in enumerate(A):
            for j, x in enumerate(row):
                B[offset + i][offset + j] = x
        offset += len(A)
    return B


def strict_escape(B, K, i):
    """Exact two K-vectors defining a strict future null ray by one radical."""
    n = len(B)
    midpoint = [-K[j][i] - K[j][-1] for j in range(n)]
    v = [-sum(row, Q()) for row in K]
    AA, BB, CC = dot(midpoint, mv(B, midpoint)), dot(midpoint, mv(B, v)), dot(v, mv(B, v))
    assert AA.sign() < 0
    eps = -AA / (4 * (1 + absq(AA) + absq(BB) + absq(CC)))
    u = [a + eps * b for a, b in zip(midpoint, v)]
    Bu = mv(B, u)
    U = dot(u, Bu)
    assert eps.sign() > 0 and U.sign() < 0 and all(x.sign() < 0 for x in Bu)
    j = next(j for j in range(n - 1) if j != i)
    p = [-K[k][j] for k in range(n)]
    P, PU = dot(p, mv(B, p)), dot(p, Bu)
    assert P.sign() > 0
    eta = P / (4 * (1 + P + absq(PU) + absq(U)))
    w = [a + eta * b for a, b in zip(p, u)]
    Bw = mv(B, w)
    W, V = dot(w, Bw), dot(u, Bw)
    assert eta.sign() > 0 and W.sign() > 0 and all(x.sign() < 0 for x in Bw)
    D = V * V - U * W
    assert D.sign() > 0 and (D - V * V).sign() > 0
    # t=(sqrt(D)-V)/W>0, z=u+t*w. Work in the formal extension R^2=D.
    z0 = [a - V * b / W for a, b in zip(u, w)]
    z1 = [b / W for b in w]
    assert dot(z0, mv(B, z0)) + D * dot(z1, mv(B, z1)) == Q()
    assert dot(z0, mv(B, z1)) == Q()
    return dict(
        common_wall=j, epsilon=eps.saved(), eta=eta.saved(),
        u=save_vector(u), w=save_vector(w), u_square=U.saved(), w_square=W.saved(),
        u_w_pairing=V.saved(), u_wall_pairings=save_vector(Bu), w_wall_pairings=save_vector(Bw),
        radical_squared=D.saved(), null_ray_constant=save_vector(z0), null_ray_radical=save_vector(z1),
        null_ray='z=u+t*w; t=(sqrt(D)-V)/W>0, V=(u,w), W=w^2, D=V^2-u^2*w^2',
        future_proof='For V<=0, (z,u)=U+tV<0. For V>0, t=-U/(sqrt(D)+V), so (z,u)=U*sqrt(D)/(sqrt(D)+V)<0.',
        strict_boundary_patch=True,
        field_scope='Real radical extension. Conjugate positive definiteness forbids a nonzero K-rational null vector.',
    )


def balanced_dual_flows(catalog):
    d = json.loads((ROOT / DUAL).read_text())
    assert d['dimension'] == 21 and d['target'] == '40'
    assert d['certificate']['kind'] == 'exact_model_obstruction'
    totals = Counter()
    lower = F(0)
    flows = []
    z = d['variables']
    for index, t in enumerate(d['certificate']['terms']):
        y = F(t['multiplier'])
        assert y <= 0
        lower += y * F(t['rhs'])
        for j, x in t['coefficients']:
            totals[j] += y * F(x)
        if t['metadata']['kind'] != 'forest':
            continue
        types = t['metadata']['types']
        ports = [[a, name, i] for a, name in enumerate(types) for i in range(len(catalog[name]))]
        assert len(ports) == 21
        vertex_mass, edge_mass = -y, -y / 2
        assert vertex_mass >= 0 and 2 * edge_mass == vertex_mass
        flows.append(dict(dual_term=index, vertex_type=types, vertex_mass=str(vertex_mass),
            self_edge_mass_per_port=str(edge_mass), ports=ports))
    assert all(F(j == z) - x >= 0 for j, x in totals.items())
    assert lower == F(d['certificate']['lower_bound']) > 40
    assert len(flows) == 506
    return dict(
        exact_original_dual_lower_bound=str(lower), literal_dual_columns=len(totals),
        original_forest_mass=str(sum((F(r['vertex_mass']) for r in flows), F(0))),
        forest_types=len(flows), balanced_port_equations=21 * len(flows), flows=flows,
        scope='The existing literal weak-duality certificate is unchanged. Every forest multiplier is equipped with nonnegative self-edge flow satisfying 2*edge_mass=vertex_mass at each labelled root port. Thus endpoint-type/root-port handshake constraints alone do not restrict these forest frequencies. This does not add consistency between face rows and edge cycles, common ambient forms, facets, or simultaneous edge stars.',
    )


def simple_port_graph():
    # Round-robin 1-factorization of K22. Each color is a root port, and
    # exactly one edge of each color meets each vertex. No loops/multiple edges.
    edges = []
    for c in range(21):
        edges.append([21, c, c])
        for k in range(1, 11):
            a, b = sorted(((c + k) % 21, (c - k) % 21))
            edges.append([a, b, c])
    assert len({tuple(sorted((a, b))) for a, b, c in edges}) == len(edges) == 231
    for v in range(22):
        assert sorted(c for a, b, c in edges if v in (a, b)) == list(range(21))
    return dict(vertices=22, edges=edges, root_port_colors=21,
        scope='A connected simple abstract graph witnessing port incidence only. All vertices may carry any fixed rank21 forest, and edge color i uses its same-port self-edge. It is not a hyperbolic polytope, a simultaneous Gram representation, or a valid two-face assignment.')


def source_manifests(paths):
    queue = [p for p in paths if p.endswith('.json')]
    seen = set()
    references = 0
    while queue:
        p = queue.pop()
        if p in seen:
            continue
        seen.add(p)
        for q, h in json.loads((ROOT / p).read_text()).get('sources', {}).items():
            assert sha256((ROOT / q).read_bytes()).hexdigest() == h, (p, q)
            references += 1
            if q.endswith('.json'):
                queue.append(q)
    return dict(manifests_checked=len(seen), source_references_checked=references)


def run():
    assert LAMBDA == Q(13, 21) and LAMBDA.conj() * LAMBDA == Q(1)
    assert (LAMBDA - 4).sign() > 0
    models = [(name, g) for name, g in C.models() if len(g) <= 21]
    assert len(models) == 68
    catalog = {name: finite_gram(g) for name, g in models}
    records = []
    largest_requirement = Q()
    maximizers = []
    all_roots = 0
    for name, g in models:
        A = catalog[name]
        J = inverse(A)
        roots = []
        for i in range(len(A)):
            B, K, rec = edge(A, J, i)
            bound = A[i][i].conj() * J[i][i].conj() / 2
            cmp = (bound - largest_requirement).sign()
            if cmp > 0:
                largest_requirement, maximizers = bound, [(name, i)]
            elif cmp == 0:
                maximizers.append((name, i))
            roots.append(rec)
            all_roots += 1
        records.append(dict(type=name, rank=len(g),
            edges=[[i, j, m] for i in g for j, m in g[i].items() if i < j],
            integral_gram=save_matrix(A), physical_positive_pivots=save_vector(pivots(A)),
            conjugate_positive_pivots=save_vector([p.conj() for p in pivots(A)]), root_witnesses=roots))
        print(name, 'roots', len(roots), flush=True)
    assert all_roots == 724 and largest_requirement == Q(30)
    assert (power(phi, 6) - largest_requirement).sign() < 0
    # Complete rank21 finite forests are multisets of the 68 distinct types.
    dp = [0] * 22
    dp[0] = 1
    for A in catalog.values():
        for rank in range(len(A), 22):
            dp[rank] += dp[rank - len(A)]
    assert dp[21] == 311574
    # Full, explicit endpoint and escape witnesses for all 21 ports of the
    # recurring exceptional forest, always in separate local edge chambers.
    A = block_diagonal([catalog['D5'], catalog['E8'], catalog['E8']])
    J = inverse(A)
    exceptional = []
    for i in range(21):
        B, K, rec = edge(A, J, i)
        assert len(B) == 22
        rec['strict_null_escape'] = strict_escape(B, K, i)
        exceptional.append(rec)
    # Diagnose the missing simultaneous constraint: two canonical copied
    # duplicates with the uniform long edge rule give two negative directions.
    # This is a check on that proposed gluing, not an exclusion of all squares.
    square_tests = []
    for i, j in combinations(range(21), 2):
        h = J[i][i] * J[j][j] - J[i][j] * J[j][i]
        assert h.sign() > 0
        alpha, beta, gamma = J[j][j] / h, -J[i][j] / h, J[i][i] / h
        a, b = 2 * alpha - A[i][i] * LAMBDA, 2 * gamma - A[j][j] * LAMBDA
        det = a * b - 4 * beta * beta
        assert a.sign() < 0 and b.sign() < 0 and det.sign() > 0
        square_tests.append(dict(roots=[i, j], symmetric_block=save_matrix([[a, 2 * beta], [2 * beta, b]]), determinant=det.saved()))
    paths = [f'scripts/{TAG}.py',
        'scripts/fresh_sqrt5_d4_two_e8_self_edge_20261016.py',
        'scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py',
        'results/fresh_sqrt5_d4_two_e8_self_edge_20261016.json',
        'results/fresh_sqrt5_forest_sensitive_weights_20261015_review.json', DUAL]
    return dict(
        status='PASS_EXACT_UNIVERSAL_SAME_TYPE_EDGES_AND_BALANCED_DUAL_PORT_FLOWS',
        sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
        recursive_source_manifests=source_manifests(paths),
        dimension=21, finite_component_types=68, rooted_component_positions=all_roots,
        all_rank21_forest_types=dp[21], forest_count_dp=dp,
        uniform_rule=dict(k=4, lambda_value=LAMBDA.saved(), lambda_conjugate=LAMBDA.conj().saved(),
            dotted_product='g=d_i*(1-phi^8)', largest_required_conjugate_bound=largest_requirement.saved(),
            maximizers=maximizers, k3_fails_for_uniform_scaled_rule=True),
        component_records=records,
        generic_argument=dict(
            assumptions='A is the positive definite integral symmetrized Gram of any finite forest; its conjugate is positive definite. Root i has norm d and residual s=1/(A^-1)_ii. Duplicate all its contacts except its mutual product g=d-d*lambda.',
            schur='The two-root Schur block is [[s,s-tau],[s-tau,s]], tau=d*lambda. Eigenvalues are tau and 2s-tau.',
            physical='lambda>2 and s<=d imply tau>2s, giving exactly one negative direction.',
            conjugate='lambda>max(d_conj*(A^-1)_ii_conj/2)=30 implies 0<tau_conj<2s_conj, giving positive definiteness.',
            endpoint='Negative Gram-dual columns at the two framing roots have wall pairings -e_i and -e_new, equal negative squares and negative mutual pairing. The interior of their positive segment has exactly the common walls active.',
            integrality='All finite reflection coefficients are integral; the new framing coefficients 2g/d=2(1-lambda) are integral. The dotted length-two cyclic product 4(1-lambda)^2 is integral and generates Q(phi).',
            spectators='Orthogonally append the other finite components. This adds only positive directions at both embeddings and leaves the endpoint and cyclic-product checks unchanged.',
            escape='For rank at least2, delete any common root j. The remaining finite-old-root Gram is positive definite and its added duplicate has residual 2tau-tau^2*c_tilde<0 because c_tilde>=1/d and tau>2d. Thus the common-root inverse diagonal is positive. A small positive perturbation of the endpoint midpoint gives a strict timelike u; a small positive u-perturbation of -G^-1 e_j gives a strict spacelike w. Their positive combination u+t*w reaches a strict future null ray. Explicit arithmetic for all21 ports of D5+2E8 is saved.',
        ),
        exceptional_vertex_type=['D5', 'E8', 'E8'], exceptional_full_finite_gram=save_matrix(A),
        exceptional_all21_edge_and_escape_witnesses=exceptional,
        latest_dual_port_flow=balanced_dual_flows(catalog),
        simple_abstract_port_incidence=simple_port_graph(),
        canonical_copied_square_diagnostic=dict(cases=len(square_tests), negative_index=2,
            cases_data=square_tests,
            scope='Two uniform duplicated roots with all cross contacts copied from the original finite diagram have a negative definite symmetric 2x2 Schur block. Hence these particular local witnesses cannot be canonically glued into a square in a Lorentz space. Other edge parameters, root identifications and two-face types remain open.'),
        scope='Universal V-arithmetic local self-edge witnesses for the frozen finite-type family, a literal old-dual check, and balanced endpoint/root-port incidence flows. No simultaneous edge-star realization, common ambient lattice, valid two-face gluing, compact or finite-volume polytope, arithmetic reflection lattice, or improved dimension bound is claimed. All rank21 displayed local edge chambers have strict ideal escape by the general argument; explicit escapes are saved for D5+2E8.',
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
