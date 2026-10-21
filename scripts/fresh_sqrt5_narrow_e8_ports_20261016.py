#!/usr/bin/env python3
"""Exact minimality14 restrictions at E8 ports and a narrow-facet incidence bound.

Copied attachments are distinguished from all same-type reattachments. The
global statement concerns vertices on a selected narrow facet, not arbitrary
whole-polytope vertex frequencies. No dimension bound is asserted.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from fresh_sqrt5_d4_two_e8_self_edge_20261016 import Q, phi, dot, mv, inverse, pivots, save_matrix, save_vector
import fresh_sqrt5_forest_sensitive_weights_20261015_review as C

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_narrow_e8_ports_20261016'
PRIOR = 'results/fresh_sqrt5_universal_self_edges_20261016.json'


def graph(A, omit=()):
    return {i: {j: 3 for j in range(len(A)) if j != i and j not in omit and A[i][j] == Q(-1)}
            for i in range(len(A)) if i not in omit}


def pattern(g, root):
    rest = {i: {j: m for j, m in row.items() if j != root} for i, row in g.items() if i != root}
    return tuple(sorted(C.family(c) for c in C.pieces(rest)))


def scalar_box():
    # If-14<g<=0 and0<g'<2, write g=a+b*phi. Then-8<b<0 and
    # -5<a<2, so the complete integer box is a=-4..1,b=-7..-1.
    assert (Q(5) - 7 / phi).sign() > 0
    assert (Q(5) - 8 / phi).sign() > 0  # used also in the sharp larger box
    assert (Q(5) - Q(F(16 * 16, 8 * 8))).sign() > 0  # sqrt5>16/8
    return [Q(a, b) for a in range(-4, 2) for b in range(-7, 0)]


def admissible_scalars(q, s):
    out = []
    assert (q - s).sign() >= 0 and (2 - q - s).sign() >= 0
    for g in scalar_box():
        if not ((g + 14).sign() > 0 and g.sign() <= 0):
            continue
        if not ((g.conj() - q + s).sign() > 0 and (q + s - g.conj()).sign() > 0):
            continue
        dotted = (g + 2).sign() < 0
        if not dotted:
            assert g.sign() < 0 and (g + 1).sign() > 0
        out.append(dict(g=g.saved(), conjugate=g.conj().saved(),
            coxeter_allowed=dotted,
            reason='DOTTED_PRODUCT' if dotted else 'REJECTED_-1<g<0_IS_NOT_A_COXETER_PRODUCT'))
    return out


def actual_edge(A, port, attachments, g):
    n = len(A)
    new = [g if j == port else Q(-int(j in attachments)) for j in range(n)]
    B = [row + [new[j]] for j, row in enumerate(A)] + [new + [Q(2)]]
    ps = pivots(B)
    assert sum(p.sign() > 0 for p in ps) == n and sum(p.sign() < 0 for p in ps) == 1
    assert all(p.conj().sign() > 0 for p in ps)
    assert all(x.a.denominator == x.b.denominator == 1 for row in B for x in row)
    K = inverse(B)
    rays = []
    for omitted in (port, n):
        v = [-K[j][omitted] for j in range(n + 1)]
        assert mv(B, v) == [-Q(j == omitted) for j in range(n + 1)]
        assert dot(v, mv(B, v)).sign() < 0
        active = [j for j in range(n + 1) if j != omitted]
        assert all(p.sign() > 0 for p in pivots([[B[j][k] for k in active] for j in active]))
        rays.append(v)
    assert dot(rays[0], mv(B, rays[1])).sign() < 0
    return dict(g=g.saved(), physical_signature=[n, 1], conjugate_signature=[n + 1, 0],
        endpoint_rays=[save_vector(v) for v in rays],
        endpoint_squares=[dot(v, mv(B, v)).saved() for v in rays],
        endpoint_mutual_pairing=dot(rays[0], mv(B, rays[1])).saved(),
        scope='Actual finite local edge with the specified endpoint diagrams. No finite-volume or completion claim.')


def copied_port_test(E):
    J = inverse(E)
    beta = 3 + 8 * phi
    assert (beta - 14).sign() > 0 and (beta + 2 - 18).sign() < 0
    rows = []
    for i in range(8):
        c = J[i][i]
        assert not c.b and c.sign() > 0
        s = 1 / c
        q = 2 - s
        candidates = admissible_scalars(q, s)
        valid = [v for v in candidates if v['coxeter_allowed']]
        row = dict(local_E8_root=i, inverse_diagonal=c.saved(), residual=s.saved(),
            conjugate_interval=save_vector([q - s, q + s]),
            complete_narrow_scalar_candidates=candidates, narrow_coxeter_products=len(valid))
        if c in (Q(14), Q(20), Q(30)):
            assert not candidates
            # A shorter absolute product has- beta<g<=0. The same embedding
            # bounds give a=-4..1,b=-8..-1. This is a proved finite box.
            shorter = []
            for a in range(-4, 2):
                for b in range(-8, 0):
                    g = Q(a, b)
                    if (g + beta).sign() > 0 and g.sign() <= 0 and (g.conj() - q + s).sign() > 0 and (q + s - g.conj()).sign() > 0:
                        shorter.append(g)
            assert not shorter
            g = -beta
            assert (g.conj() - q + s).sign() > 0 and (q + s - g.conj()).sign() > 0
            attachments = [j for j in range(8) if j != i and E[i][j] == Q(-1)]
            row['sharp_minimum_absolute_product'] = beta.saved()
            row['sharp_local_edge'] = actual_edge(E, i, attachments, g)
        rows.append(row)
    assert [r['local_E8_root'] for r in rows if not r['narrow_coxeter_products']] == [1, 2, 3]
    assert sum(r['narrow_coxeter_products'] for r in rows) == 11
    return rows


def all_same_type_reattachments(A):
    target_codes = sorted(C.code(c) for c in C.pieces(graph(A)))
    results = []
    for port in (6, 7, 8):  # local E8 positions1,2,3, with inverse diagonal14,30,20
        old = graph(A, (port,))
        nodes = sorted(old)
        S = [[A[i][j] for j in nodes] for i in nodes]
        J = inverse(S)
        u = [A[port][i] for i in nodes]
        qold = dot(u, mv(J, u))
        s = 2 - qold
        tested = 0
        accepted = []
        for degree in range(4):
            for attachments in combinations(nodes, degree):
                tested += 1
                g = {i: dict(row) for i, row in old.items()}
                g[21] = {}
                for i in attachments:
                    g[i][21] = g[21][i] = 3
                parts = C.pieces(g)
                if any(sum(map(len, c.values())) != 2 * (len(c) - 1) for c in parts):
                    continue
                if sorted(C.code(c) for c in parts) != target_codes:
                    continue
                v = [-Q(int(i in attachments)) for i in nodes]
                assert dot(v, mv(J, v)) == qold
                cross = dot(u, mv(J, v))
                candidates = admissible_scalars(cross, s)
                witnesses = [actual_edge(A, port, attachments, Q(*r['g'])) for r in candidates if r['coxeter_allowed']]
                common_neighbors = [i for i in attachments if A[port][i] == Q(-1)]
                accepted.append(dict(attachments=list(attachments),
                    copied_attachments=list(attachments) == [i for i in nodes if A[port][i] == Q(-1)],
                    cross_projection=cross.saved(), residual=s.saved(),
                    conjugate_interval=save_vector([cross - s, cross + s]),
                    narrow_scalar_candidates=candidates, narrow_actual_edge_witnesses=witnesses,
                    common_neighbors_for_integral_triangle=common_neighbors))
        assert tested == sum((1, 20, 190, 1140)) == 1351
        assert len(accepted) == (2 if port == 6 else 4)
        if port == 7:
            assert all(12 in r['common_neighbors_for_integral_triangle'] for r in accepted)
            assert all(not r['narrow_actual_edge_witnesses'] for r in accepted)
            assert sorted(Q(*r['cross_projection']).a for r in accepted) == [F(31, 30), F(41, 30), F(49, 30), F(59, 30)]
            candidates = [v for r in accepted for v in r['narrow_scalar_candidates']]
            assert len(candidates) == 1 and Q(*candidates[0]['g']) == 1 - phi and not candidates[0]['coxeter_allowed']
        else:
            assert sum(len(r['narrow_actual_edge_witnesses']) for r in accepted) == (1 if port == 6 else 2)
        results.append(dict(port=port, common_vertex_diagram=sorted(C.family(c) for c in C.pieces(old)),
            candidate_attachment_subsets=tested, complete_same_type_reattachments=accepted,
            scope='All endpoint attachments returning to D5+2E8: every new bond is3 and the new root degree is at most3, so all subsets of the20 common nodes of sizes0..3 are exhausted.'))
    return results


def endpoint_census():
    models = {C.family(g): g for name, g in C.models() if len(g) <= 21}
    assert len(models) == 68
    S = Counter(['A1', 'I2(3)', 'A4', 'D5', 'E8'])
    types = set()
    for name, g in models.items():
        for root in g:
            p = Counter(pattern(g, root))
            if not p - S:
                types.add(tuple(sorted(list((S - p).elements()) + [name])))
    rows = []
    for vertex in sorted(types):
        assert sum(len(models[name]) for name in vertex) == 21
        roots = []
        for ci, name in enumerate(vertex):
            for root in models[name]:
                remaining = Counter(list(vertex[:ci] + vertex[ci + 1:]) + list(pattern(models[name], root)))
                if remaining == S:
                    roots.append([ci, name, root])
        assert roots
        rows.append(dict(vertex_type=list(vertex), root_deletion_ports=roots, multiplicity=len(roots),
            excluded_by_central_narrow_return=vertex == ('D5', 'E8', 'E8')))
    assert len(rows) == 29
    others = [r for r in rows if not r['excluded_by_central_narrow_return']]
    assert len(others) == 28 and max(r['multiplicity'] for r in others) == 4
    assert sum(r['multiplicity'] == 4 for r in others) == 3
    return dict(common_edge_type=sorted(S.elements()), finite_component_catalogue_size=68,
        finite_completion_types=29, nonreturn_types=28, maximum_other_vertex_incidence=4,
        records=rows,
        completeness='The new root belongs to exactly one finite component. Deleting it gives a submultiset of the fixed common components; all other components stay unchanged. Exhausting all rooted finite components through rank21 therefore exhausts every completion type and every incoming port.')


def manifests(paths):
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
    prior = json.loads((ROOT / PRIOR).read_text())
    A = [[Q(*x) for x in row] for row in prior['exceptional_full_finite_gram']]
    assert len(A) == 21
    E = [[A[i][j] for j in range(5, 13)] for i in range(5, 13)]
    assert [inverse(E)[i][i] for i in range(8)] == list(map(Q, (4, 14, 30, 20, 12, 6, 2, 8)))
    copied = copied_port_test(E)
    reattachments = all_same_type_reattachments(A)
    census = endpoint_census()
    central = [7, 15]
    roles = [dict(facet_root=e, central_E8_root=e in central,
        central_outgoing_edges_inside_facet=2 - int(e in central)) for e in range(21)]
    assert sum(r['central_E8_root'] for r in roles) == 2
    assert min(r['central_outgoing_edges_inside_facet'] for r in roles) == 1
    paths = [f'scripts/{TAG}.py', PRIOR,
        'scripts/fresh_sqrt5_d4_two_e8_self_edge_20261016.py',
        'scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py']
    return dict(status='PASS_EXACT_NARROW_CENTRAL_E8_RETURN_EXCLUSION_AND_FACET_FREQUENCY_BOUND',
        sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
        recursive_sources=manifests(paths), ground_field='Q(sqrt5)',
        scalar_search=dict(ring_of_integers='Z[phi]', phi_relation='phi^2=phi+1',
            narrow_hypotheses='-14<g<=0 and0<g_conjugate<2',
            complete_coefficient_box=dict(a=list(range(-4, 2)), b=list(range(-7, 0))),
            box_proof='b=(g-g_conjugate)/sqrt5 lies strictly between-8 and0; a=g_conjugate+b/phi lies strictly between-5 and2.',
            integrality_proof='For central returns, the common A1 root has product-1 with both framing roots. Its triangular cyclic product is g. Thus arithmetic cyclic-product integrality implies g belongs to the full ring of integers Z[phi]; this is not inferred from norm-two normalization alone.'),
        copied_e8_ports=copied,
        same_type_reattachments=reattachments,
        central_edge_endpoint_census=census,
        selected_facet_incidence=dict(ambient_dimension=21, facet_dimension=20,
            exceptional_vertex_type=['D5', 'E8', 'E8'], central_root_indices=central,
            facet_root_roles=roles,
            precise_inequality='D_central+2*D_other <= sum_{U != D5+2E8} m(U)*V_U',
            coarse_inequality='D_central+2*D_other <=4*O',
            resulting_exceptional_fraction_upper_bound=str(F(4, 5)),
            proof='At a vertex on the selected narrow facet F, every central E8 port except possibly the facet normal gives an edge lying in F. Its framing pair has minimality14. The opposite endpoint cannot again have type D5+2E8. Count these edges at their nonexceptional ends; their common diagram is S, so each endpoint of type U receives at most m(U)<=4. Every edge counted has exactly one exceptional endpoint.',
            applicability='Only vertices and edges of the chosen narrow facet. No constraint on whole-polytope vertex frequencies is asserted. The theorem does not force a selected triangle/quadrilateral to meet these central ports.'),
        published_inputs=[dict(url='https://arxiv.org/pdf/0708.3991',
            location='Section2; Section4 opening, printed page10; Section4.1, printed pages11-12; Section4.2, printed page14',
            used_for='Arithmetic cyclic-product/conjugate criterion; existence of a narrow facet and minimality14 for its edge chambers; simple-facet average ensuring a triangle or quadrilateral in dimension at least6. The local interval and endpoint-type calculations are new exact deductions.')],
        scope='An exact copied-port diagnostic, complete same-type central-E8 reattachment exclusion under minimality14, and a conditional narrow-facet incidence theorem for compact arithmetic Q(sqrt5) Coxeter polytopes. No general exclusion of E8 vertices, all arithmetic edges, compact dimension21, or finite-volume examples is asserted.')


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
