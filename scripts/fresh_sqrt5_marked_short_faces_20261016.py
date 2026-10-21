#!/usr/bin/env python3
"""All-parameter short-face test for the two marked baseline38 types.

Inputs are the complete minimality14 marked-transition certificate, not one
chosen edge per port. The new-new product is free until rank and the complete
finite-label catalogue force it. The derived incidence inequality concerns
vertices of a selected narrow facet, not the whole polytope.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import isqrt
from pathlib import Path
import fresh_sqrt5_forest_sensitive_weights_20261015_review as C
from fresh_sqrt5_d4_two_e8_self_edge_20261016 import Q, phi, dot, mv, inverse, pivots, save_vector

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_marked_short_faces_20261016'
PRIOR = 'results/fresh_sqrt5_marked_maximum_transitions_20261016.json'


def integrity():
    queue, seen, count = [PRIOR], set(), 0
    while queue:
        path = queue.pop()
        if path in seen:
            continue
        seen.add(path)
        for source, digest in json.loads((ROOT / path).read_text()).get('sources', {}).items():
            assert sha256((ROOT / source).read_bytes()).hexdigest() == digest, (path, source)
            count += 1
            if source.endswith('.json'):
                queue.append(source)
    return dict(manifests_checked=len(seen), source_references_checked=count)


def rational_square(x):
    x = F(x)
    return x >= 0 and isqrt(x.numerator) ** 2 == x.numerator and isqrt(x.denominator) ** 2 == x.denominator


def field_gate():
    # A rational z is a square in Q(sqrt5) iff z or z/5 is a square in Q:
    # (r+s*sqrt5)^2=z forces r*s=0. For z=2+phi the norm is5, which
    # cannot be the rational square norm(x)^2 of a field square x^2.
    for z in (2, 3):
        assert not rational_square(z) and not rational_square(F(z, 5))
    z = 2 + phi
    assert z * z.conj() == Q(5) and not rational_square(5)
    return dict(complete_finite_labels=[2, 3, 4, 5, 6, 10],
        field_products=[dict(label=m, product=x.saved()) for m, x in ((2, Q()), (3, Q(-1)), (5, -phi))],
        excluded_square_roots=[dict(label=4, square=['2', '0']), dict(label=6, square=['3', '0']), dict(label=10, square=['2', '1'])],
        justification='The corank equation is a*d-(x-c)^2=0, with a,d,c inK and c>0. Every permitted finite label has x^2 inK, so the equation forces x=(x^2+c^2-a*d)/(2*c) inK. The products for labels4,6,10 are not inK by the saved rational-square/norm checks. Only0,-1,-phi remain. The finite-label catalogue is the same frozen Q(sqrt5) catalogue used for the complete vertex/face certificates.')


def graph(t):
    g = {i: {} for i in range(21)}
    for i, j, m in t['edges']:
        g[i][j] = g[j][i] = m
    return g


def gram(g):
    return [[Q(2 if i == j else -int(j in g[i])) for j in range(21)] for i in range(21)]


def classify(g, finite_codes):
    pieces = C.pieces(g)
    if any(sum(map(len, c.values())) != 2 * len(c) - 2 for c in pieces):
        return None
    codes = [C.code(c) for c in pieces]
    if any(c not in finite_codes for c in codes):
        return None
    return sorted(finite_codes[c] for c in codes)


def finite_corner_graph(g, i, j, r, s, m):
    out = {k: {l: value for l, value in row.items() if l not in (i, j)} for k, row in g.items() if k not in (i, j)}
    out[21], out[22] = {}, {}
    for new, choice in ((21, r), (22, s)):
        for k in choice['attachments']:
            if k in out and k < 21:
                out[k][new] = out[new][k] = 3
    if m != 2:
        out[21][22] = out[22][21] = m
    return out


def edge_basis(A, r):
    b = r['vector']
    B = [row + [b[i]] for i, row in enumerate(A)] + [b + [Q(2)]]
    ps = pivots(B)
    assert sum(z.sign() > 0 for z in ps) == 21 and sum(z.sign() < 0 for z in ps) == 1
    assert all(z.conj().sign() > 0 for z in ps)
    return B, inverse(B)


def triangle_witness(A, i, j, r):
    B, J = edge_basis(A, r)
    rays = [[-J[k][omitted] for k in range(22)] for omitted in (21, i, j)]
    for omitted, v in zip((21, i, j), rays):
        assert mv(B, v) == [-Q(k == omitted) for k in range(22)]
        assert dot(v, mv(B, v)).sign() < 0
    products = [dot(v, mv(B, w)) for v, w in combinations(rays, 2)]
    assert all(z.sign() < 0 for z in products)
    return dict(basis_root_ids=list(range(22)), omitted_roots=[21, i, j],
        future_vertex_rays=[save_vector(v) for v in rays],
        vertex_squares=[dot(v, mv(B, v)).saved() for v in rays],
        pairwise_vertex_products=save_vector(products),
        physical_signature=[21, 1], conjugate_signature=[22, 0])


def quadrilateral_witness(A, i, j, r, s, x):
    B, J = edge_basis(A, r)
    new_pairings = s['vector'] + [x]
    w = mv(J, new_pairings)
    assert dot(w, new_pairings) == Q(2)
    assert w[i].sign() > 0 and w[j].sign() < 0 and w[21].sign() > 0
    wall_data = []
    for nonzero in ({21: Q(1)}, {i: Q(1)}, {i: -w[j], j: w[i]}, {j: w[21], 21: -w[j]}):
        rhs = [nonzero.get(k, Q()) for k in range(22)]
        assert all(z.sign() >= 0 for z in rhs)
        v = [-z for z in mv(J, rhs)]
        pairings = mv(B, v) + [dot(v, new_pairings)]
        assert all(z.sign() <= 0 for z in pairings)
        assert sum(z.sign() < 0 for z in pairings) == 2
        assert dot(v, mv(B, v)).sign() < 0
        wall_data.append((v, pairings))
    rays = [v for v, p in wall_data]
    products = [dot(v, mv(B, w)) for v, w in combinations(rays, 2)]
    assert all(z.sign() < 0 for z in products)
    assert [[k for k, z in enumerate(p) if z.sign() < 0] for v, p in wall_data] == [[21, 22], [i, 22], [i, j], [j, 21]]
    return dict(basis_root_ids=list(range(22)), extra_root_id=22,
        extra_root_coefficients=save_vector(w),
        cyclic_vertex_rays=[save_vector(v) for v in rays],
        cyclic_wall_pairings=[save_vector(p) for v, p in wall_data],
        vertex_squares=[dot(v, mv(B, v)).saved() for v in rays],
        pairwise_vertex_products=save_vector(products),
        physical_signature=[21, 1, 1], conjugate_signature=[22, 0, 1],
        scope='A genuine local finite quadrilateral in the retained facet. Its two incoming edge parameters are the uniquely surviving choices from the complete narrow max-type transition catalogue. No finite-volume or full-facet completion is asserted.')


def run():
    data = json.loads((ROOT / PRIOR).read_text())
    assert data['status'] == 'PASS_EXACT_MARKED_NARROW_TRANSITIONS_RETAIN_BASELINE38_OBSTRUCTION'
    finite_codes = {C.code(g): name for name, g in C.models() if len(g) <= 21}
    assert len(finite_codes) == 68
    gate = field_gate()
    outputs = []
    for t in data['marked_types']:
        name, e, g = t['name'], t['facet_root'], graph(t)
        A, choices, records, triangles, totals = gram(g), {}, [], [], Counter()
        J = inverse(A)
        for row in data['port_records']:
            if row['source_marked_type'] != name:
                continue
            port = row['port']
            choices[port] = []
            for ai, attachment in enumerate(row['complete_marked_reattachments']):
                for gi, scalar in enumerate(attachment['all_narrow_admissible_products']):
                    b = [Q(*scalar['g']) if k == port else Q(-int(k in attachment['attachments'])) for k in range(21)]
                    jb = mv(J, b)
                    residual = 2 - dot(b, jb)
                    assert residual.sign() < 0 and residual.conj().sign() > 0
                    choices[port].append(dict(input_reattachment_index=ai, input_scalar_index=gi,
                        target_marked_type=attachment['target_marked_type'], attachments=attachment['attachments'],
                        g=scalar['g'], vector=b, inverse_vector=jb, residual=residual))
        for i, j in combinations(choices, 2):
            local = []
            for ri, r in enumerate(choices[i]):
                for si, s in enumerate(choices[j]):
                    a, d, c = r['residual'], s['residual'], dot(r['vector'], s['inverse_vector'])
                    assert c.sign() > 0
                    det = lambda x: a * d - (x - c) * (x - c)
                    dm, dz = det(Q(-2)), det(Q())
                    row = dict(choice_indices=[ri, si], schur_diagonals=save_vector([a, d]),
                        cross_projection=c.saved(), determinant_at_minus2=dm.saved(), determinant_at0=dz.saved())
                    totals['pairs_of_edge_choices'] += 1
                    if dm.sign() > 0 and dz.sign() > 0:
                        row['status'] = 'NO_REAL_FINITE_PAIRING_BY_STRICT_CONCAVITY'
                        totals['continuous_exclusions'] += 1
                    else:
                        tests = []
                        for m, x in ((2, Q()), (3, Q(-1)), (5, -phi)):
                            value = det(x)
                            test = dict(label=m, product=x.saved(), determinant=value.saved())
                            if not value:
                                totals['angular_rank_roots'] += 1
                                corner = finite_corner_graph(g, i, j, r, s, m)
                                types = classify(corner, finite_codes)
                                test['finite_fourth_vertex_type'] = types
                                if types is not None:
                                    totals['finite_quadrilaterals'] += 1
                                    test['actual_quadrilateral'] = quadrilateral_witness(A, i, j, r, s, x)
                            tests.append(test)
                        row['status'] = 'EXACT_COMPLETE_FINITE_LABEL_TEST'
                        row['finite_label_tests'] = tests
                    local.append(row)
            records.append(dict(ports=[i, j], edge_choice_pairs=local))
        # A triangle with two maximal vertices must use a nondotted edge.
        # All catalogue framings are dotted or-phi, so only the latter can
        # appear in its third finite vertex. Enumerate every possible second
        # released old root and test the complete finite component catalogue.
        for i, options in choices.items():
            for ri, r in enumerate(options):
                framing = Q(*r['g'])
                if framing != -phi:
                    assert (framing + 2).sign() < 0
                    continue
                for j in g:
                    if j in (e, i):
                        continue
                    h = {k: {l: value for l, value in rr.items() if l != j} for k, rr in g.items() if k != j}
                    h[21] = {}
                    for k in r['attachments']:
                        if k != j:
                            h[k][21] = h[21][k] = 3
                    h[i][21] = h[21][i] = 5
                    types = classify(h, finite_codes)
                    if types is None:
                        continue
                    cost = sum((C.base(c, u, v) for c in C.pieces(h) for u, v in combinations(c, 2) if e not in (u, v)), F())
                    assert cost == 35
                    triangles.append(dict(ports=sorted((i, j)), max_to_max_port=i,
                        choice_index=ri, other_maximal_type=r['target_marked_type'],
                        third_vertex_type=types, third_marked_baseline=str(cost),
                        actual_triangle=triangle_witness(A, i, j, r)))
        assert len(records) == 190
        outputs.append(dict(source_marked_type=name, facet_root=e,
            choices_by_port=[dict(port=i, choices=[{k: (v.saved() if isinstance(v, Q) else v) for k, v in c.items() if k not in ('vector', 'inverse_vector')} for c in options]) for i, options in choices.items()],
            port_pair_records=records, counts=dict(totals), triangles_with_two_maximal_vertices=triangles))
    first, second = outputs
    assert first['counts'] == dict(pairs_of_edge_choices=2542, continuous_exclusions=2527, angular_rank_roots=1, finite_quadrilaterals=1)
    assert second['counts'] == dict(pairs_of_edge_choices=486, continuous_exclusions=483)
    assert [t['ports'] for t in first['triangles_with_two_maximal_vertices']] == [[18, 19], [18, 19]]
    assert [t['ports'] for t in second['triangles_with_two_maximal_vertices']] == [[17, 18]]
    surviving_quads = [(row['ports'], case['choice_indices'], test['label'], test['finite_fourth_vertex_type'])
        for row in first['port_pair_records'] for case in row['edge_choice_pairs']
        for test in case.get('finite_label_tests', []) if 'actual_quadrilateral' in test]
    assert surviving_quads == [([18, 19], [4, 9], 3, ['A1', 'D20'])]
    assert F(155, 189) == 1 / (1 + F(34, 155))
    assert F(37) + F(155, 189) == F(7148, 189)
    paths = [f'scripts/{TAG}.py', PRIOR,
        'scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py',
        'scripts/fresh_sqrt5_d4_two_e8_self_edge_20261016.py']
    return dict(status='PASS_EXACT_ALL_PARAMETER_MARKED_SHORT_FACES_AND_FREQUENCY_BOUND',
        sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
        prior_source_integrity=integrity(), finite_component_count=68, finite_label_field_gate=gate,
        full_parameter_scope='Every admissible minimality14 marked max-to-max endpoint attachment and scalar from the frozen complete transition certificate is used, including cross transitions between A1+D20 andD21. The mutual product of the two new roots is free. The fourth corner may be any finite Q(sqrt5) type.',
        continuous_exclusion_identity='For x=-2+2t,0<=t<=1: D(x)=(1-t)D(-2)+tD(0)+4t(1-t), where D(x)=a*d-(x-c)^2. Positive endpoint values force negative index2 and rank23 for every finite fourth-corner product, contradicting an ambient22-dimensional Lorentz space.',
        all_maximal_triangle_exclusion='At a maximal source vertex, a triangle with two maximal neighbors would have one shared new wall. Its contact at either released root must simultaneously be a maximal-endpoint contact0/-1 and a framing product-phi or dotted. These sets are disjoint.',
        type_results=outputs, total_port_pairs=380, total_edge_choice_pairs=3028,
        total_continuous_exclusions=3010, remaining_exact_angular_tests=54,
        total_finite_quadrilaterals=1, total_directed_two_maximal_triangle_types=3,
        incidence_bound=dict(variables='H=#vertices of either marked baseline38 type on the selected narrow20-facet; L=#all its other vertices.',
            statements=['No triangle has3 maximal vertices.', 'A quadrilateral with3 maximal vertices has4, and every all-maximal quadrilateral uses the designated end corner at each vertex.', 'Every triangle with2 maximal vertices uses the designated end corner at each of them; its third vertex has marked baseline35.', 'Each maximal vertex has only1 designated corner, shared by these two kinds of short faces.'],
            corner_budget='2*T_2+4*Q_4<=H, where T_2 counts triangles with exactly2 maximal vertices and Q_4 counts all-maximal quadrilaterals.',
            lower_vertex_charge='Every other triangle has at least2 lower vertices and charge2; every other quadrilateral has at least2 lower vertices and charge1. Assign at most1 charge to each lower-vertex corner. There are190 such corners per lower vertex.',
            positive_face_excess='2*T+Q<=H+190*L, where T andQ count triangular and quadrilateral2-faces.',
            h_vector_comparison='35*(H+L)<5*f_2-190*(H+L)<=2*T+Q<=H+190*L; f_2 is the total number of2-faces.',
            strict_frequency_inequality='34*H<155*L', maximal_fraction_upper_bound='155/189',
            inherited_baseline_average_strict_upper_bound='7148/189',
            remaining_target='35',
            scope='Conditional geometric incidence restriction on a selected minimality14 facet of a hypothetical compact21-dimensional arithmetic Coxeter polytope overQ(sqrt5). It does not apply to arbitrary whole-polytope vertex frequencies and does not yet exclude dimension21.'),
        witness_scope='The saved finite triangles and quadrilateral are local arithmetic faces only; no complete finite-volume polytope, simultaneous full vertex star, or new example is claimed.')


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
    print('edge-choice pairs', out['total_edge_choice_pairs'], 'finite quadrilaterals', out['total_finite_quadrilaterals'])
    print(out['incidence_bound']['strict_frequency_inequality'])
