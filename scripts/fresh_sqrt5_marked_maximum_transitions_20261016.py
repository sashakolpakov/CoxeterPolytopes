#!/usr/bin/env python3
"""Exact marked max-to-max minimality14 transitions; no dimension claim.

Exhaust all endpoint attachments between the two marked baseline38 types.
The facet root is retained literally, not merely as an unmarked forest type.
Only local edge feasibility and the resulting port-flow relaxation are tested.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_forest_sensitive_weights_20261015_review as C
from fresh_sqrt5_d4_two_e8_self_edge_20261016 import Q, phi, dot, mv, inverse, pivots, save_vector

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_marked_maximum_transitions_20261016'
BASELINE = 'results/fresh_sqrt5_marked_facet_barrier_20261016.json'


def marked_code(g, marks):
    return tuple(sorted(C.code(c, marks) for c in C.pieces(g)))


def prior_integrity():
    queue, seen, references = [BASELINE], set(), 0
    while queue:
        path = queue.pop()
        if path in seen:
            continue
        seen.add(path)
        for source, digest in json.loads((ROOT / path).read_text()).get('sources', {}).items():
            assert sha256((ROOT / source).read_bytes()).hexdigest() == digest, (path, source)
            references += 1
            if source.endswith('.json'):
                queue.append(source)
    return dict(manifests_checked=len(seen), source_references_checked=references)


def scalar_box():
    # q-sqrt(st)>=0 and q+sqrt(st)<=2 are checked separately in every row.
    # Thus 0<g'<2. With -14<g<=0, b=(g-g')/sqrt5 gives -8<b<0;
    # a=g'+b/phi gives -5<a<2. No real floating point is used.
    assert (Q(5) - Q(4)).sign() > 0  # sqrt5>2
    assert (5 - 7 / phi).sign() > 0
    return [Q(a, b) for a in range(-4, 2) for b in range(-7, 0)]


def cycle_path(old, port, attachments):
    """A path from the old framing root to the new one using only3-bonds."""
    seen = {port}
    todo = [[port]]
    for path in todo:
        for j in sorted(old[path[-1]]):
            if j in seen:
                continue
            extended = path + [j]
            if j in attachments:
                return extended + [21]
            seen.add(j)
            todo.append(extended)
    raise AssertionError('No integral unit-bond cycle containing g')


def witness(gph, port, attachments, nodes, J, u, v, s, t, q, scalar):
    # Eliminate the common20 positive roots. The remaining2x2 matrix is
    # [[s,h],[h,t]], with h=g-q and determinant d<0.
    h = scalar - q
    d = s * t - h * h
    dc = s * t - h.conj() * h.conj()
    assert s.sign() > 0 and t.sign() > 0 and d.sign() < 0 and dc.sign() > 0
    A = [[Q(2 if i == j else -int(j in gph[i])) for j in range(21)] for i in range(21)]
    new = [scalar if j == port else Q(-int(j in attachments)) for j in range(21)]
    B = [row + [new[i]] for i, row in enumerate(A)] + [new + [Q(2)]]
    rays = []
    for omitted, a, b in ((port, -t / d, h / d), (21, h / d, -s / d)):
        common = [-x for x in mv(J, [a * x + b * y for x, y in zip(u, v)])]
        ray = [Q()] * 22
        for i, x in zip(nodes, common):
            ray[i] = x
        ray[port], ray[21] = a, b
        product = mv(B, ray)
        assert product == [-Q(i == omitted) for i in range(22)]
        square = dot(ray, product)
        assert square.sign() < 0
        rays.append(ray)
    mutual = dot(rays[0], mv(B, rays[1]))
    assert mutual == -h / d and mutual.sign() < 0
    mid = [a + b for a, b in zip(*rays)]
    assert dot(mid, mv(B, mid)).sign() < 0
    assert all(x.a.denominator == x.b.denominator == 1 for row in B for x in row)
    return dict(g=scalar.saved(), physical_schur_determinant=d.saved(),
        conjugate_schur_determinant=dc.saved(), physical_signature=[21, 1],
        conjugate_signature=[22, 0], endpoint_rays=[save_vector(x) for x in rays],
        endpoint_squares=[dot(x, mv(B, x)).saved() for x in rays],
        endpoint_mutual_pairing=mutual.saved(),
        edge_inside_marked_facet=True,
        scope='A genuine local finite edge with both future endpoint sectors and integral reflections in an admissible Q(sqrt5) form. No simultaneous edge star, finite volume, or polytope completion is asserted.')


def run():
    baseline = json.loads((ROOT / BASELINE).read_text())
    assert baseline['exact_maximum'] == '38'
    models = dict(C.models())
    g1 = {0: {}} | {i + 1: {j + 1: m for j, m in row.items()} for i, row in models['D20'].items()}
    types = [('A1+D20', g1, 0), ('D21', models['D21'], 19)]
    targets = {marked_code(g, [e]): name for name, g, e in types}
    scalar_candidates = scalar_box()
    rows = []
    for name, g, e in types:
        # The facet root is structurally unique in each of these two forests.
        # Consequently marking {e,port} cannot mistake an exchange of the two
        # marks for an isomorphism preserving the literal facet root.
        assert [i for i in g if marked_code(g, [i]) == marked_code(g, [e])] == [e]
        assert sum((C.base(c, i, j) for c in C.pieces(g) for i, j in combinations(c, 2) if e not in (i, j)), F()) == 38
        for port in sorted(g):
            if port == e:
                continue
            old = {i: {j: m for j, m in row.items() if j != port} for i, row in g.items() if i != port}
            nodes = sorted(old)
            S = [[Q(2 if i == j else -int(j in old[i])) for j in nodes] for i in nodes]
            assert all(x.sign() > 0 for x in pivots(S))
            J = inverse(S)
            u = [Q(-int(i in g[port])) for i in nodes]
            s = 2 - dot(u, mv(J, u))
            accepted = []
            tested = 0
            for degree in range(4):
                for attachments in combinations(nodes, degree):
                    tested += 1
                    newgraph = {i: dict(row) for i, row in old.items()}
                    newgraph[21] = {}
                    for j in attachments:
                        newgraph[j][21] = newgraph[21][j] = 3
                    if any(sum(map(len, c.values())) != 2 * len(c) - 2 for c in C.pieces(newgraph)):
                        continue
                    target = targets.get(marked_code(newgraph, [e]))
                    if target is None:
                        continue
                    v = [Q(-int(i in attachments)) for i in nodes]
                    t = 2 - dot(v, mv(J, v))
                    q = dot(u, mv(J, v))
                    assert s.sign() > 0 and t.sign() > 0 and q.sign() > 0
                    assert (q * q - s * t).sign() >= 0
                    assert (2 - q).sign() >= 0 and ((2 - q) * (2 - q) - s * t).sign() >= 0
                    path = cycle_path(g, port, attachments)
                    allowed = []
                    rejected = []
                    for scalar in scalar_candidates:
                        if (scalar + 14).sign() <= 0 or scalar.sign() > 0:
                            continue
                        h = scalar - q
                        dc = s * t - h.conj() * h.conj()
                        if dc.sign() <= 0:
                            continue
                        if (scalar + 2).sign() < 0 or scalar == -phi:
                            assert (s * t - h * h).sign() < 0
                            allowed.append(dict(g=scalar.saved(), conjugate=scalar.conj().saved(),
                                coxeter_type='DOTTED' if (scalar + 2).sign() < 0 else '5_BOND',
                                physical_schur_determinant=(s * t - h * h).saved(),
                                conjugate_schur_determinant=dc.saved()))
                        else:
                            assert scalar.sign() < 0 and (scalar + 1).sign() > 0
                            rejected.append(scalar.saved())
                    same_port = marked_code(g, [e, port]) == marked_code(newgraph, [e, 21])
                    rec = dict(target_marked_type=target, attachments=list(attachments),
                        copied=list(attachments) == sorted(g[port]),
                        same_marked_port_orbit=same_port,
                        old_residual=s.saved(), new_residual=t.saved(), cross_projection=q.saved(),
                        integral_cycle_path=path,
                        all_narrow_admissible_products=allowed,
                        rejected_noncoxeter_products=rejected)
                    if allowed:
                        rec['actual_edge'] = witness(g, port, attachments, nodes, J, u, v, s, t, q, Q(*allowed[0]['g']))
                    accepted.append(rec)
            assert tested == 1351
            rows.append(dict(source_marked_type=name, facet_root=e, port=port,
                candidate_attachment_subsets=tested, complete_marked_reattachments=accepted,
                narrow_return_exists=any(r['all_narrow_admissible_products'] for r in accepted),
                narrow_same_port_return_exists=any(r['all_narrow_admissible_products'] and r['same_marked_port_orbit'] for r in accepted)))
    first = [r for r in rows if r['source_marked_type'] == 'A1+D20']
    second = [r for r in rows if r['source_marked_type'] == 'D21']
    assert len(first) == len(second) == 20
    assert all(r['narrow_same_port_return_exists'] for r in first)
    assert sum(len(r['complete_marked_reattachments']) for r in first) == 43
    assert sum(len(r['complete_marked_reattachments']) for r in second) == 23
    flow = []
    for row in first:
        candidates = [(i, r) for i, r in enumerate(row['complete_marked_reattachments']) if r['all_narrow_admissible_products'] and r['same_marked_port_orbit']]
        i, chosen = min(candidates, key=lambda z: (not z[1]['copied'], z[0]))
        flow.append(dict(source_port=row['port'], target_same_marked_port_orbit=True,
            reattachment_index=i, g=chosen['all_narrow_admissible_products'][0]['g'],
            undirected_self_type_edge_mass='1/2', induced_port_incidence='1'))
    forbidden = [r['port'] for r in second if not r['narrow_return_exists']]
    assert forbidden == [1, 2, 4, 5, 6]
    paths = [f'scripts/{TAG}.py', BASELINE,
        'scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py',
        'scripts/fresh_sqrt5_d4_two_e8_self_edge_20261016.py']
    return dict(status='PASS_EXACT_MARKED_NARROW_TRANSITIONS_RETAIN_BASELINE38_OBSTRUCTION',
        sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
        prior_source_integrity=prior_integrity(),
        dimension=21, facet_dimension=20, facet_weight_threshold='35', baseline_maximum='38',
        marked_types=[dict(name=name, facet_root=e,
            edges=[[i, j, m] for i in g for j, m in g[i].items() if i < j]) for name, g, e in types],
        scalar_enumeration=dict(field='Q(phi), phi^2=phi+1', a_range=[-4, 1], b_range=[-7, -1], count=42,
            justification='Every row has q-sqrt(st)>=0 and q+sqrt(st)<=2. Conjugate positive definiteness gives0<gprime<2. Together with-14<g<=0, this proves-8<b<0 and-5<a<2. A saved cycle of3-bonds with framing productg proves g belongs to O_K=Z[phi].'),
        attachment_completeness='Both target finite diagrams are simply laced with maximum degree3. All1351 subsets of the20 literal common roots of sizes0..3 are exhausted, and the retained facet root is included in the exact marked tree-isomorphism code.',
        port_records=rows, total_port_cases=len(rows), total_attachment_subsets=sum(r['candidate_attachment_subsets'] for r in rows),
        total_marked_reattachments=sum(len(r['complete_marked_reattachments']) for r in rows),
        D21_ports_without_maximum_return=forbidden,
        A1_D20_all20_ports_have_same_marked_port_return=True,
        baseline38_balanced_port_flow=flow,
        surviving_relaxation='Put vertex mass1 on A1+D20 with its isolatedA1 marked as the facet root. Each of its20 port orbits has a minimality14 arithmetic self-return and undirected mass1/2, giving incidence1. The central-E8 inequality is0<=0. This finite-state marked edge-flow relaxation still permits baseline38.',
        scope='Exact local marked max-to-max edge analysis and a restricted port-flow obstruction. Separate edge witnesses need not form one simultaneous star or one global narrow facet. No finite-volume realization, global dimension bound, or obstruction to every facet-weight strategy is asserted.')


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
    print('D21 forbidden ports:', out['D21_ports_without_maximum_return'])
    print('marked reattachments:', out['total_marked_reattachments'])
