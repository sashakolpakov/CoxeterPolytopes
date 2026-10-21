#!/usr/bin/env python3
"""Independent central-E8 narrow returns and finite endpoint-type counts."""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import product
import json

import fresh_h9_beta131_h300_symmetry_20261016_review as I
import fresh_sqrt5_forest_sensitive_weights_20261015_review as C

K, ROOT = I.K, I.ROOT
TAG = 'fresh_sqrt5_narrow_e8_ports_20261016_review'
PRIMARY = 'results/fresh_sqrt5_narrow_e8_ports_20261016.json'


def chain(n):
    return [[K(2) if i == j else K(-1) if abs(i-j) == 1 else K() for j in range(n)] for i in range(n)]


def subtract(multiset, removed):
    result = multiset.copy()
    result.subtract(removed)
    if any(v < 0 for v in result.values()):
        return None
    return +result


def expanded(counts):
    return tuple(sorted(counts.elements()))


def run():
    primary = I.read(PRIMARY)
    assert primary['status'] == 'PASS_EXACT_NARROW_CENTRAL_E8_RETURN_EXCLUSION_AND_FACET_FREQUENCY_BOUND'
    inverses = {n: I.inverse(chain(n))[0] for n in (1, 2, 4)}
    # Old central E8 root attaches to a fixed end of each A1,A2,A4 arm.
    # Each new central root can choose either end of the two longer arms.
    assert sum((inverses[n][0][0] for n in (1, 2, 4)), K()) == K(F(59, 30))
    residual = K(2)-sum((inverses[n][0][0] for n in (1, 2, 4)), K())
    assert residual == K(F(1, 30))
    phip = K(0, 1).conjugate()
    rows = []
    for a2, a4 in product((0, 1), (0, 3)):
        q = inverses[1][0][0]+inverses[2][0][a2]+inverses[4][0][a4]
        lower, upper = q-residual, q+residual
        assert 1 <= lower < upper <= 2
        # If -14<g<=0 and 1<g'<2, then -16<g-g'<-1.
        # Since g-g'=b*sqrt5, the integer b lies in[-7,-1].
        assert 8*K(-1, 2) > 16 and K(-1, 2) > 1
        # g'=a+b*phi', b in[-7,-1], and -1<phi'<0 imply -6<a<2.
        assert -1 < phip < 0
        candidates = []
        for b in range(-7, 0):
            for a in range(-5, 2):
                g = K(a, b)
                if -14 < g <= 0 and lower < g.conjugate() < upper:
                    assert -1 < g < 0
                    candidates.append(dict(a=a, b=b, rejection='strictly between the right-angle and label3 products'))
        rows.append(dict(arm_ends=[a2, a4], cross_cost=I.enc(q),
                         conjugate_lower=I.enc(lower), conjugate_upper=I.enc(upper),
                         integral_candidates=candidates, compatible_returns=0))
    costs = sorted(I.dec([row['cross_cost']])[0] for row in rows)
    assert costs == [K(F(v, 30)) for v in (31, 41, 49, 59)]
    assert sum(len(row['integral_candidates']) for row in rows) == 1
    assert next(row['integral_candidates'][0] for row in rows if row['integral_candidates'])['a'] == 1
    assert next(row['integral_candidates'][0] for row in rows if row['integral_candidates'])['b'] == -1
    primary_returns = next(row for row in primary['same_type_reattachments'] if row['port'] == 7)['complete_same_type_reattachments']
    assert len(primary_returns) == 4
    for row in rows:
        cost = K.from_json(row['cross_cost'])
        match = next(r for r in primary_returns if K(F(r['cross_projection'][0]), F(r['cross_projection'][1])) == cost)
        assert 12 in match['attachments'] and 12 in match['common_neighbors_for_integral_triangle']
        assert K(*map(F, match['residual'])) == residual
        assert len(match['narrow_scalar_candidates']) == len(row['integral_candidates'])
        for candidate, checked in zip(match['narrow_scalar_candidates'], row['integral_candidates']):
            assert K(*map(F, candidate['g'])) == K(checked['a'], checked['b'])
            assert candidate['coxeter_allowed'] is False
    # Complete one-root extensions of S are found by deleting each root of
    # each connected finite component, then appending the untouched S pieces.
    models = {name: graph for name, graph in C.models() if len(graph) <= 21}
    canonical = {C.code(graph): name for name, graph in models.items()}
    assert len(models) == len(canonical) == 68
    support = Counter(['A1', 'A2', 'A4', 'D5', 'E8'])
    source_type = tuple(sorted(['D5', 'E8', 'E8']))
    deleted = {}
    for name, graph in models.items():
        for root in graph:
            rest = {i: {j:m for j,m in row.items() if j != root} for i,row in graph.items() if i != root}
            pieces = Counter(canonical[C.code(part)] for part in C.pieces(rest))
            deleted[name, root] = pieces
    completions = set()
    for (name, root), pieces in deleted.items():
        untouched = subtract(support, pieces)
        if untouched is not None:
            untouched[name] += 1
            types = expanded(untouched)
            assert sum(len(models[t]) for t in types) == 21
            completions.add(types)
    assert len(completions) == 29 and source_type in completions
    endpoint_rows = []
    for types in sorted(completions):
        ports = []
        for occurrence, name in enumerate(types):
            spectators = Counter(types)
            spectators[name] -= 1
            spectators = +spectators
            for root in models[name]:
                if spectators+deleted[name, root] == support:
                    ports.append([occurrence, name, root])
        assert ports
        endpoint_rows.append(dict(types=list(types), multiplicity=len(ports), ports=ports))
    source_row = next(row for row in endpoint_rows if tuple(row['types']) == source_type)
    assert source_row['multiplicity'] == 2
    alternatives = [row for row in endpoint_rows if tuple(row['types']) != source_type]
    assert len(alternatives) == 28 and max(row['multiplicity'] for row in alternatives) == 4
    maximizers = [row['types'] for row in alternatives if row['multiplicity'] == 4]
    assert sorted(maximizers) == sorted([
        sorted(['A2', 'A2', 'A4', 'D5', 'E8']),
        sorted(['A4', 'A4', 'D5', 'E8']),
        sorted(['A1', 'A2', 'D5', 'D5', 'E8'])])
    canonical_type = lambda names: tuple(sorted('A2' if n == 'I2(3)' else 'B2' if n == 'I2(4)' else n for n in names))
    primary_counts = {canonical_type(row['vertex_type']): row['multiplicity']
                      for row in primary['central_edge_endpoint_census']['records']}
    assert primary_counts == {tuple(row['types']): row['multiplicity'] for row in endpoint_rows}
    facet = primary['selected_facet_incidence']
    assert facet['ambient_dimension'] == 21 and facet['facet_dimension'] == 20
    assert facet['central_root_indices'] == [7, 15]
    assert [row['facet_root'] for row in facet['facet_root_roles']] == list(range(21))
    for row in facet['facet_root_roles']:
        assert row['central_outgoing_edges_inside_facet'] == (1 if row['facet_root'] in (7, 15) else 2)
    # Facet dimension20 has20 outgoing in-facet edge ports at each vertex.
    # At a source vertex, removing at most one facet-normal root leaves at
    # least one of its two central E8 ports. These edges cannot return to D.
    # Each target U has at most the counted m(U) matching deleted-root ports.
    # Therefore Dc+2Do<=sum m(U)VU<=4O, giving D<=4O andD/(D+O)<=4/5.
    paths = [f'scripts/{TAG}.py', PRIMARY,
             'scripts/fresh_h9_beta131_h300_symmetry_20261016_review.py',
             'scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py',
             'src/coxeter_polytopes/exact/golden.py']
    return dict(status='PASS_INDEPENDENT_CENTRAL_E8_NARROW_RETURN_AND_FACET_PORT_BOUND',
                sources={p: I.sha(p) for p in paths}, central_reattachments=rows,
                endpoint_types=endpoint_rows, endpoint_type_count=29,
                alternative_endpoint_types=28, maximum_alternative_ports=4,
                multiplicity_four_types=maximizers,
                selected_facet_inequality='D_central+2D_other <= sum_U m(U)V_U <=4O',
                selected_facet_source_proportion_bound='4/5',
                scope='Independent central-E8 four-reattachment narrow integral enumeration and exhaustive finite component/deletion count, plus selected-narrow-facet edge incidence argument. Other copied-port sharp bounds and the full primary scalar tables are outside this review. The narrow-facet existence theorem is a published dependency. No whole-polytope frequency bound or new dimension ceiling.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    output = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(output.read_text())
    else:
        output.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['endpoint_type_count'], result['maximum_alternative_ports'])
