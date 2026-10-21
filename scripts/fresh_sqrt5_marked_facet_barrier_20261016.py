#!/usr/bin/env python3
"""Exact marked-facet baseline barrier; no optimizer or dimension claim.

The dimension20 facet threshold is35. The inherited nonnegative baseline
already costs38 at permitted marked vertex types unaffected by the new
central-E8 incidence inequality. Therefore this proposed restricted model
cannot work without additional structural information or a changed baseline.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import comb
from pathlib import Path
import fresh_sqrt5_forest_sensitive_weights_20261015_review as C

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_marked_facet_barrier_20261016'
NARROW = 'results/fresh_sqrt5_narrow_e8_ports_20261016.json'


def h_threshold(m):
    choose2 = lambda k: comb(k, 2) if k >= 2 else 0
    paired = [F(choose2(j) + choose2(m - j), 2) for j in range(m + 1)]
    minimum = min(paired)
    assert paired[0] > minimum
    return dict(dimension=m, paired_f2_over_f0_coefficients=list(map(str, paired)),
        minimum_paired_coefficient=str(minimum), corners_per_vertex=comb(m, 2),
        strict_average_weight_threshold=str(5 * minimum - comb(m, 2)),
        strictness='h_0=h_m=1; their paired coefficient exceeds the minimum. All h_j are nonnegative.',
        identity='f_0=sum h_j; f_2=sum binom(m-j,2)h_j; total corner weight>=5f_2-binom(m,2)f_0.')


def run():
    narrow = json.loads((ROOT / NARROW).read_text())
    assert narrow['status'] == 'PASS_EXACT_NARROW_CENTRAL_E8_RETURN_EXCLUSION_AND_FACET_FREQUENCY_BOUND'
    hs = [h_threshold(m) for m in (20, 21)]
    assert [h['strict_average_weight_threshold'] for h in hs] == ['35', '40']
    models = [(name, g) for name, g in C.models() if len(g) <= 21]
    assert len(models) == 68
    data = []
    for name, g in models:
        weights = {(a, b): C.base(g, a, b) for a, b in combinations(g, 2)}
        total = sum(weights.values(), F(0))
        debits = {a: sum((w for pair, w in weights.items() if a in pair), F(0)) for a in g}
        # Direct pair summation excludes the marked facet root at the vertex.
        for a in g:
            assert total - debits[a] == sum((w for pair, w in weights.items() if a not in pair), F(0))
        data.append(dict(name=name, rank=len(g), total=total, debits=debits,
            edges=[[i, j, m] for i in g for j, m in g[i].items() if i < j], weights=weights))
    # All-rank-to21 unbounded finite-forest DP, retaining every optimal multiset.
    dp = [(F(0), {()})]
    for rank in range(1, 22):
        best = None
        forests = set()
        for c in data:
            if c['rank'] > rank:
                continue
            value, residuals = dp[rank - c['rank']]
            value += c['total']
            choices = {tuple(sorted(r + (c['name'],))) for r in residuals}
            if best is None or value > best:
                best, forests = value, choices
            elif value == best:
                forests |= choices
        assert best is not None
        dp.append((best, forests))
    marked = []
    for c in data:
        rest_value, forests = dp[21 - c['rank']]
        for root, debit in sorted(c['debits'].items()):
            value = c['total'] - debit + rest_value
            marked.append(dict(marked_component=c['name'], marked_root=root,
                marked_component_rank=c['rank'], baseline_total=str(c['total']),
                removed_pair_total=str(debit), residual_rank=21 - c['rank'],
                residual_optimum=str(rest_value), residual_optimal_forests=sorted(forests),
                full_facet_corner_total=str(value)))
    assert len(marked) == 724
    maximum = max(F(r['full_facet_corner_total']) for r in marked)
    assert maximum == 38
    maximizers = [r for r in marked if F(r['full_facet_corner_total']) == maximum]
    assert [(r['marked_component'], r['marked_root']) for r in maximizers] == [('A1', 0), ('D21', 19)]
    assert maximizers[0]['residual_optimal_forests'] == [('D20',)]
    assert maximizers[1]['residual_optimal_forests'] == [()]
    # The full refined incidence inequality is inactive at both maximizers.
    incoming = {tuple(r['vertex_type']): r['multiplicity'] for r in narrow['central_edge_endpoint_census']['records']}
    aliases = {name: C.family(g) for name, g in models}
    witnesses = []
    for r in maximizers:
        for residual in r['residual_optimal_forests']:
            types = tuple(sorted(aliases[name] for name in residual + (r['marked_component'],)))
            assert types != ('D5', 'E8', 'E8') and incoming.get(types, 0) == 0
            witnesses.append(dict(vertex_type=types, marked_component=r['marked_component'],
                marked_root=r['marked_root'], baseline_facet_total=r['full_facet_corner_total'],
                central_E8_inequality_left=0, central_E8_inequality_right=0,
                feasible_frequency='Assign frequency1 to this marked vertex type and0 to every other type in the incidence relaxation.'))
    exceptional_names = ('D5', 'E8', 'E8')
    by_name = {c['name']: c for c in data}
    exceptional_total = sum((by_name[name]['total'] for name in exceptional_names), F(0))
    exceptional_marked = []
    for ci, name in enumerate(exceptional_names):
        for root, debit in sorted(by_name[name]['debits'].items()):
            exceptional_marked.append(dict(component=ci, type=name, root=root,
                facet_corner_total=str(exceptional_total - debit)))
    assert max(F(r['facet_corner_total']) for r in exceptional_marked) == 34
    paths = [f'scripts/{TAG}.py', NARROW, 'scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py']
    return dict(status='PASS_EXACT_MARKED_FACET_BASELINE38_BARRIER_ABOVE_THRESHOLD35',
        sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
        h_vector_thresholds=hs,
        baseline_rule='One unit for connected marked roots at graph distance<=2; H4 gets one additional unit when the pair is its5-bond or has distance3. Disconnected pairs have zero baseline. Facet corners exclude the marked facet root.',
        finite_components=[dict(type=c['name'], rank=c['rank'], edges=c['edges'],
            total=str(c['total']), marked_root_debits=[dict(root=i, debit=str(v)) for i, v in sorted(c['debits'].items())],
            pair_weights=[dict(pair=list(pair), weight=str(w)) for pair, w in sorted(c['weights'].items())]) for c in data],
        unmarked_forest_dp=[dict(rank=i, optimum=str(v), all_optimal_forests=sorted(fs)) for i, (v, fs) in enumerate(dp)],
        marked_component_cases=marked, exact_maximum=str(maximum), maximizers=maximizers,
        incidence_feasible_maximizers=witnesses,
        incidence_constrained_baseline_maximum='38',
        exceptional_marked_baseline=exceptional_marked,
        nonnegative_extension_obstruction='Any rule retaining this baseline and adding only nonnegative corner bonuses has cost>=38 at each saved witness. The new central-E8 incidence inequality is0<=0 there. Hence no vertex-frequency optimization using only that incidence restriction can force the required upper bound35.',
        scope='Exact obstruction to the proposed marked-facet baseline-plus-nonnegative-bonuses model with only the new central-E8 incidence restriction. The baseline alone is not asserted to satisfy all face-payment requirements. No actual polytope realization of the extremizing frequency, no failure of every facet-weight strategy, and no dimension bound is asserted.')


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
