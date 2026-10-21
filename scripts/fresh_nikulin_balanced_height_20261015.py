#!/usr/bin/env python3
"""Exact balanced-height obstruction to every preset E8 unit-charge pair.

Self-contained Fraction matrix arithmetic. No project imports, numerical tests,
label cutoff, or inference of a new dimension upper bound.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_nikulin_balanced_height_20261015'


def star(arms):
    n = 1 + sum(arms)
    G = [[F(2 if i == j else 0) for j in range(n)] for i in range(n)]
    k = 1
    for arm in arms:
        previous = 0
        for _ in range(arm):
            G[k][previous] = G[previous][k] = F(-1)
            previous, k = k, k + 1
    return G


def inverse(A):
    n = len(A)
    B = [list(row) + [F(i == j) for j in range(n)] for i, row in enumerate(A)]
    for j in range(n):
        k = next(i for i in range(j, n) if B[i][j])
        B[j], B[k] = B[k], B[j]
        pivot = B[j][j]
        B[j] = [x / pivot for x in B[j]]
        for i in range(n):
            if i != j:
                factor = B[i][j]
                B[i] = [a - factor * b for a, b in zip(B[i], B[j])]
    assert [row[:n] for row in B] == [[F(i == j) for j in range(n)] for i in range(n)]
    J = [row[n:] for row in B]
    assert multiply(A, J) == [[F(i == j) for j in range(n)] for i in range(n)]
    return J


def multiply(A, B):
    return [[sum((a * b for a, b in zip(row, col)), F(0))
             for col in zip(*B)] for row in A]


def pivots(A):
    A = [row[:] for row in A]
    result = []
    while A:
        pivot = A[0][0]
        result.append(pivot)
        assert pivot > 0
        A = [[A[i][j] - A[i][0] * A[0][j] / pivot
              for j in range(1, len(A))] for i in range(1, len(A))]
    return result


def field_data(arms, kernel):
    G = star(arms)
    delta = list(map(F, kernel))
    assert len(G) == len(delta) and min(delta) == 1
    assert all(x > 0 for x in delta)
    assert multiply(G, [[x] for x in delta]) == [[F(0)] for _ in G]
    deletions = []
    for omitted in range(len(G)):
        kept = [i for i in range(len(G)) if i != omitted]
        B = [[G[i][j] for j in kept] for i in kept]
        diag_pivots = pivots(B)
        J = inverse(B)
        assert all(x >= 0 for row in J for x in row)
        costs = [F(0)] * len(G)
        for a, k in enumerate(kept):
            costs[k] = J[a][a]
        deletions.append(dict(omitted=omitted, kept=kept, inverse=J,
                              pivots=diag_pivots, costs=costs))
    return G, delta, deletions


def strings(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {k: strings(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [strings(x) for x in value]
    return value


def run():
    G, delta, deletions = field_data([1, 2, 5], [6, 3, 4, 2, 5, 4, 3, 2, 1])
    target = F(16, 9)
    cases = []
    legal_count = 0
    for i, j in combinations(range(9), 2):
        candidates = []
        for k in range(9):
            a, b = deletions[i]['costs'][k], deletions[j]['costs'][k]
            legal = a <= 2 and b <= 2
            legal_count += int(legal)
            # Each endpoint extension is PSD precisely when its Schur residual
            # 2-cost is nonnegative. These are independently literal unit rows.
            candidates.append(dict(node=k, row_p=[int(t == k) for t in range(9)],
                                   row_q=[int(t == k) for t in range(9)],
                                   heights=[delta[k], delta[k]], M=F(1),
                                   weighted_energy=F(0), endpoint_costs=[a, b],
                                   endpoint_schur_residuals=[2-a, 2-b],
                                   total_cost=a+b, endpoint_budget_legal=legal,
                                   violates_unit_charge=legal and a+b < target))
        allowed = [r for r in candidates if r['endpoint_budget_legal']]
        witness = min(allowed, key=lambda r: (r['total_cost'], r['node']))
        assert witness['violates_unit_charge']
        assert witness['total_cost'] > 0
        required_alpha = target / witness['total_cost']
        assert required_alpha > 1
        cases.append(dict(omissions=[i,j], unit_common_rows=candidates,
                          selected_witness=witness,
                          necessary_charge_from_common_rows=required_alpha))
    assert len(cases) == 36
    best = min(r['necessary_charge_from_common_rows'] for r in cases)
    best_pairs = [r['omissions'] for r in cases if r['necessary_charge_from_common_rows'] == best]
    assert best == F(8,7) and best_pairs == [[1,6]]
    current = next(r for r in cases if r['omissions'] == [1,6])['selected_witness']
    assert current['node'] == 8 and current['endpoint_costs'] == [F(8,9),F(2,3)]
    assert current['total_cost'] == F(14,9) and target-current['total_cost'] == F(2,9)

    G7, delta7, del7 = field_data([1,3,3], [4,2,3,2,1,3,2,1])
    e7_costs = [del7[3]['costs'][4], del7[6]['costs'][4]]
    assert e7_costs == [F(1,2),F(1)]
    assert delta7[4] == 1 and sum(e7_costs) == F(3,2) < F(14,9)
    # Algebraic interpretation only: if a matching 8/7 upper lemma were proved,
    # its accounting constant would be277/7 and the integer threshold N=166.
    assert 18*best+19 == F(277,7)
    assert (F(166)-F(277,7))/9 == F(295,21) > 14
    assert (F(165)-F(277,7))/9 < 14
    deps = ['scripts/fresh_nikulin_unit_charge_20261015.py',
            'results/fresh_nikulin_unit_charge_20261015.json',
            'results/fresh_nikulin_unit_charge_20261015.md',
            f'scripts/{TAG}.py']
    old = json.loads((ROOT/deps[1]).read_text())
    ordinary_e8 = next(r for r in old['exceptional_cases'] if r['type']=='E8' and r['forbidden'] is None)
    ordinary_e7 = next(r for r in old['exceptional_cases'] if r['type']=='E7' and r['forbidden'] is None)
    assert ordinary_e8['omissions'] == [1,6] and ordinary_e7['omissions'] == [3,6]
    assert [[F(x) for x in row] for row in ordinary_e8['gram']] == G
    assert [[F(x) for x in row] for row in ordinary_e7['gram']] == G7
    result = dict(status='PASS_BALANCED_HEIGHT_UNIT_CHARGE_OBSTRUCTION',
                  sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in deps},
                  arithmetic='Standalone exact Fraction matrices; no project imports',
                  E8=dict(gram=G,kernel=delta,target=target,deletions=deletions,
                          distinct_omission_pairs=cases,
                          identical_omissions='For i=j, p=q=e_i gives C=W=0 and M=1, violating16/9.',
                          best_necessary_charge=best,best_pairs=best_pairs),
                  current_E7=dict(gram=G7,kernel=delta7,omissions=[3,6],common_row_node=4,
                                  heights=[F(1),F(1)],M=F(1),weighted_energy=F(0),
                                  endpoint_costs=e7_costs,total_cost=F(3,2),target=F(14,9),
                                  deficit=F(1,18)),
                  counts=dict(distinct_E8_pairs=36,common_unit_rows=324,
                              endpoint_budget_legal_rows=legal_count,violated_pairs=36,
                              trivial_identical_pairs=9),
                  scope='Refutes the proposed uniform local unit-charge lemma even at M=1 '
                        'for every preset E8 omission pair. The rows obey the stated nonnegative '
                        'contact gap and each endpoint cost<=2; only magnitudes0/1 occur. '
                        'Does not prove realization of an entire finite-volume polytope, does not '
                        'exclude a global argument using cross-factor slack or different choices, '
                        'and does not establish an8/7 upper lemma or any new dimension bound.',
                  hypothetical_eight_sevenths=dict(accounting_constant=F(277,7),
                                                   integer_N_threshold=166,
                                                   conditional_dimension_bound=166,
                                                   upper_lemma_proved=False))
    return strings(result)


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    result=run()
    destination=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(destination.read_text())
    else:
        destination.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'])
