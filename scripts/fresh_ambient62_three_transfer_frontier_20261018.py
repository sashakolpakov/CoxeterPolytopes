#!/usr/bin/env python3
"""Compose three local transfers and recount the exact ambient62 frontier.

Standard library only. No geometric census, forest DP, or optimizer is run.
The rank62 forest upper certificate is inherited; its changed actual costs
are checked pointwise against its unchanged representative rule.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from math import prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_ambient62_three_transfer_frontier_20261018'
FRONT = 'fresh_ambient62_payment_frontier_20261018'
CURRENT = 'fresh_sparse_coupon_order11_13_dimension62_20261018'
BUDGET = 'fresh_sparse_coupon_rank63_forest_budget_20261018'
D13 = 'fresh_order13_d_center_pair_zero_transfer_20261018'
B16 = 'fresh_order16_affine_b_singleton_projection_transfer_20261018'
D17 = 'fresh_order17_e6_unit_coupon_zero_transfer_20261018'
STEMS = (FRONT, CURRENT, BUDGET, D13, B16, D17)
N, D = 62, 61
DROPS = {12: F(1, 2), 15: F(1), 16: F(2)}


def load_inputs():
    docs = {s: json.loads((ROOT / f'results/{s}.json').read_text()) for s in STEMS}
    for s, doc in docs.items():
        assert doc['status'].startswith('PASS'), s
        for p, h in doc['sources'].items():
            assert sha256((ROOT / p).read_bytes()).hexdigest() == h, (s, p)
    return docs


def interfaces(docs):
    a, b, c = (docs[s] for s in (D13, B16, D17))
    assert docs[CURRENT]['dimension_bound'] == 62
    assert all(x['triangle_minimum'] == 2 and x['quadrilateral_minimum'] == 1 for x in (a, b, c))
    assert a['all_family_order13_costs'] == dict(A='65', B='73', C='73', D='73')
    assert a['old_affine_total'] == '147/2' and a['new_affine_total'] == '73'
    assert a['new_normalized_envelope'] == '73/13'
    assert a['changed_pairs'] == [[0, 8]] and a['old_weight'] == '1/2' and a['new_weight'] == '0'
    assert sum((F(r['weight']) for r in a['weighted_pairs']), F()) == 73
    coupon_roles = [('whole_A2', 'either_endpoint', 2, F(1), F(2)),
                    ('whole_A3', 'unique_middle', 1, F(1, 2), F(1))]
    assert len(a['upgraded_coupons']) == 2
    for r, (small, role, multiplicity, old, new) in zip(a['upgraded_coupons'], coupon_roles):
        assert (r['small'], r['small_role'], r['large'], r['large_role']) == (
            small, role, 'whole_D10', 'fork_center')
        assert r['eligible_marked_pairs'] == multiplicity
        assert F(r['old_weight']) == F(1, 2) and F(r['new_weight']) == 1
        assert F(r['old_component_pair_total']) == old and F(r['new_component_pair_total']) == new
        assert F(r['increment']) == new - old == F(multiplicity, 2)
    assert len(a['inherited_branches_with_independent_units']) == 17
    assert all(F(r['distinct_payment']) == 1 for r in a['inherited_branches_with_independent_units'])
    assert b['new_order16_family_costs'] == dict(A='80', B='105', C='103', D='106')
    assert b['old_cost'] == 107 and b['new_cost'] == 105 and b['new_order16_maximum'] == '106'
    assert b['new_rank15_envelope'] == '53/8' and b['deleted_pairs'] == [[1, 11], [2, 12]]
    assert b['finite_weights_unchanged'] and b['added_finite_coupons'] == []
    assert len(b['retained_unit_pairs']) == 105 and b['affected_maps'] == 10
    assert c['new_order17_costs'] == dict(A=85, B=111, C=108, D=116)
    assert c['old_cost'] == '118' and c['new_cost'] == '116' and c['new_order17_envelope'] == '116/17'
    assert c['changed_pairs'] == [[3, 15], [3, 16], [9, 13], [9, 14]]
    weights = {tuple(r['pair']): F(r['weight']) for r in c['weighted_pairs']}
    assert len(weights) == 136 and sum(weights.values(), F()) == 116
    assert all(weights[tuple(p)] == 0 for p in c['changed_pairs'])
    assert all(weights[tuple(p)] == 1 for p in c['retained_short8_pairs'])
    coupon = c['finite_coupon']
    assert (coupon['small_type'], coupon['small_mark'], coupon['large_type'], coupon['large_mark']) == (
        'whole E6', 'endpoint of either length2 arm', 'whole D11', 'short fork leaf')
    assert coupon['marked_pairs_per_component_pair'] == 4
    assert [F(coupon[k]) for k in ('old_weight', 'new_weight', 'old_total', 'new_total', 'added_total')] == [F(1, 2), 1, 2, 4, 2]
    assert c['compatible_with_B16_singleton_deletions'] and c['compatible_with_D13_center_zero']
    assert len(c['changed_cases']) == 4 and all(F(r['total']) == 1 for r in c['changed_cases'])
    return dict(order13=a['all_family_order13_costs'], order16=b['new_order16_family_costs'],
                order17=c['new_order17_costs'], D13_upgraded_coupons=a['upgraded_coupons'],
                D17_upgraded_coupon=coupon, inherited_D17_priority=c['actual_corner_priority'],
                all_family_total_savings_by_rank={str(r): str(v) for r, v in DROPS.items()},
                triangle_minimum=2, quadrilateral_minimum=1)


def d13_extra(a, b):
    return sum(2 if small == ('A', 2) else 1 if small == ('A', 3) else 0
               for small, large in ((a, b), (b, a)) if large == ('D', 10))


def finite_check(docs):
    budget, front = docs[BUDGET], docs[FRONT]
    scenario = front['finite_scenarios'][1]
    assert scenario['conditional_changes'] == dict(D13_zero_coupons=True)
    assert scenario['exact_maximum'] == '441' and scenario['curvature'] == '434'
    assert scenario['actual_attaining_forests'] == [['D17', 'D17', 'D28']]
    assert budget['representative_allowance']['actual_total'] == '2'
    assert budget['representative_allowance']['upper_total'] == '4'
    checks = 0
    changed = []
    for row in budget['type_pair_checks']:
        a, b = map(tuple, row['types'])
        if max(a[1], b[1]) > N:
            continue
        ra, rb = map(tuple, row['representatives'])
        extra = 4 if {a, b} == {('E', 6), ('D', 11)} else 0
        actual = row['actual_doubled'] + d13_extra(a, b) + extra
        upper = row['upper_doubled'] + d13_extra(ra, rb)
        assert actual <= upper
        checks += 1
        if extra:
            assert actual == upper == 8
            changed.append(dict(types=row['types'],new_actual_doubled=actual,unchanged_upper_doubled=upper))
    assert checks == 35721 and len(changed) == 2
    comparisons = []
    for row in scenario['alternatives']:
        counts = Counter(row['forest'])
        actual = row['actual_doubled'] + 4 * counts['E6'] * counts['D11']
        assert actual <= row['doubled_upper'] <= 882
        comparisons.append(dict(counts=row['counts'],new_actual_doubled=actual,
                                inherited_upper_doubled=row['doubled_upper']))
    assert len(comparisons) == 1227
    witness = ['D17', 'D17', 'D28']
    # Literal costs in the frozen finite rule; no cross coupon involves these pairs.
    assert docs[CURRENT]['finite_component_costs']['D17'] == 125
    D28 = 8 * 28 - 33 - max(0, 28 - 32)
    assert D28 == 191 and 125 + 125 + D28 == 441
    assert any(r['forest'] == witness and r['actual_doubled'] == r['doubled_upper'] == 882
               for r in scenario['alternatives'])
    return dict(exact_maximum='441', curvature='434', finite_slack='-7',
                actual_attaining_forest=witness, actual_component_costs=[125, 125, D28],
                simultaneous_type_pairs=checks, changed_E6_pairs=changed,
                complete_large_vectors=1227, inherited_interaction_keys=921,
                saved_formal_attainer_comparisons=comparisons,
                inherited_finite_scenario_index=1, Bellman_recomputed=False,
                scope='The D13 increments are in frozen scenario1. The E6/D11 increment fills its already present upper allowance. Pointwise simultaneous replacement preserves that complete upper bound; the unchanged actual D-only witness attains441.')


def shapes():
    def large(left, minimum, prefix):
        if left == 0:
            yield prefix
        for r in range(minimum, left + 1):
            yield from large(left - r, r, prefix + (r,))
    for low in range(9):
        yield from large(D - low, 9, (low,) if low else ())


def old_envelope(r):
    if r <= 9:
        return F(r, 2)
    small = {10:F(54,11),11:F(21,4),12:F(147,26),13:F(6),14:F(94,15),15:F(107,16),16:F(118,17)}
    if r in small:
        return small[r]
    m = r + 1
    costs = [5*m]
    for offset, cap, first in ((-6,13,8),(-15,14,6),(3,12,10)):
        costs.append(5*m + offset + sum(min(m-k, cap) for k in range(first, first+4)))
    return F(max(costs), m)


def summary(rows):
    return dict(count=len(rows), negative=sum(F(r['normalized_margin']) < 0 for r in rows),
                zero=sum(F(r['normalized_margin']) == 0 for r in rows),
                minimum=min(rows,key=lambda r:(F(r['normalized_margin']),r['ranks'])),
                minimum_absolute=min(rows,key=lambda r:(F(r['absolute_margin']),r['ranks'])),
                negative_rows=sorted((r for r in rows if F(r['normalized_margin']) < 0),
                                     key=lambda r:(F(r['normalized_margin']),r['ranks'])))


def scalar_check(front):
    saved = {tuple(r['ranks']):r for r in front['exact_tail_rows']}
    domain = sorted(shapes())
    assert len(domain) == len(set(domain)) == len(saved) == 2641 and set(domain) == set(saved)
    current, conditional = [], []
    for rs in domain:
        row = saved[rs]
        h = [1]
        for r in rs:
            # Direct bounded convolution, separate from the producer's sliding window.
            h = [sum(h[max(0,j-r):min(j+1,len(h))]) for j in range(len(h)+r)]
        V = prod(r+1 for r in rs)
        assert len(h) == D+1 and h == h[::-1] and sum(h) == V
        G = sum(v*min(j,D-j) for j,v in enumerate(h))
        assert 2*G == D*V - sum(v*abs(2*j-D) for j,v in enumerate(h))
        curvature = F(5*G-D*V,2)
        old_cost = sum((old_envelope(r) for r in rs), F())
        old_margin = curvature/V-old_cost
        assert row['h_polynomial'] == h and row['link_vertices'] == V and row['G'] == G
        assert F(row['curvature']) == curvature and F(row['normalized_cost']) == old_cost
        assert F(row['normalized_margin']) == old_margin and F(row['absolute_margin']) == V*old_margin
        gain = sum(F(rs.count(r),r+1)*drop for r,drop in DROPS.items())
        for out, additional in ((current,F()),(conditional,F(rs.count(13),28))):
            margin = old_margin+gain+additional
            out.append(dict(ranks=list(rs),link_vertices=V,normalized_cost=str(old_cost-gain-additional),
                            normalized_margin=str(margin),absolute_margin=str(V*margin)))
    actual, hypothetical = summary(current), summary(conditional)
    expected = [([10,10,13,13,15],'-11779/379456'),([9,9,13,15,15],'-271/10240'),
                ([9,13,13,13,13],'-113/5488')]
    assert [(r['ranks'],r['normalized_margin']) for r in actual['negative_rows']] == expected
    assert actual['negative'] == 3 and actual['zero'] == 0
    assert hypothetical['negative'] == hypothetical['zero'] == 0
    assert hypothetical['minimum']['ranks'] == [10,10,11,15,15]
    assert hypothetical['minimum']['normalized_margin'] == '2741/371712'
    assert hypothetical['minimum_absolute']['ranks'] == [61]
    assert hypothetical['minimum_absolute']['absolute_margin'] == '73'
    return dict(complete_profiles=2641,original_saved_rows_compared=2641,
                domain='At most one part1..8, with arbitrary repeated parts>=9; sum61. No high-part concentration.',
                current_rows=current,current_summary=actual,
                conditional_order14=dict(assumed_total='167/2',old_total='84',
                    geometric_coverage_proved=False,additional_finite_coupons_included=False,
                    rows=conditional,summary=hypothetical),
                inherited_concentration_scope='Only low9 moves; costs below9 stay linear. Positive terminal bounds extend by decreasing curvature and vertex product. No absolute lower bound is inferred from negative terminal rows.')


def run():
    docs = load_inputs()
    interface = interfaces(docs)
    finite = finite_check(docs)
    ideal = scalar_check(docs[FRONT])
    assert docs[FRONT]['strict_global_constant'] == '4805'
    paths = [Path(__file__)] + [ROOT/f'{folder}/{s}.{ext}' for s in STEMS
                              for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    paths += [ROOT/'results/fresh_distance14_unrestricted_dimension156_20261017.md',
              ROOT/'results/fresh_cusp_link_shape_curvature_20261017.md']
    return dict(status='PASS_AMBIENT62_THREE_TRANSFER_CURRENT_FRONTIER',
                sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
                ambient_dimension=N,inherited_completed_ceiling=62,new_dimension_bound=None,
                local_interfaces=interface,finite=finite,ideal=ideal,strict_global_constant='4805',
                scope='Three proved local changes are composed under explicit actual-corner priority. Complete exact-tail scalar rows are reconstructed and compared. The finite upper/attainment certificate is inherited with pointwise absorption, not recomputed. Three ideal negatives and finite deficit7 remain; order14 half saving is conditional only.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify',action='store_true')
    ap.add_argument('--check-only',action='store_true')
    args = ap.parse_args()
    assert not (args.verify and args.check_only)
    out = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert out == json.loads(path.read_text())
    elif not args.check_only:
        assert not path.exists(), 'Existing output; use --verify.'
        path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],'finite',out['finite']['exact_maximum'],
          'profiles',out['ideal']['complete_profiles'],'negative',out['ideal']['current_summary']['negative'],
          'conditional',out['ideal']['conditional_order14']['summary']['minimum']['normalized_margin'],flush=True)
