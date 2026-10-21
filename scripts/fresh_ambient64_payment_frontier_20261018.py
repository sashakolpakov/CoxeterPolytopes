#!/usr/bin/env python3
"""Exact next frontier under the completed dimension64 payment rule.

Current finite/ideal obstructions are exact. All cost-reduction scenarios
are conditional arithmetic diagnostics, never new local payment theorems.
No project arithmetic implementation is imported.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from math import prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_ambient64_payment_frontier_20261018'
COMPLETED = 'fresh_deep_finite8_dimension64_20261018'
FINITE = 'fresh_d16_d17_d14_coupon_forest_budget_20261017'
STRIP = 'fresh_depth17_rank65_conditional_budget_20261018'
IDEAL = 'fresh_ambient65_ideal_margin_closure_20261017'
CONVEX = 'fresh_distance14_unrestricted_dimension156_20261017'
STEMS = (COMPLETED, FINITE, STRIP, IDEAL)
N, D = 64, 63


def profiles(left, threshold, minimum=1, prefix=(), low=False, high=False):
    if not left:
        yield prefix
        return
    for r in range(minimum, left + 1):
        if (r < threshold and low) or (r >= 18 and high):
            continue
        yield from profiles(left-r, threshold, r, prefix+(r,),
                            low or r < threshold, high or r >= 18)


def cost(r):
    if r <= 10:
        return F(r, 2)
    return {11: F(16, 3), 12: F(74, 13), 13: F(6),
            14: F(94, 15), 15: F(107, 16), 16: F(7)}.get(r, F(163, 23))


def scalar(rs):
    h = [1]
    for r in rs:
        new = [0] * (len(h)+r)
        for i, x in enumerate(h):
            for j in range(r+1):
                new[i+j] += x
        h = new
    V = prod(r+1 for r in rs)
    assert sum(rs) == D and sum(h) == V and h == h[::-1]
    G = sum(x*min(i, D-i) for i, x in enumerate(h))
    curvature = F(5*G-D*V, 2)
    charge = sum(map(cost, rs), F())
    margin = curvature/V-charge
    return dict(ranks=list(rs), h_polynomial=h, link_vertices=V, G=G,
                normalized_cost=str(charge), curvature=str(curvature),
                normalized_margin=str(margin), absolute_margin=str(V*margin))


def summary(rows, field='normalized_margin'):
    ordered = sorted(rows, key=lambda r: (F(r[field]), r['ranks']))
    return dict(count=len(rows), negative=sum(F(r[field]) < 0 for r in rows),
                zero=sum(F(r[field]) == 0 for r in rows), minimum=ordered[0],
                negative_rows=[r for r in ordered if F(r[field]) < 0])


def scenario(rows, drops):
    out = []
    for row in rows:
        rs = row['ranks']
        gain = sum(F(rs.count(rank))*drop/(rank+1)
                   for rank, drop in drops.items())
        out.append(dict(ranks=rs, normalized_gain=str(gain),
                        normalized_margin=str(F(row['normalized_margin'])+gain)))
    return dict(total_pair_cost_drops_by_affine_rank={str(k): str(v) for k, v in drops.items()},
                summary=summary(out), rows=out, geometric_coverage_proved=False)


def one_parameter_threshold(rows, varied, fixed):
    constraints, unfixable = [], []
    for row in rows:
        rs = row['ranks']
        margin = F(row['normalized_margin']) + sum(
            F(rs.count(rank))*drop/(rank+1) for rank, drop in fixed.items())
        if margin >= 0:
            continue
        if not rs.count(varied):
            unfixable.append(dict(ranks=rs, margin=str(margin)))
        else:
            constraints.append((-(varied+1)*margin/rs.count(varied), rs))
    bound, witness = max(constraints)
    return dict(varied_affine_rank=varied,
                fixed_total_pair_cost_drops={str(k): str(v) for k, v in fixed.items()},
                required_drop_for_nonnegative_margins=str(bound),
                active_profile=witness, unfixable_profiles=unfixable,
                strict_note='A strict increase above this threshold makes all margins positive if no unfixable profile remains.')


def current_finite(old, strip):
    full, small = [r['dimension_maxima'] for r in old['results']]
    assert small[N]['doubled_cost'] == 944
    assert small[N]['forest'] == ['D16']*4
    costs = {r['rank']: r for r in strip['type_domination']}
    bounds = []
    for r in range(18, N+1):
        c = costs[r]['new_doubled_costs']
        assert c['D'] >= c['A'] and c['D'] >= c['B']
        assert costs[r]['no_cross_coupons'] is True
        val = F(c['D']+full[N-r]['doubled_cost'], 2)
        bounds.append(dict(large_rank=r, doubled_component_cost=c['D'],
                           inherited_remainder=full[N-r], upper_bound=str(val)))
    maximum = max(F(r['upper_bound']) for r in bounds)
    assert maximum == 469 < 472
    return dict(exact_maximum='472', attaining_forest=['D16']*4,
                curvature='464', deficit='8',
                all_small_maximum=small[N], large_component_bounds=bounds,
                large_component_upper='469', large_component_attainer=['D64'],
                proof='All-small rank64 bound472 is inherited and attained unchanged. Strip any component of rank>=18: no cross coupon involves it, D dominates, and the remaining modified forest is bounded by the old unrestricted table. All47 bounds are<=469, attained by wholeD64.')


# The conditional finite optimizer uses the unchanged five-coordinate coupon
# environment of the frozen finite theorem. Only stated same-component costs
# vary; this does not establish their geometric availability.
def representative(r):
    return ('E', 6) if r == 6 else ('A' if r <= 9 else 'D', r)


def base(t, deep_zero, drop16, drop17, depth=17):
    kind, r = t
    if r == 16:
        return {'A': 232, 'B': 208, 'D': 236-2*drop16}[kind]
    if r == 17:
        return {'A': 228, 'B': 210, 'D': 250-2*drop17}[kind]
    if r <= 15:
        return r*(r-1)
    old = 16*r-(66 if kind == 'D' and r >= 20 else 64 if kind == 'D' else 68)
    changed = max(0, r-(2*depth+{'A': 8, 'B': 9, 'D': 10}[kind]))
    return old-(2 if deep_zero else 1)*changed


def cross(a, b):
    out = 0
    for small, large in ((a, b), (b, a)):
        if small[1] <= 5:
            if large == ('D', 10):
                out += 2*small[1]
            elif large[0] == 'D' and 11 <= large[1] <= 16:
                out += 3*small[1]
            elif large in (('A', 13), ('A', 15)):
                out += small[1]
            elif large == ('A', 17):
                out += 2*small[1]
        if (small, large) in ((('A', 1), ('D', 11)), (('A', 1), ('D', 12)),
                             (('A', 2), ('B', 10)), (('A', 4), ('D', 13)),
                             (('A', 1), ('A', 11)), (('A', 2), ('A', 11))):
            out += 2
        if (small, large) == (('A', 2), ('D', 11)):
            out += 4
        if (small, large) == (('E', 6), ('D', 11)):
            out += 2
        if small in (('A', 1), ('A', 2)) and large == ('D', 14):
            out += small[1]
    return out


def finite_scenario(old, deep_zero=False, drop16=0, drop17=0, depth=17):
    assert 0 <= drop16 <= 2 and -2 <= drop17 <= 1 and depth in (11, 17)
    def component(t):
        return base(t, deep_zero, drop16, drop17, depth)
    # In this range current representatives still dominate A/B16 and A/B17.
    # All other domination gains remain unchanged, except for A/B/D>=18;
    # their explicit new D dominance is checked below for every relevant rank.
    for r in range(18, N+1):
        assert component(('D', r)) >= max(component((k, r)) for k in ('A', 'B'))
    types = ([('A', r) for r in range(1, N+1)]+
             [('B', r) for r in range(2, N+1)]+
             [('D', r) for r in range(4, N+1)]+
             [('E', 6), ('E', 7), ('E', 8), ('F', 4), ('H', 3), ('H', 4), ('I2*', 2)])
    for t in types:
        rep = t if t == ('A', 17) else representative(t[1])
        assert component(rep) >= component(t)
        for other in types:
            assert cross(rep, other) >= cross(t, other)
    saved_domain = [r for r in old['large_domain'] if r['used_rank'] <= N]
    domain = []
    for size in range(N//10+1):
        for rs in combinations_with_replacement(range(10, 18), size):
            if sum(rs) > N:
                continue
            c = Counter(rs)
            counts = [c[r] for r in range(10, 18)]
            key = (2*c[10]+3*sum(c[r] for r in range(11, 17))+2*c[17],
                   c[11]+c[12], c[13], c[11], c[14])
            domain.append(dict(component_counts=counts, used_rank=sum(rs), interaction_key=list(key)))
    assert sorted(domain, key=lambda r: r['component_counts']) == sorted(
        [{k: r[k] for k in ('component_counts', 'used_rank', 'interaction_key')}
         for r in saved_domain], key=lambda r: r['component_counts'])
    keys = sorted({tuple(r['interaction_key']) for r in domain})
    tables = {}
    checks = dict(bellman_upper=0, attaining_predecessors=0, literal_interactions=0, global_alternatives=0)
    for key in keys:
        broad, upgrade, a4, d11, d14 = key
        weights = {}
        for r in tuple(range(1, 10))+tuple(range(17, N+1)):
            w = component(representative(r))
            if r <= 5:
                w += r*broad
            if r == 1:
                w += 2*upgrade+d14
            if r == 2:
                w += 4*d11+2*d14
            if r == 4:
                w += 2*a4
            if r == 6:
                w += 2*d11
            weights[r] = w
        v = [0]+[None]*N
        for r, w in weights.items():
            for n in range(r, N+1):
                if v[n-r] is not None:
                    candidate = v[n-r]+w
                    if v[n] is None or candidate > v[n]:
                        v[n] = candidate
        last = [0]
        for n in range(1, N+1):
            for r, w in weights.items():
                if r <= n:
                    assert v[n] >= v[n-r]+w
                    checks['bellman_upper'] += 1
            last.append(max(r for r, w in weights.items() if r <= n and v[n] == v[n-r]+w))
            checks['attaining_predecessors'] += 1
        tables[key] = dict(key=list(key), values=v, last=last)
    alternatives = []
    for row in domain:
        large = [('A' if r == 17 else 'D', r) for r, count in zip(range(10, 18), row['component_counts']) for _ in range(count)]
        key = tuple(row['interaction_key'])
        t = tables[key]
        for r in tuple(range(1, 10))+tuple(range(17, N+1)):
            literal = component(representative(r))+sum(cross(representative(r), u) for u in large)
            broad, upgrade, a4, d11, d14 = key
            expected = component(representative(r))
            expected += r*broad if r <= 5 else 0
            expected += 2*upgrade+d14 if r == 1 else 4*d11+2*d14 if r == 2 else 2*a4 if r == 4 else 2*d11 if r == 6 else 0
            assert literal == expected
            checks['literal_interactions'] += 1
        remaining = N-row['used_rank']
        value = sum(component(u) for u in large)+t['values'][remaining]
        forest = list(large)
        while remaining:
            r = t['last'][remaining]
            forest.append(representative(r))
            remaining -= r
        assert sum(r for _, r in forest) == N
        assert value == sum(component(u) for u in forest)+sum(cross(a, b) for i, a in enumerate(forest) for b in forest[i+1:])
        alternatives.append((value, sorted(k+str(r) for k, r in forest)))
        checks['global_alternatives'] += 1
    maximum, witness = max(alternatives)
    return dict(feature_depth=depth, deep_pair_weight='0' if deep_zero else '1/2', D16_drop=drop16,
                D17_drop=drop17, exact_maximum=str(F(maximum, 2)), attaining_forest=witness,
                finite_slack=str(F(464)-F(maximum, 2)), geometric_coverage_proved_here=False,
                uses_current_proved_rule=(depth == 17 and not deep_zero and drop16 == drop17 == 0),
                type_domination_checks=len(types)**2,
                large_vectors=len(domain), interaction_keys=len(keys), checks=checks,
                tables=[tables[k] for k in keys])


def run():
    data = {s: json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s, packet in data.items():
        assert packet['status'].startswith('PASS'), s
        for p, digest in packet['sources'].items():
            assert sha256((ROOT/p).read_bytes()).hexdigest() == digest, (s, p)
    assert data[COMPLETED]['dimension_bound'] == 64
    finite = current_finite(data[FINITE], data[STRIP])
    rows = [scalar(rs) for rs in profiles(D, 10)]
    current = summary(rows)
    assert current['count'] == 1513 and current['negative'] == 24 and current['zero'] == 0
    assert current['minimum']['normalized_margin'] == '-62665/419628'
    assert current['minimum']['ranks'] == [10, 10, 11, 16, 16]
    expanded = [scalar(rs) for rs in profiles(D, 9)]
    assert len(expanded) == 2419
    old_map = {tuple(r['ranks']): r for r in rows}
    expanded_map = {tuple(r['ranks']): r for r in expanded}
    assert all(expanded_map[k] == v for k, v in old_map.items())
    scenarios = [scenario(expanded, drops) for drops in (
        {10: F(1)}, {10: F(1, 2), 16: F(1, 2)},
        {10: F(1, 2), 16: F(1)}, {10: F(1), 11: F(1)},
        {10: F(1, 2), 11: F(1)})]
    assert scenarios[0]['summary']['negative'] == 1
    assert scenarios[1]['summary']['negative'] == scenarios[1]['summary']['zero'] == 0
    assert scenarios[1]['summary']['minimum']['normalized_margin'] == '167/419628'
    thresholds = [one_parameter_threshold(expanded, 16, {10: F(1, 2)}),
                  one_parameter_threshold(expanded, 10, {16: F(1)}),
                  one_parameter_threshold(expanded, 10, {11: F(1)})]
    assert thresholds[0]['required_drop_for_nonnegative_margins'] == '24517/49368'
    assert all(not r['unfixable_profiles'] for r in thresholds)
    finite_scenarios = [finite_scenario(data[FINITE], **kw) for kw in (
        {}, dict(drop16=2), dict(deep_zero=True, drop16=2),
        dict(deep_zero=True, drop16=2, drop17=1),
        dict(deep_zero=True, drop16=2, drop17=-2, depth=11))]
    assert [r['exact_maximum'] for r in finite_scenarios[:4]] == ['472', '469', '465', '464']
    assert finite_scenarios[-1]['exact_maximum'] == '464'
    obstruction = expanded_map[(3, 10, 10, 10, 10, 10, 10)]
    assert obstruction['normalized_margin'] == '-1927/644204'
    paths = [Path(__file__)] + [ROOT/f'{folder}/{s}.{ext}' for s in STEMS
             for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    paths.append(ROOT/f'results/{CONVEX}.md')
    return dict(status='PASS_AMBIENT64_COMPLETE_PAYMENT_FRONTIER_DIAGNOSTIC',
                sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
                ambient_dimension=N, inherited_dimension_ceiling=64,
                finite_current=finite, ideal_current_summary=current, ideal_current_rows=rows,
                expanded_ideal_domain_rows=expanded, ideal_conditional_scenarios=scenarios,
                ideal_thresholds=thresholds, finite_conditional_scenarios=finite_scenarios,
                low_rank_obstruction=obstruction,
                strengthened_global_constant='5120',
                hypothetical_joint_conclusion='Finite costs<=464 plus positive ideal margins would exclude64, using the existing higher tail; neither new local cost reduction is proved here.',
                scope='Complete current scalar frontier and exact conditional cost diagnostics. Inherited finite type/coupon domination and geometric payment coverage remain bound inputs. No new local transfer, changed payment theorem, global bound or finite-volume example.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    print(result['status'], 'finite', result['finite_current']['exact_maximum'],
          'ideal profiles', result['ideal_current_summary']['count'],
          'negative', result['ideal_current_summary']['negative'], flush=True)
