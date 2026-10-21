#!/usr/bin/env python3
"""Exact next frontier; conditional costs are not geometric payment rules.

Standard-library arithmetic only. All current finite tables are reconstructed
and compared with the frozen enlarged-coupon certificate through rank63.
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
TAG = 'fresh_ambient63_payment_frontier_20261018'
COMPLETED = 'fresh_order11_order17_dimension63_20261018'
FINITE = 'fresh_e6_b10_half_coupon_forest_budget_20261018'
IDEAL = 'fresh_ambient64_ideal_margin_closure_20261018'
STEMS = (COMPLETED, FINITE, IDEAL)
N, D = 63, 62
BETA = F(899, 2)


def profiles(left, minimum=1, prefix=(), low=False, high=False):
    if not left:
        yield prefix
        return
    for r in range(minimum, left+1):
        if (r < 9 and low) or (r >= 18 and high):
            continue
        yield from profiles(left-r, r, prefix+(r,), low or r < 9, high or r >= 18)


def ideal_cost(r):
    if r <= 9:
        return F(r, 2)
    return {10: F(109, 22), 11: F(16, 3), 12: F(74, 13), 13: F(6),
            14: F(94, 15), 15: F(107, 16), 16: F(118, 17)}.get(r, F(163, 23))


def scalar(rs):
    h = [1]
    for r in rs:
        out = [0]*(len(h)+r)
        for j, value in enumerate(h):
            for k in range(r+1):
                out[j+k] += value
        h = out
    V = prod(r+1 for r in rs)
    assert len(h) == D+1 and h == h[::-1] and sum(h) == V
    G = sum(x*min(j, D-j) for j, x in enumerate(h))
    # A separate cumulative-distribution identity checks the same G.
    cumulative = 0
    absolute = F()
    for j in range(D):
        cumulative += h[j]
        if j < D//2:
            absolute += F(cumulative, V)
        else:
            absolute += F(V-cumulative, V)
    assert F(G, V) == F(D, 2)-absolute
    curvature = F(5*G-D*V, 2)
    charge = sum(map(ideal_cost, rs), F())
    margin = curvature/V-charge
    return dict(ranks=list(rs), h_polynomial=h, link_vertices=V, G=G,
                curvature=str(curvature), normalized_cost=str(charge),
                normalized_margin=str(margin), absolute_margin=str(V*margin))


def summarize(rows):
    ordered = sorted(rows, key=lambda r: (F(r['normalized_margin']), r['ranks']))
    return dict(count=len(rows), negative=sum(F(r['normalized_margin']) < 0 for r in rows),
                zero=sum(F(r['normalized_margin']) == 0 for r in rows), minimum=ordered[0],
                negative_rows=[r for r in ordered if F(r['normalized_margin']) < 0])


def ideal_scenario(rows, drops):
    out = []
    for row in rows:
        gain = sum(F(row['ranks'].count(r))*drop/(r+1) for r, drop in drops.items())
        out.append(dict(ranks=row['ranks'], normalized_gain=str(gain),
                        normalized_margin=str(F(row['normalized_margin'])+gain)))
    return dict(total_cost_drops_by_affine_rank={str(r): str(v) for r, v in drops.items()},
                summary=summarize(out), rows=out, geometric_coverage_proved=False)


def threshold(rows, varied, fixed):
    bounds, unfixable = [], []
    for row in rows:
        rs = row['ranks']
        margin = F(row['normalized_margin'])+sum(F(rs.count(r))*v/(r+1) for r, v in fixed.items())
        if margin >= 0:
            continue
        if varied not in rs:
            unfixable.append(dict(ranks=rs, margin=str(margin)))
        else:
            bounds.append((-(varied+1)*margin/rs.count(varied), rs))
    bound, witness = max(bounds) if bounds else (F(), [])
    return dict(varied_affine_rank=varied, fixed_total_cost_drops={str(r): str(v) for r, v in fixed.items()},
                threshold_for_nonnegative=str(bound), active_profile=witness,
                unfixable_profiles=unfixable,
                strict_note='A drop strictly above the threshold gives positive margins if no unfixable row remains.')


def representative(r):
    return ('E', 6) if r == 6 else ('A' if r <= 9 else 'D', r)


def component(t, drop15=0, drop16=0, drop17=0):
    kind, r = t
    if r == 15:
        return 210-2*drop15
    if r == 16:
        return 208 if kind == 'B' else 232-2*drop16
    if r == 17:
        return {'A': 228, 'B': 210, 'D': 254-2*drop17}[kind]
    if r <= 14:
        return r*(r-1)
    old = 16*r-(66 if kind == 'D' and r >= 20 else 64 if kind == 'D' else 68)
    deep = max(0, r-(22+{'A': 8, 'B': 9, 'D': 10}[kind]))
    return old-2*deep


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
        if (small, large) in ((('A', 2), ('D', 11)), (('E', 6), ('D', 11))):
            out += 4
        if (small, large) == (('A', 1), ('B', 10)):
            out += 1
        if small in (('A', 1), ('A', 2)) and large == ('D', 14):
            out += small[1]
    return out


def domain():
    rows = []
    for size in range(N//10+1):
        for rs in combinations_with_replacement(range(10, 18), size):
            if sum(rs) > N:
                continue
            c = Counter(rs)
            key = (2*c[10]+3*sum(c[r] for r in range(11, 17))+2*c[17],
                   c[11]+c[12], c[13], c[11], c[14])
            rows.append(dict(counts=[c[r] for r in range(10, 18)], used_rank=sum(rs), key=list(key)))
    return sorted(rows, key=lambda r: r['counts'])


def finite_scenario(frozen, changes, saved=None):
    assert set(changes) <= {'drop15', 'drop16', 'drop17'}
    base = lambda t: component(t, **changes)
    types = ([('A', r) for r in range(1, N+1)]+[('B', r) for r in range(2, N+1)]+
             [('D', r) for r in range(4, N+1)]+[('E', 6), ('E', 7), ('E', 8),
             ('F', 4), ('H', 3), ('H', 4), ('I2*', 2)])
    for t in types:
        rep = t if t == ('A', 17) else representative(t[1])
        assert base(rep) >= base(t)
        assert all(cross(rep, u) >= cross(t, u) for u in types)
    dom = domain()
    assert dom == [r for r in frozen['large_domain'] if r['used_rank'] <= N]
    keys = sorted({tuple(r['key']) for r in dom})
    assert len(dom) == 700 and len(keys) == 328
    if saved is not None:
        assert [r['key'] for r in saved['tables']] == list(map(list, keys))
    tables, items = {}, {}
    counts = Counter()
    for index, key in enumerate(keys):
        broad, extra, a4, d11, d14 = key
        weights = {}
        for r in tuple(range(1, 10))+tuple(range(17, N+1)):
            w = base(representative(r))+(r*broad if r <= 5 else 0)
            w += (2*extra+d14 if r == 1 else 4*d11+2*d14 if r == 2
                  else 2*a4 if r == 4 else 4*d11 if r == 6 else 0)
            weights[r] = w
        items[key] = weights
        if saved is None:
            # Rank-first generation, unlike the frozen item-first producer.
            values, last = [0], [0]
            for n in range(1, N+1):
                value, rank = max((values[n-r]+w, r) for r, w in weights.items() if r <= n)
                values.append(value)
                last.append(rank)
        else:
            values, last = saved['tables'][index]['values'], saved['tables'][index]['last']
        assert len(values) == len(last) == N+1 and values[0] == last[0] == 0
        for n in range(1, N+1):
            for r, w in weights.items():
                if r <= n:
                    assert values[n] >= values[n-r]+w
                    counts['bellman_upper'] += 1
            r = last[n]
            assert r in weights and r <= n and values[n] == values[n-r]+weights[r]
            counts['attaining_predecessors'] += 1
        tables[key] = dict(key=list(key), values=values, last=last)
    if not changes:
        old_tables = {tuple(t['key']): t for t in frozen['tables']}
        assert all(t['values'] == old_tables[k]['values'][:N+1] for k, t in tables.items())
    alternatives = []
    for row in dom:
        large = [('A' if r == 17 else 'D', r) for r, count in zip(range(10, 18), row['counts']) for _ in range(count)]
        key = tuple(row['key'])
        for r, weight in items[key].items():
            assert weight == base(representative(r))+sum(cross(representative(r), t) for t in large)
            counts['literal_interactions'] += 1
        table = tables[key]
        rem = N-row['used_rank']
        value = sum(map(base, large))+table['values'][rem]
        forest = list(large)
        while rem:
            rank = table['last'][rem]
            forest.append(representative(rank))
            rem -= rank
        assert sum(r for k, r in forest) == N
        assert value == sum(map(base, forest))+sum(cross(a, b) for i, a in enumerate(forest) for b in forest[i+1:])
        alternatives.append(dict(doubled_upper=value, counts=row['counts'], forest=sorted(k+str(r) for k, r in forest)))
        counts['global_alternatives'] += 1
    maximum = max(r['doubled_upper'] for r in alternatives)
    best = sorted({tuple(r['forest']) for r in alternatives if r['doubled_upper'] == maximum})
    return dict(conditional_component_cost_drops=changes, exact_maximum=str(F(maximum, 2)),
                curvature=str(BETA), finite_slack=str(BETA-F(maximum, 2)),
                attaining_forests=best, largest_alternatives=sorted(alternatives, key=lambda r: (-r['doubled_upper'], r['forest']))[:20],
                complete_large_vectors=len(dom), complete_keys=len(keys), type_partner_checks=len(types)**2,
                checks=dict(counts), tables=[tables[k] for k in keys], alternatives=alternatives,
                uses_current_proved_rule=not changes, new_geometric_rule_proved=False)


def run(saved=None):
    docs = {tag: json.loads((ROOT/f'results/{tag}.json').read_text()) for tag in STEMS}
    for tag, packet in docs.items():
        assert packet['status'].startswith('PASS'), tag
        for path, digest in packet['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == digest, (tag, path)
    assert docs[COMPLETED]['dimension_bound'] == 63
    assert docs[FINITE]['finite_maximum'] == '464'
    assert all(str(ideal_cost(r['rank'])) == r['cost'] for r in docs[IDEAL]['envelopes'])
    assert BETA == F(5*sum(min(j, D-j) for j in range(D+1))-N*D, 2)
    rows = [scalar(rs) for rs in profiles(D)]
    current = summarize(rows)
    assert len(rows) == len({tuple(r['ranks']) for r in rows}) == 2204
    assert current['negative'] == 21 and current['zero'] == 0
    assert current['minimum']['normalized_margin'] == '-22543/201344'
    assert current['minimum']['ranks'] == [10, 10, 12, 15, 15]
    ideal_changes = ({15: F(1)}, {10: F(1, 2), 11: F(1)},
                     {10: F(1, 2), 11: F(1), 12: F(1, 2)},
                     {10: F(3, 4), 11: F(1)})
    ideal_targets = [ideal_scenario(rows, changes) for changes in ideal_changes]
    assert ideal_targets[2]['summary']['negative'] == ideal_targets[2]['summary']['zero'] == 0
    assert ideal_targets[2]['summary']['minimum']['normalized_margin'] == '3505/201344'
    assert ideal_targets[3]['summary']['minimum']['normalized_margin'] == '157/8788'
    finite_changes = ({}, {'drop16': 2}, {'drop17': 2},
                      {'drop15': 1, 'drop16': 1, 'drop17': 2}, {'drop16': 2, 'drop17': 2})
    finite = [finite_scenario(docs[FINITE], c, None if saved is None else saved['finite_scenarios'][i])
              for i, c in enumerate(finite_changes)]
    assert [r['exact_maximum'] for r in finite] == ['453', '453', '453', '449', '449']
    assert finite[0]['attaining_forests'] == [
        ('D15', 'D15', 'D16', 'D17'), ('D15', 'D16', 'D16', 'D16'), ('D17', 'D17', 'D29')]
    paths = [Path(__file__)]+[ROOT/f'{folder}/{tag}.{ext}' for tag in STEMS
                             for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    paths += [ROOT/'results/fresh_distance14_unrestricted_dimension156_20261017.md',
              ROOT/'results/fresh_cusp_link_shape_curvature_20261017.md']
    return dict(status='PASS_AMBIENT63_PAYMENT_FRONTIER_DIAGNOSTIC',
                sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
                ambient_dimension=N, inherited_completed_ceiling=63,
                finite_scenarios=finite, ideal_current_rows=rows, ideal_current_summary=current,
                ideal_conditional_scenarios=ideal_targets,
                ideal_thresholds=[threshold(rows, 12, {10: F(1, 2), 11: F(1)}),
                                  threshold(rows, 10, {11: F(1)}), threshold(rows, 15, {})],
                low_concentration_threshold=9, high_concentration_threshold=17,
                strict_global_constant=str(5*(N*N//4)), new_dimension_bound=None,
                scope='Exact current finite and ideal obstructions at ambient63, plus explicitly conditional same-component and all-family ideal cost targets. Every new finite scenario keeps all enlarged coupons; no cost reduction is asserted geometrically. No dimension62 theorem or new example.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    path = ROOT/f'results/{TAG}.json'
    saved = json.loads(path.read_text()) if args.verify else None
    out = json.loads(json.dumps(run(saved)))
    if saved is not None:
        assert out == saved
    else:
        path.write_text(json.dumps(out, separators=(',', ':'))+'\n')
    print(out['status'], 'finite', out['finite_scenarios'][0]['exact_maximum'],
          'ideal profiles', out['ideal_current_summary']['count'],
          'negative', out['ideal_current_summary']['negative'], flush=True)
