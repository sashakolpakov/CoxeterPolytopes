#!/usr/bin/env python3
"""Conditional exact finite-forest budget adding full A2/D11 short-leaf coupons.

Standard library only. The combinatorial certificate is conditional on the
proposed extra coupon; it does not establish an ambient dimension bound.
"""
import argparse
import json
from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_a2_d11_full_coupon_forest_budget_20261017'
PREDECESSOR = 'fresh_d10_coupon_large_d_forest_budget_20261017'
LIMIT = 75
INPUTS = (
    PREDECESSOR,
    'fresh_large_finite_d_leaf_half_transfer_20261017',
    'fresh_finite_d17_long_pair_support_20261017',
    'fresh_rank17_reduced_forest_frontier_20261017',
    'fresh_rank17_reduced_forest_frontier_20261017_review',
    'fresh_a11_cycle_coupon_domination_20261017',
    'fresh_finite8_d12_b10_coupon_budget_20261017',
)
RESIDUAL = tuple(range(1, 10)) + tuple(range(17, LIMIT + 1))


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def label(t):
    return t[0] + str(t[1])


def parse(name):
    return name[0], int(name[1:])


def base(t):
    """Doubled same-component upper cost; exact for every representative."""
    kind, r = t
    if r <= 16:
        return r * (r - 1)
    if r == 17:
        return 228 if kind == 'A' else 254
    assert kind in ('A', 'B', 'D')
    return 16 * r - (66 if kind == 'D' and r >= 20 else
                     64 if kind == 'D' else 68)


def cross(a, b):
    """Doubled actual coupon count before any type domination."""
    total = 0
    for small, large in ((a, b), (b, a)):
        if small[1] <= 5:
            if large == ('D', 10):
                total += 2 * small[1]
            elif large[0] == 'D' and 11 <= large[1] <= 16:
                total += 3 * small[1]
            elif large in (('A', 13), ('A', 15)):
                total += small[1]
            elif large == ('A', 17):
                total += 2 * small[1]
        if (small, large) in (
                (('A', 1), ('D', 11)), (('A', 1), ('D', 12)),
                (('A', 2), ('B', 10)), (('A', 4), ('D', 13)),
                (('A', 1), ('A', 11)), (('A', 2), ('A', 11))):
            total += 2
        if (small, large) == (('A', 2), ('D', 11)):
            total += 4
    return total


def representative(t):
    kind, r = t
    if r <= 9:
        return 'A', r
    if t == ('A', 17):
        return t
    return 'D', r


def forest_cost(ts):
    return sum(map(base, ts)) + sum(cross(x, y)
        for i, x in enumerate(ts) for y in ts[i + 1:])


def domination():
    types = ([('A', r) for r in range(1, 76)] +
             [('B', r) for r in range(2, 76)] +
             [('D', r) for r in range(4, 76)] +
             [('E', 6), ('E', 7), ('E', 8), ('F', 4),
              ('H', 3), ('H', 4), ('I2*', 2)])
    rows = []
    for t in types:
        rep = representative(t)
        assert rep in types and rep[1] == t[1]
        inc = base(rep) - base(t)
        assert inc >= 0
        deltas = []
        for other in types:
            delta = cross(rep, other) - cross(t, other)
            assert delta >= 0, (t, other, rep)
            if delta:
                deltas.append([label(other), delta])
        rows.append(dict(original=label(t), representative=label(rep),
                         doubled_base_increase=inc,
                         positive_cross_increases=deltas))
    assert cross(('A', 2), ('B', 10)) == 2
    assert cross(('A', 2), ('D', 10)) == 4
    return dict(types=len(types), partner_comparisons=len(types)**2,
                rows=rows, rank17_B_upper_bound=127,
                I2_star_scope='Every I2(m), m>=4; unchanged all-pairs base and no additional type-specific coupon.')


def count_domain():
    """Component-number enumeration, independent of rank-count recursion."""
    rows = []
    for number in range(LIMIT // 10 + 1):
        for rs in combinations_with_replacement(range(10, 18), number):
            used = sum(rs)
            if used > LIMIT:
                continue
            cs = Counter(rs)
            counts = [cs[r] for r in range(10, 18)]
            broad = 2 * cs[10] + 3 * sum(cs[r] for r in range(11, 17)) + 2 * cs[17]
            key = [broad, cs[11] + cs[12], cs[13], cs[11]]
            ts = [('A' if r == 17 else 'D', r) for r in rs]
            # No coupon joins two members of this large-component list.
            assert forest_cost(ts) == sum(map(base, ts))
            rows.append(dict(component_counts=counts, used_rank=used,
                             interaction_key=key,
                             doubled_base_cost=sum(map(base, ts))))
    assert len(rows) == len({tuple(r['component_counts']) for r in rows}) == 1608
    return rows


def residual_costs(key):
    broad, upgrade, a4, d11 = key
    return {r: base(('A' if r <= 9 else 'D', r)) +
            (r * broad + (2 * upgrade if r == 1 else
                          4 * d11 if r == 2 else 2 * a4 if r == 4 else 0) if r <= 5 else 0)
            for r in RESIDUAL}


def make_table(key):
    costs = residual_costs(key)
    # Item-first unbounded knapsack; the independent verifier below checks
    # every last-item recurrence and an attaining predecessor instead.
    values = [0] + [None] * LIMIT
    for r in RESIDUAL:
        for n in range(r, LIMIT + 1):
            if values[n-r] is not None:
                value = values[n-r] + costs[r]
                if values[n] is None or value > values[n]:
                    values[n] = value
    assert all(v is not None for v in values)
    last = [0] + [max(r for r in RESIDUAL if r <= n and
                     values[n] == values[n-r] + costs[r])
                  for n in range(1, LIMIT + 1)]
    return dict(interaction_key=list(key), doubled_optima=values,
                last_rank=last)


def unpack(row, rank):
    result = []
    while rank:
        r = row['last_rank'][rank]
        assert 0 < r <= rank
        result.append(('A' if r <= 9 else 'D', r))
        rank -= r
    return result


def canonical(ts):
    return sorted(ts, key=lambda t: (t[1], t[0]))


def curvature(n):
    return F(n*(n-6), 8) if n % 2 == 0 else F((n-1)*(n-5), 8)


def generate():
    domain = count_domain()
    keys = sorted({tuple(r['interaction_key']) for r in domain})
    assert len(keys) == 357
    tables = {key: make_table(key) for key in keys}
    answers = []
    for n in range(LIMIT + 1):
        best = None
        for row in domain:
            if row['used_rank'] > n:
                continue
            key = tuple(row['interaction_key'])
            value = row['doubled_base_cost'] + tables[key]['doubled_optima'][n-row['used_rank']]
            ts = [('A' if r == 17 else 'D', r)
                  for r, c in zip(range(10, 18), row['component_counts']) for _ in range(c)]
            ts = canonical(ts + unpack(tables[key], n-row['used_rank']))
            candidate = value, tuple(map(label, ts))
            if best is None or candidate > best:
                best = candidate
        value, ws = best
        answers.append(dict(dimension=n, doubled_cost=value,
            witness_components=list(ws), finite_curvature_coefficient=str(curvature(n)),
            finite_slack=str(curvature(n)-F(value, 2))))
    return dict(large_count_vectors=domain, residual_tables=[tables[k] for k in keys],
                dimension_maxima=answers)


def direct_D_counts():
    rows = []
    for n in range(18, LIMIT + 1):
        g = [set() for _ in range(n)]
        for a, b in list(zip(range(n-2), range(1, n-1))) + [(n-3, n-1)]:
            g[a].add(b); g[b].add(a)
        short_leaves = {n-2, n-1}
        short = leaf9 = half9 = 0
        for i in range(n):
            dist = {i: 0}; queue = deque([i])
            while queue:
                a = queue.popleft()
                for b in g[a]:
                    if b not in dist:
                        dist[b] = dist[a] + 1; queue.append(b)
            assert len(dist) == n
            for j in range(i+1, n):
                if dist[j] <= 8:
                    short += 1
                elif dist[j] == 9 and (len(g[i]) == 1 or len(g[j]) == 1):
                    leaf9 += 1
                    if n >= 20 and ({i, j} & short_leaves):
                        half9 += 1
        cost = 2*(short+leaf9)-half9
        assert short == 8*n-35 and leaf9 == 3
        assert half9 == (2 if n >= 20 else 0)
        assert cost == base(('D', n))
        rows.append(dict(rank=n, short_pairs=short, leaf9_pairs=leaf9,
                         halved_pairs=half9, doubled_cost=cost))
    return rows


def verify_certificate(cert):
    domain = count_domain()
    assert cert['large_count_vectors'] == domain
    old = json.loads((ROOT/f'results/{PREDECESSOR}.json').read_text())
    assert old['status'] == 'PASS_EXACT_D10_COUPON_LARGE_D_FINITE_BUDGET'
    for a, b in zip(domain, old['certificate']['large_count_vectors']):
        assert a['component_counts'] == b['component_counts']
        assert a['used_rank'] == b['used_rank']
        assert a['doubled_base_cost'] == b['doubled_base_cost']
        assert a['interaction_key'][:3] == b['interaction_key']
    tables = {tuple(r['interaction_key']): r for r in cert['residual_tables']}
    assert len(tables) == len(cert['residual_tables']) == 357
    assert set(tables) == {tuple(r['interaction_key']) for r in domain}
    # Check the compressed formula against actual type-pair coupons for every
    # large multiset and every possible residual component, not a dummy rank.
    interactions = 0
    for row in domain:
        large = [('A' if r == 17 else 'D', r)
                 for r, c in zip(range(10, 18), row['component_counts']) for _ in range(c)]
        costs = residual_costs(tuple(row['interaction_key']))
        for r in RESIDUAL:
            small = ('A' if r <= 9 else 'D', r)
            assert costs[r] == base(small) + sum(cross(small, t) for t in large)
            interactions += 1
    alternatives = 0
    for key, row in tables.items():
        values = row['doubled_optima']; last = row['last_rank']
        assert len(values) == len(last) == LIMIT+1 and values[0] == last[0] == 0
        costs = residual_costs(key)
        for n in range(1, LIMIT+1):
            assert isinstance(values[n], int)
            for r in RESIDUAL:
                if r <= n:
                    assert values[n] >= values[n-r] + costs[r]
                    alternatives += 1
            r = last[n]
            assert r in RESIDUAL and r <= n
            assert values[n] == values[n-r] + costs[r]
    global_checks = 0
    assert len(cert['dimension_maxima']) == LIMIT+1
    for n, row in enumerate(cert['dimension_maxima']):
        assert row['dimension'] == n
        value = row['doubled_cost']
        for large in domain:
            if large['used_rank'] <= n:
                table = tables[tuple(large['interaction_key'])]
                assert value >= large['doubled_base_cost'] + table['doubled_optima'][n-large['used_rank']]
                global_checks += 1
        ts = list(map(parse, row['witness_components']))
        assert sum(r for kind, r in ts) == n
        assert all(t == representative(t) for t in ts)
        assert forest_cost(ts) == value
        assert F(row['finite_curvature_coefficient']) == curvature(n)
        assert F(row['finite_slack']) == curvature(n)-F(value, 2)
        assert value == old['certificate']['dimension_maxima'][n]['doubled_cost']
    assert cert['dimension_maxima'][66]['doubled_cost'] == 990
    assert cert['dimension_maxima'][66]['witness_components'] == ['D66']
    return dict(large_count_vectors=len(domain), compressed_interaction_tables=len(tables),
        interaction_formula_checks=interactions, bellman_upper_alternatives=alternatives,
        attaining_bellman_predecessors=len(tables)*LIMIT,
        global_upper_alternatives=global_checks, typed_global_attainers=LIMIT+1)


def source_bindings():
    paths = [Path(__file__)]
    for stem in INPUTS:
        data = json.loads((ROOT/f'results/{stem}.json').read_text())
        assert data['status'].startswith('PASS'), stem
        for name, expected in data['sources'].items():
            assert digest(ROOT/name) == expected, (stem, name)
        paths.extend(ROOT/f'{folder}/{stem}.{ext}'
                     for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md')))
    return {str(p.relative_to(ROOT)): digest(p) for p in paths}


def run(verify=False):
    path = ROOT/f'results/{TAG}.json'
    if verify:
        saved = json.loads(path.read_text())
        cert = saved['certificate']
    else:
        cert = generate()
    counts = verify_certificate(cert)
    out = dict(status='PASS_CONDITIONAL_A2_D11_FULL_COUPON_FINITE_BUDGET',
        sources=source_bindings(), maximum_rank=LIMIT,
        conditional_extra_coupon='Upgrade each existing half between an A2 mark and a short fork leaf of whole D11 to a unit; total increment2 per A2/D11 component pair.',
        missing_geometric_obligation='Geometric face-payment proof for the full A2/D11 upgrade and any ideal-pair reductions using it; no such proof is established here.',
        unchanged_maxima_at_every_rank_through=LIMIT,
        large_types=['D'+str(r) for r in range(10, 17)]+['A17'],
        residual_ranks=list(RESIDUAL), interaction_key_columns=['2D10+3D11to16+2A17', 'D11+D12', 'D13', 'D11'],
        doubled_representative_cross_formula='a*(2*d10+3*tD11to16+2*a17)+2*a1*(d11+d12)+2*a4*d13+4*a2*d11',
        type_domination=domination(), direct_D_graph_counts=direct_D_counts(),
        certificate=cert, counts=counts, new_dimension_bound=None,
        scope='No project or producer arithmetic imports. Complete exact finite-forest maxima through rank75 conditional on the proposed full A2/D11 short-fork-leaf upgrade, whose geometric proof is not supplied or assumed completed here. All other payment rules are frozen bound inputs. Item-first Bellman generation is certified by all last-item upper inequalities and attaining predecessors; typed global attainers and direct D18..D75 graph counts are checked. Does not recompute ideal budgets, prove geometric face coverage, or assert an ambient dimension theorem.')
    out = json.loads(json.dumps(out))
    if verify:
        assert saved == out
    else:
        path.write_text(json.dumps(out, sort_keys=True, separators=(',', ':'))+'\n')
    print(out['status'], counts, flush=True)
    print([(r['dimension'],str(F(r['doubled_cost'],2)),r['finite_slack'],r['witness_components'])
           for r in cert['dimension_maxima'] if r['dimension'] >= 64], flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    args = p.parse_args()
    run(args.verify)
