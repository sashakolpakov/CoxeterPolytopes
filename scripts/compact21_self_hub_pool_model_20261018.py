#!/usr/bin/env python3
"""Exact complete compact21 model with same-hub cross-angle charge pools.

No optimizer. Full native cases are deduplicated only after new coefficients.
The pool cost belongs once to each eligible full-forest budget, never to a
Bellman item. The old D5/H4 exceptional forests have no eligible self hubs.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
import fresh_compact21_full_forest_epigraph_20261018 as E

C, D, P, ROOT, Z, T = E.C, E.D, E.P, E.ROOT, E.Z, E.T
NAMES = ('A1', 'I2(5)', 'H3')


def context_key(name, counts):
    # New pools may use the finer exact tuple; old item weights stay capped.
    return (name,)+tuple(counts)


def denominator(name, counts):
    n = counts[NAMES.index(name)]
    assert n >= 2
    rank = {'A1': 1, 'I2(5)': 2, 'H3': 3}[name]
    return rank*rank*n*(n-1)//2


def selected_columns(counts, pools):
    for k, name in enumerate(NAMES):
        if counts[k] >= 2:
            yield pools[context_key(name, counts)]


def add_pools(row, counts, pools):
    co = Counter(row)
    for j in selected_columns(counts, pools):
        co[j] += 1
    return dict(co)


def forest_row(rule, profiles, forest, pools, state=None):
    row, rhs = E.forest_row(rule, profiles, forest, state)
    return add_pools(row, D.exact_context(forest), pools), rhs


def self_families(row, edges):
    result = []
    for a, b in product((0, 2), (1, 3)):
        pieces = P.G.pieces(P.G.graph(row['white']+[a, b], edges))
        left = next(g for g in pieces if a in g)
        right = next(g for g in pieces if b in g)
        lf, rf = P.G.family(left), P.G.family(right)
        result.append(lf if b not in left and lf == rf and lf in NAMES else None)
    return result


def build():
    rule = D.Rule(); models, profiles = D.profiles(rule)
    assert E.check_existing_rows(rule, profiles) == 510
    grouped, auxiliary, next_column = E.groups(), {}, T+1
    for context, tuples in sorted(grouped.items()):
        for rank in range(2, max(t[-1] for t in tuples)+1):
            auxiliary[context, rank] = next_column; next_column += 1
    assert next_column == 94706
    keys = {context_key(name, t[:4]) for tuples in grouped.values() for t in tuples
            for k, name in enumerate(NAMES) if t[k] >= 2}
    pools = {key: next_column+i for i, key in enumerate(sorted(keys))}
    next_column += len(pools)
    unique, old_unique, face_count, cases = set(), set(), 0, 0
    for _, _, row in P.native_faces():
        face_count += 1
        edges = dict(zip(combinations(range(row['nodes']), 2), row['labels']))
        corners = [P.corner(row, edges, a, b) for a, b in product((0, 2), (1, 3))]
        families = self_families(row, edges)
        baseline = sum((c[0] for c in corners), F())
        for outside in P.spectators(23-row['nodes']):
            cases += 1; co = Counter()
            for family, (_, local, key, old_denominator, _) in zip(families, corners):
                counts = tuple(x+y for x, y in zip(local, outside))
                if key is not None:
                    factor = F(1)
                    if old_denominator is not None:
                        k, multiplier, offset = old_denominator
                        divisor = multiplier*(counts[k]+offset)
                        assert divisor > 0
                        factor /= divisor
                    co[rule.column(key, counts)] -= factor
                if family is not None:
                    co[pools[context_key(family, counts)]] -= F(1, denominator(family, counts))
            unique.add((tuple(sorted(co.items())), baseline-1))
            old_unique.add((tuple(sorted((j, a) for j, a in co.items() if j < Z)), baseline-1))
    assert face_count == 1979 and cases == 102917
    assert len(old_unique) == 51443
    del old_unique
    base = [(dict(co), rhs, 'face') for co, rhs in sorted(unique)]
    unique_count = len(base); del unique
    metadata = C.H.D.forest_metadata()
    checked = 0
    for i, (literal, rhs, kind) in enumerate(C.H.rows()):
        if kind == 'face':
            continue
        if kind == 'cap':
            assert literal == {T: F(1)} and rhs == 1
            base.append((literal, rhs, kind)); continue
        assert kind == 'forest'
        if i in metadata:
            forest = tuple(sorted(metadata[i]['types'])); state = metadata[i].get('state')
        else:
            name = {2805: 'A5', 2806: 'B4', 2807: 'B5'}[i]
            forest, state = E.SPECIAL[name], None
        row, expected = forest_row(rule, profiles, forest, pools, state)
        old = literal | {T: F(1)} if i == 2805 else literal
        assert row == add_pools(old, D.exact_context(forest), pools) and expected == rhs
        base.append((row, rhs, kind)); checked += 1
    assert checked == 510 and len(base) == unique_count+511
    coverage = {'native_faces': face_count, 'native_spectator_cases': cases,
                'exact_distinct_enlarged_face_rows': unique_count,
                'zero_pool_distinct_face_rows': 51443, 'all_unequal_new_rows_retained': True,
                'existing_forest_rows': checked, 'transfer_cap_rows': 1,
                'self_hub_pool_variables': len(pools), 'A5_source_capacity': 1,
                'pool_selector': 'exact_A1_I25_H3_A4_counts'}
    return base, coverage, rule, models, profiles, grouped, auxiliary, pools, next_column


def extended_rows(base, rule, profiles, grouped, auxiliary, pools):
    yield from base
    for context, tuples in sorted(grouped.items()):
        rows = E.context_rows(rule, profiles, context)
        maximum_rank = max(t[-1] for t in tuples)
        items = [name for name in profiles if name not in E.HUBS]
        for rank in range(2, maximum_rank+1):
            for name in items:
                remaining = rank-profiles[name][0]
                if remaining < 0 or remaining == 1:
                    continue
                row, baseline = rows[name]
                co = Counter(row); co[auxiliary[context, rank]] -= 1
                if remaining:
                    co[auxiliary[context, remaining]] += 1
                yield dict(co), -baseline, 'bellman'
        for a, b, h, d, rank in tuples:
            counts = a, b, h, d
            if counts in ((1, 0, 0, 0), (0, 0, 0, 0)):
                target = {'H4', 'E8'} if a else {'D5', 'E8'}
                fixed, baseline = E.combine(rows, zip(E.HUBS, counts))
                for name in items:
                    remaining = rank-profiles[name][0]
                    if name in target or remaining < 0 or remaining == 1:
                        continue
                    co = Counter(fixed); co.update(rows[name][0]); co[Z] -= 1
                    if remaining:
                        co[auxiliary[context, remaining]] += 1
                    yield add_pools(co, counts, pools), -baseline-rows[name][1], 'off_target_budget'
            else:
                co, baseline = E.combine(rows, zip(E.HUBS, counts)); co[Z] -= 1
                if rank:
                    co[auxiliary[context, rank]] += 1
                yield add_pools(co, counts, pools), -baseline, 'ordinary_budget'
    for forest in (tuple(sorted(('A1',)+('H4',)*5)), tuple(sorted(('A1',)+('H4',)*3+('E8',)))):
        row, rhs = forest_row(rule, profiles, forest, pools)
        yield row, rhs, 'pure_unchanged_budget'
    for name, forest in E.SPECIAL.items():
        assert not list(selected_columns(D.exact_context(forest), pools))
        for state in (('active', 'inactive') if name == 'D5' else (None,)):
            row, rhs = forest_row(rule, profiles, forest, pools, state)
            yield row, rhs, 'special_budget'


def exact_forest_maximum(rule, models, profiles, grouped, pools, pool_values):
    best, attainer = F(-1), None
    def offer(value, forest, kind):
        nonlocal best, attainer
        addition = sum((pool_values.get(j, F()) for j in selected_columns(D.exact_context(forest), pools)), F())
        value += addition
        if value > best:
            best, attainer = value, {'types': list(sorted(forest)), 'budget_kind': kind,
                                     'self_hub_pool_addition': str(addition)}
    for context, tuples in sorted(grouped.items()):
        values = D.component_values(rule, profiles, context)
        states = {0: (F(), ())}
        for rank in range(2, max(t[-1] for t in tuples)+1):
            choices = [(states[rank-profile[0]][0]+values[name],
                        tuple(sorted(states[rank-profile[0]][1]+(name,))))
                       for name, profile in profiles.items()
                       if name not in E.HUBS and rank-profile[0] in states]
            assert choices; states[rank] = max(choices)
        for a, b, h, d, rank in tuples:
            fixed = sum((n*values[name] for name, n in zip(E.HUBS, (a, b, h, d))), F())
            hubs = ('A1',)*a+('I2(5)',)*b+('H3',)*h+('A4',)*d
            if (a, b, h, d) in ((1, 0, 0, 0), (0, 0, 0, 0)):
                target = {'H4', 'E8'} if a else {'D5', 'E8'}
                for name, profile in profiles.items():
                    remaining = rank-profile[0]
                    if name in E.HUBS or name in target or remaining not in states:
                        continue
                    value, forest = states[remaining]
                    offer(fixed+values[name]+value, hubs+forest+(name,), 'off_target')
            else:
                value, forest = states[rank]
                offer(fixed+value, hubs+forest, 'ordinary')
    extra = [tuple(sorted(('A1',)+('H4',)*5)), tuple(sorted(('A1',)+('H4',)*3+('E8',)))]
    extra += list(E.SPECIAL.values())
    for forest in extra:
        context = D.exact_context(forest); values = D.component_values(rule, profiles, context)
        value = sum((values[name] for name in forest), F())
        if forest == E.SPECIAL['H4']:
            value -= rule.transfer
        elif forest == E.SPECIAL['D5']:
            value -= min(F(1), 2*rule.bonus(E.beta_key(models), context))
        elif forest in (E.SPECIAL['A4'], E.SPECIAL['A5'], E.SPECIAL['B4'], E.SPECIAL['B5']):
            value += rule.transfer
        offer(value, forest, 'explicit')
    return best, attainer
