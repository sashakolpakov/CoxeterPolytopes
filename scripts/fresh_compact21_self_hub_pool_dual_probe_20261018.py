#!/usr/bin/env python3
"""Exact new-column diagnostic for three previously unweighted self-hub families.

No optimization. Every representation of each old duplicate face row is
included when lower-bounding its demand on a new normalized angle pool.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import file_digest
from itertools import combinations, product
import json
import fresh_compact21_dual_face_consumers_20261018 as W

R, C, D, P, ROOT = W.R, W.C, W.D, W.P, W.ROOT
TAG = 'fresh_compact21_self_hub_pool_dual_probe_20261018'
GUARD_HELPER = f'scripts/{W.TAG}.py'
NAMES = ('A1', 'I2(5)', 'H3')


def capped(counts):
    a, b, h, d = counts
    return a, min(b, 2), min(h, 1), min(d, 2)


def pool_key(name, counts):
    return (name,)+capped(counts)


def denominator(name, counts):
    n = counts[NAMES.index(name)]
    assert n >= 2
    size = {'A1': 1, 'I2(5)': 2, 'H3': 3}[name]
    return size*size*n*(n-1)//2


def run():
    with (ROOT/W.DUAL).open() as stream:
        saved = json.load(stream)
    masses = {i: -F(y) for i, y in saved['exact_dual_rows']}
    rows, _, models, profiles, grouped, _, _ = R.model()
    metadata = R.metadata(51954, profiles, grouped)
    capacities = defaultdict(F)
    for i, info in metadata.items():
        if i not in masses:
            continue
        counts = tuple(info['types'].count(name) for name in R.E.HUBS) if 'types' in info else tuple(info['hub_tuple'])
        for k, name in enumerate(NAMES):
            if counts[k] >= 2:
                capacities[pool_key(name, counts)] += masses[i]
    row_index = {(tuple(sorted(row.items())), rhs): i
                 for i, (row, rhs, kind) in enumerate(rows) if kind == 'face' and i in masses}
    assert len(row_index) == sum(i < 51443 for i in masses)
    del rows, metadata
    rule = P.Rule()
    occurrences, positive_counts = Counter(), Counter()
    minima, witnesses = {}, {}
    faces = cases = 0
    for fi, (family, index, row) in enumerate(P.native_faces()):
        faces += 1
        edges = dict(zip(combinations(range(row['nodes']), 2), row['labels']))
        corners = [P.corner(row, edges, a, b) for a, b in product((0, 2), (1, 3))]
        self_families = []
        for a, b in product((0, 2), (1, 3)):
            pieces = P.G.pieces(P.G.graph(row['white']+[a, b], edges))
            left, right = next(g for g in pieces if a in g), next(g for g in pieces if b in g)
            lf, rf = P.G.family(left), P.G.family(right)
            self_families.append(lf if b not in left and lf == rf and lf in NAMES else None)
        baseline = sum((c[0] for c in corners), F())
        for outside in P.spectators(23-row['nodes']):
            cases += 1
            co, pools = Counter(), Counter()
            for self_family, (_, local, key, old_denominator, _) in zip(self_families, corners):
                counts = tuple(x+y for x, y in zip(local, outside))
                if key is not None:
                    factor = F(1)
                    if old_denominator is not None:
                        k, multiplier, offset = old_denominator
                        factor /= multiplier*(counts[k]+offset)
                    co[rule.column(key, counts)] -= factor
                if self_family is not None:
                    pools[pool_key(self_family, counts)] += F(1, denominator(self_family, counts))
            i = row_index.get((tuple(sorted(co.items())), baseline-1))
            if i is None:
                continue
            occurrences[i] += 1
            for key, factor in pools.items():
                identity = i, key
                positive_counts[identity] += 1
                if identity not in minima or factor < minima[identity]:
                    minima[identity] = factor
                    witnesses[identity] = {'native_face': fi, 'family': family, 'family_index': index,
                                           'outside_A1_I25_H3_A4': list(outside),
                                           'native_corner_forests': [sorted(c[4]) for c in corners]}
    assert faces == 1979 and cases == 102917
    assert all(occurrences[i] for i in masses if i < 51443)
    demand, terms = defaultdict(F), defaultdict(list)
    for (i, key), factor in minima.items():
        if positive_counts[i, key] == occurrences[i]:
            mass = masses[i]*factor
            demand[key] += mass
            terms[key].append({'global_face_row': i, 'demand_mass': str(mass),
                               'minimum_coefficient': str(factor), 'old_literal_occurrences': occurrences[i],
                               'witness': witnesses[i, key]})
    columns = []
    for key in sorted(set(demand) | set(capacities)):
        load = demand[key]-capacities[key]
        columns.append({'family': key[0], 'context_A1_I25_H3_A4': list(key[1:]),
                        'minimum_face_demand': str(demand[key]),
                        'budget_capacity': str(capacities[key]), 'forced_dual_load': str(load),
                        'strictly_breaks_any_old_row_allocation': load > 0,
                        'largest_face_contributors': sorted(terms[key], key=lambda t: -F(t['demand_mass']))[:8]})
    columns.sort(key=lambda c: (-F(c['forced_dual_load']), c['family'], c['context_A1_I25_H3_A4']))
    helper = GUARD_HELPER
    with (ROOT/helper).open('rb') as stream:
        helper_hash = file_digest(stream, 'sha256').hexdigest()
    return {'status': 'PASS_EXACT_SELF_HUB_POOL_COLUMN_PROBE', 'top_face_rows': [],
            'pool_columns': columns, 'strictly_breaking_columns': sum(c['strictly_breaks_any_old_row_allocation'] for c in columns),
            'native_faces_checked': faces, 'spectator_cases_checked': cases,
            'old_face_rows_with_positive_dual': len(row_index),
            'supplemental_sources': {helper: helper_hash},
            'numeric_optimizer_run': False, 'new_dimension_bound': None,
            'scope': 'Positive forced load disproves extension of this exact old dual for the new pool. '
                     'This does not prove that the enlarged payment model attains40. '
                     'Minima range over every native witness of each old identical coefficient row.'}


if __name__ == '__main__':
    # Reuse the completed diagnostic's bounded, locked, monitored execution shell.
    W.TAG, W.run = TAG, run
    W.main()
