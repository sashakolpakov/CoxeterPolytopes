#!/usr/bin/env python3
"""Independent saved-primal evaluator for the complete rank21 self-hub model.

Source only until the coordinating agent schedules execution. No optimizer.
Does not import the new pool engine, producer, Bellman builder, or pool helpers.
Shares the older independent angle evaluator, finite-diagram primitives, and
streaming native input reader. Reconstructs all pool selectors/denominators,
component coefficients, hub tuples, and a reverse-alphabet unbounded recurrence.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import file_digest
from itertools import combinations, product
import fcntl
import json
from pathlib import Path
import signal
import check_compact21_full_face_primal_independent_20261018 as A
from run_lean_batches_signal_safe_20261018 import resources

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = 'results/fresh_compact21_self_hub_epigraph_minimum_20261018.json'
HUBS = ('A1', 'I2(5)', 'H3', 'A4')
POOL_TYPES = HUBS[:3]
RANKS = (1, 2, 3, 4)
ARTIFACT_LIMIT = 8 * 1024**2
MARKS = ((0, 1), (0, 3), (2, 1), (2, 3))
SPECIAL = {
    'H4': tuple(sorted(('A1', 'H4', 'E8', 'E8'))),
    'D5': tuple(sorted(('D5', 'E8', 'E8'))),
    'A4': tuple(sorted(('A1', 'A4', 'E8', 'E8'))),
    'A5': tuple(sorted(('A5', 'E8', 'E8'))),
    'B4': tuple(sorted(('A1', 'B4', 'E8', 'E8'))),
    'B5': tuple(sorted(('B5', 'E8', 'E8'))),
}


def digest(path):
    with path.open('rb') as stream:
        return file_digest(stream, 'sha256').hexdigest()


def exact_tuples(rank):
    # No untracked rank-one irreducible exists: A1 is tracked explicitly.
    for a in range(rank + 1):
        for b in range((rank - a)//2 + 1):
            for h in range((rank - a - 2*b)//3 + 1):
                for d in range((rank - a - 2*b - 3*h)//4 + 1):
                    residual = rank - a - 2*b - 3*h - 4*d
                    if residual != 1:
                        yield (a, b, h, d), residual


def capped(q):
    return q[0], min(q[1], 2), min(q[2], 1), min(q[3], 2)


def pool_cost(q, pools):
    return sum((pools.get((name,) + q, F())
                for name, k in zip(POOL_TYPES, q) if k >= 2), F())


def native_faces():
    for rec in A.objects(A.BASE, 'faces'):
        if rec['row']['nodes'] <= 23:
            yield rec['row']
    for rec in A.objects(A.QUEUE, 'records'):
        if rec['result']['singular'] and rec['row']['nodes'] <= 23:
            yield rec['row']


def self_family(row, edges, marks):
    pieces = A.G.pieces(A.G.graph(row['white'] + list(marks), edges))
    left = next(g for g in pieces if marks[0] in g)
    right = next(g for g in pieces if marks[1] in g)
    if marks[1] in left:
        return None
    lf, rf = A.G.family(left), A.G.family(right)
    return lf if lf == rf and lf in POOL_TYPES else None


def read_old_values(saved):
    indexed = {j: F(v) for j, v in saved['positive_bonus_values']}
    assert len(indexed) == len(saved['positive_bonus_values'])
    assert all(isinstance(j, int) and j >= 0 and v > 0
               for j, v in indexed.items())
    values, used, records = {}, set(), 0
    for j, rec in A.M.objects('variable_order'):
        assert j == records
        records += 1
        if j in indexed:
            key = (rec['orbit'], rec['A1_count'], rec['I25_context'],
                   rec['H3_context'], rec['A4_context'])
            assert key not in values
            values[key] = indexed[j]
            used.add(j)
        if records % 4096 == 0:
            A.memory_gate()
    assert records == 93045 and used == set(indexed)
    return values, records


def read_pools(saved, old_columns, grouped):
    expected_keys = sorted({(name,) + q
                            for q, _ in exact_tuples(21)
                            for name, k in zip(POOL_TYPES, q) if k >= 2})
    auxiliary_count = sum(max(0, max(r for _, r in offers) - 1)
                          for offers in grouped.values())
    assert auxiliary_count == 1659
    first = old_columns + 2 + auxiliary_count  # old z,t retain their indices
    order = saved['pool_variable_order']
    assert len(order) == len(expected_keys)
    by_index = {}
    for offset, (rec, key) in enumerate(zip(order, expected_keys)):
        assert (rec['family'],) + tuple(rec['context']) == key
        assert rec['global_column'] == first + offset
        by_index[first + offset] = key
    positive = {j: F(v) for j, v in saved['positive_pool_values']}
    assert len(positive) == len(saved['positive_pool_values'])
    assert all(j in by_index and v > 0 for j, v in positive.items())
    return {by_index[j]: v for j, v in positive.items()}


def check_faces(values, pools):
    count = cases = 0
    minimum = None
    witness = None
    for fi, row in enumerate(native_faces()):
        count += 1
        edges = dict(zip(combinations(range(row['nodes']), 2), row['labels']))
        corners = [A.angle(row, edges, marks) for marks in MARKS]
        same = [self_family(row, edges, marks) for marks in MARKS]
        for outside, _ in exact_tuples(23 - row['nodes']):
            payment = F()
            for (base, local, orbit, divisor_rule), family in zip(corners, same):
                q = tuple(x + y for x, y in zip(local, outside))
                payment += base
                if orbit is not None:
                    bonus = values.get((orbit,) + capped(q), F())
                    if divisor_rule is not None:
                        slot, multiple, offset = divisor_rule
                        divisor = multiple * (q[slot] + offset)
                        assert divisor > 0
                        bonus /= divisor
                    payment += bonus
                if family is not None:
                    slot = POOL_TYPES.index(family)
                    k, rank = q[slot], RANKS[slot]
                    assert k >= 2
                    # Independent unordered component-pair / facet-pair count.
                    pairs = rank**2 * (k*(k-1)//2)
                    payment += pools.get((family,) + q, F()) / pairs
            assert payment >= 1, (fi, outside, str(payment))
            cases += 1
            if minimum is None or payment < minimum:
                minimum, witness = payment, (fi, outside)
        if fi % 32 == 0:
            A.memory_gate()
            resources(start=False)
    assert (count, cases) == (1979, 102917)
    return count, cases, minimum, witness


def component_profiles(models):
    """Literal pair/root multiplicities, independent of the model's profiles."""
    result = {}
    source = models['A4']
    for name, graph in models.items():
        baseline = F()
        tables = [Counter() for _ in range(6)]
        for u, v in combinations(graph, 2):
            base = F(A.G.base(graph, u, v))
            baseline += base
            if not base:
                tables[0][A.encode(A.G.code(graph, (u, v)))] += 1
        family = A.G.family(graph)
        for root in graph:
            rooted = A.encode(A.G.code(graph, (root,)))
            for slot, (hub, prefix) in enumerate(zip(POOL_TYPES, ('A1x', 'I25x', 'H3x'))):
                if name not in POOL_TYPES[:slot+1]:
                    tables[slot+1][prefix + family + ':' + rooted] += 1
            if name not in HUBS:
                # Sum all four literal source roots; source multiplicity d
                # cancels the old 4d normalization.
                for sr in source:
                    key = 'A4sx' + family + ':' + A.encode([
                        A.encode(A.G.code(source, (sr,))), rooted])
                    tables[4][key] += F(1, 4)
        if name == 'A4':
            for u, v in product(graph, repeat=2):
                codes = sorted((A.encode(A.G.code(graph, (u,))),
                                A.encode(A.G.code(graph, (v,)))))
                tables[5]['A4pair:' + A.encode(codes)] += F(1, 8)
        result[name] = len(graph), baseline, tables
    return result


def component_costs(profiles, c, values):
    costs = {}
    enabled = (True, c[0] > 0, c[1] > 0, c[2] > 0, c[3] > 0, c[3] >= 2)
    for name, (_, baseline, tables) in profiles.items():
        total = baseline
        for use, table in zip(enabled, tables):
            if use:
                total += sum((multiplicity * values.get((key,) + c, F())
                              for key, multiplicity in table.items()), F())
        costs[name] = total
    return costs


def extreme_e8_key(graph):
    fork = next(i for i in graph if len(graph[i]) == 3)
    ends = {}
    for start in graph[fork]:
        previous, node, length = fork, start, 1
        while len(graph[node]) == 2:
            nxt = next(i for i in graph[node] if i != previous)
            previous, node, length = node, nxt, length + 1
        ends[length] = node
    assert set(ends) == {1, 2, 4}
    return A.encode(A.G.code(graph, (ends[2], ends[4])))


def check_forests(saved, values, pools, grouped):
    models = {name: graph for name, graph in A.G.models() if len(graph) <= 21}
    assert len(models) == 68
    profiles = component_profiles(models)
    transfer = F(saved['transfer_value'])
    assert 0 <= transfer <= 1
    beta = values.get((extreme_e8_key(models['E8']), 0, 0, 0, 0), F())
    maximum, attaining = F(-1), None
    offers = 0

    def offer(value, forest, kind):
        nonlocal maximum, attaining, offers
        q = tuple(forest.count(name) for name in HUBS)
        assert sum(profiles[name][0] for name in forest) == 21
        value += pool_cost(q, pools)
        offers += 1
        if value > maximum:
            maximum, attaining = value, (tuple(sorted(forest)), kind)

    for ci, (c, tuples) in enumerate(sorted(grouped.items())):
        costs = component_costs(profiles, c, values)
        cap = max(r for _, r in tuples)
        # Reverse-alphabet, type-outer unbounded recurrence; repeated copies
        # of each item are admitted by increasing rank in the inner loop.
        dp = {0: (F(), ())}
        residual_types = sorted(set(models) - set(HUBS), reverse=True)
        for name in residual_types:
            size = profiles[name][0]
            for rank in range(size, cap + 1):
                if rank - size not in dp:
                    continue
                old, forest = dp[rank - size]
                proposal = old + costs[name], tuple(sorted(forest + (name,)))
                if rank not in dp or proposal > dp[rank]:
                    dp[rank] = proposal
        assert all(rank in dp for rank in range(2, cap + 1))
        for q, rank in tuples:
            hubs = tuple(name for name, multiplicity in zip(HUBS, q)
                         for _ in range(multiplicity))
            fixed = sum((costs[name] for name in hubs), F())
            if q in ((0, 0, 0, 0), (1, 0, 0, 0)):
                targets = {'D5', 'E8'} if q[0] == 0 else {'H4', 'E8'}
                for name in residual_types:
                    remaining = rank - profiles[name][0]
                    if name in targets or remaining not in dp:
                        continue
                    tail, forest = dp[remaining]
                    offer(fixed + costs[name] + tail, hubs + (name,) + forest, 'off_target')
            else:
                tail, forest = dp[rank]
                offer(fixed + tail, hubs + forest, 'ordinary')
        if ci % 8 == 0:
            A.memory_gate()
            resources(start=False)
    explicit = [tuple(sorted(('A1',) + ('H4',)*5)),
                tuple(sorted(('A1',) + ('H4',)*3 + ('E8',)))]
    explicit += list(SPECIAL.values())
    for forest in explicit:
        q = tuple(forest.count(name) for name in HUBS)
        costs = component_costs(profiles, capped(q), values)
        value = sum((costs[name] for name in forest), F())
        if forest == SPECIAL['H4']:
            value -= transfer
        elif forest == SPECIAL['D5']:
            value -= min(F(1), 2*beta)
        elif forest in (SPECIAL['A4'], SPECIAL['A5'], SPECIAL['B4'], SPECIAL['B5']):
            value += transfer
        assert pool_cost(q, pools) == 0  # all inherited exceptions
        offer(value, forest, 'explicit')
    assert maximum == F(saved['exact_primal_upper']), (str(maximum), saved['exact_primal_upper'])
    return maximum, attaining, offers


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', default=DEFAULT)
    args = parser.parse_args()
    path = (ROOT / args.candidate).resolve()
    assert path.is_relative_to(ROOT) and path.stat().st_size < ARTIFACT_LIMIT
    candidate_hash = digest(path)
    with path.open() as stream:
        saved = json.load(stream)
    assert saved['status'] in {
        'PASS_EXACT_SELF_HUB_COMPLETE_MODEL_PRIMAL',
        'PASS_EXACT_SELF_HUB_COMPLETE_MODEL_AT_MOST_40',
    }
    pins = saved['sources']
    for relative, expected in pins.items():
        pinned = (ROOT / relative).resolve()
        assert pinned.is_relative_to(ROOT) and digest(pinned) == expected
    grouped = defaultdict(list)
    for q, rank in exact_tuples(21):
        grouped[capped(q)].append((q, rank))
    assert sum(map(len, grouped.values())) == 729
    values, columns = read_old_values(saved)
    pools = read_pools(saved, columns, grouped)
    count, cases, payment, witness = check_faces(values, pools)
    maximum, attainer, offers = check_forests(saved, values, pools, grouped)
    assert digest(path) == candidate_hash
    for relative, expected in pins.items():
        assert digest(ROOT / relative) == expected
    A.memory_gate()
    print(json.dumps({
        'status': 'PASS_INDEPENDENT_SELF_HUB_FULL_PRIMAL',
        'candidate_sha256': candidate_hash,
        'native_faces': count, 'native_spectator_cases': cases,
        'minimum_payment': str(payment), 'minimum_witness': witness,
        'old_variable_records': columns, 'positive_old_bonuses': len(values),
        'positive_self_hub_pools': len(pools), 'exact_hub_tuples': 729,
        'forest_offers': offers, 'exact_primal_upper': str(maximum),
        'independent_attainer': attainer,
        'primal_at_most_40': maximum <= 40,
        'new_model_helpers_imported': False, 'optimizer_run': False,
        'geometric_coverage_reproved': False, 'new_dimension_bound': None,
    }, separators=(',', ':')), flush=True)


if __name__ == '__main__':
    with (ROOT / 'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        signal.alarm(59)
        resources()
        main()
        resources(start=False)
