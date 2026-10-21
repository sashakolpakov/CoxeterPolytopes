#!/usr/bin/env python3
"""Separate graph-count and Bellman verifier; imports no project arithmetic."""
import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = 'fresh_a17_endpoint_d11_upgrade_forest_budget_20261017'
TAG = BASE+'_review'
LIMIT = 75


def rank(typ):
    return int(typ[1:])


def graph(typ):
    n = rank(typ)
    if typ[0] in 'AB':
        edges = [(i, i+1) for i in range(n-1)]
    else:
        assert typ[0] == 'D' and n >= 4
        edges = [(i, i+1) for i in range(n-3)]+[(n-3, n-2), (n-3, n-1)]
    adj = [set() for _ in range(n)]
    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)
    return adj


@lru_cache(None)
def internal(typ):
    adj = graph(typ)
    n = len(adj)
    total = 0
    for i in range(n):
        distance = {i: 0}
        queue = [i]
        for u in queue:
            for v in adj[u]:
                if v not in distance:
                    distance[v] = distance[u]+1
                    queue.append(v)
        assert len(distance) == n
        for j in range(i+1, n):
            paid = n <= 17 or distance[j] <= 8 or (distance[j] == 9 and (len(adj[i]) == 1 or len(adj[j]) == 1))
            if typ == 'A17' and (i, j) == (0, 16):
                paid = False
            total += 2*paid
    return total


def broad_positions(typ):
    n = rank(typ)
    if typ[0] == 'D' and 11 <= n <= 16:
        return {n-3, n-2, n-1}
    if typ in ('A13', 'A15'):
        return {n//2}
    if typ == 'A17':
        return {5, 11}
    return set()


@lru_cache(None)
def cross(a, b, upgrade):
    ra, rb = rank(a), rank(b)
    if ra > 5 and rb > 5 or ra <= 5 and rb <= 5:
        return 0
    if ra > rb:
        return cross(b, a, upgrade)
    # a is the unique small component. Set a weight per actual pair, taking
    # a maximum when an upgrade replaces an inherited half.
    total = 0
    for i in range(ra):
        for j in range(rb):
            value = int(j in broad_positions(b))
            if a == 'A1' and (b == 'D12' or upgrade and b == 'D11') and j in (rb-2, rb-1):
                value = max(value, 2)
            if a == 'A2' and b == 'B10' and j == 9:
                value = max(value, 1)
            if a == 'A4' and i in (0, 3) and b == 'D13' and j == 9:
                value = max(value, 1)
            total += value
    return total


def direct_forest(types, upgrade):
    return sum(internal(t) for t in types)+sum(cross(a, b, upgrade) for a, b in combinations(types, 2))


def domains():
    # Enumerate by component number, as nondecreasing multisets, independently
    # of the producer's successive-rank count-vector recursion.
    vectors = set()
    for count in range(8):
        def walk(left, minimum, total, chosen):
            if not left:
                hist = Counter(chosen)
                vectors.add(tuple(hist[r] for r in range(10, 18)))
                return
            for r in range(minimum, 18):
                if total+r*left <= LIMIT:
                    walk(left-1, r, total+r, chosen+(r,))
        walk(count, 10, 0, ())
    assert len(vectors) == 1608
    return vectors


def domination(upgrade):
    types = ['A'+str(n) for n in range(1, LIMIT+1)]
    types += ['B'+str(n) for n in range(2, LIMIT+1)]
    types += ['D'+str(n) for n in range(4, LIMIT+1)]
    def representative(t):
        n = rank(t)
        if n <= 9:
            return 'A'+str(n)
        if n == 10:
            return 'B10'
        if n == 17 and t == 'A17':
            return t
        return 'D'+str(n)
    comparisons = 0
    for t in types:
        rep = representative(t)
        assert internal(rep) >= internal(t)
        for partner in types:
            assert cross(rep, partner, upgrade) >= cross(t, partner, upgrade)
            comparisons += 1
    return dict(classical_type_profiles=len(types), ordered_partner_comparisons=comparisons,
        rank17_representatives=['A17', 'D17'],
        exceptional_input='All exceptional types and arbitrary I2 lie at rank<=8, where the complete base and broad small eligibility are unchanged; special coupons only favor the retained A representatives.')


def verify_variant(v, vectors):
    upgrade = v['D11_upgrade']
    residual_types = ['A'+str(n) for n in range(1, 10)]+['D'+str(n) for n in range(17, LIMIT+1)]
    large_types = ['B10']+['D'+str(n) for n in range(11, 17)]+['A17']
    assert v['large_types'] == large_types
    assert v['residual_ranks'] == [rank(t) for t in residual_types]
    rows = v['large_count_vectors']
    assert len(rows) == len(vectors) and {tuple(r['component_counts']) for r in rows} == vectors
    tables = {tuple(t['coupon_counts']): t for t in v['residual_tables']}
    assert len(tables) == len(v['residual_tables'])
    expected_keys = set()
    for r in rows:
        c = r['component_counts']
        key = (sum(c[1:7]), c[7], c[0], c[2], c[3], c[1] if upgrade else 0)
        assert list(key) == r['coupon_counts']
        expected_keys.add(key)
        assert r['used_rank'] == sum(n*k for n, k in zip(range(10, 18), c))
        assert r['doubled_base_cost'] == sum(internal(t)*k for t, k in zip(large_types, c))
    assert set(tables) == expected_keys
    inequalities = predecessors = 0
    for key, table in tables.items():
        t, u, b10, b12, b13, b11 = key
        # Build an equivalent coupon environment from independently defined
        # actual graph types. For the non-upgrade variant generic D14 carries
        # the indistinguishable D11/D14/D15/D16 broad coupons.
        environment = (['A17']*u+['B10']*b10+['D12']*b12+['D13']*b13+
            ['D11']*b11+['D14']*(t-b12-b13-b11))
        assert t-b12-b13-b11 >= 0
        weights = {rank(r): internal(r)+sum(cross(r, p, upgrade) for p in environment) for r in residual_types}
        values, last = table['doubled_optima'], table['last_rank']
        assert len(values) == len(last) == LIMIT+1 and values[0] == last[0] == 0
        for n in range(1, LIMIT+1):
            for r, w in weights.items():
                if r <= n:
                    assert values[n] >= values[n-r]+w
                    inequalities += 1
            r = last[n]
            assert r in weights and 0 < r <= n
            assert values[n] == values[n-r]+weights[r]
            predecessors += 1
    maxima = v['dimension_maxima']
    assert [r['dimension'] for r in maxima] == list(range(LIMIT+1))
    global_checks = 0
    for row in rows:
        table = tables[tuple(row['coupon_counts'])]
        for n in range(row['used_rank'], LIMIT+1):
            assert maxima[n]['doubled_cost'] >= row['doubled_base_cost']+table['doubled_optima'][n-row['used_rank']]
            global_checks += 1
    for n, m in enumerate(maxima):
        witness = m['witness_components']
        assert sum(rank(t) for t in witness) == n
        assert direct_forest(witness, upgrade) == m['doubled_cost']
        beta2 = None if n < 3 else (n*(n-6) if n % 2 == 0 else (n-1)*(n-5))//4
        assert m['doubled_curvature_slack'] == (None if beta2 is None else beta2-m['doubled_cost'])
    assert maxima[69]['doubled_cost'] == 1088
    return dict(D11_upgrade=upgrade, type_domination=domination(upgrade),
        independent_large_vectors=len(vectors), residual_tables=len(tables),
        Bellman_upper_bound_checks=inequalities, attaining_predecessor_checks=predecessors,
        global_upper_bound_checks=global_checks, graph_counted_witnesses=len(maxima),
        rank69_cost='544', rank69_finite_curvature_slack='0')


def run():
    path = ROOT/f'results/{BASE}.json'
    d = json.loads(path.read_text())
    assert d['status'] == 'PASS_EXACT_A17_ENDPOINT_AND_CONDITIONAL_D11_FOREST_BUDGETS'
    for name, expected in d['sources'].items():
        assert sha256((ROOT/name).read_bytes()).hexdigest() == expected, name
    vectors = domains()
    out = [verify_variant(v, vectors) for v in d['variants']]
    paths = [Path(__file__), ROOT/f'scripts/{BASE}.py', path, ROOT/f'results/{BASE}.md']
    return dict(status='PASS_SEPARATE_A17_D11_FOREST_CERTIFICATE_REVIEW',
        sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
        producer_arithmetic_imported=False, variants=out,
        scope='Independent multiset domain, graph-based marked-pair costs, every Bellman bound and attaining predecessor, every global bound and all typed maximizing forests. Geometric A17/D11 face proofs and finite classification remain separate inputs.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    args = p.parse_args()
    out = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], [(v['D11_upgrade'], v['Bellman_upper_bound_checks'],
        v['global_upper_bound_checks']) for v in out['variants']], flush=True)
