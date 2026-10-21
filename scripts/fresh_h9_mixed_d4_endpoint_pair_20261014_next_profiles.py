#!/usr/bin/env python3
"""Complete local spherical endpoint profiles on71 saved rank-eight supports.

No cap search: physical acuteness and actual endpoint location remain separate.
The all-label completeness proof uses only cycle/star/path principal minors.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
from fresh_h9_mixed_d4_coupling_20261014_review import Q, E, dm, dv, pair

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_mixed_d4_endpoint_pair_20261014_next_profiles'
BASE = 'fresh_h9_mixed_d4_endpoint_pair_20261014'
PROFILE_BASE = 'fresh_h9_mixed_d4_endpoint_profiles_20261014'
ZERO, ONE, PHI = Q(), Q(1), Q(F(1, 2), F(1, 2))
SQUARES = {3: ONE, 4: Q(2), 5: Q(F(3, 2), F(1, 2))}


def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()


def encode(x):
    if isinstance(x, (Q, E)):
        return x.enc()
    if isinstance(x, F):
        return str(x)
    if isinstance(x, (list, tuple)):
        return [encode(y) for y in x]
    if isinstance(x, dict):
        return {str(k): encode(v) for k, v in x.items()}
    return x


def inverse(M):
    n = len(M)
    A = [row[:] + [Q(i == j) for j in range(n)] for i, row in enumerate(M)]
    for j in range(n):
        p = next(i for i in range(j, n) if A[i][j])
        A[j], A[p] = A[p], A[j]
        d = A[j][j]
        A[j] = [x / d for x in A[j]]
        for i in range(n):
            if i != j and A[i][j]:
                c = A[i][j]
                A[i] = [x - c * y for x, y in zip(A[i], A[j])]
    J = [row[n:] for row in A]
    assert all(sum((M[i][k] * J[k][j] for k in range(n)), ZERO) == Q(i == j)
               for i in range(n) for j in range(n))
    return J


def pd(M):
    A = [row[:] for row in M]
    for i in range(len(A)):
        d = A[i][i]
        assert d.sign() > 0
        for j in range(i + 1, len(A)):
            for k in range(i + 1, len(A)):
                A[j][k] = A[j][k] - A[j][i] * A[i][k] / d


def components(nodes, edges):
    adj = {i: set() for i in nodes}
    for (i, j), m in edges.items():
        if m > 2:
            adj[i].add(j)
            adj[j].add(i)
    left = set(nodes)
    out = []
    while left:
        seen = {min(left)}
        todo = list(seen)
        for i in todo:
            for j in sorted(adj[i] - seen):
                seen.add(j)
                todo.append(j)
        out.append(sorted(seen))
        left -= seen
    return out, adj


def component_type(nodes, edges):
    n = len(nodes)
    local = {p: m for p, m in edges.items() if set(p) <= set(nodes) and m > 2}
    if n == 1:
        return 'A1'
    _, adj = components(nodes, local)
    assert len(local) == n - 1
    if n == 2:
        m = next(iter(local.values()))
        return {3: 'A2', 4: 'B2'}.get(m, f'I2({m})')
    forks = [i for i in nodes if len(adj[i]) == 3]
    if not forks:
        assert max(map(len, adj.values())) == 2
        start = next(i for i in nodes if len(adj[i]) == 1)
        ordered = []
        prev, cur = None, start
        while True:
            nxt = adj[cur] - ({prev} if prev is not None else set())
            if not nxt:
                break
            j = next(iter(nxt))
            ordered.append(local[tuple(sorted((cur, j)))])
            prev, cur = cur, j
        if all(m == 3 for m in ordered):
            return f'A{n}'
        heavy = [(i, m) for i, m in enumerate(ordered) if m != 3]
        assert len(heavy) == 1
        pos, m = heavy[0]
        if m == 4 and pos in [0, n - 2]:
            return f'B{n}'
        if m == 4 and n == 4 and pos == 1:
            return 'F4'
        if m == 5 and n in [3, 4] and pos in [0, n - 2]:
            return f'H{n}'
        raise AssertionError(('unclassified positive path', nodes, ordered))
    assert len(forks) == 1 and all(m == 3 for m in local.values())
    center = forks[0]
    arms = []
    for start in adj[center]:
        prev, cur, length = center, start, 1
        while adj[cur] - {prev}:
            nxt = next(iter(adj[cur] - {prev}))
            prev, cur, length = cur, nxt, length + 1
        arms.append(length)
    arms.sort()
    if arms[:2] == [1, 1]:
        return f'D{n}'
    assert (n, tuple(arms)) in [(6, (1, 2, 2)), (7, (1, 2, 3)), (8, (1, 2, 4))]
    return f'E{n}'


def profiles(S, labels):
    edges = {tuple(sorted((i, j))): labels[tuple(sorted((i, j)))]
             for i, j in combinations(S, 2)}
    assert all(m in [2, 3, 5] for m in edges.values())
    old_components, _ = components(S, edges)
    old_types = [component_type(c, edges) for c in old_components]
    where = {i: k for k, comp in enumerate(old_components) for i in comp}
    G = [[Q(2) if i == j else {2: ZERO, 3: -ONE, 5: -PHI}[edges[tuple(sorted((i, j)))]]
          for j in S] for i in S]
    pd(G)
    J = inverse(G)
    assert all(x.sign() >= 0 for row in J for x in row)
    assert all((J[i][i] - F(1, 2)).sign() >= 0 for i in range(8))
    assert all(not J[a][b] for a, i in enumerate(S) for b, j in enumerate(S)
               if where[i] != where[j])
    rows, rejected, tested = [], Counter(), 0
    for count in range(4):
        for touched in combinations(S, count):
            if len({where[i] for i in touched}) != count:
                continue
            choices = product([3, 4, 5], repeat=count) if count < 3 else [(3, 3, 3)]
            for contact_labels in choices:
                tested += 1
                cost = sum((SQUARES[m] * J[S.index(i)][S.index(i)]
                            for i, m in zip(touched, contact_labels)), ZERO)
                rho = Q(2) - cost
                if rho.sign() <= 0:
                    rejected['affine' if not rho else 'indefinite'] += 1
                    continue
                new_edges = dict(edges)
                new_edges.update({(i, 14): m for i, m in zip(touched, contact_labels)})
                comps, _ = components(S + [14], new_edges)
                types = [component_type(c, new_edges) for c in comps]
                attached = next(k for k, comp in enumerate(comps) if 14 in comp)
                rows.append(dict(contacts=[[i, m] for i, m in zip(touched, contact_labels)],
                                 projection_square=cost, residual_norm=rho,
                                 endpoint_components=types,
                                 new_node_component=types[attached]))
    tails = []
    for comp in old_components:
        if len(comp) != 1:
            continue
        i = comp[0]
        assert J[S.index(i)][S.index(i)] == F(1, 2)
        tails.append(dict(wall=i, integer_labels='every integer m >= 6',
                          other_contacts='all right', projection_square='2 cos(pi/m)^2',
                          residual_norm='2 sin(pi/m)^2', new_node_component='I2(m)',
                          other_components=[component_type(c, edges) for c in old_components if c != comp]))
    assert tested == len(rows) + sum(rejected.values())
    return dict(support=S, support_components=old_components, support_component_types=old_types,
                normalized_Gram=G, normalized_inverse=J,
                finite_profiles=rows, dihedral_tail_profiles=tails,
                finite_trial_count=tested, rejected_trial_counts=dict(rejected))


def run():
    primary = f'results/{BASE}.json'
    data = json.loads((ROOT / primary).read_text())
    for p, h in data['sources'].items():
        assert sha(p) == h, p
    old_profiles = json.loads((ROOT / f'results/{PROFILE_BASE}.json').read_text())
    for p, h in old_profiles['sources'].items():
        assert sha(p) == h, p
    M, N = dm(data['form']), dm(data['normals'])
    assert len(N) == 14
    norms = [pair(n, n, M) for n in N]
    assert all(n.sign() > 0 for n in norms)
    labels, types = {}, Counter()
    targets = {2: Q(), 3: Q(F(1, 4)), 5: Q(F(3, 8), F(1, 8))}
    for i, j in combinations(range(14), 2):
        g = pair(N[i], N[j], M)
        assert g.sign() <= 0
        if (g * g - norms[i] * norms[j]).sign() > 0:
            labels[i, j] = 0
            types['dotted'] += 1
        else:
            ms = [m for m, target in targets.items() if g * g == E(target) * norms[i] * norms[j]]
            assert len(ms) == 1, (i, j)
            labels[i, j] = ms[0]
            types[str(ms[0])] += 1
    assert len(labels) == 91
    queue = data['finite_corner_directions']
    assert len(queue) == 135
    rays = data['complete_ray_census']['rays']
    output, priority, repeated = {}, [], []
    old_by_support = {tuple(r['support']): r for r in old_profiles['profiles'].values()}
    for index, old in enumerate(queue):
        if old['kind'] != 'open_to_light_cone':
            continue
        start, other = rays[old['start_ray']], rays[old['other_cone_ray']]
        S = old['support']
        assert len(S) == 8 and len(start['incident']) == 9
        assert start['kind'] == 'finite' and other['kind'] == 'spacelike'
        assert set(S) == set(start['incident']) - {old['released_wall']}
        assert set(S) <= set(other['incident'])
        row = profiles(S, labels)
        row.update(direction_index=index, edge_id=index, released_wall=old['released_wall'],
                   start_ray=old['start_ray'], other_cone_ray=old['other_cone_ray'])
        row['from_new_finite_vertex'] = any(i >= 12 for i in start['incident'])
        if row['from_new_finite_vertex']:
            priority.append(index)
        else:
            repeated.append(index)
            assert tuple(S) in old_by_support
            original = old_by_support[tuple(S)]
            for k in ('finite_profiles', 'dihedral_tail_profiles', 'normalized_Gram', 'normalized_inverse'):
                assert encode(row[k]) == original[k]
        output[str(index)] = row
    assert len(output) == len({tuple(r['support']) for r in output.values()}) == 71
    assert len(priority) == 21 and len(repeated) == 50
    assert any(12 in r['support'] for r in output.values())
    assert any(13 in r['support'] for r in output.values())
    accepted = sum(len(r['finite_profiles']) for r in output.values())
    tails = sum(len(r['dihedral_tail_profiles']) for r in output.values())
    trials = sum(r['finite_trial_count'] for r in output.values())
    paths = [primary, f'scripts/{BASE}.py', f'scripts/{TAG}.py',
        f'scripts/{PROFILE_BASE}.py', f'results/{PROFILE_BASE}.json',
        f'scripts/{PROFILE_BASE}_scope_review.py', f'results/{PROFILE_BASE}_scope_review.json',
        'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py']
    return encode(dict(status='PASS_COMPLETE_LOCAL_SPHERICAL_PROFILES_71_MIXED_D4_OPEN_EDGES',
        sources={p: sha(p) for p in paths}, original_direction_queue_count=135,
        normalized_wall_norm=2, proposed_new_node_id=14, open_edge_count=71,
        new_vertex_open_direction_count=21, repeated_old_open_direction_count=50,
        priority_direction_indices=priority, repeated_old_direction_indices=repeated,
        finite_trial_count=trials, finite_profile_count=accepted, unbounded_tail_count=tails,
        priority_finite_profile_count=sum(len(output[str(i)]['finite_profiles']) for i in priority),
        priority_unbounded_tail_count=sum(len(output[str(i)]['dihedral_tail_profiles']) for i in priority),
        all_old_pair_counts=dict(types), profiles=output,
        completeness=old_profiles['completeness'],
        residual_parameterization=old_profiles['residual_parameterization'],
        scope='Complete necessary local finite-endpoint profiles for all71 saved open finite-corner directions of the14-wall chamber. Original135 queue indices preserved. Includes arbitrary isolated I2(m),m>=6. No global cap feasibility, affine endpoint classification, first-new-mirror enumeration or finite-volume claim. Physical frames are inherited from the bound full cone and are not rebuilt here.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    out = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['finite_profile_count'], 'finite profiles',
          result['unbounded_tail_count'], 'tails;', result['priority_finite_profile_count'],
          'profiles on21 new directions')


