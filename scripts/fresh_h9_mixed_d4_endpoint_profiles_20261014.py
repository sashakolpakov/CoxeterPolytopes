#!/usr/bin/env python3
"""Complete local spherical endpoint profiles on 52 saved rank-eight supports.

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
TAG = 'fresh_h9_mixed_d4_endpoint_profiles_20261014'
BASE = 'fresh_h9_mixed_d4_coupling_20261014'
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
                new_edges.update({(i, 12): m for i, m in zip(touched, contact_labels)})
                comps, _ = components(S + [12], new_edges)
                types = [component_type(c, new_edges) for c in comps]
                attached = next(k for k, comp in enumerate(comps) if 12 in comp)
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
    paths = [f'results/{BASE}.json', f'results/{BASE}_review.json', f'results/{BASE}_frontier.json']
    data, review, frontier = [json.loads((ROOT / p).read_text()) for p in paths]
    for doc in [data, review, frontier]:
        for p, h in doc['sources'].items():
            assert sha(p) == h, p
    M, N = dm(data['form']), dm(data['final_normals'])
    norms = [pair(n, n, M) for n in N]
    labels = {}
    for row in data['final_pairs']:
        i, j = row['pair']
        g = pair(N[i], N[j], M)
        assert g.sign() <= 0
        if row['type'] == 'dotted':
            assert (g * g - norms[i] * norms[j]).sign() > 0
            labels[i, j] = 0
        else:
            m = row['label']
            target = {2: Q(), 3: Q(F(1, 4)), 5: Q(F(3, 8), F(1, 8))}[m]
            assert g * g == E(target) * norms[i] * norms[j]
            labels[i, j] = m
    assert len(labels) == 66
    ids = frontier['open_direction_ids']
    priority = frontier['source_wall_releases_at_new_polar_vertices']
    assert len(ids) == 52 and priority == [18, 19, 20, 21, 45, 46, 47, 48]
    output = {}
    for edge_id in priority + [i for i in ids if i not in priority]:
        old = frontier['queue'][edge_id]
        assert old['id'] == edge_id and old['kind'] == 'open_to_light_cone'
        row = profiles(old['support'], labels)
        row.update(edge_id=edge_id, released_wall=old['released_wall'], start_ray=old['start_ray'])
        # Exact unnormalized oriented Lorentz frame for the old edge line.
        # U=X/sqrt(A), V=Y/sqrt(Y^2), with V pointing into the outgoing edge.
        X = dv(old['start_point'])
        q = dv(data['final_complete_ray_census']['rays'][old['other_cone_ray']]['vector'])
        A = -pair(X, X, M)
        Xq = pair(X, q, M)
        Y = [A * v + Xq * u for u, v in zip(X, q)]
        Yn = pair(Y, Y, M)
        assert A.sign() > 0 and Yn.sign() > 0 and not pair(X, Y, M)
        assert all(not pair(v, N[i], M) for v in [X, Y] for i in old['support'])
        assert not pair(X, N[old['released_wall']], M)
        assert pair(Y, N[old['released_wall']], M).sign() < 0
        assert pair(X, dv(data['center']), M).sign() < 0
        row['oriented_residual_frame'] = dict(X=X, X_negative_norm=A, Y=Y, Y_norm=Yn,
            normalization='U=X/sqrt(X_negative_norm), V=Y/sqrt(Y_norm)',
            forward_released_wall_pairing=pair(Y, N[old['released_wall']], M))
        output[edge_id] = row
    accepted = sum(len(r['finite_profiles']) for r in output.values())
    tails = sum(len(r['dihedral_tail_profiles']) for r in output.values())
    components_count = Counter(p['new_node_component'] for r in output.values() for p in r['finite_profiles'])
    priority_summary = {i: dict(support=output[i]['support'], support_types=output[i]['support_component_types'],
                               finite_profiles=len(output[i]['finite_profiles']), tails=len(output[i]['dihedral_tail_profiles']),
                               new_component_counts=dict(Counter(p['new_node_component'] for p in output[i]['finite_profiles'])))
                        for i in priority}
    paths += [f'scripts/{TAG}.py', f'scripts/{BASE}_review.py', f'scripts/{BASE}_frontier.py']
    return encode(dict(status='PASS_COMPLETE_LOCAL_SPHERICAL_PROFILES_52_MIXED_D4_OPEN_EDGES',
                       sources={p: sha(p) for p in paths}, normalized_wall_norm=2,
                       profile_scalar_field='Q(sqrt5), [a,b] means a+b*sqrt5',
                       physical_frame_field=data['field'],
                       open_edge_count=52, finite_profile_count=accepted, unbounded_tail_count=tails,
                       new_node_component_counts=dict(sorted(components_count.items())), priority_edges=priority_summary,
                       profiles=output,
                       completeness='Cycles have a nonpositive all-ones quadratic form. Four contacts give a nonpositive four-leaf star. A three-contact star can only have labels3/3/3. A label>=6 with any other adjacent edge gives a nonpositive three-node path, so every such label is a single attachment to an isolated old A1. All remaining3/4/5 cases are tested by exact positive Schur complement. Contacts in distinct old components have zero inverse cross terms.',
                       residual_parameterization='For normalized support roots s_i and contact magnitudes c_i=2cos(pi/m_i), put p=-sum_i (G^-1 c)_i s_i and rho=2-c^T G^-1 c>0. With any orthonormal residual frame U^2=-1,V^2=1,U.V=0 and U future, every norm2 root with those support contacts is p+sqrt(rho)*(sinh(t)*U+epsilon*cosh(t)*V), t real, epsilon=+-1. Additional old-wall Coxeter gates, future endpoint location and firstness are not imposed here.',
                       oriented_edge_reduction='In the saved forward frame, an endpoint cap retaining the finite start must use epsilon=+1 and t>0. Its future endpoint is cosh(t)*U+sinh(t)*V. For T=exp(t)>1, acuteness to an old root w becomes (sqrt(rho)/2)*(U+V).w*T^2+p.w*T+(sqrt(rho)/2)*(V-U).w<=0. This is a necessary continuous gate, not an assertion that a parameter passes all old angles or inequalities.',
                       scope='Complete necessary local profiles for a new finite endpoint on each saved open edge. Includes H3/H4 and every isolated I2(m) tail. Not a cap existence or compatibility result, not an affine endpoint census, not an enumeration of vertex-free edges, and not a finite-volume completion.'))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = run()
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['finite_profile_count'], result['unbounded_tail_count'], flush=True)
