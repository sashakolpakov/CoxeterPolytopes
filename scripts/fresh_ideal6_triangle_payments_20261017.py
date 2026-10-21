#!/usr/bin/env python3
"""Exact triangle census for finite9/ideal6 with feature and small-factor payments.

The universal supporting-white argument is written in the accompanying report.
No quadrilateral theorem or new ambient dimension bound is asserted here.
"""
import argparse
import json
from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_ideal6_triangle_payments_20261017'
CENSUS = 'fresh_high_rank_quasilanner_census_20261017'


def paths(adj, nodes, start):
    nodes = set(nodes)
    out = {start: [start]}
    todo = deque([start])
    while todo:
        i = todo.popleft()
        for j in sorted((adj[i] & nodes) - out.keys()):
            out[j] = out[i] + [j]
            todo.append(j)
    return out


def positive_pivots(matrix):
    a = [list(map(F, row)) for row in matrix]
    pivots = []
    while a:
        q = a[0][0]
        assert q > 0
        pivots.append(q)
        a = [[a[i][j] - a[i][0]*a[0][j]/q
              for j in range(1, len(a))] for i in range(1, len(a))]
    return list(map(str, pivots))


def positive_kernel(matrix):
    a = [list(map(F, row)) for row in matrix]
    n = len(a)
    pivot_columns = []
    row = 0
    for col in range(n):
        p = next((i for i in range(row, n) if a[i][col]), None)
        if p is None:
            continue
        a[row], a[p] = a[p], a[row]
        q = a[row][col]
        a[row] = [x/q for x in a[row]]
        for i in range(n):
            if i != row:
                q = a[i][col]
                a[i] = [x-q*y for x, y in zip(a[i], a[row])]
        pivot_columns.append(col)
        row += 1
    free = sorted(set(range(n)) - set(pivot_columns))
    assert len(free) == 1
    z = [F(int(i == free[0])) for i in range(n)]
    for row, col in enumerate(pivot_columns):
        z[col] = -a[row][free[0]]
    assert all(x > 0 for x in z)
    assert all(sum(F(x)*y for x, y in zip(r, z)) == 0 for r in matrix)
    return z


def graph(diagram):
    adj = {i: set() for i in range(diagram['order'])}
    labels = {}
    for i, j, m in diagram['edges']:
        adj[i].add(j)
        adj[j].add(i)
        labels[min(i,j), max(i,j)] = m
    return adj, labels


def corner_record(adj, labels, marks, a, b):
    nodes = sorted(set(adj) - (set(marks) - {a,b}))
    reach = paths(adj, nodes, a)
    p = reach.get(b)
    if p is None:
        return dict(pair=[a,b], path=None, length=None,
                    feature_witness=None, elementary_payment=False)
    component = set(reach)
    features = {i for i in component if len(adj[i] & component) >= 3}
    for (i,j), m in labels.items():
        if m > 3 and i in component and j in component:
            features.update((i,j))
    witness = None
    for start in (a,b):
        routes = paths(adj, component, start)
        for f in sorted(features):
            route = routes[f]
            if len(route)-1 <= 6:
                witness = dict(endpoint=start, feature=f, path=route,
                               feature_neighbors=sorted(adj[f] & component),
                               incident_labels=[[j,labels[min(f,j),max(f,j)]]
                                                for j in sorted(adj[f] & component)])
                break
        if witness:
            break
    length = len(p)-1
    return dict(pair=[a,b], path=p, length=length,
                feature_witness=witness,
                elementary_payment=length <= 6 or (length <= 9 and witness is not None))


def exceptional(diagram):
    assert diagram['index'] == 25 and diagram['order'] == 10
    adj, labels = graph(diagram)
    assert set(labels.values()) == {3}
    g = [[2 if i == j else -int(j in adj[i]) for j in range(10)] for i in range(10)]
    affine = list(range(9))
    f = [[g[i][j] for j in affine] for i in affine]
    z = positive_kernel(f)
    # The positive-definite eight-node principal matrix and the one-dimensional
    # kernel prove that the full nine-node matrix is positive semidefinite.
    pivots = positive_pivots([r[:-1] for r in f[:-1]])
    assert diagram['deletions'][9] == [dict(kind='affine',name='affine_E8',rank=8)]
    path = [1,0,4,5,6,7,8,9]
    nodes = sorted(set(range(10)) - {2})
    assert [i for i in nodes if i not in path] == [3]
    assert all(j in adj[i] and labels[min(i,j),max(i,j)] == 3 for i,j in zip(path,path[1:]))
    assert sum(len(adj[i] & set(path)) for i in path) == 14
    assert not (adj[3] & set(nodes))
    assert adj[1] & set(nodes) == {0}
    assert adj[0] & set(nodes) == {1,4}
    assert {1,0} <= set(affine)
    corner_pivots = positive_pivots([[g[i][j] for j in nodes] for i in nodes])
    return dict(diagram=25, marks=[1,2,9], affine_nodes=affine,
                affine_type='affine_E8', affine_kernel=list(map(str,z)),
                affine_rank=8, positive_principal_pivots=pivots,
                short_paid_pair=[1,2], short_path=[1,0,2],
                replacement_pair=[1,9], replacement_path=path,
                finite_core_type='A8+A1', finite_core_pivots=corner_pivots,
                fixed_simple_leaf=1, fixed_degree_two_neighbor=0,
                maximum_exceptional_affine_order=9,
                universal_argument='All extra supporting whites are orthogonal to affine nodes0..8. Thus leaf1 and its degree-two neighbor0 persist at corner(1,9). No classical irreducible affine diagram has this simple-leaf configuration; exceptional affine factors have at most9 nodes and receive the whole-factor payment. A finite corner receives the distance7 payment.')


def run():
    data = json.loads((ROOT/f'results/{CENSUS}.json').read_text())
    rows = []
    exceptions = []
    short_failures = []
    for diagram in data['diagrams']:
        if diagram['order'] < 8:
            continue
        adj, labels = graph(diagram)
        for marks in combinations(range(diagram['order']), 3):
            corners = [corner_record(adj, labels, marks, a, b)
                       for a,b in combinations(marks,2)]
            count = sum(c['elementary_payment'] for c in corners)
            row = dict(diagram=diagram['index'], order=diagram['order'],
                       marks=list(marks), corners=corners,
                       elementary_payment_count=count)
            if sum(c['length'] is not None and c['length'] <= 6 for c in corners) < 2:
                short_failures.append([diagram['index'],list(marks)])
            if count < 2:
                exceptions.append([diagram['index'],list(marks)])
                row['universal_exception_resolution'] = 'affine_E8_kernel_and_fixed_leaf'
            rows.append(row)
    assert len(rows) == 920
    assert Counter(r['order'] for r in rows) == {8:224,9:336,10:360}
    assert len(short_failures) == 9
    assert exceptions == [[25,[1,2,9]]]
    special = exceptional(next(d for d in data['diagrams'] if d['index'] == 25))
    stems = [CENSUS,'fresh_ideal_distance7_mixed_replacement_20261017',
             'fresh_cusp_factor_payment_framework_20261017']
    sources = [Path(__file__)]
    sources += [ROOT/f'{folder}/{stem}.{ext}' for stem in stems
                for folder,ext in [('scripts','py'),('results','json'),('results','md')]]
    return dict(status='PASS_FINITE9_IDEAL6_TRIANGLE_PAYMENTS',
                sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sources},
                rule=dict(finite_distance_at_most=9,ideal_distance_at_most=6,
                          ideal_feature_distance_at_most=6,ideal_feature_pair_distance_at_most=9,
                          ideal_whole_factor_order_at_most=9),
                rows=rows, counts_by_order=dict(Counter(r['order'] for r in rows)),
                pure_distance6_failures=short_failures,
                elementary_rule_exceptions=exceptions, exceptional_resolution=special,
                lower_order_argument='For a connected core of order at most7 choose a Steiner tree through the three marks. At least two marked pairs have paths avoiding the third mark, and every such path has at most6 edges.',
                scope='Universal triangle lower bound2 using the stated finite and ideal rule, the inherited complete minimal-indefinite classification and actual corner geometry. No quadrilateral theorem, vertex budget or ambient dimension conclusion.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    data = json.loads(json.dumps(run()))
    target = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(target.read_text())
    else:
        target.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'], '920 triples, one affine-kernel resolution', flush=True)
