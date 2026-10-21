#!/usr/bin/env python3
"""Conservative Q(sqrt5) face-weight LP, including all finite extensions.

SciPy proposes weights; rational arithmetic checks certificates. This is a
method diagnostic until its geometric coverage receives independent review.
"""
import argparse
import json
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

import fresh_sqrt5_sealed_weights_20261015 as T

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_uniform_face_weights_20261015'


def models(bound):
    out = []
    def add(name, n, edges):
        H = {i: {} for i in range(n)}
        for a, b, m in edges:
            H[a][b] = H[b][a] = m
        out.append((name, H))
    for n in range(1, bound + 1):
        add(f'A{n}', n, [(i, i+1, 3) for i in range(n-1)])
    for n in range(2, bound + 1):
        add(f'B{n}', n, [(i, i+1, 4 if i == 0 else 3) for i in range(n-1)])
    for n in range(4, bound + 1):
        add(f'D{n}', n, [(i, i+1, 3) for i in range(n-2)] + [(1, n-1, 3)])
    for n in (6, 7, 8):
        if n <= bound:
            add(f'E{n}', n, [(i, i+1, 3) for i in range(n-2)] + [(2, n-1, 3)])
    if bound >= 4:
        add('F4', 4, [(0, 1, 3), (1, 2, 4), (2, 3, 3)])
    for n in (3, 4):
        if n <= bound:
            add(f'H{n}', n, [(i, i+1, 5 if i == 0 else 3) for i in range(n-1)])
    for m in (5, 6, 10):
        add(f'I2({m})', 2, [(0, 1, m)])
    return out


@lru_cache(maxsize=100000)
def small_key(labels, n, marks):
    g = dict(zip(combinations(range(n), 2), labels))
    H = T.graph(list(range(n)), g)
    return T.keystr(T.tree_key(H, marks))


def induced_key(H, nodes, a, b):
    nodes = tuple(sorted(nodes))
    return small_key(tuple(H[i].get(j, 2) for i, j in combinations(nodes, 2)),
                     len(nodes), tuple(sorted((nodes.index(a), nodes.index(b)))))


def build(bound):
    unions = json.loads((ROOT/'results/fresh_sqrt5_all_quad_unions_20261015.json').read_text())
    quotient = json.loads((ROOT/'results/fresh_sqrt5_all_quad_unions_20261015_classes.json').read_text())
    triangles = json.loads((ROOT/'results/fresh_sqrt5_triangle_weights_20261015.json').read_text())
    faces = []
    partial = set()
    def face(kind, ident, white, cornerpairs, g):
        coef = Counter()
        disconnected = 0
        for a, b in cornerpairs:
            H = T.graph(white + [a, b], g)
            C = T.component(H, a)
            if b not in C:
                disconnected += 1
                continue
            pk = induced_key(H, C, a, b)
            coef[pk] += 1
            partial.add(pk)
        faces.append(dict(kind=kind, index=ident, required=2 if kind == 'triangle' else 1,
                          coefficients=sorted(coef.items()), disconnected=disconnected))
    for ci, cls in enumerate(quotient['marked_isomorphism_classes']):
        row = unions['labelled_union_supports'][cls['representative']]
        face('quadrilateral', ci, row['white'], product((0, 2), (1, 3)), T.edges(row))
    for i, row in enumerate(triangles['markings']):
        face('triangle', i, row['white'], combinations(row['black'], 2),
             dict(zip(combinations(range(row['order']), 2), row['labels'])))
    print('faces', len(faces), 'partial orbits', len(partial), 'empty',
          [(f['kind'], f['index']) for f in faces if not f['coefficients']], flush=True)
    full = []
    component_rows = []
    extensions = set()
    subset_counts = Counter()
    for name, H in models(bound):
        orbits = defaultdict(list)
        pair_keys = {}
        for a, b in combinations(H, 2):
            pk = T.keystr(T.tree_key(H, (a, b)))
            orbits[pk].append((a, b))
            pair_keys[a, b] = pk
        component_rows.append(dict(type=name, rank=len(H),
                                   labels=[H[a].get(b, 2) for a, b in combinations(H, 2)],
                                   orbits=[dict(key=k, multiplicity=len(v)) for k, v in sorted(orbits.items())]))
        full.extend(orbits)
        subsets = {frozenset((i,)) for i in H}
        # Every connected subset is obtainable by deleting a leaf, inductively.
        # Hence these frontiers exhaust all connected induced subtrees <=8 nodes.
        for size in range(2, min(8, len(H)) + 1):
            subsets = {S | {j} for S in subsets for i in S for j in H[i] if j not in S}
            subset_counts[size] += len(subsets)
            for S in sorted(subsets, key=lambda z: tuple(sorted(z))):
                for a, b in combinations(sorted(S), 2):
                    pk = induced_key(H, S, a, b)
                    if pk in partial:
                        extensions.add((pk, pair_keys[a, b]))
        if name in ('A29', 'B29', 'D29', 'I2(10)'):
            print('extensions through', name, len(extensions), flush=True)
    assert len(full) == len(set(full))
    full = sorted(full)
    partial = sorted(partial)
    full_idx = {k: i for i, k in enumerate(full)}
    partial_idx = {k: i + len(full) for i, k in enumerate(partial)}
    z_index = len(full) + len(partial)
    assert {p for p, _ in extensions} == set(partial)
    rows, rhs, meta = [], [], []
    def addrow(coef, b, info):
        rows.append(dict(sorted(coef.items())))
        rhs.append(b)
        meta.append(info)
    # Literal and marked-class multiplicities do not alter a weight inequality.
    unique_faces = {}
    for f in faces:
        key = (f['required'], tuple(tuple(t) for t in f['coefficients']))
        unique_faces.setdefault(key, []).append((f['kind'], f['index']))
    for (required, coef), witnesses in sorted(unique_faces.items()):
        addrow({partial_idx[k]: -v for k, v in coef}, -required,
               dict(kind='face', required=required, witnesses=witnesses))
    for pk, fk in sorted(extensions):
        addrow({partial_idx[pk]: 1, full_idx[fk]: -1}, 0,
               dict(kind='extension', partial=pk, full=fk))
    component_coefs = []
    for row in component_rows:
        component_coefs.append({full_idx[o['key']]: o['multiplicity'] for o in row['orbits']})
    return dict(full=full, partial=partial, components=component_rows, component_coefs=component_coefs,
                rows=rows, rhs=rhs, meta=meta, faces=faces, z_index=z_index,
                connected_subset_counts={str(k): v for k, v in sorted(subset_counts.items())},
                extension_constraints=len(extensions), distinct_face_constraints=len(unique_faces))


def max_forest(components, totals, bound):
    values = [F(0)] + [None] * bound
    forests = [()] + [None] * bound
    for n in range(1, bound + 1):
        candidates = [(values[n-r['rank']] + v, tuple(sorted(forests[n-r['rank']] + (i,))))
                      for i, (r, v) in enumerate(zip(components, totals)) if r['rank'] <= n]
        values[n], forests[n] = max(candidates)
    return values, forests


@lru_cache(None)
def marked_distance(key):
    paths = []
    def visit(node, path):
        marked, children = node
        if marked:
            paths.append(path)
        for i, (_, child) in enumerate(children):
            visit(child, path + (i,))
    visit(json.loads(key), ())
    assert len(paths) == 2
    common = 0
    for a, b in zip(*paths):
        if a != b:
            break
        common += 1
    return len(paths[0]) + len(paths[1]) - 2*common


def distance_six_certificate(data, bound):
    weights = [int(marked_distance(k) <= 6) for k in data['full'] + data['partial']]
    totals = [sum(v*weights[j] for j, v in row.items()) for row in data['component_coefs']]
    values, forests = max_forest(data['components'], totals, bound)
    x = weights + [values[bound]]
    assert all(sum(v*x[j] for j, v in row.items()) <= b
               for row, b in zip(data['rows'], data['rhs']))
    return dict(rule='Weight1 for connected marked pairs at graph distance at most6; otherwise0.',
                exact_maximum=str(values[bound]),
                maximizing_forest=[data['components'][i]['type'] for i in forests[bound]],
                component_totals={r['type']: str(v) for r, v in zip(data['components'], totals)},
                rank_maxima=list(map(str, values)))


def solve(data, bound):
    rows, rhs, meta = data['rows'], data['rhs'], data['meta']
    comps, component_coefs = data['components'], data['component_coefs']
    zi = data['z_index']
    forests = set()
    def cut(forest):
        forest = tuple(sorted(forest))
        assert forest not in forests
        assert sum(comps[i]['rank'] for i in forest) == bound
        forests.add(forest)
        coef = Counter({zi: -1})
        for i in forest:
            coef.update(component_coefs[i])
        rows.append(dict(sorted(coef.items())))
        rhs.append(0)
        meta.append(dict(kind='forest', components=[comps[i]['type'] for i in forest]))
    a1 = next(i for i, r in enumerate(comps) if r['type'] == 'A1')
    for i, r in enumerate(comps):
        q, rem = divmod(bound, r['rank'])
        forest = tuple(sorted((i,) * q + (a1,) * rem))
        if forest not in forests:
            cut(forest)
    objective = np.zeros(zi + 1)
    objective[zi] = 1
    iterations = []
    for iteration in range(1000):
        rr, cc, vv = [], [], []
        for i, row in enumerate(rows):
            for j, v in row.items():
                rr.append(i); cc.append(j); vv.append(v)
        A = coo_matrix((vv, (rr, cc)), shape=(len(rows), zi+1)).tocsr()
        result = linprog(objective, A_ub=A, b_ub=np.array(rhs), bounds=(0, None), method='highs')
        if not result.success:
            return dict(status='NUMERICAL_LP_FAILURE_NO_THEOREM', message=result.message,
                        iterations=iterations, constraints=len(rows))
        x = [F(float(v)).limit_denominator(10**7) for v in result.x]
        assert all(v >= 0 for v in x)
        slack = [F(b) - sum(v*x[j] for j, v in row.items()) for row, b in zip(rows, rhs)]
        assert min(slack) >= 0, ('primal rationalization', min(slack))
        totals = [sum(v*x[j] for j, v in coef.items()) for coef in component_coefs]
        values, best_forests = max_forest(comps, totals, bound)
        print('LP', iteration, 'objective', x[zi], 'forest maximum', values[bound],
              [comps[i]['type'] for i in best_forests[bound]], flush=True)
        iterations.append(dict(objective=str(x[zi]), forest_maximum=str(values[bound]),
                               forest=[comps[i]['type'] for i in best_forests[bound]]))
        if values[bound] > x[zi]:
            cut(best_forests[bound])
            continue
        y = [F(float(v)).limit_denominator(10**7) for v in result.ineqlin.marginals]
        assert all(v <= 0 for v in y)
        dual = [F(0)] * (zi+1)
        for weight, row in zip(y, rows):
            if weight:
                for j, v in row.items():
                    dual[j] += weight*v
        assert all(v <= (1 if j == zi else 0) for j, v in enumerate(dual)), 'dual rationalization'
        assert sum(v*b for v, b in zip(y, rhs)) == x[zi]
        return dict(status='PASS_EXACT_CONSERVATIVE_WEIGHT_OPTIMUM', optimum=str(x[zi]),
                    target=3*(bound-1), meets_target=x[zi] <= 3*(bound-1),
                    exact_primal=list(map(str, x)),
                    nonzero_dual=[dict(row=i, multiplier=str(v), constraint=meta[i])
                                  for i, v in enumerate(y) if v],
                    component_totals={r['type']: str(v) for r, v in zip(comps, totals)},
                    rank_maxima=[str(v) for v in values],
                    maximizing_forest=[comps[i]['type'] for i in best_forests[bound]],
                    minimum_slack=str(min(slack)), iterations=iterations,
                    forest_constraints=[m for m in meta if m['kind']=='forest'],
                    constraints=len(rows))
    raise AssertionError('cutting-plane iteration guard')


def run(bound):
    data = build(bound)
    print('models', len(data['components']), 'full orbits', len(data['full']),
          'extensions', data['extension_constraints'], 'face rows', data['distinct_face_constraints'], flush=True)
    solution = solve(data, bound)
    simple = distance_six_certificate(data, bound)
    if solution['status'] == 'PASS_EXACT_CONSERVATIVE_WEIGHT_OPTIMUM':
        assert F(simple['exact_maximum']) == F(solution['optimum'])
    paths = [f'scripts/{TAG}.py', 'scripts/fresh_sqrt5_sealed_weights_20261015.py',
             'results/fresh_sqrt5_all_quad_unions_20261015.json',
             'results/fresh_sqrt5_all_quad_unions_20261015_classes.json',
             'results/fresh_sqrt5_triangle_weights_20261015.json']
    return dict(status=solution['status'], rank=bound,
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                face_count=len(data['faces']), distinct_face_constraints=data['distinct_face_constraints'],
                full_pair_orbits=len(data['full']), partial_pair_orbits=len(data['partial']),
                extension_constraints=data['extension_constraints'],
                connected_subset_counts=data['connected_subset_counts'],
                components=data['components'], full_orbit_order=data['full'], partial_orbit_order=data['partial'],
                face_rows=data['faces'], solution=solution, distance_six_certificate=simple,
                coverage_argument='Every actual triangle contains a minimal Lanner core with three black nodes; every actual quadrilateral contains two opposite minimal indefinite pieces. The complete catalogues supply their marked cores. At each corner the connected marked-pair component of a core is an induced connected subtree of the actual finite component, of rank at most8. All such induced embeddings into every allowed irreducible finite component of rank at most the specified bound are included. Auxiliary weights are no larger than every corresponding entire-component weight. Disconnected core marked pairs are charged zero. Nonnegative weights make that conservative even if later white vertices connect them. The rank-partition dynamic program includes every vertex forest.',
                scope='Exact optimization of a conservative component-local weighting method, zero on disconnected entire-component pairs. Symbolic dotted and principal-union conditions may include unrealizable or unextendible faces; extension inequalities also include arithmetically impossible completions. An optimum above target obstructs only this conservative system, not Coxeter polytopes or unrestricted weight methods. An optimum at or below target requires independent geometric coverage review before any dimension claim.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--rank', type=int, default=29)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    assert 8 <= args.rank <= 29
    # Normalize tuple-valued canonical metadata to its JSON representation.
    out = json.loads(json.dumps(run(args.rank)))
    path = ROOT/'results'/f'{TAG}_rank{args.rank}.json'
    if args.verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out, indent=2) + '\n')
    print(out['status'], out['solution'].get('optimum'), flush=True)
