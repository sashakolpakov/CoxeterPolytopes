#!/usr/bin/env python3
"""Exact finite-vertex counterexample to multi-ideal triangle-port transitivity.

Standalone rational/radical SymPy arithmetic; no project arithmetic imports.
The twenty triple intersections give the entire bounded projective section.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_finite_vertex_t2_port_counterexample_20261017'
S = s.simplify


def sign(x):
    x = S(x)
    ans = s.sign(x)
    assert ans in (-1, 0, 1), ('unresolved exact sign', x, ans)
    return int(ans)


def dot(a, b):
    return S(-a[0]*b[0] + sum(a[i]*b[i] for i in range(1, 4)))


def enc(v):
    return [str(S(x)) for x in v]


def normchart(v):
    assert sign(v[0]) > 0
    return tuple(s.expand(s.radsimp(S(s.sympify(x)/s.sympify(v[0])))) for x in v)


def run():
    r2, r3 = s.sqrt(2), s.sqrt(3)
    names = ['B', 'L', 'D', 'R', 'U', 'C']
    normals = [
        (0, 0, 0, -r2),
        (r2/2, -r2, 0, r2/2),
        (r2/2, 0, -r2, r2/2),
        (s.sqrt(6)/2, r2, 0, s.sqrt(6)/2),
        (s.sqrt(6)/2, 0, r2, s.sqrt(6)/2),
        (2, r3, r3, 0),
    ]
    G = s.Matrix([[dot(a, b) for b in normals] for a in normals])
    assert all(G[i, i] == 2 for i in range(6))
    assert s.Matrix(normals).rank() == 4
    relations = []
    finite = {s.Integer(0): 2, s.Integer(-1): 3, -r3: 6}
    for i, j in itertools.combinations(range(6), 2):
        g = S(G[i, j])
        if g in finite:
            relation = dict(kind='finite', label=finite[g])
        elif g == -2:
            relation = dict(kind='parallel')
        else:
            assert sign(g+2) < 0
            relation = dict(kind='dotted')
        relations.append(dict(pair=[names[i], names[j]], product=str(g), **relation))

    center = tuple(map(s.Rational, ['5/4', '0', '0', '3/4']))
    assert dot(center, center) == -1
    assert all(sign(dot(center, n)) < 0 for n in normals)
    dual = tuple(S(sum(n[i] for n in normals)) for i in range(4))
    assert sign(dual[0]) > 0 and sign(dot(dual, dual)) < 0
    # h=-<sum normals, .> is strictly positive on the nonzero feasible cone:
    # equality would force all six wall products zero, hence the vector zero.
    # Therefore h=1 is bounded and all its vertices occur in the 20 triples.
    chart = [dual[0], -dual[1], -dual[2], -dual[3]]
    vertices, trials = {}, []
    for ids in itertools.combinations(range(6), 3):
        A = s.Matrix([[-normals[i][0], *normals[i][1:]] for i in ids] + [chart])
        det = S(A.det())
        if not det:
            trials.append(dict(support=list(ids), status='singular_section_frame'))
            continue
        v = tuple(S(x) for x in A.inv()*s.Matrix([0, 0, 0, 1]))
        products = [dot(v, n) for n in normals]
        bad = [i for i, g in enumerate(products) if sign(g) > 0]
        if bad:
            trials.append(dict(support=list(ids), status='infeasible', blockers=bad))
            continue
        w = normchart(v)
        active = tuple(i for i, n in enumerate(normals) if not dot(w, n))
        vertices[w] = active
        trials.append(dict(support=list(ids), status='feasible', active=list(active), ray=enc(w)))

    expected = {
        'A': (1, -s.Rational(1, 2), -s.Rational(1, 2), 0),
        'I1': (1, -s.Rational(1, 2), r3/2, 0),
        'I2': (1, r3/2, -s.Rational(1, 2), 0),
        'Iinf': (1, 0, 0, 1),
        'F': (s.Rational(3, 2), r3/2, r3/2, s.Rational(1, 2)),
        'G': (1, r3/2, r3/6, 0),
        'H': (1, r3/6, r3/2, 0),
    }
    expected = {name: normchart(v) for name, v in expected.items()}
    assert set(vertices) == set(expected.values()) and len(vertices) == 7
    active = {name: vertices[v] for name, v in expected.items()}
    finite_names = [name for name, v in expected.items() if sign(dot(v, v)) < 0]
    ideal_names = [name for name, v in expected.items() if not dot(v, v)]
    assert sorted(finite_names) == ['A', 'F', 'G', 'H']
    assert sorted(ideal_names) == ['I1', 'I2', 'Iinf']
    faces = {names[i]: sorted(name for name, a in active.items() if i in a) for i in range(6)}
    assert faces == {
        'B': ['A', 'G', 'H', 'I1', 'I2'],
        'L': ['A', 'I1', 'Iinf'], 'D': ['A', 'I2', 'Iinf'],
        'R': ['F', 'G', 'I2', 'Iinf'], 'U': ['F', 'H', 'I1', 'Iinf'],
        'C': ['F', 'G', 'H'],
    }
    # Positive sums give strict relative-interior facet witnesses.
    facet_witnesses = {}
    for i, name in enumerate(names):
        w = tuple(S(sum(expected[v][j] for v in faces[name])) for j in range(4))
        assert sign(dot(w, w)) < 0
        assert not dot(w, normals[i])
        assert all(sign(dot(w, n)) < 0 for j, n in enumerate(normals) if j != i)
        facet_witnesses[name] = enc(w)

    edges = set()
    for i, j in itertools.combinations(range(6), 2):
        ends = tuple(sorted(name for name, a in active.items() if i in a and j in a))
        assert len(ends) <= 2
        if len(ends) == 2:
            edges.add(ends)
    assert len(edges) == 11 and 7-len(edges)+6 == 2
    assert sorted(v for e in edges if 'A' in e for v in e if v != 'A') == ['I1', 'I2', 'Iinf']
    mixed_edges = [e for e in sorted(edges) if sum(v in finite_names for v in e) == 1]
    t2 = [f for f, vs in faces.items() if len(vs) == 3 and sum(v in ideal_names for v in vs) == 2]
    assert t2 == ['L', 'D'] and len(mixed_edges) == 6

    cusps = []
    for name in ideal_names:
        ids = list(active[name]); H = G.extract(ids, ids)
        for k in range(1, len(ids)+1):
            for jj in itertools.combinations(range(len(ids)), k):
                assert sign(H.extract(jj, jj).det()) >= 0
        assert H.rank() == 2
        if name != 'Iinf':
            ker = H.nullspace(); assert len(ker) == 1
            z = ker[0]
            if sign(z[0]) < 0: z = -z
            assert all(sign(x) > 0 for x in z)
            relation = tuple(S(sum(z[k]*normals[i][j] for k, i in enumerate(ids))) for j in range(4))
            assert sign(relation[0]) > 0 and normchart(relation) == expected[name]
            links = 'bounded Euclidean (2,3,6) triangle'
        else:
            assert ids == [1, 2, 3, 4]
            assert G[1, 3] == G[2, 4] == -2
            assert all(G[i, j] == 0 for i in (1, 3) for j in (2, 4))
            for ii in ((1, 3), (2, 4)):
                relation = tuple(S(sum(normals[i][j] for i in ii)) for j in range(4))
                assert normchart(relation) == expected[name]
            links = 'bounded Euclidean rectangle, affine A1 x affine A1'
        cusps.append(dict(vertex=name, active=[names[i] for i in ids], rank=2, link=links))

    paths = [f'scripts/{TAG}.py', 'results/fresh_multi_ideal_triangle_scope_20261017.md']
    return dict(
        status='PASS_EXACT_FINITE_VERTEX_T2_PORT_COUNTEREXAMPLE',
        scope='Finite-volume H3 Coxeter example refutes unique simplex-port blocks at a finite vertex for T2; neither proposed numerical T2 inequality is refuted.',
        sources={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        field='Q(sqrt2,sqrt3)', metric=[-1, 1, 1, 1], normal_labels=names,
        normals=[enc(v) for v in normals], gram=[enc(list(G.row(i))) for i in range(6)],
        pair_relations=relations, strict_center=enc(center), strict_dual=enc(dual),
        triple_trials=trials, vertices={k: enc(v) for k, v in expected.items()},
        vertex_norms={k: str(dot(v, v)) for k, v in expected.items()},
        active_sets={k: [names[i] for i in a] for k, a in active.items()},
        faces=faces, facet_witnesses=facet_witnesses, edges=[list(e) for e in sorted(edges)],
        cusps=cusps, counts=dict(facets=6, vertices=7, finite=4, ideal=3, exterior=0,
                               edges=11, T2=2, E_fi=6),
        obstruction=dict(finite_vertex='A', triangles=t2, shared_edge=['A', 'Iinf'],
                         third_two_face='B', third_two_face_sides=5,
                         missing_port_pair=['I1', 'I2']),
        distinction_from_prior='The prior T2 square pyramid fails at a simple ideal vertex along an ideal-ideal edge; this example fails at a finite vertex along a finite-ideal edge.',
    )


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--verify', action='store_true')
    args = ap.parse_args(); data = run(); path = ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert json.loads(path.read_text()) == data
    else:
        path.write_text(json.dumps(data, sort_keys=True, indent=2)+'\n')
    print(data['status'], 'saved_replay=', args.verify, data['counts'], flush=True)
