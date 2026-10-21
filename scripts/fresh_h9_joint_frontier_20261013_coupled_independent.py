#!/usr/bin/env python3
"""Separate Fraction/interval audit of the two PD7-face completion caps.

No cap-producer arithmetic is imported.  Positive normalization radicals are
eliminated by using projective normals and their exact, variable squared norms.
The earlier independent frontier engine supplies only polynomial arithmetic,
outward rational intervals, and a freshly replayed 27-wall census.
"""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import fresh_h9_joint_frontier_20261013_independent as I
import fresh_h9_joint_frontier_20261013_all_edges_independent as A

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_joint_frontier_20261013_coupled_independent'
SOURCE = ROOT/'results/fresh_h9_makarov_joint_completion_20261013_coupled_exact.json'
FIRST = ROOT/'results/fresh_h9_makarov_joint_completion_20261013_face_bridge.json'
TARGET = ROOT/'results/fresh_h9_makarov_joint_completion_20261013_old_bank_gate.json'
OUT = ROOT/'results'/f'{TAG}.json'
SUPPORT = [0, 2, 3, 5, 7, 8, 9]
TWO = (F(2), F(0))

def enc(q):
    return list(map(str, q))

def kvalue(q):
    assert not(set(q)-{0})
    return q.get(0, I.Z)

def finite_cliques(adj, size):
    stack = [((), set(range(len(adj))))]
    while stack:
        chosen, rest = stack.pop()
        if len(chosen) == size:
            yield chosen
            continue
        if len(chosen)+len(rest) < size:
            continue
        for j in sorted(rest):
            stack.append((chosen+(j,), {i for i in rest & adj[j] if i > j}))

def run():
    previous = A.run()
    d = json.loads(SOURCE.read_text())
    for p, h in d['sources'].items():
        assert I.sha(ROOT/p) == h
    base = json.loads(I.SOURCE.read_text())
    physical = json.loads(I.WALLS.read_text())
    all_edges = json.loads(A.SOURCE.read_text())
    Q = [[I.K(x) for x in row] for row in physical['form']]
    W = [[I.dec(x) for x in row] for row in base['normalized_normals']]
    c = [I.con(I.K(x)) for x in physical['strict_center']]
    old_count = len(W)
    assert old_count == 27 and len(d['new_caps']) == 2
    assert I.sign(I.pair(c, c, Q)) < 0
    E = [I.con((F(x), F(0))) for x in json.loads(TARGET.read_text())['escape']['point']]
    assert I.sign(I.pair(E, E, Q)) > 0
    assert I.sign(I.pair(E, c, Q)) < 0
    assert all(I.sign(I.pair(E, w, Q)) < 0 for w in W)
    norms = [TWO]*27
    centers = []
    cuts = []
    for row in d['new_caps']:
        v = [I.con(I.K(x)) for x in row['projective_normal']]
        q = kvalue(I.pair(v, v, Q))
        assert q == I.K(row['projective_norm']) and I.ks(q) > 0
        radical = I.K(row['radicand'])
        assert I.ks(radical) > 0 and I.km(radical, q) == TWO
        assert row['normal'] == [[['0', '0'], enc(kvalue(x))] for x in v]
        center = kvalue(I.pair(v, c, Q))
        cut = kvalue(I.pair(v, E, Q))
        assert I.ks(center) < 0 and I.ks(cut) > 0
        centers.append(enc(center)); cuts.append(enc(cut))
        W.append(v); norms.append(q)
    assert cuts == d['target_projective_pairings']
    first = json.loads(FIRST.read_text())
    # The saved first cap is bound independently of the second producer.
    first_row = first['exact_cap'] if 'exact_cap' in first else first
    if 'projective_normal' not in first_row:
        first_row = first['cap']
    assert first_row['projective_normal'] == d['new_caps'][0]['projective_normal']
    CV = [I.cov(w, Q) for w in W]
    G = [[None]*29 for _ in W]
    adj = [set() for _ in W]
    for i, q in enumerate(norms):
        G[i][i] = q
    phi = (F(1, 2), F(1, 2))
    for row in physical['all_pair_gates']:
        if row['type'] == 'dotted':
            continue
        i, j = row['pair']
        G[i][j] = G[j][i] = {'m2': I.Z, 'm3': I.kn(I.ONE), 'm5': I.kn(phi)}[row['type']]
        adj[i].add(j); adj[j].add(i)
    counts = Counter()
    pair_types = {}
    for j in [27, 28]:
        for i in range(j):
            p = I.dot(W[j], CV[i])
            if not p:
                tag = 'm2'
                G[i][j] = G[j][i] = I.Z
                adj[i].add(j); adj[j].add(i)
            else:
                assert I.sign(p) < 0
                assert I.sign(I.sub(I.mul(p, p), I.con(I.km(norms[i], norms[j])))) > 0
                tag = 'dotted'
            pair_types[i, j] = tag
            counts[tag] += 1
    assert dict(counts) == d['new_pair_counts'] == {'m2': 15, 'dotted': 40}
    assert {(r['walls'][0], r['walls'][1]): r['type'] for r in d['new_pairs']} == pair_types
    assert adj[27] == set(SUPPORT+[28]) and adj[28] == set(SUPPORT+[27])
    assert all(I.sign(I.dot(c, x)) < 0 for x in CV)

    vertex = d['finite_vertex']
    p = [I.con(I.K(x)) for x in vertex['point']]
    pn = kvalue(I.pair(p, p, Q))
    assert pn == I.K(vertex['norm']) and I.ks(pn) < 0
    assert kvalue(I.pair(p, c, Q)) == I.K(vertex['center_pairing'])
    assert I.sign(I.pair(p, c, Q)) < 0
    incident = SUPPORT+[27, 28]
    assert incident == vertex['incident']
    signs = [I.sign(I.dot(p, x)) for x in CV]
    assert signs == vertex['all29_wall_signs']
    assert [i for i, s in enumerate(signs) if s == 0] == incident and max(signs) == 0
    M = [[G[i][j] for j in incident] for i in incident]
    assert I.positive(M)
    B = I.inverse(M)
    # Rebuild the corner as the orthogonal projection of the common center.
    cp = [kvalue(I.dot(c, CV[i])) for i in incident]
    coeff = [I.ksum(I.km(x, y) for x, y in zip(row, cp)) for row in B]
    rebuilt = [I.sub(c[a], I.sm(I.scale(t, W[i][a]) for t, i in zip(coeff, incident))) for a in range(10)]
    assert rebuilt == p
    kept = []
    for v in base['finite_vertices']:
        old = [I.dec(x) for x in v['point']]
        assert all(I.sign(I.dot(old, CV[j])) < 0 for j in [27, 28])
        kept.append(tuple(v['active']))
    nines = list(finite_cliques(adj, 9))
    pd9 = [ids for ids in nines if I.positive([[G[i][j] for j in ids] for i in ids])]
    assert len(nines) == 352 and set(pd9) == set(kept+[tuple(incident)])

    eights = list(finite_cliques(adj, 8))
    pd8 = [ids for ids in eights if I.positive([[G[i][j] for j in ids] for i in ids])]
    old_edges = {tuple(e['support']): e for e in all_edges['all_edge_supports']}
    added_edges = {tuple(i for i in incident if i != j) for j in incident}
    assert len(eights) == 1114 and set(pd8) == set(old_edges) | added_edges
    edge_classes = Counter()
    affine = 0
    for ids in pd8:
        inverse = I.inverse([[G[i][j] for j in ids] for i in ids])
        assert all(I.ks(x) >= 0 for row in inverse for x in row)
        ends = []
        for j in range(29):
            if j in ids or any(j not in adj[i] for i in ids):
                continue
            b = [G[i][j] for i in ids]
            proj = I.ksum(I.km(I.km(b[a], inverse[a][k]), b[k]) for a in range(8) for k in range(8))
            schur = I.ka(norms[j], I.kn(proj))
            if schur == I.Z:
                affine += 1
            elif I.ks(schur) > 0:
                assert tuple(sorted(ids+(j,))) in pd9
                ends.append(j)
        assert len(ends) in [0, 1, 2]
        tag = ['no_finite_ends', 'one_finite_end', 'two_finite_ends'][len(ends)]
        if ids in old_edges:
            assert tag == old_edges[ids]['kind']
        else:
            assert tag == 'one_finite_end'
        edge_classes[tag] += 1
    assert not affine
    assert edge_classes == {'two_finite_ends': 144, 'one_finite_end': 306, 'no_finite_ends': 5}
    # Exact points beyond the light cone on every newly created edge.
    continuations = []
    for j, released in enumerate(incident):
        u = [I.sm(I.scale(B[j][k], W[i][a]) for k, i in enumerate(incident)) for a in range(10)]
        assert not I.pair(p, u, Q)
        un = kvalue(I.pair(u, u, Q)); assert I.ks(un) > 0
        bound = I.km(I.kn(pn), I.kr(un))
        _, upper = I.sqrt_bounds(*I.interval_k(bound))
        for exponent in range(1, 100):
            t = upper*(1+F(1, 2**exponent))
            y = [I.sub(a, I.scale((t, F(0)), b)) for a, b in zip(p, u)]
            if I.sign(I.pair(y, c, Q)) >= 0:
                continue
            sy = [I.sign(I.dot(y, x)) for x in CV]
            if max(sy) <= 0:
                break
        else:
            raise AssertionError('no strict open-end witness')
        ids = [i for i in incident if i != released]
        assert [i for i, s in enumerate(sy) if s == 0] == ids
        yn = kvalue(I.pair(y, y, Q)); assert I.ks(yn) > 0
        continuations.append({'support': ids, 'released_wall': released, 'parameter': str(t), 'point': [enc(kvalue(x)) for x in y], 'norm': enc(yn), 'all29_wall_signs': sy})

    esc = d['escape']
    y = [I.con((F(x), F(0))) for x in esc['point']]
    yn = kvalue(I.pair(y, y, Q)); yc = kvalue(I.pair(y, c, Q))
    assert yn == I.K(esc['norm']) and I.ks(yn) > 0
    assert yc == I.K(esc['center_pairing']) and I.ks(yc) < 0
    assert [I.sign(I.dot(y, x)) for x in CV] == esc['all29_wall_signs'] == [-1]*29
    # Scale-invariant Cartan two-cycle.  Nonintegral trace over Q(sqrt5)
    # rules out algebraic integrality, independently of normalization.
    wpair = kvalue(I.dot(W[27], CV[1]))
    cycle = I.km((F(4), F(0)), I.km(I.km(wpair, wpair), I.kr(I.km(norms[27], norms[1]))))
    assert cycle == (F(1349500093, 480261426), F(94041207, 160087142))
    trace = 2*cycle[0]
    assert trace.denominator != 1
    source_paths = [SOURCE, FIRST, TARGET, I.SOURCE, I.WALLS, A.SOURCE, Path(I.__file__), Path(A.__file__), Path(__file__)]
    return {
        'status': 'PASS_SEPARATE_FRACTION_INTERVAL_29_WALL_GEOMETRY_AND_COMPLETE_EDGES',
        'sources': {str(f.relative_to(ROOT)): I.sha(f) for f in source_paths},
        'coupled_producer_dependencies_rechecked': d['sources'],
        'arithmetic': 'Independent stdlib Fraction polynomial identities and 160-bit outward rational intervals. Positive normalization radicals eliminated through variable-norm projective vectors. No cap producer imports.',
        'old_complete_frontier_replay': previous['status'],
        'walls': 29, 'old_pairs_replayed': 351, 'new_pair_counts': dict(counts),
        'new_projective_norms': [enc(q) for q in norms[27:]],
        'new_center_pairings': centers, 'old_bank_escape_cut_pairings': cuts,
        'all_original_vertices_strictly_retained': 65, 'finite_pair_nine_cliques': len(nines), 'finite_vertices': len(pd9),
        'new_finite_vertex_incident_walls': incident, 'new_vertex_norm': enc(pn),
        'finite_pair_eight_cliques': len(eights), 'actual_geodesic_edges': len(pd8),
        'edge_classification': dict(edge_classes), 'ideal_endpoint_candidates': affine,
        'finite_edge_graph_component_sizes': [65, 1], 'new_open_edge_continuations': continuations,
        'strict_spacelike_escape': esc,
        'cartan_two_cycle': {'walls': [1, 27], 'value': enc(cycle), 'trace': str(trace), 'trace_is_nonintegral': True},
        'scope': 'Exact 29-wall acute Coxeter cone, complete finite vertices and geodesic edges; still infinite volume. Any finite-volume Coxeter completion retaining walls 1 and 27 is nonarithmetic by the necessary integrality of Cartan cyclic products. No new finite-volume polytope is claimed.'
    }

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.verify:
        assert result == json.loads(OUT.read_text())
    else:
        OUT.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['finite_vertices'], result['edge_classification'], flush=True)
