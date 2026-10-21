#!/usr/bin/env python3
"""Independent quartic arithmetic and all-support census of the two H9 caps.

Uses the earlier independent Fraction implementation, never the primary
producer. Complete support enumeration replaces the primary cone clipping.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A

ROOT = Path(__file__).resolve().parents[1]
BASE = 'fresh_h9_mixed_d4_endpoint_pair_20261014'
TAG = BASE + '_review'
OLD = 'fresh_h9_mixed_d4_coupling_20261014'


def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()


def run():
    paths = [f'results/{BASE}.json', f'results/{OLD}.json',
             f'results/{OLD}_review.json', f'results/{OLD}_frontier.json']
    data, old, previous, frontier = [json.loads((ROOT / p).read_text()) for p in paths]
    for doc in [data, old, previous, frontier]:
        for p, h in doc['sources'].items():
            assert sha(p) == h, p
    assert previous['status'] == 'PASS_INDEPENDENT_QUARTIC_MIXED_D4_SUPPORT_CENSUS'
    M, oldN = A.dm(old['form']), A.dm(old['final_normals'])
    assert M == A.dm(data['form']) and A.inertia(M) == [9, 1]
    # Literal coefficients, independent of decoding the two output vectors.
    literal = [
        [('65/4', '29/4'), ('20', '9'), ('57/4', '27/4'), ('19/2', '9/2'), ('19/4', '9/4'),
         ('2', '1'), ('-9/2', '-3/2'), ('-4', '-3/2'), ('-11/4', '-3/4'), ('-9/4', '-3/4')],
        [('33/4', '15/4'), ('21/2', '9/2'), ('15/2', '3'), ('5', '2'), ('5/2', '1'),
         ('1/2', '1/2'), ('-3', '-1'), ('-2', '-1'), ('-3/2', '-1/2'), ('-3/2', '-1/2')]]
    caps = [[A.E(A.Q(a, b)) for a, b in row] for row in literal]
    assert caps == A.dm(data['added_normals'])
    N = oldN + caps
    assert N == A.dm(data['normals']) and len(N) == 14
    G = [[A.pair(v, w, M) for w in N] for v in N]
    assert G[12][12] == G[13][13] == 2 and not G[12][13]
    assert A.dec(data['mutual_product']) == G[12][13]
    assert all(not x.b for row in caps for x in row)
    rows = {tuple(r['pair']): r for r in data['new_pairs']}
    counts = Counter()
    for i, j in combinations(range(14), 2):
        kind, label = A.relation(G[i][j], G[i][i], G[j][j])
        assert kind in ['right', 'finite', 'parallel', 'dotted']
        if j < 12:
            continue
        rec = rows[i, j]
        assert rec['type'] == kind
        if label is not None:
            assert rec['label'] == label
        if 'squared_weight' in rec:
            assert A.dec(rec['squared_weight']) == G[i][j] * G[i][j] / (G[i][i] * G[j][j])
        counts[kind + (str(label) if label is not None else '')] += 1
    assert len(rows) == 25 and dict(counts) == data['new_pair_counts']
    oldE = A.dv(old['strict_escape']['vector'])
    oldcuts = [A.pair(oldE, cap, M) for cap in caps]
    assert all(g.sign() > 0 for g in oldcuts) and oldcuts == A.dv(data['old_escape_cut_pairings'])
    oldfinite = [A.dv(r['vector']) for r in old['final_complete_ray_census']['rays'] if r['kind'] == 'finite']
    assert len(oldfinite) == data['preserved_old_finite_vertex_count'] == 12
    assert all(A.pair(v, n, M).sign() < 0 for v in oldfinite for n in caps)
    h = A.dv(data['center'])
    expected_center = [sum((v[j] for v in oldfinite), A.Z) + A.dv(old['center'])[j] / 1000 for j in range(10)]
    assert h == expected_center
    assert A.pair(h, h, M) == A.dec(data['center_norm']) and A.pair(h, h, M).sign() < 0
    assert A.pair(h, A.dv(old['center']), M).sign() < 0
    assert all(A.pair(h, n, M).sign() < 0 for n in N)
    for i, n in enumerate(N):
        projection = [x - A.pair(h, n, M) * y / G[i][i] for x, y in zip(h, n)]
        assert projection == A.dv(data['actual_facet_points'][i])
        assert A.pair(projection, projection, M).sign() < 0 and A.pair(projection, h, M).sign() < 0
        assert all(A.pair(projection, w, M).sign() == (0 if i == j else -1) for j, w in enumerate(N))
    for index, cap in enumerate(caps):
        rec = data['targeted_edge_endpoints'][index]
        oldedge = frontier['queue'][rec['edge_id']]
        assert oldedge['kind'] == 'open_to_light_cone' and rec['support'] == oldedge['support']
        assert rec['start'] == oldedge['start_point']
        contact = dict(rec['contacts'])
        contact[rec['anchor']] = rec['anchor_label']
        for j in rec['support'] + [rec['anchor']]:
            g = G[12 + index][j]
            expect = {2: A.Z, 3: -A.ONE, 5: -A.PHI}[contact.get(j, 2)]
            assert g == expect
        startg = A.pair(A.dv(rec['start']), cap, M)
        assert startg.sign() < 0 and startg == A.dec(rec['start_pairing'])
    print('PASS_LITERAL_CAPS_ALL_PAIRS_CENTER_FACETS_OLD_VERTICES', flush=True)

    # Independent Gaussian nullspace enumeration of every14-choose9 support.
    # A strictly positive timelike sum of normals and spanning old roots certify
    # a bounded affine slice, so this exhausts the complete projective cone.
    A.inverse([A.mv(M, N[i]) for i in list(range(6)) + list(range(7, 11))])
    census = A.census(M, N, data['complete_ray_census'])
    assert census['tested_supports'] == 2002
    assert census['counts'] == {'finite': 15, 'spacelike': 91}
    print('PASS_INDEPENDENT_2002_SUPPORT_CENSUS', census['counts'], flush=True)
    allrows = data['complete_ray_census']['rays']
    vectors = [A.dv(row['vector']) for row in allrows]
    finite = [i for i, row in enumerate(allrows) if row['kind'] == 'finite']
    new = [i for i in finite if any(j >= 12 for j in allrows[i]['incident'])]
    assert len(new) == 3
    expected_new = {A.normkey(A.dv(row['vector'])) for row in data['new_finite_vertices']}
    assert {A.normkey(vectors[i]) for i in new} == expected_new
    assert {A.normkey(v) for v in oldfinite} <= {A.normkey(vectors[i]) for i in finite}
    for i in finite:
        ids = allrows[i]['incident']
        assert len(ids) == 9 and A.inertia([[G[a][b] for b in ids] for a in ids]) == [9, 0]
        assert A.pair(vectors[i], h, M).sign() < 0
    for rec in data['targeted_edge_endpoints']:
        hits = [i for i in finite if set(rec['support']) <= set(allrows[i]['incident'])
                and any(j >= 12 for j in allrows[i]['incident'])]
        assert len(hits) == 1
        assert allrows[hits[0]] == rec['finite_endpoint']

    # Recover actual adjacency from the complete ray supports, rather than
    # trusting that a shared rank-eight label set names an edge.
    rebuilt = []
    edge_counts = Counter()
    closed = set()
    for i in finite:
        active = allrows[i]['incident']
        for released in active:
            S = sorted(set(active) - {released})
            assert A.inertia([[G[a][b] for b in S] for a in S]) == [8, 0]
            others = [j for j, row in enumerate(allrows) if j != i and set(S) <= set(row['incident'])]
            assert len(others) == 1
            j = others[0]
            middle = [x + y for x, y in zip(vectors[i], vectors[j])]
            assert all(A.pair(middle, n, M).sign() == (0 if k in S else -1) for k, n in enumerate(N))
            kind = allrows[j]['kind']
            edge_counts[kind] += 1
            if kind == 'finite':
                assert A.pair(middle, middle, M).sign() < 0
                closed.add(tuple(sorted((i, j))))
            else:
                assert kind == 'spacelike'
                # Strict nonsupport point beyond the light crossing.
                for k in range(1, 80):
                    eps = F(1, 2**k)
                    z = [eps * x + (1 - eps) * y for x, y in zip(vectors[i], vectors[j])]
                    if A.pair(z, z, M).sign() > 0:
                        break
                else:
                    raise AssertionError('No open continuation')
                assert all(A.pair(z, n, M).sign() == (0 if k in S else -1) for k, n in enumerate(N))
            rebuilt.append(dict(start_ray=i, released_wall=released, support=S, other_cone_ray=j,
                                kind='finite_endpoint' if kind == 'finite' else 'open_to_light_cone'))
    assert rebuilt == data['finite_corner_directions']
    assert dict(edge_counts) == data['finite_corner_counts'] == {'spacelike': 71, 'finite': 64}
    assert len(closed) == 32
    adjacency = {i: set() for i in finite}
    for i, j in closed:
        adjacency[i].add(j)
        adjacency[j].add(i)
    seen = {min(finite)}
    todo = list(seen)
    for i in todo:
        for j in adjacency[i] - seen:
            seen.add(j)
            todo.append(j)
    assert seen == set(finite)
    new_closed = sorted(e for e in closed if set(e) & set(new))
    assert len(new_closed) == 4
    er = data['strict_escape']
    escape = A.dv(er['vector'])
    assert all(not x.b and not x.a.b for x in escape)
    q = A.pair(escape, escape, M)
    signs = [A.pair(escape, n, M) for n in N]
    assert q.sign() > 0 and q == A.dec(er['norm'])
    assert all(x.sign() < 0 for x in signs) and signs == A.dv(er['all_pairings'])
    assert A.pair(escape, h, M).sign() < 0 and A.pair(escape, h, M) == A.dec(er['center_pairing'])
    paths += [f'scripts/{BASE}.py', f'scripts/{TAG}.py', f'scripts/{OLD}_review.py']
    return A.enc(dict(status='PASS_INDEPENDENT_14_WALL_ENDPOINT_PAIR_AND_COMPLETE_CONE',
        sources={p: sha(p) for p in paths}, method='Independent Fraction quartic arithmetic; all2002 Gaussian support nullspaces instead of primary cone clipping. No primary producer imports.',
        new_pair_counts=dict(counts), literal_caps=caps, old_escape_cut_pairings=oldcuts,
        preserved_old_finite_vertices=12, actual_facets=14, complete_census=census,
        new_finite_ray_ids=new, new_finite_active_sets=[allrows[i]['incident'] for i in new],
        finite_corner_direction_counts=dict(edge_counts), unique_finite_edges=32,
        finite_edge_graph_connected=True, finite_edge_graph_cycle_rank=32 - 15 + 1,
        four_new_finite_edges=new_closed,
        strict_escape=dict(vector=escape, norm=q, center_pairing=A.pair(escape, h, M), all_pairings=signs),
        scope='Complete14-wall projective ray census, all finite vertices and finite-corner edges; no vertex-free-edge census or arbitrary completion theorem. Both caps cut the old escape, but the new rational strict spacelike escape proves infinite volume.'))


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
    print(result['status'], flush=True)
