#!/usr/bin/env python3
"""Independent baseline optimization and twenty explicit marked edge witnesses."""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from math import comb

import fresh_h9_beta131_h300_symmetry_20261016_review as I
import fresh_sqrt5_forest_sensitive_weights_20261015_review as C

ROOT, K = I.ROOT, I.K
TAG = 'fresh_sqrt5_marked_flow_20261016_review'
BASE = 'results/fresh_sqrt5_marked_facet_barrier_20261016.json'
FLOW = 'results/fresh_sqrt5_marked_maximum_transitions_20261016.json'


def distances(g, start):
    d, todo = {start: 0}, [start]
    for i in todo:
        for j in g[i]:
            if j not in d:
                d[j] = d[i]+1
                todo.append(j)
    return d


def rooted_tree(g, node, parent, facet, port):
    return (int(node == facet), int(node == port),
            tuple(sorted((label, rooted_tree(g, j, node, facet, port))
                         for j, label in g[node].items() if j != parent)))


def forest_code(g, facet, port):
    remaining, pieces = set(g), []
    while remaining:
        start = min(remaining)
        part = set(distances(g, start))
        remaining -= part
        assert sum(len(g[i]) for i in part) == 2*len(part)-2
        pieces.append(min(rooted_tree(g, i, None, facet, port) for i in part))
    return tuple(sorted(pieces))


def decode(v):
    return K(*map(F, v))


def run():
    base, flow = map(I.read, (BASE, FLOW))
    models = {name: g for name, g in C.models() if len(g) <= 21}
    assert len(models) == 68
    components = []
    for name, g in models.items():
        dist = {i: distances(g, i) for i in g}
        weights = {(i, j): F(dist[i][j] <= 2)+
                   F(name == 'H4' and (g[i].get(j) == 5 or dist[i][j] == 3))
                   for i, j in combinations(g, 2)}
        total = sum(weights.values(), F())
        debits = {i: sum((w for pair, w in weights.items() if i in pair), F()) for i in g}
        components.append((name, len(g), total, debits))
        old = next(row for row in base['finite_components'] if row['type'] == name)
        assert total == F(old['total'])
        assert debits == {row['root']: F(row['debit']) for row in old['marked_root_debits']}
    # A two-state knapsack: exactly zero or one component contains the mark.
    best = [[None, None] for _ in range(22)]
    best[0][0] = F()
    for rank in range(1, 22):
        for name, size, total, debits in components:
            if size > rank:
                continue
            for state in (0, 1):
                previous = best[rank-size][state]
                if previous is not None:
                    value = previous+total
                    if best[rank][state] is None or value > best[rank][state]:
                        best[rank][state] = value
            previous = best[rank-size][0]
            if previous is not None:
                value = previous+total-min(debits.values())
                if best[rank][1] is None or value > best[rank][1]:
                    best[rank][1] = value
    assert best[21][1] == 38
    maximizers = []
    for name, size, total, debits in components:
        for root, debit in debits.items():
            if total-debit+best[21-size][0] == 38:
                maximizers.append([name, root])
    assert maximizers == [['A1', 0], ['D21', 19]]
    choose2 = lambda n: comb(n, 2) if n >= 2 else 0
    minimum = min(F(choose2(j)+choose2(20-j), 2) for j in range(21))
    assert minimum == 45 and 5*minimum-comb(20, 2) == 35
    assert F(choose2(20), 2) > minimum
    # Reconstruct the chosen local witnesses directly from literal22x22 Grams.
    g = {0: {}} | {i+1: {j+1:m for j,m in row.items()} for i,row in models['D20'].items()}
    A = [[K(2 if i == j else -int(j in g[i])) for j in range(21)] for i in range(21)]
    assert len(I.pivots(A)) == 21
    J, determinant = I.inverse(A)
    assert determinant > 0
    records = {row['port']: row for row in flow['port_records'] if row['source_marked_type'] == 'A1+D20'}
    assert set(records) == set(range(1, 21))
    verified = []
    for selected in flow['baseline38_balanced_port_flow']:
        port = selected['source_port']
        rec = records[port]['complete_marked_reattachments'][selected['reattachment_index']]
        scalar = decode(selected['g'])
        assert -14 < scalar <= 0 and (scalar < -2 or scalar == -K(0, 1))
        attachments = rec['attachments']
        assert 0 not in attachments and port not in attachments
        end = {i: {j:m for j,m in row.items() if j != port} for i,row in g.items() if i != port}
        end[21] = {}
        for j in attachments:
            end[j][21] = end[21][j] = 3
        assert forest_code(g, 0, port) == forest_code(end, 0, 21)
        b = [scalar if i == port else K(-int(i in attachments)) for i in range(21)]
        schur = 2-I.dot(b, I.mv(J, b))
        assert schur < 0 and schur.conjugate() > 0
        B = [row+[b[i]] for i,row in enumerate(A)]+[b+[K(2)]]
        assert all(x.is_integral for row in B for x in row)
        assert all(B[i][j] <= 0 for i in range(22) for j in range(i))
        rays = [[decode(x) for x in ray] for ray in rec['actual_edge']['endpoint_rays']]
        assert len(rays) == 2
        squares = []
        for omitted, ray in zip((port, 21), rays):
            products = I.mv(B, ray)
            assert products == [-K(int(i == omitted)) for i in range(22)]
            square = I.dot(ray, products)
            assert square < 0 and products[0] == 0
            squares.append(square)
        mutual = I.dot(rays[0], I.mv(B, rays[1]))
        assert mutual < 0
        assert F(selected['undirected_self_type_edge_mass'])*2 == F(selected['induced_port_incidence']) == 1
        verified.append(dict(port=port, framing_product=I.enc(scalar),
                             physical_last_schur=I.enc(schur),
                             conjugate_last_schur=I.enc(schur.conjugate()),
                             endpoint_squares=I.enc(squares), future_mutual_pairing=I.enc(mutual)))
    assert sorted(row['port'] for row in verified) == list(range(1, 21))
    paths = [f'scripts/{TAG}.py', BASE, FLOW,
             'scripts/fresh_h9_beta131_h300_symmetry_20261016_review.py',
             'scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py',
             'src/coxeter_polytopes/exact/golden.py']
    return dict(status='PASS_INDEPENDENT_MARKED_BASELINE_AND_TWENTY_LOCAL_FLOW_WITNESSES',
                sources={p:I.sha(p) for p in paths}, marked_baseline_maximum='38',
                facet_required_threshold='35', maximizing_marked_components=maximizers,
                local_edge_witnesses=verified, flow_mass_per_port='1/2',
                scope='Independent BFS baseline, two-state forest optimization, distinct-color tree codes, direct Gram signatures and forty finite endpoint rays for twenty selected narrow same-port returns. Does not independently enumerate all other reattachments or D21 forbidden ports. Local witnesses and flow need not have a common geometric realization. No dimension or finite-volume claim.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    target = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(target.read_text())
    else:
        target.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'])
