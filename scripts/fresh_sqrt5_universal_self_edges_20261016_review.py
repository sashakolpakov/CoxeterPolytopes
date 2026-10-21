#!/usr/bin/env python3
"""Independent project-field checks of universal local edge witnesses."""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import json
from math import comb

import fresh_h9_beta131_h300_symmetry_20261016_review as I
import fresh_sqrt5_forest_sensitive_weights_20261015_review as C

K, ROOT = I.K, I.ROOT
TAG = 'fresh_sqrt5_universal_self_edges_20261016_review'
PRIMARY = 'results/fresh_sqrt5_universal_self_edges_20261016.json'
DUAL = 'results/fresh_sqrt5_dimension21_a4_source_orbits_20261016_cap2.json'
PHI, LAMBDA = K(0, 1), K(13, 21)


def dec(x):
    return K(F(x[0]), F(x[1]))


def vector(x):
    return list(map(dec, x))


def matrix(x):
    return list(map(vector, x))


def copied_edge(a, i):
    n = len(a)
    g = a[i][i]*(1-LAMBDA)
    b = [row+[a[j][i] if i != j else g] for j, row in enumerate(a)]
    b.append([a[i][j] if i != j else g for j in range(n)]+[a[i][i]])
    return b


def endpoints(a, j, i, rec):
    b = copied_edge(a, i)
    n, d, s = len(a), a[i][i], 1/j[i][i]
    tau = d*LAMBDA
    assert tau > 2*s > 0 and 0 < tau.conjugate() < 2*s.conjugate()
    assert b[i][-1] < -d
    assert all((2*b[u][v]/b[u][u]).is_integral for u in range(n+1) for v in range(n+1))
    rays = list(map(vector, rec['endpoint_rays']))
    squares = []
    for ray, omitted in zip(rays, (i, n)):
        values = I.mv(b, ray)
        assert values == [-K(int(k == omitted)) for k in range(n+1)]
        square = I.dot(ray, values)
        assert square < 0
        squares.append(square)
    assert squares[0] == squares[1] == dec(rec['endpoint_square'])
    mutual = I.dot(rays[0], I.mv(b, rays[1]))
    assert mutual < 0 and mutual == dec(rec['endpoint_mutual_pairing'])
    assert mutual/squares[0] == dec(rec['edge_length_cosh']) > 1
    assert dec(rec['finite_inverse_diagonal']) == j[i][i]
    assert dec(rec['finite_residual_diagonal']) == s
    cycle = 4*b[i][-1]**2/d**2
    assert cycle == dec(rec['length_two_cyclic_product']) and cycle.is_integral and cycle.b
    return b


def run():
    data = I.read(PRIMARY)
    assert data['status'] == 'PASS_EXACT_UNIVERSAL_SAME_TYPE_EDGES_AND_BALANCED_DUAL_PORT_FLOWS'
    assert LAMBDA*LAMBDA.conjugate() == 1 and LAMBDA > 30
    expected = {name: graph for name, graph in C.models() if len(graph) <= 21}
    records = {r['type']: r for r in data['component_records']}
    assert set(records) == set(expected) and len(records) == 68
    cycles = {3: K(1), 4: K(2), 5: 1+PHI, 6: K(3), 10: 2+PHI}
    catalogue, maximum, root_count = {}, K(), 0
    for name, graph in expected.items():
        rec = records[name]
        a = matrix(rec['integral_gram'])
        n = len(graph)
        assert len(a) == n and all(len(row) == n for row in a)
        assert all(a[i][i] > 0 and a[i][i].conjugate() > 0 for i in range(n))
        assert all(x.is_integral for row in a for x in row)
        for i, j in combinations(range(n), 2):
            assert a[i][j] == a[j][i]
            if j not in graph[i]:
                assert a[i][j] == 0
            else:
                assert a[i][j] < 0
                assert 4*a[i][j]**2/(a[i][i]*a[j][j]) == cycles[graph[i][j]]
        I.pivots(a)
        I.pivots([[x.conjugate() for x in row] for row in a])
        inv, determinant = I.inverse(a)
        assert determinant > 0
        roots = {r['root']: r for r in rec['root_witnesses']}
        assert set(roots) == set(range(n))
        for i in range(n):
            endpoints(a, inv, i, roots[i])
            maximum = max(maximum, (a[i][i]*inv[i][i]/2).conjugate())
            root_count += 1
        catalogue[name] = a
        print('independent self-edge type', name, n, flush=True)
    assert root_count == 724 and maximum == 30
    # Count multisets by rank multiplicities using (1-z^r)^(-m_r).
    multiplicities = Counter(map(len, catalogue.values()))
    polynomial = [1]+[0]*21
    for rank, multiplicity in sorted(multiplicities.items()):
        nxt = [0]*22
        for old in range(22):
            for count in range((21-old)//rank+1):
                nxt[old+count*rank] += polynomial[old]*comb(multiplicity+count-1, count)
        polynomial = nxt
    assert polynomial == data['forest_count_dp'] and polynomial[21] == 311574
    a = [[K() for _ in range(21)] for _ in range(21)]
    offset = 0
    for name in ('D5', 'E8', 'E8'):
        block = catalogue[name]
        for i, row in enumerate(block):
            for j, value in enumerate(row):
                a[offset+i][offset+j] = value
        offset += len(block)
    assert a == matrix(data['exceptional_full_finite_gram'])
    inv, _ = I.inverse(a)
    exceptional = data['exceptional_all21_edge_and_escape_witnesses']
    assert [r['root'] for r in exceptional] == list(range(21))
    for i, rec in enumerate(exceptional):
        b = endpoints(a, inv, i, rec)
        escape = rec['strict_null_escape']
        u, w = vector(escape['u']), vector(escape['w'])
        bu, bw = I.mv(b, u), I.mv(b, w)
        U, W, V = I.dot(u, bu), I.dot(w, bw), I.dot(u, bw)
        assert U < 0 < W and all(x < 0 for x in bu+bw)
        D = V*V-U*W
        assert D > V*V >= 0 and D == dec(escape['radical_squared'])
        z0, z1 = vector(escape['null_ray_constant']), vector(escape['null_ray_radical'])
        assert z0 == [x-V*y/W for x, y in zip(u, w)]
        assert z1 == [y/W for y in w]
        assert I.dot(z0, I.mv(b, z0))+D*I.dot(z1, I.mv(b, z1)) == 0
        assert I.dot(z0, I.mv(b, z1)) == 0
    square_rows = data['canonical_copied_square_diagnostic']['cases_data']
    assert [r['roots'] for r in square_rows] == [list(x) for x in combinations(range(21), 2)]
    for row in square_rows:
        i, j = row['roots']
        determinant = inv[i][i]*inv[j][j]-inv[i][j]**2
        alpha, beta, gamma = inv[j][j]/determinant, -inv[i][j]/determinant, inv[i][i]/determinant
        block = [[2*alpha-a[i][i]*LAMBDA, 2*beta], [2*beta, 2*gamma-a[j][j]*LAMBDA]]
        assert block == matrix(row['symmetric_block'])
        assert block[0][0] < 0 and block[0][0]*block[1][1]-block[0][1]**2 > 0
    dual = I.read(DUAL)
    total, lower, flows = Counter(), F(0), []
    for index, term in enumerate(dual['certificate']['terms']):
        y = F(term['multiplier'])
        assert y <= 0
        lower += y*F(term['rhs'])
        for column, value in term['coefficients']:
            total[column] += y*F(value)
        if term['metadata']['kind'] == 'forest':
            types = term['metadata']['types']
            ports = [[c, name, i] for c, name in enumerate(types) for i in range(len(catalogue[name]))]
            assert len(ports) == 21
            flows.append((index, types, -y, -y/2, ports))
    assert all(F(i == dual['variables']) >= value for i, value in total.items())
    assert lower == F(dual['certificate']['lower_bound']) > 40
    saved_flows = data['latest_dual_port_flow']['flows']
    assert len(flows) == len(saved_flows) == 506
    for (index, types, mass, edge_mass, ports), saved in zip(flows, saved_flows):
        assert index == saved['dual_term'] and types == saved['vertex_type']
        assert F(saved['vertex_mass']) == mass == 2*edge_mass
        assert F(saved['self_edge_mass_per_port']) == edge_mass and ports == saved['ports']
    graph = data['simple_abstract_port_incidence']['edges']
    assert len(graph) == 231 and len({tuple(sorted((a,b))) for a,b,c in graph}) == 231
    for v in range(22):
        assert sorted(c for a,b,c in graph if v in (a,b)) == list(range(21))
    paths = [f'scripts/{TAG}.py', PRIMARY, DUAL,
             'scripts/fresh_h9_beta131_h300_symmetry_20261016_review.py',
             'scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py',
             'src/coxeter_polytopes/exact/golden.py']
    return dict(status='PASS_INDEPENDENT_UNIVERSAL_LOCAL_EDGES_AND_PORT_FLOW_REVIEW',
                sources={p:I.sha(p) for p in paths}, component_types=68,
                rooted_positions=root_count, rank21_forests=polynomial[21],
                explicit_rank21_strict_null_escapes=21, copied_square_failures=210,
                port_flow_equations=21*len(flows), dual_lower_bound=str(lower),
                scope='Independent Golden arithmetic, literal finite-Gram and endpoint checks, all-rank21 forest count,21 strict real null escapes,210 canonical copied-square failures and literal old-dual/port-flow checks. Does not reconstruct LP coefficient generation, realize simultaneous stars or arbitrary squares, or prove a new dimension bound.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    output = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(output.read_text())
    else:
        output.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'])
