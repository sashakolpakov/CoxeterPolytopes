#!/usr/bin/env python3
"""Fixed28-wall algebra and bounded finite-vertex walking; no root search.

Coordinator execution only, under the parent-lock bounded wrapper.
No producer/checker body is imported or called. Default walk cap64;
--walk-limit0 skips it. A closed queue is diagnostic, not a volume claim.
"""
import argparse
from collections import Counter, deque
from fractions import Fraction as F
import hashlib
from itertools import combinations
import json
from pathlib import Path
import signal
from search_h19_omega2_finite_profiles_20261019 import (
    Q, native as N, decode, dot, mv, digest, inverse, pair_kind)
from check_arithmetic_dn_opposed_fork_models_20261019 import inertia

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT/'results/h19_twenty_eight_wall_seed_20261019.json'
PINS = {
    'data/candidates/h19_vinberg_prefix.json':
    '91c16b581e5763d6a283b6818f0e087e941fbb987b6b05a7e283febf1e32d0f1',
    'results/h19_twenty_six_wall_seed_20261019.json':
    'e310723817ce78bce5c8ed1c11f2a210925d6edc3c6fcb8faf443d5458fc011b',
    'results/h19_twenty_seven_wall_seed_20261019.json':
    'b5b529b485a01b9315be6f31bbb122c170652273d3932b1eecece85680b5ce7d',
    'scripts/check_h19_twenty_seven_wall_seed_20261019.py':
    '4fed900c65b747dbf99226d16fd0e99b7e1f3106fd3269e0ebce9882455f9cf1',
    'results/h19_omega3_finite_profiles_20261019.json':
    'aa85c4c7c5373f66f124702563a55f9c88a7b7a28fe962f78db5bbe01af54c1d',
    'scripts/search_h19_omega3_finite_profiles_20261019.py':
    'da6be4a933eba17b74701b1aac0699123410427cab4440045e4822c285b636e0',
    'scripts/search_h19_omega2_finite_profiles_20261019.py':
    '02f82e5ddfb72486593a761d8a1ca86fef8e5ae3ffaaa72e901dccf129967b3e',
    'scripts/check_arithmetic_dn_opposed_fork_models_20261019.py':
    '185d20648368696429b675c1b6c023d30bb1be262ff83d3622e4f48801e1bae9',
    'results/h19_omega3_finite_profiles_completion_20261019.md':
    'e267d45c70bdcc0d5a9c054fe7d3933108b1568eaa2aa9f9235aca8a573c3445',
}


def add(a, b, scale=1):
    return [x+scale*y for x, y in zip(a, b)]


def encoded(vector):
    return [v.encoded() for v in vector]


def solve_pivoted(matrix, rhs):
    """Exact Gauss-Jordan with row pivoting; active covectors are not PD."""
    n = len(rhs)
    assert len(matrix) == n and all(len(row) == n for row in matrix)
    rows = [list(row)+[value] for row, value in zip(matrix, rhs)]
    for i in range(n):
        pivot = next((j for j in range(i, n) if rows[j][i] != 0), None)
        assert pivot is not None, 'Singular active affine system'
        rows[i], rows[pivot] = rows[pivot], rows[i]
        value = rows[i][i]
        rows[i] = [v/value for v in rows[i]]
        for j in range(n):
            if j != i:
                scale = rows[j][i]
                rows[j] = [a-scale*b for a, b in zip(rows[j], rows[i])]
    answer = [row[-1] for row in rows]
    assert all(dot(row, answer) == value for row, value in zip(matrix, rhs))
    return answer


def inputs():
    for name, expected in PINS.items():
        assert digest(ROOT/name) == expected, name
    prior = json.loads((ROOT/'results/h19_twenty_six_wall_seed_20261019.json').read_text())
    current = json.loads((ROOT/'results/h19_twenty_seven_wall_seed_20261019.json').read_text())
    found = json.loads((ROOT/'results/h19_omega3_finite_profiles_20261019.json').read_text())
    assert prior['status'] == 'PASS_EXACT_26_WALL_LOCAL_ALGEBRA'
    assert current['status'] == 'PASS_EXACT_27_WALL_LOCAL_ALGEBRA'
    assert found['status'] == 'PASS_BOUNDED_OMEGA3_FINITE_PROFILE_SEARCH'
    for record in (prior, current, found):
        for name, expected in record['sources'].items():
            assert digest(ROOT/name) == expected, name
    phi = N(0, 1)
    G = [[Q() for _ in range(10)] for _ in range(10)]
    G[0][0], G[0][1], G[1][0], G[1][1] = -2*phi, phi, phi, Q(2)
    for i in range(8):
        G[i+2][i+2] = Q(2)
    for i, j in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)]:
        G[i+2][j+2] = G[j+2][i+2] = Q(-1)
    assert inertia(G) == [9, 1, 0] and inertia(G, True) == [10, 0, 0]
    gi = inverse(G, Q)
    ci = [row[2:] for row in gi[2:]]
    old = json.loads((ROOT/'data/candidates/h19_vinberg_prefix.json').read_text())
    roots = {f'old{i}': [N(v['a'], v['b']) for v in row]
             for i, row in enumerate(old['output']['roots_in_vinberg_order'])
             if i not in (10, 20, 21)}
    for name, covector in prior['new_covectors'].items():
        roots[name] = mv(gi, [decode(v) for v in covector])
    roots['q'] = [decode(v) for v in current['new_root']]
    assert len(found['survivors']) == 1
    h = found['survivors'][0]
    assert h['x_phi_coefficients'] == [22, 36]
    roots['h'] = [decode(v) for v in h['coordinates']]
    assert len(roots) == 28
    cov = {n: mv(G, v) for n, v in roots.items()}
    norms = {n: dot(v, cov[n]) for n, v in roots.items()}
    assert cov['h'] == [N(-101, -163), Q(), Q(), Q(), Q(), N(10, 17), Q(1), Q(), Q(), Q()]
    assert cov['h'] == [decode(v) for v in h['covector']]
    for name, root in roots.items():
        expected = Q(2) if name in ('q', 'h') else decode(prior['norms'][name])
        assert norms[name] == expected and expected.sign() > 0
        assert all(v.integral() for v in root)
        assert all((2*v/expected).integral() for v in cov[name])
    p = [Q(2), -phi]+[Q()]*8
    c = add(p, [Q(), Q(F(1, 10))]+[sum(row, Q())/10 for row in ci])
    return phi, G, ci, roots, cov, norms, p, c, h


def bounded_walk(G, roots, cov, p, limit):
    if not limit:
        return {'stop': 'not_requested', 'global_finite_volume_proved': False}
    names = list(roots)
    evaluate = lambda v: {n: dot(v, cov[n]) for n in names}
    square = lambda v: dot(v, mv(G, v))
    active_names = lambda values: tuple(n for n in names if values[n] == 0)
    start = [v/p[0] for v in p]
    initial = active_names(evaluate(start))
    assert len(initial) == 9
    queue, seen = deque([(initial, start, 0)]), {initial: 0}
    vertices, edges, retained = [], [], 0
    stop, witness = None, None

    def vertex_record(active, point, ident):
        assert point[0] == 1 and square(point).sign() < 0
        assert len(active) == 9
        assert inertia([[dot(roots[a], cov[b]) for b in active] for a in active]) == [9, 0, 0]
        solved = [Q(1)]+solve_pivoted([cov[n][1:] for n in active], [-cov[n][0] for n in active])
        assert solved == point
        return {'id': ident, 'active': list(active), 'point': encoded(point), 'square': square(point).encoded()}

    first = vertex_record(initial, start, 0)
    vertices.append(first)
    retained += len(json.dumps(first).encode())
    while queue and stop is None:
        active, point, ident = queue.popleft()
        values = evaluate(point)
        assert all(v.sign() <= 0 for v in values.values())
        matrix = [cov[n][1:] for n in active]
        for freed_index, freed in enumerate(active):
            direction = [Q()]+solve_pivoted(matrix, [Q(-int(j == freed_index)) for j in range(9)])
            slopes = evaluate(direction)
            assert slopes[freed] == -1
            assert all(slopes[n] == 0 for n in active if n != freed)
            assert all(not (values[n] == 0 and slopes[n].sign() > 0) for n in names)
            hits = [(-values[n]/slopes[n], n) for n in names
                    if values[n].sign() < 0 and slopes[n].sign() > 0]
            if not hits:
                assert direction[0] == 0 and square(direction).sign() > 0
                assert all(v.sign() <= 0 for v in slopes.values())
                stop = 'unbounded_affine_ray'
                witness = {'from_vertex': ident, 'freed': freed, 'direction': encoded(direction),
                           'direction_square': square(direction).encoded(),
                           'all_slopes': {n: v.encoded() for n, v in slopes.items()}}
                break
            tau = hits[0][0]
            for value, _ in hits[1:]:
                if (value-tau).sign() < 0:
                    tau = value
            assert tau.sign() > 0
            endpoint = add(point, direction, tau)
            ev, sq = evaluate(endpoint), square(endpoint)
            assert endpoint[0] == 1 and all(v.sign() <= 0 for v in ev.values())
            hit = [n for value, n in hits if value == tau]
            record = {'from': ident, 'freed': freed, 'parameter': tau.encoded(), 'first_hit': hit}
            if sq.sign() >= 0:
                stop = 'spacelike_endpoint' if sq.sign() > 0 else 'ideal_endpoint'
                witness = record | {'endpoint': encoded(endpoint), 'square': sq.encoded(),
                                    'products': {n: v.encoded() for n, v in ev.items()}}
                break
            endpoint_active = active_names(ev)
            assert len(endpoint_active) == 9
            if endpoint_active not in seen:
                target = len(seen)
                new = vertex_record(endpoint_active, endpoint, target)
                retained += len(json.dumps(new).encode())
                assert retained < 512*1024, 'Bounded walk output limit reached'
                seen[endpoint_active] = target
                vertices.append(new)
                queue.append((endpoint_active, endpoint, target))
            else:
                target = seen[endpoint_active]
            record['to'] = target
            retained += len(json.dumps(record).encode())
            assert retained < 512*1024, 'Bounded walk output limit reached'
            edges.append(record)
            if len(seen) >= limit:
                stop = 'finite_vertex_cap'
                break
    if stop is None:
        assert not queue
        stop = 'queue_exhausted'
    return {'stop': stop, 'finite_vertex_limit': limit, 'finite_vertices_discovered': len(seen),
            'pending_queue': len(queue), 'vertices': vertices, 'directed_edges_checked': edges,
            'witness': witness, 'global_finite_volume_proved': False,
            'scope': 'Only the finite-vertex walk from p until the first listed stopping condition. '
                     'Queue closure alone is not promoted to compactness by this script.'}


def check(walk_limit):
    phi, G, ci, roots, cov, norms, p, c, h = inputs()
    evaluate = lambda v: {n: dot(v, cov[n]) for n in roots}
    square = lambda v: dot(v, mv(G, v))
    cp = evaluate(c)
    assert square(c).sign() < 0 and dot(p, mv(G, c)).sign() < 0
    assert all(v.sign() < 0 for v in cp.values())
    assert cp['h'] == -(418+695*phi)/5
    categories, pairs, hproducts = Counter(), [], {}
    for left, right in combinations(roots, 2):
        value = dot(roots[left], cov[right])
        kind = pair_kind(value, norms[right], phi) if norms[left] == 2 else None
        if norms[left] != 2:
            if value == 0:
                kind = 'right'
            elif norms[right] == 2:
                kind = pair_kind(value, norms[left], phi)
            else:
                assert value.sign() < 0 and (value*value-norms[left]*norms[right]).sign() > 0
                kind = 'ultraparallel'
        assert not kind.startswith('reject'), (left, right, kind)
        categories[kind] += 1
        pairs.append([left, right, value.encoded(), kind])
        if right == 'h':
            hproducts[left] = [value.encoded(), kind]
    assert len(pairs) == 378 and hproducts == h['wall_products']
    assert categories == {'angle3': 36, 'angle5': 26, 'right': 190, 'ultraparallel': 126}
    for name, root in roots.items():
        witness = add(c, root, -cp[name]/norms[name])
        values = evaluate(witness)
        assert square(witness).sign() < 0 and dot(c, mv(G, witness)).sign() < 0
        assert values[name] == 0 and all(v.sign() < 0 for n, v in values.items() if n != name)
    ps = add(p, roots['s'], -cov['s'][0])
    pq = add(p, [Q(), Q()]+ci[2], (11-2*phi)/10)
    ph = [decode(v) for v in h['endpoint']]
    assert dot(pq, cov['h']) == (1-52*phi)/5
    vertices = {}
    for name, point in [('p', p), ('p_s', ps), ('p_q', pq), ('p_h', ph)]:
        values, sq = evaluate(point), square(point)
        active = [n for n, v in values.items() if v == 0]
        assert len(active) == 9 and all(v.sign() <= 0 for v in values.values())
        assert sq.sign() < 0 and dot(p, mv(G, point)).sign() < 0
        assert inertia([[dot(roots[a], cov[b]) for b in active] for a in active]) == [9, 0, 0]
        assert values['r'].sign() < 0 and values['w'].sign() < 0
        vertices[name] = {'point': encoded(point), 'active': active, 'square': sq.encoded(),
                          'h_product': values['h'].encoded()}
    assert len({tuple(tuple(v) for v in item['point']) for item in vertices.values()}) == 4
    directions = [[Q(), Q(1)]+[Q()]*8]+[[Q(), Q()]+row for row in ci]
    initial_edges, pvalues = [], evaluate(p)
    for i, direction in enumerate(directions):
        slopes = evaluate(direction)
        assert slopes[f'old{i}'].sign() < 0
        assert all(slopes[f'old{j}'] == 0 for j in range(9) if j != i)
        hits = [(-pvalues[n]/slope, n) for n, slope in slopes.items()
                if pvalues[n].sign() < 0 and slope.sign() > 0]
        assert hits
        tau = hits[0][0]
        for value, _ in hits[1:]:
            if (value-tau).sign() < 0:
                tau = value
        point = add(p, direction, tau)
        values, sq = evaluate(point), square(point)
        active = [n for n, v in values.items() if v == 0]
        assert tau.sign() > 0 and sq.sign() < 0 and all(v.sign() <= 0 for v in values.values())
        assert len(active) == 9
        assert inertia([[dot(roots[a], cov[b]) for b in active] for a in active]) == [9, 0, 0]
        initial_edges.append({'freed_wall': f'old{i}', 'parameter': tau.encoded(),
                              'first_hit': [n for v, n in hits if v == tau], 'square': sq.encoded()})
    # Exact premises for the coordinator's separate bounded-slice proof.
    assert cov['old0'] == [-phi, Q(-2)]+[Q()]*8
    assert cov['r'][0] == 0 and cov['r'][1].sign() > 0 and all(v == 0 for v in cov['r'][2:])
    for j in range(8):
        assert roots[f'old{j+1}'] == [Q(), Q()]+[Q(-int(k == j)) for k in range(8)]
    assert roots['old9'][0] == 1 and roots['old9'][1] == 0
    assert all(v.sign() > 0 for v in roots['old9'][2:])
    walk = bounded_walk(G, roots, cov, p, walk_limit)
    pair_hash = hashlib.sha256(json.dumps(pairs, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    return {'status': 'PASS_EXACT_28_WALL_AND_BOUNDED_EDGE_DIAGNOSTIC',
            'checker_sha256': digest(Path(__file__)), 'sources': PINS,
            'walls': 28, 'wall_pairs': 378, 'wall_pair_categories': dict(categories),
            'wall_pairs_sha256': pair_hash, 'new_wall_products': hproducts,
            'actual_facet_projection_witnesses': 28, 'finite_vertex_witnesses': vertices,
            'initial_p_edge_first_hits': initial_edges, 'initial_p_finite_edges': 9,
            'slice_x0_only_zero_premises_checked': True,
            'bounded_edge_walk': walk, 'global_finite_volume_proved': False,
            'inherited_only': 'The two full compact disjoint41-vertex H8 facets. '
                              'Their82 vertices plus four checked off-cap vertices give >=86; no census.',
            'scope': 'Fixed28-wall algebra,28 actual facets,four finite vertices,nine initial edges, '
                     'and the explicitly capped walk. No complete root census, volume or dimension claim.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--verify', action='store_true')
    parser.add_argument('--walk-limit', type=int, default=64)
    args = parser.parse_args()
    assert args.walk_limit == 0 or 2 <= args.walk_limit <= 64
    signal.alarm(45)
    result = check(args.walk_limit)
    payload = json.dumps(result, sort_keys=True, indent=2)+'\n'
    assert len(payload.encode()) < 1024**2
    if args.verify:
        assert json.loads(CERT.read_text()) == result, 'Saved certificate differs'
    else:
        with CERT.open('x') as stream:
            stream.write(payload)
    signal.alarm(0)
    print(json.dumps({'status': result['status'], 'walls': 28, 'initial_finite_edges': 9,
                      'walk_stop': result['bounded_edge_walk']['stop']}))


if __name__ == '__main__':
    main()
