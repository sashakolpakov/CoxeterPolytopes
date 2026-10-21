#!/usr/bin/env python3
"""Fixed27-wall exact extension; inherited26-wall PASS, no root search.

Run through run_bounded_research_check_20261018.py --parent-lock -- THIS
--write or --verify. Imported26-wall pure helpers never invoke check().
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
from check_h19_twenty_six_wall_seed_20261019 import (
    Q, phi_number as N, dot, mv, add, digest, canonical_digest, inertia, inverse)

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT/'results/h19_twenty_seven_wall_seed_20261019.json'
PINS = {
    'results/h19_twenty_six_wall_seed_20261019.json':
    'e310723817ce78bce5c8ed1c11f2a210925d6edc3c6fcb8faf443d5458fc011b',
    'scripts/check_h19_twenty_six_wall_seed_20261019.py':
    'f0ac1d28225a37be17e9b61f2e436ea2fc3b1511efb32b02212b33539fecb363',
    'scripts/check_arithmetic_dn_opposed_fork_models_20261019.py':
    '185d20648368696429b675c1b6c023d30bb1be262ff83d3622e4f48801e1bae9',
    'data/candidates/h19_vinberg_prefix.json':
    '91c16b581e5763d6a283b6818f0e087e941fbb987b6b05a7e283febf1e32d0f1',
    'results/h19_omega2_finite_profiles_20261019.json':
    '46ab5c5d83415941362c8ab22df5a42d7c96ed595069f01a8f5e3597ae445813',
    'scripts/search_h19_omega2_finite_profiles_20261019.py':
    '02f82e5ddfb72486593a761d8a1ca86fef8e5ae3ffaaa72e901dccf129967b3e',
    'results/h19_omega2_finite_profiles_completion_20261019.md':
    '5568426ae65b447a9aef5ea97b943499e6ac908de0bcdcc037efb7ffccd8043a',
}


def decode(pair):
    return Q(F(pair[0]), F(pair[1]))


def check():
    for name, expected in PINS.items():
        assert digest(ROOT/name) == expected, name
    prior = json.loads((ROOT/'results/h19_twenty_six_wall_seed_20261019.json').read_text())
    found = json.loads((ROOT/'results/h19_omega2_finite_profiles_20261019.json').read_text())
    assert prior['status'] == 'PASS_EXACT_26_WALL_LOCAL_ALGEBRA'
    assert prior['walls'] == 26 and prior['wall_pairs'] == 325
    assert found['status'] == 'PASS_BOUNDED_OMEGA2_FINITE_PROFILE_SEARCH'
    for record in (prior, found):
        for name, expected in record['sources'].items():
            assert digest(ROOT/name) == expected, name
    phi, D = N(0, 1), N(1, 5)
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
    for name, encoded in prior['new_covectors'].items():
        roots[name] = mv(gi, [decode(v) for v in encoded])
    assert len(roots) == 26 and set(roots) == set(prior['norms'])
    candidates = found['survivors']
    small = next(s for s in candidates if s['x_phi_coefficients'] == [1, 3])
    large = next(s for s in candidates if s['x_phi_coefficients'] == [44, 72])
    assert len(candidates) == 2 and small['profile'] == 'profile_44'
    roots['q'] = [decode(v) for v in small['coordinates']]
    assert roots['q'] == [N(1, 3), N(-1, -2)]+[N(a, b) for a, b in [
        (7, 10), (14, 20), (20, 30), (16, 24), (12, 18), (8, 12), (4, 6), (10, 15)]]
    cov = {name: mv(G, root) for name, root in roots.items()}
    norms = {name: dot(root, cov[name]) for name, root in roots.items()}
    assert cov['q'] == [N(-8, -11), Q(1), Q(), Q(1), phi]+[Q()]*5
    assert cov['q'] == [decode(v) for v in small['covector']]
    mutual = dot([decode(v) for v in large['coordinates']], cov['q'])
    assert mutual == 4+5*phi and mutual.sign() > 0
    for name, root in roots.items():
        expected = decode(prior['norms'][name]) if name != 'q' else Q(2)
        assert norms[name] == expected and expected.sign() > 0
        assert all(v.integral() for v in root)
        assert all((2*v/norms[name]).integral() for v in cov[name])
    categories, pairs, qproducts = Counter(), [], {}
    for left, right in combinations(roots, 2):
        value = dot(roots[left], cov[right])
        if value == 0:
            kind = 'right'
        elif norms[left] == 2 and norms[right] == 2:
            if value in (Q(-1), -phi):
                kind = 'angle3' if value == -1 else 'angle5'
            else:
                assert (value+2).sign() <= 0, (left, right)
                kind = 'parallel' if value == -2 else 'ultraparallel'
        else:
            assert value.sign() < 0 and (value*value-norms[left]*norms[right]).sign() > 0
            kind = 'ultraparallel'
        categories[kind] += 1
        pairs.append([left, right, value.encoded(), kind])
        if right == 'q':
            qproducts[left] = [value.encoded(), kind]
    assert len(pairs) == 351 and qproducts == small['wall_products']
    assert categories == {'angle3': 34, 'angle5': 26, 'right': 183, 'ultraparallel': 108}
    p = [Q(2), -phi]+[Q()]*8
    center = add(p, [Q(), Q(F(1, 10))]+[sum(row, Q())/10 for row in ci])
    evaluate = lambda point: {name: dot(point, v) for name, v in cov.items()}
    cp = evaluate(center)
    assert dot(center, mv(G, center)).sign() < 0 and dot(p, mv(G, center)).sign() < 0
    assert all(v.sign() < 0 for v in cp.values())
    assert cp['q'] == -(68+95*phi)/10
    for name, root in roots.items():
        witness = add(center, root, -cp[name]/norms[name])
        values = evaluate(witness)
        assert dot(witness, mv(G, witness)).sign() < 0 and dot(center, mv(G, witness)).sign() < 0
        assert values[name] == 0 and all(v.sign() < 0 for n, v in values.items() if n != name)
    ps = add(p, roots['s'], -cov['s'][0])
    tau = (11-2*phi)/10
    pq = add(p, [Q(), Q()]+ci[2], tau)
    assert tau == decode(small['endpoint_parameter'])
    assert dot(pq, mv(G, pq)) == F(71, 2)-22*phi
    vertices = {}
    white = ['old0', 'old1', 'old2']+[f'old{i}' for i in range(4, 9)]
    for name, point, expected in [
        ('p', p, {f'old{i}' for i in range(9)}),
        ('p_s', ps, {'old0', 'old1', 's'} | {f'old{i}' for i in range(3, 9)}),
        ('p_q', pq, set(white+['q'])),
    ]:
        values = evaluate(point)
        active = sorted(n for n, v in values.items() if v == 0)
        assert set(active) == expected and all(v.sign() <= 0 for v in values.values())
        square = dot(point, mv(G, point))
        assert square.sign() < 0 and dot(p, mv(G, point)).sign() < 0
        assert inertia([[dot(roots[a], cov[b]) for b in active] for a in active]) == [9, 0, 0]
        assert values['r'].sign() < 0 and values['w'].sign() < 0
        vertices[name] = {'active': active, 'square': square.encoded(), 'q_product': values['q'].encoded()}
    # The nine original p-rays are local diagnostics, never a final-cone census.
    directions = [[Q(), Q(1)]+[Q()]*8]+[[Q(), Q()]+row for row in ci]
    first_hits, edge_kinds, pvalues = [], Counter(), evaluate(p)
    for i, direction in enumerate(directions):
        slopes = evaluate(direction)
        assert slopes[f'old{i}'].sign() < 0
        assert all(slopes[f'old{j}'] == 0 for j in range(9) if j != i)
        hits = [(-pvalues[n]/slope, n) for n, slope in slopes.items()
                if pvalues[n].sign() < 0 and slope.sign() > 0]
        assert hits
        first = hits[0][0]
        for value, _ in hits[1:]:
            if (value-first).sign() < 0:
                first = value
        endpoint = add(p, direction, first)
        values, square = evaluate(endpoint), dot(endpoint, mv(G, endpoint))
        assert first.sign() > 0 and all(v.sign() <= 0 for v in values.values())
        kind = {-1: 'finite', 0: 'ideal', 1: 'spacelike'}[square.sign()]
        active = [n for n, v in values.items() if v == 0]
        if kind == 'finite':
            assert len(active) == 9
            assert inertia([[dot(roots[a], cov[b]) for b in active] for a in active]) == [9, 0, 0]
        edge_kinds[kind] += 1
        first_hits.append({'freed_wall': f'old{i}', 'parameter': first.encoded(),
                           'first_hit_walls': [n for v, n in hits if v == first],
                           'square': square.encoded(), 'kind': kind})
    escape = add(p, directions[4], (8+21*phi)/44)
    ep = evaluate(escape)
    square = dot(escape, mv(G, escape))
    assert all(v.sign() <= 0 for v in ep.values())
    assert ep['q'] == (5*phi-18)/11 and ep['q'].sign() < 0
    assert square == (1557-955*phi)/484 and square.sign() > 0
    assert dot(p, mv(G, escape)).sign() < 0
    return {'status': 'PASS_EXACT_27_WALL_LOCAL_ALGEBRA', 'walls': 27, 'wall_pairs': 351,
        'checker_sha256': digest(Path(__file__)), 'sources': PINS,
        'new_root': [v.encoded() for v in roots['q']], 'new_covector': [v.encoded() for v in cov['q']],
        'reflection_integrality_checked': True, 'wall_pair_categories': dict(categories),
        'wall_pairs_sha256': canonical_digest(pairs), 'new_wall_products': qproducts,
        'strict_center_products': {n: v.encoded() for n, v in cp.items()},
        'actual_facet_projection_witnesses': 27, 'finite_vertex_witnesses': vertices,
        'alternate_survivor_product': mutual.encoded(),
        'initial_p_edge_first_hits': first_hits, 'initial_p_edge_kinds': dict(edge_kinds),
        'spacelike_escape_square': square.encoded(),
        'spacelike_escape_products': {n: v.encoded() for n, v in ep.items()},
        'inherited_only': 'The26-wall PASS and the two complete disjoint41-vertex H8 facets. '
                          'The >=85 count uses those82 vertices plus p,p_s,p_q, not a new census.',
        'scope': 'Fixed27-wall algebra/projections/three finite vertices/nine initial p-ray diagnostics. '
                 'The feasible spacelike witness has the separate written infinite-volume interpretation. '
                 'No root or complete vertex census, maximal completion, nonreflectivity or dimension claim.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = check()
    encoded = json.dumps(result, indent=2, sort_keys=True)+'\n'
    assert len(encoded.encode()) < 1024**2
    if args.verify:
        assert json.loads(CERT.read_text()) == result, 'Saved certificate differs'
    else:
        with CERT.open('x') as stream:
            stream.write(encoded)
    print(json.dumps({n: result[n] for n in ('status', 'walls', 'wall_pairs', 'initial_p_edge_kinds')}))


if __name__ == '__main__':
    main()
