#!/usr/bin/env python3
"""Exact local algebra for26 walls; no root search or global volume census.

Run ONLY through run_bounded_research_check_20261018.py --parent-lock --
scripts/check_h19_twenty_six_wall_seed_20261019.py --write (or --verify).
The wrapper owns the shared lock and all resource/deadline/RSS gates.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import combinations
import json
from pathlib import Path
from check_arithmetic_dn_opposed_fork_models_20261019 import field, inertia, inverse

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / 'results/h19_twenty_six_wall_seed_20261019.json'
Q = field(5)
PINS = {
    'scripts/check_arithmetic_dn_opposed_fork_models_20261019.py':
    '185d20648368696429b675c1b6c023d30bb1be262ff83d3622e4f48801e1bae9',
    'data/candidates/h19_vinberg_prefix.json':
    '91c16b581e5763d6a283b6818f0e087e941fbb987b6b05a7e283febf1e32d0f1',
    'data/known_examples/bugaenko_compact_h8.json':
    'c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c',
    'results/compact_h9_discriminant19_mixed_norm_seed_20261019.md':
    '70fa010732d837d5c044701b4c18ed5e68cb36c18cc0954dc7d7f369250a5c87',
    'results/compact_h9_discriminant19_polar_caps_research_20261019.md':
    '792634e935866ed78d95d0a93c04fdd7a8484f3605b5c50a66c1d376c9599f1b',
    'results/compact_h9_discriminant19_twentyfive_escape_20261019.md':
    '1cc74516c9e22b7150c64a95b4c4251dec9d3c02f238dbb4c946ef937e37c096',
    'results/compact_h9_discriminant19_nonpolar_pell_mirror_20261019.md':
    'b817f3a362144053d8d0ff8ca4f71de2c56faecf8b357ed1751fb802b4b130c4',
}


def phi_number(a=0, b=0):
    """Native a+b*phi represented in the helper's rational sqrt5 basis."""
    return Q(F(a) + F(b)/2, F(b)/2)


def power(value, exponent):
    answer = Q(1)
    for _ in range(exponent):
        answer *= value
    return answer


def dot(left, right):
    return sum((x*y for x, y in zip(left, right)), Q())


def mv(matrix, vector):
    return [dot(row, vector) for row in matrix]


def add(left, right, factor=Q(1)):
    return [x + factor*y for x, y in zip(left, right)]


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(65536), b''):
            result.update(chunk)
    return result.hexdigest()


def canonical_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True,
                                    separators=(',', ':')).encode()).hexdigest()


def check():
    for name, expected in PINS.items():
        assert digest(ROOT/name) == expected, name
    phi = phi_number(0, 1)
    b, D = 4 + phi, 1 + 5*phi
    powers = [power(phi, i) for i in range(10)]
    C = [[Q(2 if i == j else 0) for j in range(8)] for i in range(8)]
    for i, j in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)]:
        C[i][j] = C[j][i] = Q(-1)
    G = [[Q() for _ in range(10)] for _ in range(10)]
    G[0][0], G[0][1], G[1][0], G[1][1] = -2*phi, phi, phi, Q(2)
    for i in range(8):
        for j in range(8):
            G[i+2][j+2] = C[i][j]
    assert inertia(G) == [9, 1, 0] and inertia(G, True) == [10, 0, 0]
    with (ROOT/'data/candidates/h19_vinberg_prefix.json').open() as stream:
        native = json.load(stream)['output']['roots_in_vinberg_order']
    native = [[phi_number(v['a'], v['b']) for v in row] for row in native]
    assert len(native) == 22 and all(len(v) == 10 for v in native)
    assert all(dot(v, mv(G, v)) == 2 for v in native)
    roots = {f'old{i}': v for i, v in enumerate(native) if i not in (10, 20, 21)}
    zero = lambda: [Q() for _ in range(8)]
    coordinates = lambda pairs: [phi_number(a, b0) for a, b0 in pairs]
    omega0 = list(map(Q, [4, 7, 10, 8, 6, 4, 2, 5]))
    omega1 = list(map(Q, [7, 14, 20, 16, 12, 8, 4, 10]))
    assert mv(C, omega0) == [Q(1)] + [Q()]*7
    assert mv(C, omega1) == [Q(), Q(1)] + [Q()]*6
    assert dot(omega0, mv(C, omega0)) == 4
    assert dot(omega1, mv(C, omega1)) == 14
    u, t, M = 10 + 17*phi, 12 + 19*phi, 8*phi - 5
    assert 7*t*t-D*u*u == 1 and D*u-M*t == 3-phi
    assert 7*t-M*u == 2*(phi-1)
    roots['r'] = [Q(1), Q(2)] + zero()
    roots['k'] = [2+3*phi, Q()] + coordinates([
        (6, 10), (12, 20), (18, 30), (15, 24),
        (12, 18), (8, 12), (4, 6), (9, 16)])
    roots['w'] = [2*powers[4], -powers[5]] + [b*v for v in [
        2*powers[3], 4*powers[3], 6*powers[3], 5*powers[3],
        4*powers[3], 3*powers[3], powers[4], 3*powers[3]]]
    roots['a'] = [2*powers[2], -powers[3]] + coordinates([
        (6, 8), (12, 16), (18, 24), (15, 20),
        (12, 16), (9, 11), (5, 6), (9, 12)])
    roots['d'] = [2*powers[5], -powers[6]] + coordinates([
        (22, 36), (43, 72), (64, 108), (53, 89),
        (42, 70), (31, 51), (17, 28), (32, 55)])
    roots['z'] = [2*powers[4], -powers[5]] + [b*powers[2]*v for v in omega0]
    roots['s'] = [2*u, -phi*u] + [t*v for v in omega1]
    assert len(roots) == 26
    cov = {name: mv(G, v) for name, v in roots.items()}
    norm = {name: dot(v, cov[name]) for name, v in roots.items()}
    for name, vector in roots.items():
        assert norm[name] == (2*b if name in ('r', 'w', 'z') else 2), name
        assert all(v.integral() for v in vector), name
        assert all((2*v/norm[name]).integral() for v in cov[name]), name
    expected = {
        'r': [Q(), b] + zero(),
        'k': [-6-10*phi, 3+5*phi] + zero(),
        'w': [-b*powers[5], Q()] + zero(),
        'a': [-b*powers[3], Q()] + zero(),
        'd': [-b*powers[6], Q()] + zero(),
        'z': [-b*powers[5], Q()] + zero(),
        's': [-D*u, Q()] + zero(),
    }
    for name, entries in {
        'k': {4: Q(1), 7: 2*phi}, 'w': {5: b*phi, 6: b},
        'a': {4: phi, 5: Q(1), 6: powers[2]},
        'd': {0: Q(1), 5: 3+4*phi, 6: powers[5], 7: 2*phi},
        'z': {0: b*powers[2]}, 's': {1: t},
    }.items():
        for i, value in entries.items():
            expected[name][i+2] = value
    assert all(cov[name] == vector for name, vector in expected.items())
    pairs, categories = [], Counter()
    for left, right in combinations(roots, 2):
        value = dot(roots[left], cov[right])
        if value == 0:
            kind = 'right'
        elif norm[left] == 2 and norm[right] == 2:
            if value in (Q(-1), -phi):
                kind = 'angle3' if value == -1 else 'angle5'
            else:
                assert (value+2).sign() <= 0, (left, right)
                kind = 'parallel' if value == -2 else 'ultraparallel'
        else:
            assert value.sign() < 0
            assert (value*value-norm[left]*norm[right]).sign() > 0, (left, right)
            kind = 'ultraparallel'
        categories[kind] += 1
        pairs.append([left, right, value.encoded(), kind])
    assert len(pairs) == 325
    ci = inverse(C, Q)
    p = [Q(2), -phi] + zero()
    center = add(p, [Q(), Q(F(1, 10))] +
                 [sum(row, Q())/10 for row in ci])
    evaluate = lambda point: {name: dot(point, g) for name, g in cov.items()}
    cp = evaluate(center)
    assert dot(center, mv(G, center)).sign() < 0
    assert dot(p, mv(G, center)).sign() < 0
    assert all(v.sign() < 0 for v in cp.values())
    assert cp['s'] == -(808+1311*phi)/10
    for name, root in roots.items():
        witness = add(center, root, -cp[name]/norm[name])
        values = evaluate(witness)
        assert dot(witness, mv(G, witness)).sign() < 0
        assert dot(center, mv(G, witness)).sign() < 0
        assert values[name] == 0
        assert all(v.sign() < 0 for other, v in values.items() if other != name)
    ps = add(p, roots['s'], D*u)
    vertices = {}
    for name, point, expected_active in [
        ('p', p, [f'old{i}' for i in range(9)]),
        ('p_s', ps, ['old0', 'old1']+[f'old{i}' for i in range(3, 9)]+['s']),
    ]:
        values = evaluate(point)
        active = [wall for wall, v in values.items() if v == 0]
        assert active == expected_active and all(v.sign() <= 0 for v in values.values())
        square = dot(point, mv(G, point))
        assert square.sign() < 0 and dot(p, mv(G, point)).sign() < 0
        matrix = [[dot(roots[a0], cov[b0]) for b0 in active] for a0 in active]
        assert inertia(matrix) == [9, 0, 0]
        vertices[name] = {'active': active, 'square': square.encoded()}
    assert dot(ps, mv(G, ps)) == -2*D-2*D*D*u*u
    lam = D/M
    X = [Q(2), -phi] + [lam*v for v in omega1]
    xp = evaluate(X)
    assert all(value.sign() <= 0 for name, value in xp.items() if name != 's')
    assert xp['s'] == 4*lam*(phi-1) and xp['s'].sign() > 0
    xsquare = dot(X, mv(G, X))
    assert xsquare == 2*(33+32*phi)/power(11+14*phi, 2)
    assert xsquare.sign() > 0
    assert xsquare.a*xsquare.a-5*xsquare.b*xsquare.b == F(4*19*59, 79**2)
    # Nine fixed rays from p: diagnostic first intersections, not a cone census.
    directions = [[Q(), Q(1)]+zero()] + [[Q(), Q()]+list(row) for row in ci]
    edge_records, edge_kinds = [], Counter()
    pvalues = evaluate(p)
    for i, direction in enumerate(directions):
        slopes = evaluate(direction)
        assert slopes[f'old{i}'].sign() < 0
        assert all(slopes[f'old{j}'] == 0 for j in range(9) if j != i)
        hits = [(-pvalues[name]/slope, name) for name, slope in slopes.items()
                if pvalues[name].sign() < 0 and slope.sign() > 0]
        assert hits
        first = hits[0][0]
        for value, _ in hits[1:]:
            if (value-first).sign() < 0:
                first = value
        endpoint = add(p, direction, first)
        ep = evaluate(endpoint)
        assert first.sign() > 0 and all(v.sign() <= 0 for v in ep.values())
        square = dot(endpoint, mv(G, endpoint))
        kind = {-1: 'finite', 0: 'ideal', 1: 'spacelike'}[square.sign()]
        edge_kinds[kind] += 1
        active = [name for name, v in ep.items() if v == 0]
        if kind == 'finite':
            assert len(active) == 9
            assert inertia([[dot(roots[a0], cov[b0]) for b0 in active]
                            for a0 in active]) == [9, 0, 0]
        edge_records.append({'freed_wall': f'old{i}', 'parameter': first.encoded(),
                             'first_hit_walls': [name for value, name in hits if value == first],
                             'square': square.encoded(), 'kind': kind})
    return {
        'status': 'PASS_EXACT_26_WALL_LOCAL_ALGEBRA', 'walls': 26, 'wall_pairs': 325,
        'sources': PINS, 'checker_sha256': digest(Path(__file__)),
        'norms': {name: value.encoded() for name, value in norm.items()},
        'new_covectors': {name: [v.encoded() for v in cov[name]] for name in expected},
        'reflection_integrality_checked': True, 'wall_pair_categories': dict(categories),
        'wall_pairs_sha256': canonical_digest(pairs),
        'strict_center_products': {name: value.encoded() for name, value in cp.items()},
        'actual_facet_projection_witnesses': 26, 'finite_vertex_witnesses': vertices,
        'old25_escape_square': xsquare.encoded(), 'escape_product_with_s': xp['s'].encoded(),
        'initial_p_edge_first_hits': edge_records, 'initial_p_edge_kinds': dict(edge_kinds),
        'inherited_only': 'Two full41-vertex Bugaenko H8 facets and their completeness; '
                          'the83/84 vertex count uses inherited disjoint facet sets, not a new census.',
        'scope': 'Exact26-wall algebra, actual-facet projection and listed vertex witnesses. '
                 'Nine p-edge diagnostics only; no complete root, vertex, cusp or final-cone census. '
                 'No finite-volume, compactness, nonreflectivity or dimension conclusion.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = check()
    if args.verify:
        with CERT.open() as stream:
            assert json.load(stream) == result, 'Saved certificate differs'
    else:
        with CERT.open('x') as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write('\n')
    print(json.dumps({key: result[key] for key in
                      ('status', 'walls', 'wall_pairs', 'initial_p_edge_kinds')}))


if __name__ == '__main__':
    main()
