#!/usr/bin/env python3
"""Exact fixed seed and universal profile algebra; ROOT SCHEDULES ALL RUNS.

No height loop, root census, cone census, optimizer, or compactness claim.
The wrapper owns resource gates and the shared lock. Import runs no check.
"""
import argparse
import hashlib
import importlib.util
import json
import signal
from fractions import Fraction as F
from math import isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / 'results/compact_h9_golden_discriminant59_edges_20261020.json'
PINS = {
    'results/compact_h9_golden_discriminant59_edge_route_20261020.md': '725d8775cebec51d8c18b4b0a0d3be0255e8b1472621adb4d3cefa656adb429e',
    'data/known_examples/bugaenko_compact_h8.json': 'c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c',
    'data/known_examples/bugaenko_vinberg_transcript.json': 'c4178fc4678bacb48c430b2cec18e16577078824fdd6dea39375e63ed6ae015a',
    'scripts/check_arithmetic_dn_opposed_fork_models_20261019.py': '185d20648368696429b675c1b6c023d30bb1be262ff83d3622e4f48801e1bae9',
}


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def check():
    paths = list(PINS) + [str(Path(__file__).resolve().relative_to(ROOT))]
    sources = {p: digest(ROOT/p) for p in paths}
    assert all(sources[p] == value for p, value in PINS.items())
    helper_path = ROOT/'scripts/check_arithmetic_dn_opposed_fork_models_20261019.py'
    spec = importlib.util.spec_from_file_location('pinned_dn_pure_algebra', helper_path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)  # definitions only; never helper.main/model/run
    Q = helper.field(5)
    g = lambda a=0, b=0: Q(F(a)+F(b)/2, F(b)/2)
    phi, zero = g(0, 1), g()
    sq = lambda x: x*x
    enc = lambda x: [str(x.a-x.b), str(2*x.b)]  # a+b*phi
    vec = lambda values: [g(x) for x in values]
    add = lambda a, b: [x+y for x, y in zip(a, b)]
    scale = lambda t, a: [t*x for x in a]
    mv = lambda m, a: [sum(x*y for x, y in zip(row, a)) for row in m]
    transpose = lambda m: list(map(list, zip(*m)))
    norm = lambda x: (x*x.conjugate()).a
    data = json.loads((ROOT/'data/known_examples/bugaenko_compact_h8.json').read_text())
    decode = lambda value: g(value['a'], value['b'])
    C = [[decode(x) for x in row[1:]] for row in data['lattice_gram'][1:]]
    Ci = helper.inverse(C, Q)
    omega = [[Ci[i][j] for i in range(8)] for j in range(8)]
    b, delta = g(7, 5), g(5, 12)
    assert norm(b) == 59 and delta == phi*b
    assert all(59 % p for p in (2, 3, 5, 7))
    assert (34*34-34-1) % 59 == 0 and (7+5*34) % 59 == 0
    assert (b.conjugate()).sign() > 0 and (4-b.conjugate()).sign() > 0
    assert (59*(4-b.conjugate())-4).sign() > 0
    assert (-delta.conjugate()-2).sign() > 0 and (25-delta).sign() > 0
    G = [[zero for _ in range(10)] for _ in range(10)]
    G[0][0], G[0][1], G[1][0], G[1][1] = -2*phi, g(-1), g(-1), g(4, 2)
    for i in range(8):
        for j in range(8):
            G[i+2][j+2] = C[i][j]
    dot = lambda a, b: helper.bilinear(a, G, b)
    assert helper.inertia(G) == [9, 1, 0]
    assert helper.inertia(G, True) == [10, 0, 0]
    embed = lambda root: [decode(root[0]), zero] + list(map(decode, root[1:]))
    old = list(map(embed, data['simple_roots']))
    r, c = [phi-1, g(-2)]+vec([0]*8), vec([1, 1]+[0]*8)
    h = [g(3, 2), g(1, 2)]+vec([0]*8)
    e0 = vec([1]+[0]*9)
    u = [sq(phi), phi]+scale(sq(phi), omega[7])
    assert omega[7] == vec([2, 3, 4, 6, 5, 4, 3, 2])
    assert omega[1] == vec([5, 8, 10, 15, 12, 9, 6, 3])
    assert dot(r, r) == 2*b and dot(c, c) == 2 and dot(u, u) == 2
    assert dot(r, c) == -b and dot(h, h) == -2*delta and dot(h, c) == 0
    assert r[0]*c[1]-r[1]*c[0] == sq(phi)
    assert mv(G, r) == [zero, -b]+vec([0]*8)
    sx = [g(9, 16), zero]+scale(g(6, 10), omega[1])
    sy = scale(F(1, 2), add(scale(phi-1, sx), scale(-1, r)))
    permutation = {0: 7, 2: 6, 3: 5, 4: 4, 5: 3, 6: 2, 7: 0}
    columns = [sx, sy]
    for j in range(8):
        columns.append(scale(-1, old[10]) if j == 1 else vec([0, 0]+[int(i == permutation[j]) for i in range(8)]))
    S = transpose(columns)
    identity = [[g(int(i == j)) for j in range(10)] for i in range(10)]
    assert all(x.integral() for row in S for x in row)
    assert [mv(S, column) for column in columns] == transpose(identity)
    assert [[dot(a, b) for b in columns] for a in columns] == G
    reversal = [8, 7, 6, 5, 4, 3, 2, 1, 0, 10, 9]
    assert [mv(S, x) for x in old] == [old[i] for i in reversal]
    assert mv(S, r) == r
    sc, su = mv(S, c), mv(S, u)
    roots = old + [r, c, u, sc, su]
    assert len(roots) == 16 and dot(c, su) == g(-62, -98)
    assert dot(u, su) == g(-72, -116)
    native = [(g(1), 3), (g(2), 4), (sq(phi), 5), (g(3), 6), (phi+2, 10)]
    # Fundamental-wall labels stay native; arbitrary intersecting mirrors
    # also allow angles 2*pi/5 and 3*pi/10, with conjugate cosine squares.
    arrangement_finite = [value for value, _ in native]+[2-phi, 3-phi]
    for value in arrangement_finite:
        assert norm(value).denominator == 1 and int(norm(value)) % 59
    for value in [zero]+arrangement_finite:
        assert int(norm(value-4)) % 59
    gram, pairs = [[dot(x, y) for y in roots] for x in roots], []
    for i, root in enumerate(roots):
        q, cov = gram[i][i], mv(G, root)
        assert q == 2 or q == 2*b
        assert all((2*x*y/q).integral() for x in root for y in cov)
        for j in range(i):
            product = gram[i][j]
            assert product.sign() <= 0
            if product == 0:
                label = 2
            else:
                cosine4 = 4*sq(product)/(q*gram[j][j])
                matches = [m for x, m in native if x == cosine4]
                assert matches or (cosine4-4).sign() > 0
                label = matches[0] if matches else 'ultraparallel'
            pairs.append([j, i, label, enc(product)])
    direction = add(scale(-1, c), [zero, zero]+[sum(row) for row in Ci])
    center = add(h, scale(F(1, 1000000), direction))
    assert dot(center, center).sign() < 0 and center[0].sign() > 0
    assert all(dot(center, x).sign() < 0 for x in roots)
    center = add(center, mv(S, center))
    assert all(dot(center, x).sign() < 0 for x in roots)
    assert dot(center, center).sign() < 0 and mv(S, center) == center

    def polynomial2(fun):
        a, px, nx, pz, nz, both = [fun(g(x), g(z)) for x, z in ((0, 0), (1, 0), (-1, 0), (0, 1), (0, -1), (1, 1))]
        return [a, (px-nx)/2, (pz-nz)/2, (px+nx)/2-a, (pz+nz)/2-a, both-px-pz+a]

    profiles = []
    rows = [('E8_extension', 1, 0, True, 2), ('isolated_short', 4, 0, False, 2),
            ('A2', 3, 1, False, 2), ('I2_5', 3-phi, phi, False, 2), ('isolated_long', 4*b, 0, False, 2*b)]
    tau = delta/(2*sq(phi))
    assert 4*sq(sq(phi))-delta == 3 and (delta-sq(tau)).sign() > 0
    assert helper.bilinear(omega[6], C, omega[6]) == 6
    assert helper.bilinear(omega[6], C, omega[7]) == 3
    for name, R, k, extension, q in rows:
        R, k, q = [x if isinstance(x, Q) else g(x) for x in (R, k, q)]
        def root_at(X, Z):
            spatial = scale((Z-3)/2 if extension else Z/2, omega[7])
            if extension:
                spatial = add(spatial, omega[6])
            return [(h[0]*X-k)/2, (h[1]*X-k)/2]+spatial
        assert polynomial2(lambda X, Z: dot(root_at(X, Z), root_at(X, Z))-q) == [-R/2, zero, zero, -delta/2, g(F(1, 2)), zero]
        endpoint = add(h, [zero, zero]+scale(tau, omega[7]))
        assert polynomial2(lambda X, Z: dot(endpoint, root_at(X, Z))) == [zero, -delta, tau, zero, zero, zero]
        norm_bound = g(norm(R))/(3*(-delta.conjugate()))
        if name in ('E8_extension', 'I2_5'):
            assert (1-norm_bound).sign() > 0
            outcome = 'no_nonzero_integral_X_by_norm'
        elif name in ('A2', 'isolated_short'):
            assert (4-norm_bound).sign() > 0
            assert (sq(phi)-R/3).sign() > 0 and (sq(phi)-R.conjugate()/(-delta.conjugate())).sign() > 0
            if name == 'A2':
                assert root_at(g(1), 2*sq(phi)) == u
                assert sq(2*sq(phi)) == delta+3
                outcome = 'unique_X_1_Z_2phi_squared'
            else:
                assert norm(delta+4) == 45 and isqrt(45)**2 != 45
                outcome = 'X_1_but_Z_not_field_square'
        else:
            assert ((1-(3+delta)/(3*59))).sign() > 0
            assert mv(G, root_at(zero, g(1)))[9] == g(F(1, 2))
            outcome = 'Z_over_2b_integral_nonzero_but_norm_below_1'
        profiles.append({'name': name, 'R': enc(R), 'outcome': outcome})

    def vertex(indices):
        small = [[gram[i][j] for j in indices] for i in indices]
        assert helper.inertia(small) == [9, 0, 0]
        inv = helper.inverse(small, Q)
        coefficients = mv(inv, [dot(center, roots[i]) for i in indices])
        point = list(center)
        for i, coefficient in zip(indices, coefficients):
            point = add(point, scale(-coefficient, roots[i]))
        assert dot(point, point).sign() < 0 and point[0].sign() > 0
        assert [i for i, root in enumerate(roots) if dot(point, root) == 0] == sorted(indices)
        assert all(dot(point, root).sign() <= 0 for root in roots)
        return [enc(x/point[0]) for x in point]

    white = [1, 9, 2, 3, 4, 5, 6]
    vertices = [vertex(sorted(indices+[11])) for indices in data['geometry']['vertices']]
    vertices += [vertex(sorted(white+[7, 12])), vertex(sorted(white+[12, 13])), vertex(sorted(white+[8, 13]))]
    assert len(vertices) == len({json.dumps(x) for x in vertices}) == 44
    counter = [g(-2, -1), -sq(phi)]+scale(sq(phi), vec([1]+[0]*7))
    assert dot(counter, counter) == 2 and dot(counter, c) == -1
    assert sources == {p: digest(ROOT/p) for p in paths}, 'input changed during check'
    return {'status': 'PASS_FIXED_L59_SEED_AND_EXACT_PROFILE_ALGEBRA', 'sources': sources,
            'scope': 'written norm-pool/contact/classification/unit and chamber inputs; no scientific census or compactness claim',
            'owns_lock': False, 'field_encoding': 'a+b*phi, coefficients rational strings',
            'b': enc(b), 'delta': enc(delta), 'norm_pool': [enc(g(2)), enc(2*b)],
            'roots': [[enc(x) for x in root] for root in roots], 'pairs': pairs,
            'S_columns': [[enc(x) for x in column] for column in columns], 'strict_center': list(map(enc, center)),
            'five_profiles': profiles, 'first_parameter': enc(tau), 'finite_vertices_verified_not_exhaustive': vertices,
            'finite_vertex_lower_bound': 44, 'full_compactness': 'OPEN; fixed seed has uncapped c edges'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    for option in ('write', 'verify', 'stdout'):
        mode.add_argument('--'+option, action='store_true')
    args = parser.parse_args()
    if args.write and CERT.exists():
        raise SystemExit('refusing to overwrite certificate')
    signal.alarm(45)
    result = json.loads(json.dumps(check(), sort_keys=True))
    signal.alarm(0)
    if args.verify:
        assert result == json.loads(CERT.read_text()), 'saved certificate differs'
    elif args.write:
        with CERT.open('x') as stream:
            json.dump(result, stream, sort_keys=True, indent=2)
            stream.write('\n')
    else:
        print(json.dumps(result, sort_keys=True))
        return
    print(json.dumps({'status': result['status'], 'roots': 16, 'pairs': len(result['pairs']),
                      'finite_vertex_lower_bound': 44, 'mode': 'replay' if args.verify else 'generation'}))


if __name__ == '__main__':
    main()
