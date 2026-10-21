#!/usr/bin/env python3
"""Complete bounded D7+c edge census and fixed period algebra; root schedules runs.

The wrapper owns all resource gates and the shared lock. This module executes
no check on import. No prefix loop, optimizer, cone census or unit search.
"""
import argparse
import hashlib
import importlib.util
import json
import signal
from fractions import Fraction as F
from math import gcd, isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT/'results/compact_h9_golden59_d7_first_edge_20261020.json'
PINS = {
    'results/compact_h9_golden59_d7_edge_interval_20261020.md': 'e4715a78f22628bf7a27051ae6588adae7803f191a2a08f476d0a02fa2223f52',
    'results/compact_h9_golden59_long_facet_period_20261020.md': '699e9922fc0cc3ca6f16d5619e02ab345fbcb39ed25ddb355d32a4120ef14f9d',
    'results/compact_h9_golden_discriminant59_edges_20261020.json': '23c650dd70cb297af4923bbc5dbd30e19f672d6329c7c014ea90567a85831a99',
    'scripts/check_compact_h9_golden_discriminant59_edges_20261020.py': '01a569b66b77c38056255a94dce765e4784d809dfe86735ac5bef869755edf3d',
    'scripts/check_arithmetic_dn_opposed_fork_models_20261019.py': '185d20648368696429b675c1b6c023d30bb1be262ff83d3622e4f48801e1bae9',
    'data/known_examples/bugaenko_compact_h8.json': 'c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c',
}


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def check():
    paths = list(PINS)+[str(Path(__file__).resolve().relative_to(ROOT))]
    sources = {p: digest(ROOT/p) for p in paths}
    assert all(sources[p] == value for p, value in PINS.items())
    spec = importlib.util.spec_from_file_location('pinned_quadratic_definitions', ROOT/'scripts/check_arithmetic_dn_opposed_fork_models_20261019.py')
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)  # pure definitions only; no model/check/main
    Q = helper.field(5)
    g = lambda a=0, b=0: Q(F(a)+F(b)/2, F(b)/2)
    enc = lambda x: [str(x.a-x.b), str(2*x.b)]
    dec = lambda x: g(*x)
    zero, phi = g(), g(0, 1)
    sq = lambda x: x*x
    add = lambda a, b: [x+y for x, y in zip(a, b)]
    scale = lambda t, a: [t*x for x in a]
    mv = lambda m, a: [sum(x*y for x, y in zip(row, a)) for row in m]
    tr = lambda m: list(map(list, zip(*m)))
    mm = lambda a, b: [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]
    norm = lambda x: (x*x.conjugate()).a
    integer_norm = lambda x: int(norm(x)) if norm(x).denominator == 1 else None
    data = json.loads((ROOT/'data/known_examples/bugaenko_compact_h8.json').read_text())
    C = [[g(v['a'], v['b']) for v in row[1:]] for row in data['lattice_gram'][1:]]
    Ci = helper.inverse(C, Q)
    omega = tr(Ci)
    b, delta = g(7, 5), g(5, 12)
    assert delta == phi*b and norm(b) == 59
    G = [[zero for _ in range(10)] for _ in range(10)]
    G[0][0], G[0][1], G[1][0], G[1][1] = -2*phi, g(-1), g(-1), g(4, 2)
    for i in range(8):
        for j in range(8):
            G[i+2][j+2] = C[i][j]
    dot = lambda x, y: helper.bilinear(x, G, y)
    basis = [[g(int(i == j)) for i in range(10)] for j in range(10)]
    spatial = lambda z: [zero, zero]+z
    h = [g(3, 2), g(1, 2)]+[zero]*8
    c = add(basis[0], basis[1])
    w0 = spatial(omega[0])
    assert omega[0] == list(map(g, [4, 5, 7, 10, 8, 6, 4, 2]))
    assert dot(h, h) == -2*delta and dot(c, c) == 2 and dot(h, c) == 0
    assert dot(w0, w0) == 4
    saved = json.loads((ROOT/'results/compact_h9_golden_discriminant59_edges_20261020.json').read_text())
    roots = [[dec(x) for x in row] for row in saved['roots']]
    assert len(roots) == 16 and roots[12] == c
    center = list(map(dec, saved['strict_center']))
    assert dot(center, center).sign() < 0 and all(dot(center, r).sign() < 0 for r in roots)

    # Exact universal norm coefficients for all eight written spherical profiles.
    profiles = [('isolated_short', 0, None, g(8), g(2)),
                ('A2_on_c', 1, None, g(6), g(2)),
                ('I2_5_on_c', phi, None, g(6)-2*phi, g(2)),
                ('isolated_long', 0, None, 8*b, 2*b),
                ('D8', 0, 7, g(4), g(2)), ('E8_at_1', 0, 1, g(1), g(2)),
                ('E8_at_2', 0, 2, g(1), g(2)), ('D9_bridge', 1, 7, g(2), g(2))]
    profile_rows = []
    for name, k, j, R, q in profiles:
        const = scale(-k/2 if isinstance(k, Q) else F(-k, 2), c)
        if j is not None:
            m = helper.bilinear(omega[0], C, omega[j])
            const = add(const, spatial(add(omega[j], scale(-m/4, omega[0]))))
            assert m == {7: 2, 1: 5, 2: 7}[j]
        ux, uz = scale(F(1, 2), h), scale(F(1, 4), w0)
        coefficients = [dot(const, const)-q, 2*dot(const, ux), 2*dot(const, uz),
                        dot(ux, ux), dot(uz, uz), 2*dot(ux, uz)]
        assert coefficients == [-R/4, zero, zero, -delta/2, g(F(1, 4)), zero]
        profile_rows.append({'name': name, 'R': enc(R), 'norm_coefficients_1_X_Z_X2_Z2_XZ': list(map(enc, coefficients))})
    residues = {i*i % 59 for i in range(59)}
    exclusions = [(8, 0), (6, 0), (6, -2), (2, 0)]
    assert [(a+34*d) % 59 for a, d in exclusions] == [8, 6, 56, 2]
    assert all((a+34*d) % 59 not in residues for a, d in exclusions)
    assert (34*34-34-1) % 59 == (7+5*34) % 59 == 0

    alpha, gamma = g(38, 61), g(9, 14)
    beta = b*gamma
    v = add(scale(alpha, h), scale(beta, w0))
    Xv, Zv = 2*alpha, 4*beta
    assert beta == g(133, 213) and sq(Zv)-2*delta*sq(Xv) == 8*b
    assert dot(v, v) == 2*b and gcd(integer_norm(alpha), integer_norm(beta)) == 1
    assert (5000*b-sq(Xv)).sign() > 0 and (274-Xv).sign() > 0
    assert (-2*delta.conjugate()-4).sign() > 0 and (32-gamma).sign() > 0
    assert (2*b.conjugate()*g(F(81, 625))-1).sign() > 0

    # These are complete integer strips derived in the note, not a height prefix.
    Y_list, t_list = [], []
    for d in range(13):
        m = (isqrt(5*d*d)-d)//2
        for a in range(m-1, m+3):
            y = g(a, d)
            if y.sign() > 0 and (25-y).sign() > 0 and (1+2*delta.conjugate()*sq(y.conjugate())).sign() >= 0:
                Y_list.append((a, d))
    expected_Y = [(1, 1), (1, 2), (2, 3), (3, 5), (4, 6), (4, 7), (5, 8), (6, 9), (6, 10), (7, 11)]
    assert Y_list == expected_Y
    Y_norms = [integer_norm(1+2*delta*sq(g(*x))) for x in Y_list]
    assert Y_norms == [99, 641, 2059, 5771, 5401, 5571, 15489, 1531, 20249, 24131]
    assert all(n > 0 and isqrt(n)**2 != n for n in Y_norms)
    # Both E8 placements use the SAME larger list, with X even as an extra filter.
    E8_domain = [x for x in Y_list if not (x[0] % 2 or x[1] % 2)]
    assert E8_domain == [(4, 6), (6, 10)]
    for d in range(17):
        m = (isqrt(5*d*d)-d)//2
        for a in range(m-1, m+3):
            t = g(a, d)
            if t.sign() > 0 and (gamma-t).sign() >= 0 and (1-2*b.conjugate()*sq(t.conjugate())).sign() >= 0:
                t_list.append((a, d))
    assert t_list == [(1, 2), (2, 3), (3, 5), (4, 6), (4, 7), (5, 8), (6, 10), (7, 11), (8, 13), (9, 14)]
    long_radicands = [(g(-4, 14)*sq(g(*x))+1-phi) for x in t_list]
    t_norms = list(map(integer_norm, long_radicands))
    assert t_norms == [305, 1181, 3475, 1895, 1189, 9481, 11071, 12659, 25205, 1681]
    square_indices = [i for i, n in enumerate(t_norms) if n >= 0 and isqrt(n)**2 == n]
    assert square_indices == [9] and long_radicands[9] == sq(alpha) == g(5165, 8357)
    # In a field, these two roots are exactly +/-alpha. Positive contact chooses +.
    assert alpha.sign() > 0 and gamma.sign() > 0
    assert all(x.integral() for x in v)
    cov = mv(G, v)
    assert all((x/b).integral() for x in cov)
    assert all((2*x*y/(2*b)).integral() for x in v for y in cov)

    # Native compatibility is checked AFTER the complete census, never as a filter.
    new_pairs = []
    for i, old in enumerate(roots):
        product, q = dot(v, old), dot(old, old)
        assert q == 2 or q == 2*b
        assert product.sign() <= 0
        ratio4 = 4*sq(product)/(2*b*q)
        assert product == 0 or (ratio4-4).sign() > 0
        new_pairs.append([i, enc(product), 2 if product == 0 else 'ultraparallel'])
    A = g(10, 18)
    H = add(scale(A, h), scale(alpha, w0))
    tau = delta*Xv/Zv
    assert alpha/A == tau and (delta/2-sq(tau)).sign() > 0
    assert dot(H, H) == g(-4)/phi and dot(h, H).sign() < 0
    all_roots = roots+[v]
    white = [2, 3, 4, 5, 6, 7, 9, 12]
    assert [i for i, old in enumerate(all_roots) if dot(H, old) == 0] == white+[16]
    assert all(dot(H, old).sign() <= 0 for old in all_roots)
    assert helper.inertia([[dot(all_roots[i], all_roots[j]) for j in white+[16]] for i in white+[16]]) == [9, 0, 0]
    cv = dot(center, v)
    eps = g(1) if cv.sign() <= 0 else (-dot(h, v))/(2*(1+cv))
    if (eps-1).sign() > 0:
        eps = g(1)
    new_center = add(h, scale(eps, center))
    assert dot(new_center, new_center).sign() < 0 and dot(h, new_center).sign() < 0
    assert all(dot(new_center, old).sign() < 0 for old in all_roots)

    # Fixed period columns in the ORIGINAL full lattice basis, not a sublattice.
    columns = [[g(175, 285), g(119, 193)]+scale(g(99, 160), omega[0]),
               [g(-63, -103), g(-43, -70)]+scale(g(-36, -57), omega[0]),
               add([g(-44, -71), g(-30, -48)]+scale(g(-25, -40), omega[0]), basis[2])]+basis[3:]
    M = tr(columns)
    assert all(x.integral() for column in columns for x in column)
    assert [[dot(x, y) for y in columns] for x in columns] == G
    I = [[g(int(i == j)) for j in range(10)] for i in range(10)]
    Mi = mm(mm(helper.inverse(G, Q), tr(M)), G)
    assert mm(M, Mi) == I and all(x.integral() for row in Mi for x in row)
    P = add(scale(2, basis[0]), scale(phi, w0))
    assert dot(P, P) == g(-4)/phi and mv(M, roots[11]) == v
    assert mv(M, roots[0]) == c and mv(M, P) == H
    assert all(mv(M, basis[j+2]) == basis[j+2] for j in range(1, 8))
    old_vertex = [0, 2, 3, 4, 5, 6, 7, 9]
    assert old_vertex in data['geometry']['vertices']
    assert [i for i, old in enumerate(roots) if dot(P, old) == 0] == old_vertex+[11]
    assert all(dot(P, old).sign() <= 0 for old in roots)
    trace = sum(M[i][i] for i in range(10))
    assert trace == g(40, 55) and (trace-10).sign() > 0
    assert sources == {p: digest(ROOT/p) for p in paths}, 'input changed during check'
    return {'status': 'PASS_COMPLETE_D7_EDGE_AND_FIXED_PERIOD_ALGEBRA', 'sources': sources,
            'owns_lock': False, 'field_encoding': 'a+b*phi, rational coefficient strings',
            'scope': 'complete finite arithmetic domain; written norm-pool/spherical/actual-edge/chamber-germ premises remain explicit',
            'profile_norms': profile_rows, 'excluded_mod59': [8, 6, 56, 2],
            'integer_pairs_tested': {'Y': 52, 't': 68}, 'Y_domain': Y_list, 'Y_norms': Y_norms,
            'E8_even_subdomain_both_placements': E8_domain, 't_domain': t_list, 't_norms': t_norms,
            'unique_positive_root': list(map(enc, v)), 'new_pairs': new_pairs,
            'endpoint': list(map(enc, H)), 'first_parameter': enc(tau), 'finite_corner': 'D7+A1_short+A1_long',
            'strict_center': list(map(enc, new_center)), 'period_columns': [list(map(enc, x)) for x in columns],
            'period_inverse': [list(map(enc, row)) for row in Mi], 'period_trace': enc(trace),
            'old_actual_vertex': old_vertex+[11], 'new_actual_vertex': white+[16],
            'no_new_example': True, 'global_chamber_and_coset_argument': 'separate pinned written period proof'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    for name in ('write', 'verify', 'stdout'):
        mode.add_argument('--'+name, action='store_true')
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
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write('\n')
    else:
        print(json.dumps(result, sort_keys=True))
        return
    print(json.dumps({'status': result['status'], 'pairs_tested': 120,
                      'unique_edge_wall': True, 'period_trace': result['period_trace'],
                      'mode': 'replay' if args.verify else 'generation'}))


if __name__ == '__main__':
    main()
