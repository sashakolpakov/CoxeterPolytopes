#!/usr/bin/env python3
"""Fixed literal-H19 full-edge/period algebra; root schedules execution.

Outer --parent-lock owns resources. No body runs on import. No height,
Pell or cone census; full calibration/germ/index arguments remain written.
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
CERT = ROOT/'results/h19_d7_edge_period_20261020.json'
HELPER = 'scripts/check_arithmetic_dn_opposed_fork_models_20261019.py'
FIXTURE = 'data/known_examples/bugaenko_compact_h8.json'
TRANSCRIPT = 'data/known_examples/bugaenko_vinberg_transcript.json'
PINS = {
    HELPER: '185d20648368696429b675c1b6c023d30bb1be262ff83d3622e4f48801e1bae9',
    FIXTURE: 'c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c',
    TRANSCRIPT: 'c4178fc4678bacb48c430b2cec18e16577078824fdd6dea39375e63ed6ae015a',
    'results/bugaenko_vinberg_transcript.md': '113e54694c6e2c05a2de4b8bdd7a5a072aabd2572aba760ecb8d4e039ebda719',
    'results/compact_h9_discriminant19_d7_edge_period_20261020.md': 'ecc61a2dcff6df1678cc94533aa201c70cd0b262e26deab5ef85f5d64ebaf824',
    'results/review_h19_d7_edge_period_resource_20261020.md': '55428c652f61f2b59134f69f54f34ac336bf8084d1846328bce54e2557275602',
    'results/review_h19_d7_edge_period_root_20261020.md': 'd873ea0465023dafc7029fa0011230f6cbe0f38bef67fa33b3d0c22c78a9ff44',
    'results/review_h19_d7_edge_period_batch_guard_20261020.md': '16fe2efb19c6d607ed766eeeb825551ade90e7c3288d318b631df86b6f92a336',
    'results/compact_h9_discriminant19_mixed_norm_seed_20261019.md': '70fa010732d837d5c044701b4c18ed5e68cb36c18cc0954dc7d7f369250a5c87',
    'results/compact_h9_discriminant19_polar_caps_research_20261019.md': '792634e935866ed78d95d0a93c04fdd7a8484f3605b5c50a66c1d376c9599f1b',
    'results/compact_h9_golden_discriminant59_edge_route_20261020.md': '725d8775cebec51d8c18b4b0a0d3be0255e8b1472621adb4d3cefa656adb429e',
    'results/og10_golden_h9_infinite_stabilizer.md': '1d9a71240d499e498395f5c1c7eade7c2f307f2fed5e44c2578e02b1357392a1',
    'scripts/run_bounded_research_check_20261018.py': '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
    'scripts/build_illustrated_ra_manuscript_20261018.py': '27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15',
    'scripts/run_lean_batches_signal_safe_20261018.py': 'c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d',
    'scripts/run_lean_batches_20261018.py': '97822ecf4d9aef2bf04aa429b582c4cc36cf6d1b6f0364579251278f66fae72c',
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
    assert all(sources[p] == h for p, h in PINS.items())
    spec = importlib.util.spec_from_file_location('pinned_h19_golden_definitions', ROOT/HELPER)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)  # definitions only; never model/main
    Q = helper.field(5)
    g = lambda a=0, b=0: Q(F(a)+F(b)/2, F(b)/2)
    enc = lambda x: [str(x.a-x.b), str(2*x.b)]
    decode = lambda x: g(x['a'], x['b'])
    phi, zero, one = g(0, 1), g(), g(1)
    sq = lambda x: x*x
    def power(x, n):
        value = one
        for _ in range(n):
            value = value*x
        return value
    norm = lambda x: (x*x.conjugate()).a
    add = lambda a, b: [x+y for x, y in zip(a, b)]
    scale = lambda t, a: [t*x for x in a]
    mv = lambda m, a: [sum(x*y for x, y in zip(row, a)) for row in m]
    tr = lambda m: list(map(list, zip(*m)))
    mm = lambda a, b: [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]
    trace = lambda m: sum(m[i][i] for i in range(len(m)))
    assert sq(phi) == phi+1 and phi.integral() and not (phi/2).integral()
    assert phi.conjugate() == 1-phi and norm(phi) == -1
    fixture = json.loads((ROOT/FIXTURE).read_text())
    transcript = json.loads((ROOT/TRANSCRIPT).read_text())
    assert fixture['dimension'] == 8 and fixture['certificate_kind'] == 'compact_reflective'
    assert len(fixture['geometry']['vertices']) == 41 and fixture['geometry']['cusps'] == []
    assert transcript['root_norm_normalization']['remaining_representative'] == '2'
    assert transcript['height_zero']['all_root_count'] == 240
    assert [x['solution_count'] for x in transcript['levels']] == [1, 1, 0, 0, 0, 0, 0, 1]
    assert transcript['output']['wall_count'] == 11
    C = [[g(2*int(i == j)) for j in range(8)] for i in range(8)]
    for i, j in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)]:
        C[i][j] = C[j][i] = g(-1)
    pi = [0, 7, 1, 2, 3, 4, 5, 6]
    Cs = [[decode(x) for x in row[1:]] for row in fixture['lattice_gram'][1:]]
    assert [[C[pi[i]][pi[j]] for j in range(8)] for i in range(8)] == Cs
    Ci = helper.inverse(C, Q)
    weights = tr(Ci)
    assert weights[0] == list(map(g, [4, 7, 10, 8, 6, 4, 2, 5]))
    assert weights[7] == list(map(g, [5, 10, 15, 12, 9, 6, 3, 8]))
    assert all(x.integral() for row in Ci for x in row)
    gram = [[zero for _ in range(10)] for _ in range(10)]
    gram[0][0], gram[0][1], gram[1][0], gram[1][1] = -2*phi, phi, phi, g(2)
    for i in range(8):
        for j in range(8):
            gram[i+2][j+2] = C[i][j]
    dot = lambda x, y: helper.bilinear(x, gram, y)
    b, D = 4+phi, 1+5*phi
    assert D == phi*b and norm(b) == 19 and (-2*phi)*2-sq(phi) == -D
    assert helper.inertia(gram) == [9, 1, 0] and helper.inertia(gram, True) == [10, 0, 0]
    basis = [[g(int(i == j)) for i in range(10)] for j in range(10)]
    spatial = lambda x: [zero, zero]+x
    weight0 = spatial(weights[0])
    assert dot(weight0, weight0) == 4
    p, pr = [g(2), -phi]+[zero]*8, basis[0]
    r, c = [one, g(2)]+[zero]*8, scale(-1, basis[1])
    assert dot(r, r) == 2*b and dot(c, c) == 2 and dot(r, c) == -b
    assert r[0]*c[1]-r[1]*c[0] == -1
    assert mv(gram, r) == [zero, b]+[zero]*8
    assert dot(p, p) == -2*D and mv(gram, p) == [-D]+[zero]*9
    assert dot(p, c) == dot(pr, r) == zero and dot(pr, c).sign() < 0 and dot(p, r).sign() < 0

    # Whole-arrangement finite-angle and binary nonintervention scalars.
    values = [one, g(2), sq(phi), g(3), 2+phi, 2-phi, 3-phi]
    assert list(map(norm, values)) == list(map(F, [1, 4, 1, 9, 5, 1, 5]))
    assert [norm(x-4) for x in values] == list(map(F, [9, 4, 5, 1, 1, 5, 1]))
    assert all(int(norm(x)) % 19 and int(norm(x-4)) % 19 for x in values)
    assert norm(g(-4)) == 16 and (15*15-15-1) % 19 == (4+15) % 19 == 0
    assert (b.conjugate()).sign() > 0 and (4-b.conjugate()).sign() > 0
    assert (4-b.conjugate()-g(F(3, 5))).sign() > 0
    integer_norm_upper = g(4)/(19*(4-b.conjugate()))
    assert integer_norm_upper.sign() > 0 and (1-integer_norm_upper).sign() > 0
    # Native M/permutation calibration: equality of all output roots as sets.
    def embed(raw):
        z = [zero]*8
        for i, value in enumerate(raw[1:]):
            z[pi[i]] = decode(value)
        return [decode(raw[0]), zero]+z
    wallsF = list(map(embed, fixture['simple_roots']))
    vorder = list(map(embed, transcript['output']['roots_in_vinberg_order']))
    assert len(wallsF) == len(vorder) == 11
    assert all(vorder.count(x) == wallsF.count(x) == 1 for x in wallsF)
    assert all(dot(x, x) == 2 and dot(x, r) == 0 for x in wallsF)
    assert [[dot(x, y) for y in wallsF] for x in wallsF] == [[decode(x) for x in row] for row in fixture['wall_gram']]
    gcap = [phi, zero]+scale(phi, weights[0])
    old9 = [one, zero]+scale(phi, list(map(g, [2, 4, 6, 5, 4, 3, 2, 3])))
    kcap = [g(2, 3), zero]+list(map(lambda t: g(*t), [(6, 10), (12, 20), (18, 30), (15, 24),
                                                        (12, 18), (8, 12), (4, 6), (9, 16)]))
    assert wallsF[0] == gcap and wallsF[8] == old9 and wallsF[10] == kcap
    assert all(dot(pr, x).sign() <= 0 and dot(p, x).sign() <= 0 for x in wallsF)

    # Eight literal finite profiles, checked as full quadratic identities.
    profiles = [('isolated_short', g(8), zero, None, g(2)), ('A2', g(6), one, None, g(2)),
                ('I2_5', 6-2*phi, phi, None, g(2)), ('isolated_long', 8*b, zero, None, 2*b),
                ('D8', g(4), zero, 6, g(2)), ('E8_at_1', one, zero, 1, g(2)),
                ('E8_at_7', one, zero, 7, g(2)), ('D9', g(2), one, 6, g(2))]
    profile_rows = []
    for name, R, kval, j, q in profiles:
        const = scale(-kval/2, c)
        if j is not None:
            m = helper.bilinear(weights[0], C, weights[j])
            assert m == {6: g(2), 1: g(7), 7: g(5)}[j]
            residual = helper.bilinear(weights[j], C, weights[j])-sq(m)/4
            assert residual == (one if j == 6 else g(F(7, 4)))
            const = add(const, spatial(add(weights[j], scale(-m/4, weights[0]))))
        ux, uz = scale(F(1, 2), p), scale(F(1, 4), weight0)
        co = [dot(const, const)-q, 2*dot(const, ux), 2*dot(const, uz),
              dot(ux, ux), dot(uz, uz), 2*dot(ux, uz)]
        assert co == [-R/4, zero, zero, -D/2, g(F(1, 4)), zero]
        profile_rows.append({'name': name, 'R': enc(R), 'k': enc(kval), 'j': j,
                             'norm_coefficients': list(map(enc, co))})
    residues = {i*i % 19 for i in range(19)}
    assert [8 % 19, (6-2*15) % 19, 2 % 19] == [8, 14, 2]
    assert all(x not in residues for x in (8, 14, 2))
    z = [2*power(phi, 4), -power(phi, 5)]+scale(b*sq(phi), weights[0])
    covz = mv(gram, z)
    assert dot(z, z) == 2*b and norm(z[1]) == -1
    assert covz == scale(b, [-power(phi, 5), zero]+[sq(phi)]+[zero]*7)
    assert all(x.integral() for x in z) and all((x*y/b).integral() for x in z for y in covz)
    Xz, Zz = 2*power(phi, 4), 4*b*sq(phi)
    assert sq(Zz)-2*D*sq(Xz) == 8*b
    tau = power(phi, 3)/2
    assert D*Xz/Zz == tau
    V = add(p, scale(tau, weight0))
    assert dot(V, V) == 3-2*phi and dot(V, V).sign() < 0
    K, conjugate_a = power(phi, 8)/(2*b), -2*D.conjugate()
    assert K == g(44, 71)/38 and (g(F(17, 4))-K).sign() > 0
    assert conjugate_a == 10*phi-12 and (conjugate_a-4).sign() > 0
    assert (phi-g(F(8, 5))).sign() > 0 and (g(F(13, 8))-phi).sign() > 0
    assert (2*phi-1-g(F(11, 5))).sign() > 0
    assert F(17, 16) < F(441, 400) and F(17, 4) < F(441, 100)
    assert F(51, 2) < F(2601, 100) and F(3, 2) < F(25, 16)
    # Complete coarse strips give B=0,1,2, A=-3,...,5 in all three rows.
    assert F(51, 10)+F(5, 4) < 3*F(11, 5) and 2*F(13, 8) < 4
    short_domains = {}
    for name, R, divisor, parity, expected in [('E8_both', one, 2, False, []),
                                              ('D8', g(4), 2, False, []),
                                              ('A2', g(6), 1, True, [(1, 1)])]:
        found = []
        for Bc in range(3):
            for Ac in range(-3, 6):
                value = g(Ac, Bc)
                X = divisor*value
                if value.sign() <= 0 or (R*K-sq(X)).sign() < 0 or (
                        R/conjugate_a-sq(X.conjugate())).sign() < 0:
                    continue
                if parity and (Ac % 2 != 1 or Bc % 2 != 1):
                    continue
                found.append((Ac, Bc))
        assert found == expected
        short_domains[name] = found
    for Ac in range(2):
        for Bc in range(2):
            assert ((1-phi*g(Ac, Bc))/2).integral() == (Ac == Bc == 1)
            assert ((-phi*g(Ac, Bc))/2).integral() == (Ac == Bc == 0)
    radicand = 2*D*power(phi, 4)+6
    assert radicand == g(40, 56) and norm(radicand) == 704 and isqrt(704)**2 != 704
    assert (b.conjugate()-g(F(27, 8))).sign() > 0
    assert F(4, 27) < F(4, 25)
    long_domain = []
    for Bc in range(2):
        for Ac in range(-1, 3):
            t = g(Ac, Bc)
            if t.sign() > 0 and (sq(phi)-t).sign() >= 0 and (
                    1-2*b.conjugate()*sq(t.conjugate())).sign() >= 0:
                long_domain.append((Ac, Bc))
    assert long_domain == [(1, 1)]
    assert 2*b*power(phi, 4)-phi*power(phi, 8) == 1
    assert z == add(scale(power(phi, 4), p), scale(b*sq(phi), weight0))

    # Actual calibrated target vertex and positive corner Grams.
    whites = [scale(-1, basis[i+2]) for i in range(1, 8)]
    B = [one, zero]+scale(phi/2, weights[0])
    activeF = [i for i, x in enumerate(wallsF) if dot(B, x) == 0]
    assert activeF == [0, 2, 3, 4, 5, 6, 7, 9] and activeF in fixture['geometry']['vertices']
    assert dot(B, B) == 1-phi and dot(B, r) == zero
    known = wallsF+[r, c, z]
    assert all(dot(B, x).sign() <= 0 and dot(V, x).sign() <= 0 for x in known)
    assert [i for i, x in enumerate(known) if dot(V, x) == 0] == [2, 3, 4, 5, 6, 7, 9, 12, 13]
    for corner in (whites+[c, z], whites+[r, gcap]):
        assert helper.inertia([[dot(x, y) for y in corner] for x in corner]) == [9, 0, 0]
    assert all(dot(z, x) == 0 for x in whites+[c])

    # Native full-lattice matrices; no restriction to a finite cap subgroup.
    P = [sq(phi), -phi]+scale(sq(phi), weights[0])
    assert dot(P, P) == -2*phi and dot(P, r) == -D and dot(P, c) == -1
    assert mv(gram, P) == [-power(phi, 5), one]+[sq(phi)]+[zero]*7
    wcols = basis[:]
    wcols = [x[:] for x in wcols]
    wcols[2] = [zero, zero]+list(map(g, [1, 3, 5, 4, 3, 2, 1, 3]))
    wcols[3], wcols[9] = scale(-1, basis[9]), scale(-1, basis[3])
    for i in range(2, 7):
        wcols[i+2] = scale(-1, basis[i+2])
    w0 = tr(wcols)
    Jcols = [add(scale(-1, e), scale(-dot(P, e)/phi, P)) for e in basis]
    J = tr(Jcols)
    H = mm(J, w0)
    Acols = [[g(9, 16), zero]+scale(g(6, 10), weights[7]),
             [g(-4, -8), one]+scale(g(-3, -5), weights[7])]
    reversal = {0: 6, 1: 5, 2: 4, 3: 3, 4: 2, 5: 1, 6: 0}
    for i in range(8):
        Acols.append(scale(-1, kcap) if i == 7 else basis[reversal[i]+2])
    A = tr(Acols)
    I = [[g(int(i == j)) for j in range(10)] for i in range(10)]
    for matrix in (w0, J, H, A):
        columns = tr(matrix)
        assert all(x.integral() for col in columns for x in col)
        assert mm(matrix, matrix) == I
        assert [[dot(x, y) for y in columns] for x in columns] == gram
    assert mv(w0, P) == P and mv(w0, V) == V and mm(w0, J) == H
    assert mv(J, r) == z and mv(J, c) == gcap
    assert [mv(H, x) for x in whites] == [whites[6]]+whites[1:6]+[whites[0]]
    assert mv(H, z) == r and mv(H, c) == gcap and dot(P, V) == -phi
    assert mv(H, V) == add(P, scale(-1, V)) == scale(phi-1, B)
    assert mv(A, r) == r and A[1] == basis[1]
    fixture_reversal = [8, 7, 6, 5, 4, 3, 2, 1, 0, 10, 9]
    assert [mv(A, x) for x in wallsF] == [wallsF[i] for i in fixture_reversal]
    AP = mv(A, P)
    assert AP[0] == g(8, 13) and AP[1] == -phi and AP[2] == g(27, 42)
    assert dot(P, AP) == g(-20, -34) and trace(mm(A, w0)) == g(4, 6)
    AH = mm(A, H)
    assert mm(AH, mm(H, A)) == I
    assert trace(AH) == -trace(mm(A, w0))-dot(P, AP)/phi == g(10, 14)
    assert (trace(AH)-10).sign() > 0
    assert P[0].sign() > 0 and mv(A, pr)[0].sign() > 0 and mv(AH, P)[0].sign() > 0
    direction = [zero, one]+[sum(row) for row in Ci]
    center = add(p, scale(F(1, 1000), direction))
    assert dot(center, center).sign() < 0 and all(dot(center, x).sign() < 0 for x in known)
    for x in known:
        q, cov = dot(x, x), mv(gram, x)
        assert q == 2 or q == 2*b
        assert all((2*a*y/q).integral() for a in x for y in cov)
    assert sources == {p: digest(ROOT/p) for p in paths}, 'source changed during check'
    return {'status': 'PASS_FIXED_H19_COMPLETE_EDGE_AND_INFINITE_PERIOD_ALGEBRA',
            'sources': sources, 'owns_lock': False, 'field_encoding': 'a+b*phi; rational strings',
            'scope': 'full primitive and H8 calibration, first-contact/germ and arithmetic-index arguments remain written inputs; no general dimension claim',
            'native_E8_permutation': pi, 'arrangement_angle_norms': list(map(lambda x: int(norm(x)), values)),
            'binary_integer_norm_upper': enc(integer_norm_upper), 'profiles': profile_rows,
            'all_height_mod19_exclusions': [8, 14, 2], 'short_complete_domains': short_domains,
            'A2_radicand': enc(radicand), 'A2_radicand_norm': 704, 'long_complete_domain': long_domain,
            'first_wall': list(map(enc, z)), 'actual_endpoint': list(map(enc, V)),
            'calibrated_target': list(map(enc, B)), 'target_fixture_vertex': activeF,
            'P': list(map(enc, P)), 'A_columns': [list(map(enc, x)) for x in Acols],
            'H_columns': [list(map(enc, x)) for x in tr(H)], 'AH_trace': enc(trace(AH)),
            'strict_common_center': list(map(enc, center)), 'infinite_order_period': True,
            'no_new_example': True, 'literal_host': 'E8+native[[-2phi,phi],[phi,2]]'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.write and CERT.exists():
        raise SystemExit('refusing to overwrite certificate')
    signal.alarm(45)
    result = json.loads(json.dumps(check(), sort_keys=True))
    payload = json.dumps(result, sort_keys=True, separators=(',', ':'))+'\n'
    assert len(payload.encode()) <= 61440
    assert result['sources'] == {p: digest(ROOT/p) for p in result['sources']}
    if args.verify:
        assert CERT.stat().st_size <= 61440
        assert result == json.loads(CERT.read_text()), 'saved certificate differs'
    else:
        with CERT.open('x') as stream:
            stream.write(payload)
    signal.alarm(0)
    print(json.dumps({'status': result['status'], 'mode': 'replay' if args.verify else 'generation',
                      'profiles': 8, 'unique_first_wall': True, 'period_trace': result['AH_trace']}))


if __name__ == '__main__':
    main()
