#!/usr/bin/env python3
"""Fixed R5 interval and integral chamber period; no height/root search.

Run with -B through the parent-lock bounded wrapper. Geometry, completeness
of the six finite profiles, and actual-germ provenance are written inputs.
"""
import argparse
import hashlib
import json
import math
import os
import re
import signal
from fractions import Fraction as F
from pathlib import Path
from check_arithmetic_dn_opposed_fork_models_20261019 import field, inverse

ROOT = Path(__file__).resolve().parents[1]
STEM = 'golden_dyadic_r5_chamber_period_20261020'
CERT = ROOT/('results/'+STEM+'.json')
PAPER = 'results/golden_dyadic_r5_chamber_period_candidate_20261020.md'
FREEZE = 'results/golden_dyadic_d4_research_20261020_freeze.json'
FREEZE_SHA = '22dd341c2d79651742aa93ebdc8675e8fc55025fdd0cbcb678b128f7e59aa37b'
PINS = {
    'results/golden_dyadic_r5_chamber_period_preparation_20261020.md':
        '27db896fd103ff7e0397596beacfb1a45d33994223768336d2bac627ce124470',
    PAPER: 'e2aeff93b42f87673e55399804df28a9eee3713ac2e7cf064c6f885e5bfce671',
    'results/golden_dyadic_unit_vertex_roofs_20261020.md':
        'e0c6b78e7c9c9ece0b77ef75cdfdcdc54d25c5177b1b474a3bcbead9515653e6',
    'results/golden_dyadic_a3_b2_b3_actual_edge_20261020.md':
        '400ba9e7fb1369816a9cdfa10d51a22a75a5babec0a3b5aa7f128e2fbb9327dc',
    'results/review_golden_dyadic_r5_chamber_period_root_20261020.md':
        '74de10e5edb7a0348e4493fa3d2500d68f5f362afa274048ebfcf15738e963a4',
    'results/review_golden_dyadic_r5_chamber_period_batch_guard_20261020.md':
        '81bea9b58ddc9a2cab61146adbb2308b8eea71dd0b6c60ca42ac3fb06df32a14',
    'scripts/check_arithmetic_dn_opposed_fork_models_20261019.py':
        '185d20648368696429b675c1b6c023d30bb1be262ff83d3622e4f48801e1bae9',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
    'scripts/build_illustrated_ra_manuscript_20261018.py':
        '27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15',
    'scripts/run_lean_batches_signal_safe_20261018.py':
        'c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d',
}


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def check():
    sources = {p: digest(ROOT/p) for p in PINS}
    assert sources == PINS, 'pinned proof/control source changed'
    # The freeze's own hash is literal; all 25 enclosed source hashes are
    # independently checked, including its original generation/replay ledgers.
    assert digest(ROOT/FREEZE) == FREEZE_SHA
    sources[FREEZE] = FREEZE_SHA
    packet = json.loads((ROOT/FREEZE).read_text())
    assert packet['file_count'] == len(packet['files']) == 25
    assert len({row['path'] for row in packet['files']}) == 25
    for row in packet['files']:
        path = ROOT/row['path']
        assert path.stat().st_size == row['bytes']
        assert digest(path) == row['sha256']
        if row['path'] in sources:
            assert sources[row['path']] == row['sha256']
        sources[row['path']] = row['sha256']
    own = Path(__file__).relative_to(ROOT).as_posix()
    sources[own] = digest(ROOT/own)
    Q = field(5)
    z, phi = Q(0), Q(F(1, 2), F(1, 2))
    p = [Q(1)]
    for _ in range(10):
        p.append(p[-1]*phi)
    ip = [Q(1)/v for v in p]
    assert phi*phi == phi+1 and phi.conjugate().sign() < 0
    sqrt5 = 2*phi-1

    def unit(i):
        return [Q(int(i == j)) for j in range(10)]

    def add(*vectors):
        return [sum((v[j] for v in vectors), z) for j in range(10)]

    def mul(a, v):
        return [a*x for x in v]

    def transpose(a):
        return [list(row) for row in zip(*a)]

    def mm(a, b):
        return [[sum((x*y for x, y in zip(row, col)), z)
                 for col in zip(*b)] for row in a]

    def mv(a, v):
        return [sum((x*y for x, y in zip(row, v)), z) for row in a]

    def ident(n):
        return [[Q(int(i == j)) for j in range(n)] for i in range(n)]

    def enc(v):
        return [x.encoded() for x in v]

    def em(a):
        return [enc(row) for row in a]

    g3 = [[Q(4), Q(-2), Q(-2)], [Q(-2), Q(4), z],
          [Q(-2), z, phi+2]]
    G = [[z for _ in range(10)] for _ in range(10)]
    G[0][0] = -phi
    for i in range(1, 7):
        G[i][i] = Q(1)
    for i in range(3):
        for j in range(3):
            G[i+7][j+7] = g3[i][j]
    Gi = inverse(G, Q)  # Nonzero block pivots; no indefinite pivot search.

    def dot(v, w):
        return sum((v[i]*sum((G[i][j]*w[j] for j in range(10)), z)
                    for i in range(10)), z)

    def pd(roots):
        a = [[dot(v, w) for w in roots] for v in roots]
        pivots = []
        for i in range(len(a)):
            pivot = a[i][i]
            assert pivot.sign() > 0
            pivots.append(pivot)
            for j in range(i+1, len(a)):
                for k in range(i+1, len(a)):
                    a[j][k] = a[j][k]-a[j][i]*a[i][k]/pivot
        return pivots

    # Exact equality with <-phi>+I5+D4(O), checked in both directions.
    embedding = [
        [p[3]]+[z]*5+[-p[2], -phi, z, -p[3]],
        [-phi]+[z]*5+[phi, z, z, phi],
    ]+[unit(i) for i in range(1, 6)]+[
        mul(Q(2), unit(8)), [z]*6+[Q(-1), Q(-1), Q(-1), Q(1)],
        [-p[3]]+[z]*5+[p[2], p[2], Q(-1), p[3]]]
    E = transpose(embedding)
    original_metric = [[Q(0) for _ in range(10)] for _ in range(10)]
    original_metric[0][0] = -phi
    for i in range(1, 10):
        original_metric[i][i] = Q(1)
    assert mm(mm(transpose(E), original_metric), E) == G
    Ei = mm(mm(Gi, transpose(E)), original_metric)
    assert mm(Ei, E) == mm(E, Ei) == ident(10)

    def integral_original(v):
        return all(x.integral() for x in v) and (sum(v[6:], z)/2).integral()

    standard = [unit(i) for i in range(6)]
    for i, j, sign in ((6, 7, -1), (7, 8, -1), (8, 9, -1), (8, 9, 1)):
        standard.append(add(unit(i), mul(Q(sign), unit(j))))
    assert all(integral_original(v) for v in embedding+standard)
    assert all(x.integral() for v in standard for x in mv(Ei, v))
    assert (-phi.conjugate()).sign() > 0

    V0, a3, s, wp = unit(0), unit(7), unit(8), unit(9)
    es = [unit(i) for i in range(1, 7)]
    beta = [add(mul(Q(-1), es[i]), es[i+1]) for i in range(5)]+[mul(Q(-1), es[5])]
    # Recover omega3 in the OLD (a3,s,r) basis, independently of its source use.
    r = add(mul(-(phi+2)*ip[3], a3), mul(-2*ip[1], s), mul(-2*ip[3], wp))
    old3 = [a3, s, r]
    old3g = [[dot(v, w) for w in old3] for v in old3]
    assert old3g == [[Q(4), Q(-2), z], [Q(-2), Q(4), -2*phi], [z, -2*phi, Q(4)]]
    old3i = inverse(old3g, Q)
    omega3 = add(*(mul(old3i[2][i], old3[i]) for i in range(3)))
    assert [dot(omega3, v) for v in old3] == [z, z, Q(1)]
    p0 = [p[3]]+[p[2]]*4+[z]*5
    u = [ip[1], phi]+[z]*8
    a2, b0 = add(V0, wp), add(p0, wp)
    N = add([p[3]]+[phi]*4+[z]*5, mul(-2*phi, omega3))
    assert dot(N, N) == 2 and integral_original(mv(E, N))
    old_h4, new_h4 = [u]+beta[:3], [N, beta[3], beta[2], beta[1]]
    C = [[dot(v, w) for w in old_h4] for v in old_h4]
    assert C == [[Q(2), -phi, z, z], [-phi, Q(2), Q(-1), z],
                 [z, Q(-1), Q(2), Q(-1)], [z, z, Q(-1), Q(2)]]
    assert [[dot(v, w) for w in new_h4] for v in new_h4] == C
    pd(old_h4)
    Ci = inverse(C, Q)
    assert all(x.integral() for row in Ci for x in row)
    assert Ci == [[4*p[4], 3*p[5], 2*p[5], p[5]],
                  [3*p[5], 6*p[4], 4*p[4], 2*p[4]],
                  [2*p[5], 4*p[4], 2*(3-phi)*p[4], (3-phi)*p[4]],
                  [p[5], 2*p[4], (3-phi)*p[4], 2*p[2]]]

    def projection(v, roots):
        cg = inverse([[dot(a, b) for b in roots] for a in roots], Q)
        coeff = mv(cg, [dot(v, a) for a in roots])
        return add(*(mul(coeff[i], a) for i, a in enumerate(roots)))

    zeta = add(es[0], mul(Q(-1), projection(es[0], new_h4)))
    eta = add(es[4], mul(Q(-1), projection(es[4], new_h4)))
    assert all(x.integral() for x in zeta+eta)
    P1, g = add(eta, mul(Q(-1), zeta)), add(mul(2*ip[1], eta), mul(Q(-1), zeta))
    assert [dot(zeta, zeta), dot(eta, eta), dot(zeta, eta)] == [-32*phi-19, -18*phi-11, -24*phi-15]
    assert [dot(P1, P1), dot(g, g), dot(P1, g)] == [-2*phi, Q(1), z]
    L = add(mul(Q(3), beta[3]), mul(Q(2), beta[2]), beta[1])
    assert P1 == add(mul(Q(-1), es[0]), es[4], mul(p[5], N), mul(p[3], L))
    assert g == add([2*p[4]]+[p[3]]*5+[z]*4, mul(-4*p[2], omega3))
    assert mv(E, g) == [5*phi+2]+[p[3]]*4+[z, -p[4], z, z, -p[4]]
    Qroot = add(g, es[5])
    assert dot(Qroot, Qroot) == 2 and g == add(Qroot, beta[5])
    assert dot(V0, P1) == -p[9] and mv(E, P1)[0].sign() > 0
    B3 = [a2, a3, s]
    assert [[dot(v, w) for w in B3] for v in B3] == [
        [Q(2), Q(-2), z], [Q(-2), Q(4), Q(-2)], [z, Q(-2), Q(4)]]
    W = new_h4+B3+[beta[5]]
    pd(W)
    assert all(dot(P1, v) == dot(g, v) == 0 for v in new_h4+B3)
    assert dot(P1, es[5]) == dot(g, es[5]) == 0
    S = add(V0, mul(Q(-1), projection(V0, W)))
    assert S == add(mul(p[8]/2, P1), mul(-2*p[5], g))
    ts = 4*ip[3]
    assert S == mul(p[8]/2, add(P1, mul(-ts, g)))
    assert dot(S, S).sign() < 0 and mv(E, S)[0].sign() > 0

    # Six fixed finite profiles; positive white projections and all contact
    # squares checked independently against the retained B3+A1 Gram.
    profiles = [('isolated1', 1, 0, 0, 0), ('isolated2', 2, 0, 0, 0),
                ('isolated4', 4, 0, 0, 0), ('B2', 2, 1, 0, 0),
                ('F4', 2, 0, 1, 0), ('B4', 4, 0, 0, 2)]
    B3i = inverse([[dot(v, w) for w in B3] for v in B3], Q)
    white_data = []
    for name, q, c, f, h in profiles:
        v = add(*(mul(-f*B3i[0][i]-h*B3i[2][i], B3[i]) for i in range(3)))
        kappa = dot(v, v)
        assert kappa == F(3, 2)*f*f+F(1, 2)*h*h+f*h
        assert (q-c*c-kappa).sign() > 0
        assert [dot(v, v0) for v0 in B3] == [Q(-f), z, Q(-h)]
        assert dot(v, P1) == dot(v, g) == dot(v, es[5]) == 0
        # Cross term vanishes, so q=-2phi a²+b²+c²+kappa identically.
        white_data.append({'role': name, 'q': q, 'c': c, 'kappa': kappa.encoded()})
    assert p[8]/4-2*phi == F(13, 4)*p[2]
    assert ts*p[4]/2 == 2*phi
    assert (Q(F(121, 64))-Q(64)/(13*p[2])).sign() > 0
    assert (Q(F(324, 25))-8*phi).sign() > 0
    assert (3*sqrt5-Q(F(11, 8))-Q(F(18, 5))).sign() > 0
    assert (2*sqrt5-Q(F(18, 5))).sign() > 0
    assert (4*sqrt5-phi*F(18, 5)+phi.conjugate()*F(11, 8)).sign() > 0
    strips, strip_checks = {}, 0
    for q in (1, 2, 4):
        found = []
        for m in range(-4, 5):
            for n in range(-2, 2):
                A = m+n*phi
                strip_checks += 1
                if A.sign() >= 0:
                    continue
                if (A*A-Q(16*q)/(13*p[2])).sign() > 0:
                    continue
                if (A.conjugate()*A.conjugate()-2*q*phi).sign() > 0:
                    continue
                found.append(A)
        expected = [-ip[1]] if q < 4 else [Q(-1), -ip[1], -ip[2], -2*ip[1]]
        assert sorted(a.encoded() for a in found) == sorted(a.encoded() for a in expected)
        strips[q] = found
    assert strip_checks == 108
    earlier = []
    for name, q, c, f, h in profiles:
        kappa = Q(F(3, 2)*f*f+F(1, 2)*h*h+f*h)
        for A in strips[q]:
            a = A/2
            if q == 4 and not a.integral():
                continue
            b2 = q+2*phi*a*a-c*c-kappa
            if q < 4:
                assert not b2.integral()
                reason = 'b_square_not_integral'
            else:
                assert a == -ip[1] and not (b2/4).integral()
                reason = 'b_square_over_four_not_integral'
            earlier.append({'role': name, 'A': A.encoded(), 'b_square': b2.encoded(), 'exclusion': reason})

    def rational_square(a):
        return a >= 0 and math.isqrt(a.numerator)**2 == a.numerator and math.isqrt(a.denominator)**2 == a.denominator

    assert all(not rational_square(a) for a in (F(2), F(2, 5), F(1, 2), F(1, 10)))
    at_zero = []
    for name, q, c, f, h in profiles:
        b2 = q-c*c-Q(F(3, 2)*f*f+F(1, 2)*h*h+f*h)
        at_zero.append({'role': name, 'b_square': b2.encoded()})
    assert [row['b_square'] for row in at_zero] == [Q(v).encoded() for v in (1, 2, 4, 1, F(1, 2), 2)]
    assert dot(mul(Q(2), g), mul(Q(2), g)) == 4
    assert g[6] == 0 and Qroot[6] == 1
    assert all(x.integral() for x in g+Qroot)
    for normal in [g, Qroot]:
        q = dot(normal, normal)
        assert dot(normal, V0).sign() < 0 and dot(normal, S).sign() < 0 and dot(normal, P1) == 0
        assert all((2*dot(normal, unit(i))/q).integral() for i in range(10))

    Pold = add(mul(p[2], p0), mul(Q(-2), omega3))
    old_corner = old_h4+[b0, a3, s, beta[4], beta[5]]
    new_corner = new_h4+[a2, a3, s, Qroot, beta[5]]
    old_pivots, new_pivots = pd(old_corner), pd(new_corner)
    assert [[dot(v, w) for w in old_corner] for v in old_corner] == [
        [dot(v, w) for w in new_corner] for v in new_corner]
    assert dot(Pold, Pold) == dot(P1, P1) == -2*phi
    assert all(dot(Pold, v) == 0 for v in old_corner)
    assert all(dot(P1, v) == 0 for v in new_corner)
    assert mv(E, Pold)[0].sign() > 0
    assert p0 == add(mul(p[2]/2, Pold), mul(-phi/2, add(mul(Q(3), b0), mul(Q(2), a3), s)))
    Dp0 = add(mul(p[2]/2, P1), mul(-phi/2, add(mul(Q(3), a2), mul(Q(2), a3), s)))
    assert Dp0 == [26*phi+17, 2*p[5]]+[p[6]]*4+[z, 4*p[4], 2*p[4], 6*p[4]]
    assert p0 == add(mul(ip[3], V0), mul(4*phi, u), mul(p[2], add(mul(Q(3), beta[0]), mul(Q(2), beta[1]), beta[2])))
    DV0 = mul(p[3], add(Dp0, mul(-4*phi, N), mul(-p[2], L)))
    assert DV0 == [60*phi+37, 22*phi+14]+[19*phi+12]*4+[z, 28*phi+16, 14*phi+8, 42*phi+24]
    fvec = add(mul(ip[1], N), mul(-ip[2], DV0))
    columns = [DV0, fvec, add(fvec, beta[3]), add(fvec, beta[3], beta[2]),
               add(fvec, beta[3], beta[2], beta[1]), mul(Q(-1), g), es[5],
               a3, s, add(a2, mul(Q(-1), Dp0))]
    D = transpose(columns)

    def parse(token):
        if 'p' not in token:
            return Q(int(token))
        match = re.fullmatch(r'([+-]?\d*)p([+-]\d+)?', token)
        assert match is not None, 'unsafe/unexpected matrix token'
        coefficient, constant = match.groups()
        coefficient = {'': 1, '+': 1, '-': -1}.get(coefficient, coefficient)
        return int(coefficient)*phi+int(constant or 0)

    text = (ROOT/PAPER).read_text()
    blocks = re.findall(r'```text\n(.*?)\n```', text, flags=re.S)
    assert len(blocks) == 1
    parsed = [[parse(v) for v in row.split()] for row in blocks[0].strip().splitlines()]
    assert len(parsed) == 10 and all(len(row) == 10 for row in parsed)
    assert D == parsed and all(x.integral() for row in D for x in row)
    assert mm(mm(transpose(D), G), D) == G
    Di = mm(mm(Gi, transpose(D)), G)
    assert all(x.integral() for row in Di for x in row)
    assert mm(Di, D) == mm(D, Di) == ident(10)
    assert all(mv(D, v) == w for v, w in zip(old_corner, new_corner))
    assert mv(D, Pold) == P1 and mv(D, p0) == Dp0
    assert mv(D, es[4]) == mul(Q(-1), g) and mv(D, es[5]) == es[5]
    for normal in old_corner+new_corner:
        q = dot(normal, normal)
        assert q in (1, 2, 4) and all(x.integral() for x in normal)
        assert all((2*dot(normal, unit(i))/q).integral() for i in range(10))
    trace = sum((D[i][i] for i in range(10)), z)
    assert trace == 10+11*phi and (trace-10).sign() > 0
    assert all(digest(ROOT/path) == pin for path, pin in sources.items()), 'source changed during check'
    return {
        'status': 'FIXED_R5_INTERVAL_AND_CHAMBER_PERIOD_PASS', 'owns_lock': False,
        'sources': sources, 'frozen_packet_file_count': 25,
        'field_encoding': '[a,b] means a+b*sqrt5; rational strings',
        'basis': ['V0', 'e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'a3', 's', 'wP'],
        'basis_in_original_coordinates': em(embedding), 'basis_gram': em(G),
        'standard_original_basis_in_B': [enc(mv(Ei, v)) for v in standard],
        'H4_gram': em(C), 'H4_gram_inverse': em(Ci),
        'zeta': enc(zeta), 'eta': enc(eta), 'P1': enc(P1), 'g': enc(g), 'Q': enc(Qroot),
        'S': enc(S), 't_s': ts.encoded(), 'profiles': white_data,
        'coefficient_pairs_checked': strip_checks,
        'negative_A_strips': {str(q): enc(v) for q, v in strips.items()},
        'earlier_profile_exclusions': earlier, 'a_zero_profile_squares': at_zero,
        'simultaneous_identity': 'g=Q+beta6; g is nonsimple, Q unique simple',
        'old_corner': [enc(v) for v in old_corner], 'new_corner': [enc(v) for v in new_corner],
        'old_corner_pd_pivots': enc(old_pivots), 'new_corner_pd_pivots': enc(new_pivots),
        'old_apex': enc(Pold), 'new_apex': enc(P1),
        'D': em(D), 'D_inverse': em(Di), 'D_original': em(mm(mm(E, D), Ei)),
        'trace': trace.encoded(), 'future_orientation': 'D(P_old)=P1, both future timelike',
        'written_geometry_inputs': 'complete profiles, actual-wall/full-germ provenance and chamber matching remain pinned paper proofs',
        'written_implication': 'integral future chamber stabilizer with trace>10 has infinite order; literal R5 chamber is noncompact/infinite volume',
        'height_searches': 0,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--write', action='store_true')
    modes.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.write and CERT.exists():
        parser.error('refusing to overwrite certificate')
    assert int(os.environ.get('COXETER_CHECK_OWNER', '-1')) == os.getppid(), 'bounded parent wrapper required'
    signal.alarm(45)
    data = check()
    payload = json.dumps(data, sort_keys=True, indent=2)+'\n'
    assert len(payload.encode()) <= 131072, 'certificate exceeds128KiB'
    if args.write:
        with CERT.open('x') as stream:
            stream.write(payload)
    else:
        assert CERT.stat().st_size <= 131072
        assert data == json.loads(CERT.read_text()), 'fresh result differs from certificate'
    print(json.dumps({'status': 'PASS', 'mode': 'write' if args.write else 'verify',
                      'coefficient_pairs': 108, 'trace': '10+11phi',
                      'certificate_bytes': len(payload.encode()),
                      'certificate_sha256': hashlib.sha256(payload.encode()).hexdigest()}))


if __name__ == '__main__':
    main()
