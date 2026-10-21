#!/usr/bin/env python3
"""Fixed Q13 far-edge/glue/period algebra; root schedules all executions.

Use outer --parent-lock. No work on import, height search or cone census.
Actual germs and the chamber/index argument remain written premises.
"""
import argparse
import hashlib
import importlib.util
import json
import signal
import sys
from fractions import Fraction as F
from math import isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT/'results/compact_h9_sqrt13_far_end_period_20261020.json'
HELPER = 'scripts/check_compact_h9_sqrt13_two_edges_20261020.py'
FIRST = 'results/compact_h9_sqrt13_two_edges_20261020.json'
SECOND = 'results/compact_h9_sqrt13_a8_edge_20261020.json'
PINS = {
    HELPER: '6839ac3b07c8bc87056712692e526deacc27bdcc13d33f7fb07a72fe8e64af31',
    'scripts/check_compact_h9_sqrt13_a8_edge_20261020.py': '63dd91e189e3534024b6e2a2223adfacd5633c4a65cd87ef56177dea0a340bda',
    FIRST: '7608882aae2a4d0ed050cf10753552d61cb02cb7b7ead8774461ddcb04b9a753',
    SECOND: '80fe12cdb1f4818d6107333385eef9a70c81ba45a40501b5a9ad910a58fc26c9',
    'results/compact_h9_sqrt13_b9_two_edge_route_20261020.md': 'b922647902d32beccd98fbe5a2781ddcbb9cafaeb10a12902eb3ebbce809547f',
    'results/review_compact_h9_sqrt13_two_edges_resource_20261020.md': '20b3381d2d20ff7efb0d4df7b0642d1b8dc0eecb1924005bc0502878467f96de',
    'results/compact_h9_sqrt13_a8_edge_interval_20261020.md': 'e3e62a6ec4d834e2e0c7d4abc4d28961839f1a467ac06966ef08854757814228',
    'results/review_compact_h9_sqrt13_a8_edge_resource_20261020.md': 'ad0325a20866f248edca7ac77864f434adb85723af2761ae61dd5bf01d051128',
    'results/compact_h9_sqrt13_a9_far_end_and_period_20261020.md': '86ffab56f4696298fac4ae34b8bc8772dd1e7a6daf2e95e48d56291e810bbcf2',
    'results/review_compact_h9_sqrt13_far_end_period_resource_20261020.md': '10bbac3e674717ec24698e8c1b7ce76326980025e3ef246a15e3c0ae32d5c934',
    'results/review_compact_h9_sqrt13_far_edge_root_20261020.md': 'fdf2f9b40e5700c368489cf912c27f85f376700a091b4a3d7407c68526c42447',
    'results/review_compact_h9_sqrt13_far_end_period_batch_guard_20261020.md': '9a402457ad5688594ec864c9105aeed7c6da4072c114dfe042bbd795f9eb359a',
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
    spec = importlib.util.spec_from_file_location('pinned_q13_period_definitions', ROOT/HELPER)
    helper = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = helper
    spec.loader.exec_module(helper)  # pure definitions; never helper.check/main
    Q, enc = helper.Q, helper.enc
    add, scale, mv, mm, tr = helper.add, helper.scale, helper.mv, helper.mm, helper.transpose
    omega, zero = Q(0, 1), Q()
    eta, one = omega+1, Q(1)
    sq = lambda x: x*x
    basis = [[Q(int(i == j)) for i in range(10)] for j in range(10)]
    gram = [[zero for _ in range(10)] for _ in range(10)]
    for i in range(10):
        gram[i][i] = -eta if i == 0 else one
    dot = lambda x, y: sum(a*b for a, b in zip(x, mv(gram, y)))
    alpha = [add(scale(-1, basis[i]), basis[i+1]) for i in range(1, 9)]
    alpha.append(scale(-1, basis[9]))
    p = basis[0]
    v = add(p, scale(omega, basis[1]))
    u = add(scale(eta, p), scale(eta+1, add(basis[1], basis[2])))
    d = add(u, basis[3])
    a = [2*eta-1, eta+1]+[eta]*8
    apex = [9*eta+1]+[5*eta+2]*9
    beta = [a]+alpha[:7]
    P0 = [5*eta+1]+[sq(eta)]*8+[2*eta+1]
    R0 = [7*eta+1]+[4*eta+2]*8+[3*eta-1]
    P = scale(1/eta, P0)
    assert dot(P0, P0) == -eta and dot(R0, R0) == 9 and dot(P0, R0) == 0
    assert dot(P, P) == -1/eta and all(x.integral() for x in P)
    assert P0[0]*R0[1]-P0[1]*R0[0] == eta
    assert eta.norm() == -1 and all(dot(x, y) == zero for x in (P0, R0) for y in beta)
    assert apex == add(scale(eta+3, P0), scale(-eta, R0))
    endpoint = add(scale(eta+3, P0), scale(eta, R0))
    assert endpoint == [53*eta+15]+[33*eta+10]*8+[21*eta+8]
    assert dot(apex, apex) == dot(endpoint, endpoint) == -10*eta

    # Native two-variable plane and all six attached placements.
    xd = [one]+[zero]*8+[2*eta-1]
    zd = [zero]+[one]*8+[-eta-5]
    assert xd == add(scale(2*omega, P0), scale(-eta, R0))
    assert zd == add(scale(-eta-4, P0), scale(eta+2, R0))
    profile_rows = []
    attached = [(1, 1, 1), (1, 8, 1), (2, 2, 4), (2, 7, 4), (2, 1, 10), (2, 8, 10)]
    for q, k, R in attached:
        const = [zero]+[Q(int(i < k-1)) for i in range(8)]+[-k-eta+4]
        Ac, Zc = 5-eta-k, (eta+3)*k-9
        Ax, Az, Zx, Zz = 2*omega, -eta-4, -9*eta, 9*(eta+2)
        assert -dot(P, const) == Ac and dot(R0, const) == Zc
        assert -dot(P, xd) == Ax and -dot(P, zd) == Az
        assert dot(R0, xd) == Zx and dot(R0, zd) == Zz
        assert [dot(const, x) for x in beta] == [Q(-int(i+1 == k)) for i in range(8)]
        white = add(const, add(scale(-Ac, P0), scale(-Zc/9, R0)))
        assert dot(white, P0) == dot(white, R0) == zero
        assert dot(white, white) == Q(F(k*(9-k), 9))
        co = [dot(const, const)-q, 2*dot(const, xd), 2*dot(const, zd),
              dot(xd, xd), dot(zd, zd), 2*dot(xd, zd)]
        pell = [sq(Zc)-9*eta*sq(Ac)-R, 2*Zc*Zx-18*eta*Ac*Ax,
                2*Zc*Zz-18*eta*Ac*Az, sq(Zx)-9*eta*sq(Ax),
                sq(Zz)-9*eta*sq(Az), 2*Zx*Zz-18*eta*Ax*Az]
        assert co == [x/9 for x in pell] and R == 9*q-k*(9-k)
        assert Ac+k+eta-5 == zero and (Zc-(eta+3)*k)/9+1 == zero
        profile_rows.append({'norm': q, 'k': k, 'R': R, 'A_constant': enc(Ac),
                             'Z_constant': enc(Zc), 'norm_coefficients': list(map(enc, co))})
    mod3 = []
    for q, w in [(1, 0), (2, 1)]:
        assert (w*w-w-3) % 3 == 0
        assert all((-(w+1)*A*A-q) % 3 != 0 for A in range(3))
        mod3.append([q, w])
    # Signed complete interval: |A|<=1, |A'|<2 => omega coefficient0.
    assert sq(eta+3)-9*eta == 10
    assert (4-10*eta/9).sign() > 0 and (2*omega-4).sign() > 0
    assert (10*eta).sign() > 0
    # sqrt10 outsideK: rational-square alternatives10 and130 both fail.
    assert all(isqrt(n)**2 != n for n in (10, 130))

    def reconstruct(A, Z, k):
        const = [zero]+[Q(int(i < k-1)) for i in range(8)]+[-k-eta+4]
        return add(add(scale(A+k+eta-5, P0),
                       scale((Z-(eta+3)*k)/9+1, R0)), const)

    b0 = add(a, [sum(x[i] for x in alpha[:8]) for i in range(10)])
    b = [8*eta+3]+[5*eta+2]*8+[sq(eta)]
    cases, outgoing = [], []
    for q, k, R in attached:
        for A in (-1, 0, 1):
            if 10*A*A > R:
                continue
            row = {'norm': q, 'k': k, 'R': R, 'A': A}
            if R == 10 and A == 0:
                row['outcome'] = 'sqrt10 outsideK'
            else:
                Z = Q(1) if R == 1 else Q(2) if R == 4 else eta+3
                assert sq(Z)-9*eta*A*A == R and Z.sign() > 0
                row['Z'] = enc(Z)
                if not ((Z-(eta+3)*k)/9).integral():
                    row['outcome'] = 'fails native glue'
                else:
                    r = reconstruct(Q(A), Z, k)
                    assert all(x.integral() for x in r) and dot(r, r) == q
                    assert -dot(P, r) == A and dot(R0, r) == Z
                    assert dot(apex, r) == -eta*((eta+3)*A+Z)
                    assert dot(endpoint, r) == -eta*((eta+3)*A-Z)
                    if dot(apex, r).sign() >= 0:
                        assert A == -1 and r == scale(-1, b0) and dot(apex, r) == zero
                        row['outcome'] = 'known root through startingQ'
                    else:
                        assert A == 1 and k == 1 and r == b
                        outgoing.append(r)
                        row['outcome'] = 'unique first candidate'
            cases.append(row)
    assert len(outgoing) == 1 and outgoing[0] == b
    assert all(x.integral() for x in b) and any(x.norm() in (1, -1) for x in b)
    assert dot(b, b) == 2
    assert all((x*y).integral() for x in b for y in mv(gram, b))
    old = alpha+[v, d, a]
    assert [dot(r, r) for r in old] == [Q(2)]*8+[one, Q(2), Q(2), Q(2)]
    products = [dot(b, r) for r in old]
    assert products == [zero]*7+[-2*eta-1, -sq(eta), -15*eta-5, -40*eta-11, Q(-1)]
    for i, r in enumerate(old):
        value = 4*sq(products[i])/(2*dot(r, r))
        if i <= 6 or i == 11:
            assert value == (one if i == 11 else zero)
        else:
            assert (value-4).sign() > 0
    assert all(dot(endpoint, r).sign() <= 0 for r in old+[b])
    assert [i for i, r in enumerate(old+[b]) if dot(endpoint, r) == zero] == list(range(7))+[11, 12]
    corner = list(reversed(beta))+[b]
    cg = [[dot(x, y) for y in corner] for x in corner]
    assert cg == [[Q(2 if i == j else -1 if abs(i-j) == 1 else 0) for j in range(9)] for i in range(9)]
    assert [[dot(x, y) for y in beta+[alpha[7]]] for x in beta+[alpha[7]]] == cg
    assert dot(apex, alpha[7]) == zero and all(dot(apex, r).sign() <= 0 for r in old)
    helper.positive_definite(cg)

    # Bind actual earlier certificates without running their scientific bodies.
    prior = []
    for path in (FIRST, SECOND):
        assert (ROOT/path).stat().st_size <= 15360
        data = json.loads((ROOT/path).read_text())
        assert all(PINS.get(p) == h for p, h in data['sources'].items())
        prior.append(data)
    assert prior[0]['status'] == 'PASS_FIXED_TWO_EDGES_AND_CHAMBER_INVOLUTION_ALGEBRA'
    assert prior[1]['status'] == 'PASS_FIXED_COMPLETE_A8_INTERVAL_AND_A9_ENDPOINT_ALGEBRA'
    assert prior[0]['roots'] == [list(map(enc, r)) for r in old[:11]]
    assert prior[1]['unique_wall'] == list(map(enc, a)) and prior[1]['actual_endpoint'] == list(map(enc, apex))
    Gcols = [[3*eta, 4*eta+1, 3*eta+1]+[zero]*7,
             [-eta-1, -2*eta, -eta-1]+[zero]*7,
             scale(-1, u)]+basis[3:]
    assert prior[0]['involution_columns'] == [list(map(enc, x)) for x in Gcols]
    G = tr(Gcols)
    I = [[Q(int(i == j)) for j in range(10)] for i in range(10)]
    assert mm(G, G) == I and [[dot(x, y) for y in Gcols] for x in Gcols] == gram
    gammas = [add(a, [sum(x[i] for x in alpha[:7]) for i in range(10)])]
    for lo, hi in [(0, 6), (1, 5), (2, 4)]:
        gammas.append([sum(x[i] for x in alpha[lo:hi]) for i in range(10)])
    assert [[dot(x, y) for y in gammas] for x in gammas] == [[Q(2*int(i == j)) for j in range(4)] for i in range(4)]
    assert all(dot(P0, r) == zero for r in gammas)
    wcols, jcols = [], []
    for e in basis:
        wcols.append(add(e, [sum(-dot(r, e)*r[i] for r in gammas) for i in range(10)]))
        jcols.append(add(scale(-1, e), scale(-2*dot(P0, e)/eta, P0)))
    w0, JP = tr(wcols), tr(jcols)
    H = mm(JP, w0)
    assert mm(w0, w0) == mm(JP, JP) == mm(H, H) == I and mm(w0, JP) == H
    Hcols = tr(H)
    assert all(x.integral() for col in Hcols for x in col)
    assert [[dot(x, y) for y in Hcols] for x in Hcols] == gram
    assert mv(w0, P0) == mv(H, P0) == P0 and mv(H, R0) == scale(-1, R0)
    assert mv(w0, alpha[7]) == b0 and mv(H, alpha[7]) == b
    assert [mv(H, r) for r in beta] == list(reversed(beta))
    assert mv(H, apex) == endpoint and mv(JP, P0) == P0 and P0[0].sign() > 0
    HG = mm(H, G)
    assert all(x.integral() for row in HG for x in row)
    assert mm(HG, mm(G, H)) == I
    assert [[dot(x, y) for y in tr(HG)] for x in tr(HG)] == gram
    gamma_products = [dot(r, mv(G, r)) for r in gammas]
    assert gamma_products == [-63*eta-19, 1-2*eta, -eta, Q(2)]
    assert dot(P, mv(G, P)) == -51*eta-14 and dot(P0, mv(G, P0)) == -552*eta-167
    trace = lambda m: sum(m[i][i] for i in range(10))
    assert trace(G) == 6 and trace(w0) == 2 and trace(H) == 0
    assert trace(mm(w0, G)) == 66*eta+22
    assert trace(HG) == -trace(mm(w0, G))-2*eta*dot(P, mv(G, P)) == 268*eta+80
    assert (trace(HG)-10).sign() > 0 and mv(HG, p)[0].sign() > 0
    center = [one]+[Q(F(i, 100)) for i in range(9, 0, -1)]
    assert dot(center, center).sign() < 0 and all(dot(center, r).sign() < 0 for r in old+[b])
    assert sources == {p: digest(ROOT/p) for p in paths}, 'source changed during check'
    return {'status': 'PASS_FIXED_FAR_A9_EDGE_AND_INFINITE_PERIOD_ALGEBRA',
            'sources': sources, 'owns_lock': False, 'field_encoding': 'a+b*omega; rational strings',
            'scope': 'written primitive alphabet, spherical inventory, actual germs and arithmetic chamber/index proof remain explicit inputs',
            'orthogonal_basis': [list(map(enc, P0)), list(map(enc, R0))],
            'isolated_mod3': mod3, 'attached_profiles': profile_rows, 'complete_signed_cases': cases,
            'unique_wall': list(map(enc, b)), 'old_wall_products': list(map(enc, products)),
            'actual_endpoint': list(map(enc, endpoint)), 'H_columns': [list(map(enc, r)) for r in Hcols],
            'gamma_G_products': list(map(enc, gamma_products)), 'HG_trace': enc(trace(HG)),
            'strict_common_center': list(map(enc, center)), 'infinite_order_period': True,
            'no_new_example': True, 'literal_lattice_scope': '<-eta>+I9 over Z[(1+sqrt13)/2]'}


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
    assert len(payload.encode()) <= 15360
    assert result['sources'] == {p: digest(ROOT/p) for p in result['sources']}
    if args.verify:
        assert CERT.stat().st_size <= 15360
        assert result == json.loads(CERT.read_text()), 'saved certificate differs'
    else:
        with CERT.open('x') as stream:
            stream.write(payload)
    signal.alarm(0)
    print(json.dumps({'status': result['status'], 'mode': 'replay' if args.verify else 'generation',
                      'unique_first_wall': True, 'period_trace': result['HG_trace']}))


if __name__ == '__main__':
    main()
