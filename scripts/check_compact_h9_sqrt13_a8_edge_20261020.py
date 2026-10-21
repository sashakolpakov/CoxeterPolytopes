#!/usr/bin/env python3
"""Fixed quadratic13 A8-edge algebra; outer wrapper owns lock/resources.

No work runs on import. This checks the tiny complete written interval,
not a height prefix. Root schedules generation and a separate replay.
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
CERT = ROOT/'results/compact_h9_sqrt13_a8_edge_20261020.json'
HELPER = 'scripts/check_compact_h9_sqrt13_two_edges_20261020.py'
PINS = {
    HELPER: '6839ac3b07c8bc87056712692e526deacc27bdcc13d33f7fb07a72fe8e64af31',
    'results/compact_h9_sqrt13_b9_two_edge_route_20261020.md': 'b922647902d32beccd98fbe5a2781ddcbb9cafaeb10a12902eb3ebbce809547f',
    'results/review_compact_h9_sqrt13_two_edges_resource_20261020.md': '20b3381d2d20ff7efb0d4df7b0642d1b8dc0eecb1924005bc0502878467f96de',
    'results/compact_h9_sqrt13_a8_edge_interval_20261020.md': 'e3e62a6ec4d834e2e0c7d4abc4d28961839f1a467ac06966ef08854757814228',
    'results/review_compact_h9_sqrt13_a8_edge_resource_20261020.md': 'ad0325a20866f248edca7ac77864f434adb85723af2761ae61dd5bf01d051128',
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
    spec = importlib.util.spec_from_file_location('pinned_q13_definitions', ROOT/HELPER)
    helper = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = helper
    spec.loader.exec_module(helper)  # definitions only; never helper.check/main
    Q, enc = helper.Q, helper.enc
    add, scale, mv = helper.add, helper.scale, helper.mv
    omega, one = Q(0, 1), Q(1)
    eta, zero = omega+1, Q()
    sq = lambda x: x*x
    assert sq(omega) == omega+3 and sq(eta) == 3*eta+1
    assert eta.conjugate() == -1/eta and eta.norm() == -1
    assert omega.integral() and not (omega/2).integral()
    basis = [[Q(int(i == j)) for i in range(10)] for j in range(10)]
    gram = [[zero for _ in range(10)] for _ in range(10)]
    for i in range(10):
        gram[i][i] = -eta if i == 0 else one
    dot = lambda x, y: sum(a*b for a, b in zip(x, mv(gram, y)))
    p = basis[0]
    alpha = [add(scale(-1, basis[i]), basis[i+1]) for i in range(1, 9)]
    alpha.append(scale(-1, basis[9]))
    v = add(p, scale(omega, basis[1]))
    u = add(scale(eta, p), scale(eta+1, add(basis[1], basis[2])))
    d = add(u, basis[3])
    old = alpha+[v, d]
    e = [zero]+[one]*9
    xa, Za = 2*eta-1, 9*eta+1
    a = add(scale(xa, p), add(scale(eta, e), basis[1]))
    assert dot(a, a) == 2 and all(x.integral() for x in a)
    assert eta in a and eta.norm() == -1  # primitive coordinate ideal
    cov = mv(gram, a)
    assert all((x*y).integral() for x in a for y in cov)  # norm2 reflection
    assert sq(Za)-9*eta*sq(xa) == 10
    tau = eta*xa/Za
    assert xa.sign() > 0 and Za.sign() > 0 and (eta/9-sq(tau)).sign() > 0

    # Eight literal placements. The displayed coefficients are those of
    # (r,r)-q in1,x,Z,x^2,Z^2,xZ; no finite classification is inferred here.
    profiles = [('isolated_short', 1, 0, 9), ('isolated_long', 2, 0, 18),
                ('B9_at_1', 1, 1, 1), ('B9_at_8', 1, 8, 1),
                ('D9_at_2', 2, 2, 4), ('D9_at_7', 2, 7, 4),
                ('A9_at_1', 2, 1, 10), ('A9_at_8', 2, 8, 10)]
    profile_rows = []
    for name, norm, k, R in profiles:
        const = [zero]+[Q(F(9-k if i < k else -k, 9)) for i in range(9)]
        zdir = scale(F(1, 9), e)
        co = [dot(const, const)-norm, 2*dot(const, p), 2*dot(const, zdir),
              dot(p, p), dot(zdir, zdir), 2*dot(p, zdir)]
        assert co == [Q(F(-R, 9)), zero, zero, -eta, Q(F(1, 9)), zero]
        assert R == 9*norm-k*(9-k)
        assert [dot(const, r) for r in alpha[:8]] == [Q(-int(j+1 == k)) for j in range(8)]
        profile_rows.append({'name': name, 'norm': norm, 'k': k, 'R': R,
                             'norm_coefficients': list(map(enc, co))})

    # Isolated roots have native9z^2-eta*x^2=q, not just the Z equation.
    mod3 = []
    for q, w in [(1, 0), (2, 1)]:
        assert (w*w-w-3) % 3 == 0
        eta_residue = (w+1) % 3
        assert all((-eta_residue*x*x-q) % 3 != 0 for x in range(3))
        mod3.append({'norm': q, 'omega_residue': w, 'eta_residue': eta_residue})

    # Paper coarse bounds and complete coefficient box, before square tests.
    assert (eta-3).sign() > 0 and (Q(F(10, 3))-eta).sign() > 0
    s = 2*omega-1
    assert (s-Q(F(18, 5))).sign() > 0 and sq(s) == 13
    assert sq(xa) == 5+8*eta
    chains = [(sq(xa)/10, Q(F(19, 6)), Q(F(81, 25))),
              (eta/9, Q(F(10, 27)), Q(F(25, 64))),
              (2*sq(xa)/5, Q(F(38, 3)), Q(F(324, 25))),
              (4*eta/9, Q(F(40, 27)), Q(F(25, 16)))]
    assert all((middle-left).sign() > 0 and (right-middle).sign() > 0
               for left, middle, right in chains)
    assert (6-xa).sign() > 0 and (4-10*eta/9).sign() > 0
    # Thus0<x<6,|x'|<2 gives B=0,1,2; -4<=A<=5 since2omega<5.
    assert (s-2).sign() > 0 and (3*s-8).sign() > 0 and (5-2*omega).sign() > 0
    expected = {1: [], 4: [(1, 0), (1, 1)],
                10: [(1, 0), (0, 1), (1, 1), (2, 1), (3, 1), (1, 2)]}
    domains = {}
    for R, values in expected.items():
        found = []
        for B in range(3):
            for A in range(-4, 6):
                x = Q(A, B)
                if x.sign() > 0 and (R*sq(xa)-10*sq(x)).sign() >= 0 and (
                        R*eta/9-sq(x.conjugate())).sign() >= 0:
                    found.append((A, B))
        assert found == values
        domains[str(R)] = found

    expected_norms = {4: [43, 1231], 10: [289, 901, 3259, 4861, 1819, 2809]}
    square_rows, candidates = [], []
    for R in (4, 10):
        for index, pair in enumerate(expected[R]):
            x = Q(*pair)
            rad = 9*eta*sq(x)+R
            norm = rad.norm()
            assert norm.denominator == 1
            n = int(norm)
            assert n == expected_norms[R][index]
            row = {'R': R, 'x': enc(x), 'radicand': enc(rad), 'rational_norm': n}
            if isqrt(n)**2 != n:
                row['outcome'] = 'nonsquare rational norm'
            else:
                assert R == 10 and index in (0, 5)
                Z = eta+3 if index == 0 else Za
                assert Z.sign() > 0 and sq(Z) == rad
                row['positive_square_root'] = enc(Z)
                row['integral_placements'] = []
                for k in (1, 8):
                    z = (Z-k)/9
                    if not z.integral():
                        continue
                    root = [x]+[z+int(i < k) for i in range(9)]
                    assert all(y.integral() for y in root) and dot(root, root) == 2
                    row['integral_placements'].append(k)
                    candidates.append((k, root))
            square_rows.append(row)
    assert len(candidates) == 1 and candidates[0] == (1, a)
    assert sq(eta+3) == 9*eta+10 and sq(Za) == 82+261*eta
    assert all(not ((eta+3-k)/9).integral() for k in (1, 8))
    assert ((Za-1)/9).integral() and not ((Za-8)/9).integral()

    # Postselection native/cap products, never used as interval filters.
    products = [dot(a, r) for r in old]
    assert products == [Q(-1)]+[zero]*7+[-eta, -2*eta-2, -7*eta-2]
    pair_rows = []
    for i, r in enumerate(old):
        product = products[i]
        q = dot(r, r)
        assert q == (one if i == 8 else Q(2))
        ratio4 = 4*sq(product)/(2*q)
        if i == 0:
            assert ratio4 == one
            label = 3
        elif 1 <= i <= 7:
            assert ratio4 == zero
            label = 2
        else:
            assert product.sign() < 0 and (ratio4-4).sign() > 0
            label = 'ultraparallel'
        pair_rows.append([i, enc(product), label])

    apex = [Za]+[5*eta+2]*9
    assert eta*xa == 5*eta+2 and apex == add(scale(Za, p), scale(eta*xa, e))
    assert dot(apex, apex) == -10*eta and dot(p, apex).sign() < 0
    all_roots = old+[a]
    assert all(dot(apex, r).sign() <= 0 for r in all_roots)
    assert [i for i, r in enumerate(all_roots) if dot(apex, r) == zero] == list(range(8))+[11]
    assert [dot(apex, r) for r in (alpha[8], v, d)] == [-5*eta-2, -16*eta-6, -44*eta-12]
    corner = [a]+alpha[:8]
    corner_gram = [[dot(x, y) for y in corner] for x in corner]
    assert corner_gram == [[Q(2 if i == j else -1 if abs(i-j) == 1 else 0)
                           for j in range(9)] for i in range(9)]
    helper.positive_definite(corner_gram)
    # A9 principal determinant recurrence D_n=2D_(n-1)-D_(n-2).
    prev, det = 1, 2
    for _ in range(2, 10):
        prev, det = det, 2*det-prev
    assert det == 10
    center = [one]+[Q(F(i, 100)) for i in range(9, 0, -1)]
    assert dot(center, center).sign() < 0 and dot(center, p).sign() < 0
    assert all(dot(center, r).sign() < 0 for r in all_roots)
    assert sources == {p: digest(ROOT/p) for p in paths}, 'source changed during check'
    return {'status': 'PASS_FIXED_COMPLETE_A8_INTERVAL_AND_A9_ENDPOINT_ALGEBRA',
            'sources': sources, 'owns_lock': False,
            'scope': 'primitive alphabet, spherical inventory and full actual-germ interpretation remain written inputs; no volume/reflectivity claim',
            'field_encoding': 'a+b*omega; rational strings', 'profiles': profile_rows,
            'isolated_native_mod3_obstructions': mod3, 'complete_domains': domains,
            'square_table': square_rows, 'unique_wall': list(map(enc, a)),
            'new_wall_products': pair_rows, 'first_parameter': enc(tau),
            'actual_endpoint': list(map(enc, apex)), 'corner_gram': [[enc(x) for x in row] for row in corner_gram],
            'corner_determinant': det, 'strict_common_center': list(map(enc, center)),
            'new_example': False, 'finite_volume_asserted': False}


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
                      'profiles': 8, 'unique_first_wall': True, 'corner': 'A9'}))


if __name__ == '__main__':
    main()
