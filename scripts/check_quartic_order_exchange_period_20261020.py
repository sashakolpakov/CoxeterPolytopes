#!/usr/bin/env python3
"""Fixed quartic five-coordinate identity certificate; no root search.

Run only through run_bounded_research_check_20261018.py --parent-lock.
Actualness, full mirror arrangements, finite Coxeter classification and
the arithmetic lattice theorem are explicit written inputs, not proved
by this matrix check. --write refuses replacement; --verify recomputes.
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

from fresh_sqrt2sqrt5_lanner_catalogue_20261016_review import K

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT/'results/quartic_order_exchange_period_20261020.json'
PINS = {
    'results/compact_h9_biquadratic_unit_diagonal_first_edge_20261020.md':
        '05b6473eeadb4befa89c82af2f7002c7441ada3d1878c93d4630355a42685efe',
    'results/compact_h9_biquadratic_unit_second_edge_20261020.md':
        'b3a7c8f4e2afffe13fb49dd80ee2ed03095608003fabea85e2e6b74c9df2a52b',
    'results/compact_biquadratic_unit_order_exchange_obstruction_20261020.md':
        '043f1f1f6605f6e0b80c1051fd6d288f119e6fbb65487af45d092e1669320299',
    'results/review_compact_biquadratic_order_exchange_root_20261020.md':
        '05d0d24ad21c3bc9ce586b9cec4603aa5dcd1ca4236fe0f6698be9ff9e59b140',
    'results/review_biquadratic_unit_order_exchange_batch_guard_20261020.md':
        '3ebff7f5d68b30b5f71a7ccaa082801dd9c5e8c3b1f7b6350b2a009849071c39',
    'scripts/fresh_sqrt2sqrt5_lanner_catalogue_20261016_review.py':
        '122cfeea5190517d9203a61d933fbe964485652be04c5fdca157e492eca2c811',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
}


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def check():
    sources = {name: digest(ROOT/name) for name in PINS}
    assert sources == PINS, 'pinned source changed'
    own = Path(__file__).relative_to(ROOT).as_posix()
    sources[own] = digest(ROOT/own)
    zero, one = K(0), K(1)
    s, phi = K(0, 0, 1), K(F(1, 2), F(1, 2))
    a, alpha, k, d = phi-1+s, 1+s, s-1, 2-s
    kval = a*alpha*phi
    h, ii, jj = alpha*kval, kval*(phi+s), kval*phi
    cc, uu, vv = s*alpha*phi, alpha*phi*phi, phi*(phi+s)
    assert s*s == 2 and phi*phi == phi+1
    # Avoid exponentiation: the arithmetic helper has deliberately few methods.
    assert a*(a-2)*(alpha*phi)*(alpha*phi) == 1
    assert ii+jj == a*(h-1) and ii-jj == s*kval
    assert cc == kval-alpha and kval*kval == a/(a-2)
    places = ((1, 1), (-1, 1), (1, -1), (-1, -1))
    assert [a.emb(*p).sign() for p in places] == [1, -1, -1, -1]

    def unit(i):
        return [K(int(i == j)) for j in range(5)]

    def add(*vectors):
        return [sum((v[i] for v in vectors), zero) for i in range(5)]

    def scale(c, v):
        return [c*x for x in v]

    def transpose(m):
        return [list(row) for row in zip(*m)]

    def mm(m, n):
        return [[sum((x*y for x, y in zip(row, col)), zero)
                 for col in zip(*n)] for row in m]

    def mv(m, v):
        return [sum((x*y for x, y in zip(row, v)), zero) for row in m]

    def identity():
        return [unit(i) for i in range(5)]

    def gram(columns):
        return mm(mm(transpose(columns), metric), columns)

    def integral(matrix):
        return all(x.integral() for row in matrix for x in row)

    def encode(x):
        aa, bb, c, dd = x.v
        return list(map(str, (aa-bb, 2*bb, c-dd, 2*dd)))

    def ev(v):
        return [encode(x) for x in v]

    def em(m):
        return [ev(row) for row in m]

    o, e1, e2, e3, e4 = [unit(i) for i in range(5)]
    metric = [[zero for _ in range(5)] for _ in range(5)]
    for i in range(5):
        metric[i][i] = -a if i == 0 else one
    tp = add(scale(h, o), scale(ii, e1), scale(jj, e2))
    tq = scale(kval, add(o, e1, e2))
    r = add(o, scale(phi, e1))
    c = add(scale(cc, o), scale(uu, e1), scale(vv, e2))
    beta = add(scale(-one, e1), e2)
    qe, qf = add(o, scale(a, e2)), add(e1, scale(-one, e2))
    ac = [add(scale(kval, tp), scale(kval, e3), scale(kval, e4)),
          r, c, add(tp, scale(a, e4)), add(e3, scale(-one, e4))]
    bc = [add(scale(h, tq), scale(ii, e3), scale(jj, e4)),
          add(tq, scale(phi, e3)),
          add(scale(cc, tq), scale(uu, e3), scale(vv, e4)), qe, qf]
    amat, bmat = transpose(ac), transpose(bc)
    g = [[zero for _ in range(5)] for _ in range(5)]
    g[0][0] = -a
    g[1][1], g[1][2], g[2][1], g[2][2] = d, -one, -one, K(2)
    g[3][3], g[3][4], g[4][3], g[4][4] = a*(a-1), -a, -a, K(2)
    assert gram(amat) == gram(bmat) == g
    assert integral(amat) and integral(bmat)
    # Explicit block inverse: no elimination on an indefinite matrix.
    gi = [[zero for _ in range(5)] for _ in range(5)]
    gi[0][0] = -one/a
    for i, row in enumerate(((K(2), one), (one, d)), 1):
        for j, x in enumerate(row, 1):
            gi[i][j] = x/(k*k)
    for i, row in enumerate(((K(2), a), (a, a*(a-1))), 3):
        for j, x in enumerate(row, 3):
            gi[i][j] = x/(a*(a-2))
    assert integral(gi) and mm(g, gi) == mm(gi, g) == identity()
    ai, bi = mm(mm(gi, transpose(amat)), metric), mm(mm(gi, transpose(bmat)), metric)
    assert integral(ai) and integral(bi)
    assert mm(ai, amat) == mm(amat, ai) == identity()
    assert mm(bi, bmat) == mm(bmat, bi) == identity()
    period, inverse = mm(bmat, ai), mm(amat, bi)
    assert integral(period) and integral(inverse)
    assert mm(period, inverse) == mm(inverse, period) == identity()
    assert mm(mm(transpose(period), metric), period) == metric
    assert mm(period, amat) == bmat
    # Full actual-germ status is the written input; verify every oriented image.
    simple_a = [ac[1], ac[2], scale(-one, ac[4]),
                add(scale(cc, ac[3]), scale(uu, ac[4]))]
    simple_b = [bc[1], bc[2], scale(-one, bc[4]),
                add(scale(cc, bc[3]), scale(uu, bc[4]))]
    assert [mv(period, v) for v in simple_a] == simple_b
    assert ac[0][0].sign() > 0 and bc[0][0].sign() > 0
    assert mv(period, ac[0]) == bc[0]
    expected = 196*phi+121+s*(138*phi+86)
    trace = sum((period[i][i] for i in range(5)), zero)
    assert trace == expected and (trace-5).sign() > 0
    assert alpha*kval*kval*kval-alpha*alpha*kval*kval+(3+s)*kval == trace
    # The two forbidden B3 orbit projections really are nonintegral.
    forbidden = [add(scale(-one, beta), scale(-one/s, c)),
                 add(scale(-one, c), scale(-one/s, beta))]
    assert all(not all(x.integral() for x in v) for v in forbidden)
    assert {name: digest(ROOT/name) for name in sources} == sources, 'source changed during check'
    return {
        'status': 'PASS_FIXED_QUARTIC_ORDER_EXCHANGE_IDENTITIES',
        'basis': ['1', 'phi', 'sqrt2', 'sqrt2*phi'],
        'source_sha256': sources, 'active_dimension': 5,
        'frames_A_B': [em(amat), em(bmat)], 'block_gram': em(g),
        'period': em(period), 'inverse': em(inverse),
        'active_trace': encode(trace), 'full_N9_trace': encode(trace+5),
        'oriented_simple_images': [list(map(ev, simple_a)), list(map(ev, simple_b))],
        'future_point_images': [ev(ac[0]), ev(bc[0])],
        'nonintegral_B3_projections': list(map(ev, forbidden)),
        'written_inputs': [
            'Arithmetic preflight and complete first and second actual edges',
            'Actual PD-face completion and both FULL mirror-arrangement identities',
            'Finite Coxeter classification including the nonparabolic B3 orbit',
            'Actual full germs identify one chamber; free action of W on chambers',
            'Borel--Harish-Chandra arithmetic lattice theorem',
            'Uniform n>=8 inventory threshold and fixed spectator extension'],
        'scope': 'Exact five-coordinate arithmetic only; no geometric enumeration or new dimension ceiling.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.write:
        assert not CERT.exists(), 'refusing to replace saved evidence'
    result = check()
    if args.write:
        with CERT.open('x') as stream:
            json.dump(result, stream, indent=2)
            stream.write('\n')
    else:
        with CERT.open() as stream:
            assert json.load(stream) == result, 'saved certificate differs'
    print(result['status'], 'generation' if args.write else 'fresh replay')


if __name__ == '__main__':
    main()
