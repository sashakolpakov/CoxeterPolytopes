#!/usr/bin/env python3
"""All-real, all-label obstruction to cutting E while right to T+31."""
import argparse
import json
from pathlib import Path

import fresh_nonarithmetic_two_caps_20261016_review as V

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_right31_all_real_escape_20261016'
R, K, Z, O = V.R, V.K, V.Z, V.O


def run():
    base = V.read(V.BASE)
    primary = V.read(V.PRIMARY)
    review_path = 'results/fresh_nonarithmetic_two_caps_20261016_review.json'
    review = V.read(review_path)
    assert review['status'] == 'PASS_INDEPENDENT_BOUNDED_TWO_CAP_BRANCH_REVIEW'
    assert review['sources'][V.PRIMARY] == V.digest(V.PRIMARY)
    Q, N, E = R.m4(base['form']), R.m4(base['normals']), R.v4(base['strict_escape']['point'])
    pair = lambda x, y: V.dot(x, V.mv(Q, y))
    support = [0, 1, 2, 5, 6, 10, 30, 31]
    frame = next(f for f in primary['linear_frames'] if f['anchor'] == 3)
    p, k = list(map(V.dec, frame['particular'])), list(map(V.dec, frame['kernel']))
    rows = [V.mv(Q, N[i]) for i in support+[3]]
    assert V.rank(rows) == 9
    assert V.mv(rows, p) == [Z]*8+[-O]
    assert V.mv(rows, k) == [Z]*9
    assert len(p) == len(k) == 10
    pp, kk, pk = pair(p, p), pair(k, k), pair(p, k)
    assert R.sg(kk) > 0 and pk == Z
    assert V.rank(rows+[V.mv(Q, k)]) == 10
    pE, kE, k7 = pair(p, E), pair(k, E), pair(k, N[7])
    assert R.sg(kE) < 0 and R.sg(k7) < 0
    tau = -pE/kE
    assert R.sg(tau) < 0
    boundary = [a+tau*b for a, b in zip(p, k)]
    assert pair(boundary, E) == Z
    qboundary = pp+tau*tau*kk
    assert pair(boundary, boundary) == qboundary
    assert R.sg(qboundary-5*O) > 0
    assert pair(N[3], N[3]) == 2*O
    # u=a(p+t k), a=-u.N3>0. Escape cutting forces t<tau<0.
    # Therefore q(t)>q(tau)>5 and the invariant at N3 is2/q(t)<2/5.
    assert R.sg(O-2*O/5) > 0
    paths = [f'scripts/{TAG}.py', V.BASE, V.PRIMARY, review_path,
             'scripts/fresh_nonarithmetic_two_caps_20261016_review.py',
             'scripts/fresh_h9_joint_corner_20261015_bridge_review.py',
             'scripts/fresh_h9_h8_all_high_20261014_review.py',
             'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
    return dict(
        status='PASS_ALL_REAL_RIGHT31_ESCAPE_CUTTER_OBSTRUCTION',
        sources={p: V.digest(p) for p in paths},
        support=support, old_angle_blocker=3, kernel_blocker=7,
        ambient_dimension=9, ambient_vector_dimension=10,
        linear_rank=9, particular=R.enc(p), kernel=R.enc(k),
        particular_norm=R.enc(pp), kernel_norm=R.enc(kk), cross_pairing=R.enc(pk),
        escape_point=R.enc(E), particular_escape_pairing=R.enc(pE),
        kernel_escape_pairing=R.enc(kE), kernel_wall7_pairing=R.enc(k7),
        parameter_upper_bound=R.enc(tau), boundary_vector=R.enc(boundary),
        boundary_norm=R.enc(qboundary), strict_rational_norm_lower_bound='5',
        strict_Cartan_cycle_upper_bound='2/5',
        case_split=[
            'u.N3>0: obtuse to a retained wall.',
            'u.N3=0: u=b*k; u.E>0 and k.E<0 force b<0, hence u.N7>0.',
            'u.N3<0: a=-u.N3>0 and u=a*(p+t*k). Escape cutting implies t<tau<0. The square q(t)=pp+kk*t*t exceeds q(tau)>5. Since N3 has square2, the Cartan cycle is2/q(t), strictly between0 and2/5, below the first nonzero Coxeter value1.'
        ],
        conclusion='No real spacelike normal orthogonal to the literal eight support walls can both cut the saved escape and satisfy Coxeter compatibility with retained walls3and7.',
        scope='All real coefficients and all integer Coxeter labels, without a norm, field or height cutoff. The support and escape are fixed. This does not exclude non-cutting caps, caps with different support contacts, later cutters in a chain, or arbitrary finite-volume completions.'
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    output = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(output.read_text())
    else:
        output.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], 'all real normals; Cartan cycle <2/5')
