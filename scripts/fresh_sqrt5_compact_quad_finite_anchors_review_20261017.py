#!/usr/bin/env python3
"""Independent integer-lattice reconstruction of the compact quad anchors."""
import argparse
import json
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_compact_quad_finite_anchors_review_20261017'
PRIMARY = 'fresh_sqrt5_compact_quad_finite_anchors_20261017'


def sign(p, q):
    """Sign of p+q*sqrt(5), using integer comparisons only."""
    if q == 0:
        return (p > 0)-(p < 0)
    if p == 0:
        return (q > 0)-(q < 0)
    if p > 0 and q > 0:
        return 1
    if p < 0 and q < 0:
        return -1
    assert p*p != 5*q*q
    if p > 0:
        return 1 if p*p > 5*q*q else -1
    return 1 if 5*q*q > p*p else -1


def run():
    primary = json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    for p, h in primary['sources'].items():
        assert sha256((ROOT/p).read_bytes()).hexdigest() == h
    direct, squares, rational_direct, rational_squares = [], [], [], []
    for q in range(-3, 1):
        for p in range(-15, 13):
            if (p-q) % 2:
                continue
            if not (sign(p+12, q) > 0 and sign(p, q) <= 0
                    and sign(p+4, -q) > 0 and sign(4-p, q) > 0):
                continue
            if sign(p+4, q) >= 0 and (p, q) not in ((0, 0), (-2, 0), (-1, -1)):
                continue
            (rational_direct if q == 0 else direct).append([p, q])
    finite_squares = {(0, 0), (2, 0), (4, 0), (3, 1), (6, 0), (5, 1)}
    for q in range(-1, 18):
        for p in range(-69, 94):
            if (p-q) % 2:
                continue
            if not (sign(p, q) >= 0 and sign(72-p, -q) > 0
                    and sign(p, -q) >= 0 and sign(8-p, q) > 0):
                continue
            if sign(p-8, q) <= 0 and (p, q) not in finite_squares:
                continue
            (rational_squares if q == 0 else squares).append([p, q])
    def convert(rows):
        return sorted([2*a+b, b] for a, b in (r['phi_coefficients'] for r in rows))
    for ours, field in ((direct, 'direct_product_anchors'), (squares, 'squared_product_anchors'),
                        (rational_direct, 'rational_direct_excluded'),
                        (rational_squares, 'rational_squares_excluded')):
        assert sorted(ours) == convert(primary[field]), field
    assert list(map(len, (direct, squares, rational_direct, rational_squares))) == [8, 59, 2, 4]
    paths = [f'scripts/{TAG}.py', f'scripts/{PRIMARY}.py', f'results/{PRIMARY}.json',
             f'results/{PRIMARY}.md']
    return dict(status='PASS_INDEPENDENT_INTEGER_ANCHOR_CENSUS_AND_GEOMETRIC_REVIEW',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                basis='(p+q*sqrt5)/2, p congruent q modulo2',
                direct_anchors=sorted(direct), squared_anchors=sorted(squares),
                rational_direct=sorted(rational_direct), rational_squares=sorted(rational_squares),
                lattice_box_proof='Subtracting the two embeddings bounds q to -3..0 for a direct product and -1..17 for its square. The primary safe a-b boxes, transformed by p=2a+b,q=b, lie in p=-15..12 and -69..93 respectively. All lattice points of these rectangles with the required parity are tested, using integer-square sign comparisons only.',
                physical_review='Adjacent residual2-by2 blocks are positive definite and the white inverse is nonnegative. Thus diagonal<=2, each adjacent product is strictly greater than-2 and at most0, and each opposite residual product is no larger than its full product. If both full products<=-6, sums of opposite basis vectors have a negative definite2-by2 Gram: both diagonal<=-8 and offdiagonal magnitude<8. This contradicts index1.',
                identity_review='The v-squared coefficient is minus the opposite residual determinant, by the two determinant entries containing v. An identically-zero specialization forces a singular opposite block with positive diagonals and nonpositive offdiagonal. Its positive kernel lifts through -W_inverse*B to a nonnegative full combination. The strict center makes this null vector nonzero and future; acuteness puts it in every retained halfspace, contradicting compactness.',
                integrality_review='Squared normalized products are arithmetic cycle invariants in O_K. If a product itself lies in K, its square integrality makes it integral by integral closure. Conjugate equality to+/-2 would force physical square4, impossible in the compact acute setting. Rational exclusions require the stated ADE fixed-contact edge hypothesis.',
                scope='Independent reconstruction of every finite anchor with no producer or project arithmetic imports, plus a written mathematical review of the physical and positive-null arguments. No diagram-by-diagram quadratic solving, dimension theorem, formal proof-assistant verification or external review.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    data = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(data, indent=2)+'\n')
    print(data['status'], len(data['direct_anchors']), len(data['squared_anchors']))
