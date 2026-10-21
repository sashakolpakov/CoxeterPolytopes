#!/usr/bin/env python3
"""Independent 384-case reconstruction after changing the 4--5 block.

Only prior independent Fraction polynomial arithmetic is imported. Both
block inverses and all third-wall projections are rebuilt from their Gram.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h9_t8_middle_20261014_review_algebra as K

P = K.P
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_t8_all_blocks_20261014_review_algebra'
Z, O = K.Z, K.O
add, neg, mul, scale, sign, const = K.add, K.neg, K.mul, K.scale, K.sign, K.const

def run():
    rows = []
    count = 0
    minimum = None
    cos5 = (F(1,4), F(0), F(1,4), F(0))
    for second in (4, 5):
        J = [[Z for _ in range(4)] for _ in range(4)]
        for offset, label in ((0, 5), (2, second)):
            diag = scale(O, 2) if label == 4 else (F(2), F(0), F(2,5), F(0))
            off = P.SQ2 if label == 4 else (F(1), F(0), F(3,5), F(0))
            cosine = P.COS4 if label == 4 else cos5
            for i in range(2):
                J[offset+i][offset+i] = diag
                J[offset+i][offset+(i^1)] = off
            for i in range(2):
                for j in range(2):
                    entry = add(J[offset+i][offset+j], neg(mul(cosine, J[offset+(i^1)][offset+j])))
                    assert entry == (O if i == j else Z)
        def foreign(w): return [neg(P.HALF) if j == w else Z for j in (4,5,7,8)]
        def dot(x, y):
            value = Z
            for i in range(4):
                for j in range(4): value = add(value, mul(x[i], mul(J[i][j], y[j])))
            return value
        fz = [neg(P.HALF), neg(P.HALF), Z, Z]
        qvalue = add((F(25,19)-1, F(0), F(0), F(0)), dot(fz,fz))
        assert qvalue == (F(25,19)+F(1,2), F(0), F(1,2), F(0))
        for block, ia, ib, ma, mb, own, lower, delta in P.cases():
            if ia not in (4,5) and ib not in (4,5): continue
            if second == 5 and ia not in (7,8) and ib not in (7,8): continue
            e,a,b,d = [own[n] for n in ('e','a','b','d')]
            c = P.var('c')
            B = P.minus(const(1), P.sq(c))
            fx, fy = foreign(ia), foreign(ib)
            costs = [dot(fx,fx), dot(fy,fy), dot(fx,fy)]
            third = [add(Z if m == 2 else neg(P.HALF), neg(dot(f,fz)))
                     for m,f in ((ma,fx),(mb,fy))]
            u,v,p = map(const,costs)
            ac,bc = map(const,third)
            Q = const(qvalue)
            h = P.plus(delta,p)
            pc = P.plus(P.sq(e), P.sc(P.times(c,e,a),2), P.sq(a))
            cc = P.plus(P.sq(b), P.sc(P.times(c,b,d),2), P.sq(d))
            cross = P.plus(P.times(e,b), P.times(a,d), P.times(c,P.plus(P.times(e,d),P.times(a,b))))
            qa = P.minus(P.times(P.minus(const(1),u),B), pc)
            qb = P.minus(P.times(P.minus(const(1),v),B), cc)
            ab = P.sc(P.plus(P.times(h,B),cross), -1)
            M = [[qa,ab,P.times(ac,B)], [ab,qb,P.times(bc,B)],
                 [P.times(ac,B),P.times(bc,B),P.sc(P.times(Q,B),-1)]]
            aa = P.plus(P.times(P.minus(P.times(P.minus(const(1),u),P.minus(const(1),v)),P.sq(h)),B),
                        P.sc(P.times(P.minus(const(1),u),cc),-1), P.sc(P.times(P.minus(const(1),v),pc),-1),
                        P.sc(P.times(h,cross),-2), P.sq(P.minus(P.times(e,d),P.times(a,b))))
            nn = P.plus(P.sc(P.times(Q,aa),-1),
                        P.sc(P.times(B,P.plus(P.times(P.minus(const(1),u),P.sq(bc)),
                                             P.times(P.minus(const(1),v),P.sq(ac)),P.sc(P.times(h,ac,bc),2))),-1),
                        P.times(cc,P.sq(ac)), P.times(pc,P.sq(bc)), P.sc(P.times(cross,ac,bc),-2))
            assert P.determinant(M) == P.times(P.sq(B),nn)
            assert P.determinant([r[:2] for r in M[:2]]) == P.times(B,aa)
            cert = []
            for name, pol in (('N',nn),('minus_A',P.sc(aa,-1))):
                degree, terms, coefficients = P.bernstein(pol,lower)
                count += len(coefficients)
                for _,value in coefficients:
                    if minimum is None or sign(add(value,neg(minimum))) < 0: minimum = value
                cert.append(dict(polynomial=name,degree=degree,
                                 unit_power_terms=[dict(index=k,value=v) for k,v in sorted(terms.items())],
                                 positive_Bernstein_coefficients=[dict(index=k,value=v) for k,v in coefficients]))
            rows.append(dict(block=block,foreign=[ia,ib],third_labels=[ma,mb],
                             second_middle_label=second,lower_bounds=lower,
                             foreign_costs=costs,third_pairings=third,Q=qvalue,certificates=cert))
    assert len(rows) == 384 and count == 5472
    assert minimum == (F(261,2000),F(0),F(-31,1000),F(0))
    paths = [f'scripts/{TAG}.py', 'scripts/fresh_h9_t8_middle_20261014_review_algebra.py',
             'scripts/fresh_h9_t8_block_polynomials_20261014_review.py']
    return P.encode(dict(status='PASS_INDEPENDENT_384_ALL_BLOCKS_QUARTIC_CASES',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        field_basis=[1,2,5,10],case_count=len(rows),positive_coefficient_count=count,
        minimum_positive_coefficient=minimum,cases=rows,
        scope='Independent polynomial algebra for the 384 changed cases at internal4-5 label5. Literal domains, inherited cases, and unrestricted-label extension are a separate proof obligation.'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify: assert result == json.loads(path.read_text())
    else: path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['case_count'],result['positive_coefficient_count'],flush=True)
