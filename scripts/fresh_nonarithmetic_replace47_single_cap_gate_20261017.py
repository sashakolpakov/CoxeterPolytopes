#!/usr/bin/env python3
"""Continuous norm obstruction to closing the whole remaining B region with one cap."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_nonarithmetic_replace47_close12_20261017 as C
P,U,I,D=C.P,C.U,C.I,C.D
ROOT=C.ROOT;TAG='fresh_nonarithmetic_replace47_single_cap_gate_20261017'


def run():
    z=C.geometry();ns=z['normals'];pair=z['pair']
    c=U.rnd(pair(ns[7],ns[52]));d=U.rnd(pair(ns[30],ns[47]))
    gram=[[I(F(4,3)),-D.theta(5),I(-1),c],
          [-D.theta(5),I(2),d,-D.theta(4)],
          [I(-1),d,I(2),I()],
          [c,-D.theta(4),I(),I(2)]]
    inverse,pivots=P.C53.inverse(gram)
    first=[row[:2]for row in gram[:2]];first_pd=P.pd(first)
    fi,fp=P.C53.inverse(first)
    assert all(fi[i][i].l>2 for i in range(2)) and fi[0][1].l>0 and fi[1][0].l>0
    endpoint747=[[I(F(4,3)),I(-1)],[I(-1),I(2)]]
    endpoint3052=[[I(2),-D.theta(4)],[-D.theta(4),I(2)]]
    e1,_=P.C53.inverse(endpoint747);e2,_=P.C53.inverse(endpoint3052)
    assert e1[1][1].l<=F(4,5)<=e1[1][1].h
    assert e2[1][1].l<=1<=e2[1][1].h
    # Exact simple Gram inverses give a²<5/2 and b²<2; the closed
    # rectangle is a relaxation. Separate convexity reduces its maximum
    # to four corners without needing definiteness of the2x2 form.
    q=[[inverse[i][j]for j in [2,3]]for i in [2,3]]
    assert q[0][0].l>0 and q[1][1].l>0
    bounds=[U.sqrt(I(F(5,2))),D.theta(4)];corners=[]
    for i,j in [(0,0),(0,1),(1,0),(1,1)]:
        v=[bounds[0]if i else I(),bounds[1]if j else I()]
        value=U.rnd(U.idot(v,U.imv(q,v)))
        assert value.h<2
        corners.append(dict(upper_choices=[i,j],contacts=P.enc(v),forced_square=value.enc(),strict_deficit=U.rnd(2-value).enc()))
    assert max(F(c['forced_square'][1])for c in corners)<F(219,125)
    return dict(status='PASS_ALL_LABEL_SINGLE_CAP_OBSTRUCTION_FOR_REMAINING_B_COMPONENT',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in[Path(__file__),Path(C.__file__),ROOT/'results'/(C.TAG+'.json'),ROOT/'results/fresh_nonarithmetic_replace47_close12_boundary_20261017.json']},
        support=C.S,residual_frame=[7,30,47,52],frame_gram=[P.enc(r)for r in gram],frame_inverse_pivots=pivots,
        edge730_gram=[P.enc(r)for r in first],edge730_PD_pivots=first_pd,edge730_inverse=[P.enc(r)for r in fi],
        edge730_inverse_pivots=fp,endpoint747_inverse=[P.enc(r)for r in e1],endpoint3052_inverse=[P.enc(r)for r in e2],
        endpoint_contact_bounds_squared=['5/2','2'],forced_norm_form=[P.enc(r)for r in q],rectangle_corners=corners,
        maximum_forced_square_upper=str(max(F(c['forced_square'][1])for c in corners)),
        simple_rational_norm_upper='219/125',
        proof='A norm2 cap rightS6 terminating allfour open edges7/30,7/47,30/52,47/52 must be finite to their sides. Gram(7,30) inverse has positive entries and both diagonals>2. Since any nonzero Coxeter contact has magnitude>=1, positivity of that endpoint forces contacts7=30=0. The other endpoint PD conditions then imply a²<5/2,b²<2 for a=-cap.47,b=-cap.52. The invertible4frame uniquely determines the cap and its norm q(a,b). Both diagonal coefficients are positive, hence the maximum on the relaxed closed rectangle lies at a corner. Allfour exact interval corner bounds are strictly below2, contradiction.',
        automatic_right_support='The finite endpoint on S6+7+30 has diagramH4+H4. Neither H4 admits a nonorthogonal one-root extension in a finite Coxeter diagram. Thus any Coxeter wall making this endpoint finite is automatically right to all ofS6,7,30; the rightS6 hypothesis loses no possible single-wall closure.',
        scope='Excludes a single actual-S6-facet Coxeter cap making allfour remaining B-boundary endpoints finite, for every label and coefficient field. RightS6 follows automatically from theH4² endpoint. No common-center assumption is used in the norm obstruction. It does not exclude multiplecaps, replacements or finite-volume completion by another mechanism.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    path=ROOT/'results'/(TAG+'.json')
    if a.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['maximum_forced_square_upper'],flush=True)

if __name__=='__main__':main()
