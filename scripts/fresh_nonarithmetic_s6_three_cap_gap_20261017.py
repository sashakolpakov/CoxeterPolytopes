#!/usr/bin/env python3
"""Exact rejection of the lone numerical near-hit in the three-cap scout."""
import argparse
import json
from hashlib import sha256
from pathlib import Path

import fresh_nonarithmetic_escape62_20261017_cap as C

P,U,I,D=C.P,C.U,C.I,C.D
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nonarithmetic_s6_three_cap_gap_20261017'
S=[0,1,2,3,5,6]


def run():
    z=C.geometry();normals=z['normals'];pair=z['pair']
    zero,two=I(0),I(2)
    b=U.rnd(pair(normals[8],normals[30]))
    x=U.rnd(pair(normals[44],normals[52]))
    theta4,theta15=D.theta(4),D.theta(15)
    # All four old side roots are perpendicular to S6. Literal right/finite
    # entries avoid expanding identities through ill-conditioned coordinates.
    gram3=[[two,zero,zero],[zero,two,x],[zero,x,two]]
    inverse,pivots=P.C53.inverse(gram3)
    contacts=[zero,-theta15,-theta4]
    old30=[b,zero,-theta4]
    square=U.rnd(2-U.idot(old30,U.imv(inverse,old30)))
    residual=U.rnd(2-U.idot(contacts,U.imv(inverse,contacts)))
    assert square.l>0 and residual.l>0
    particular30=U.rnd(U.idot(contacts,U.imv(inverse,old30)))
    c30=U.rnd(particular30-U.sqrt(square*residual))
    assert c30.h<-2
    gram4=[[two,zero,-theta4,c30],
           [zero,two,x,-theta15],
           [-theta4,x,two,-theta4],
           [c30,-theta15,-theta4,two]]
    dual,dual_pivots=P.C53.inverse(gram4)
    assert dual[1][1].l>0 and dual[2][2].l>0
    # A and B are the negative normalized duals to52 and44 respectively.
    value=U.rnd(2*dual[1][2]/U.sqrt(dual[1][1]*dual[2][2]))
    lower=U.rnd(value+D.theta(28));upper=U.rnd(value+D.theta(27))
    assert lower.l>0 and upper.h<0, ('gap unresolved',value.enc(),lower.enc(),upper.enc())
    paths=[f'scripts/{TAG}.py',f'scripts/{C.TAG}.py',f'results/{C.TAG}.json']
    return dict(status='PASS_EXACT_THREE_CAP_NEAR_HIT_IS_A_27_28_ANGLE_GAP',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                retained_walls=63,C_frame=S+[8,44,52],C_labels=[2]*7+[15,4],
                C_inverse_pivots=pivots,C_kernel_square=square.enc(),C_residual=residual.enc(),
                old_dotted_products=P.enc([b,x]),C_old30_product=c30.enc(),
                four_root_frame=[30,44,52,'C'],four_root_gram=[P.enc(r)for r in gram4],
                dual_gram=[P.enc(r)for r in dual],dual_inverse_pivots=dual_pivots,mutual_product=value.enc(),
                consecutive_labels=[27,28],above_negative_theta28=lower.enc(),below_negative_theta27=upper.enc(),
                exact_recipe='C has the stated nine finite contacts and the spacelike-kernel branch directed negatively along the projection of old30. A and B are the negative unit polars obtained by projecting opposite old52/44 off S6+[30,44,C] and S6+[30,52,C].',
                scope='Rejection of one prescribed three-cap configuration. The mutual product lies strictly between consecutive Coxeter values; the numerical near-hit is not an exact wall cluster. No all-label or all-topology exclusion.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)
