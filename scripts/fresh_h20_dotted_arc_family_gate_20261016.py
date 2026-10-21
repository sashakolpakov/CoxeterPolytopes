#!/usr/bin/env python3
"""Exact all-label escape obstruction for symmetric single-support arc caps.

Uses the reviewed nonabelian chain constant and standalone rational scalar
inequalities. No cap-bank enumeration or floating point input is used.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_dotted_arc_family_gate_20261016'
BASE='results/fresh_h20_s61_multicap_20261014_exact.json'
CAP='results/fresh_h20_dotted_arc_cap_20261016.json'
ABEL='results/fresh_h20_chain_abelian_gate_20261016.json'
ABEL_REVIEW='results/fresh_h20_chain_abelian_gate_20261016_review.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def run():
    base,cap,abel,rev=[json.loads((ROOT/p).read_text())for p in [BASE,CAP,ABEL,ABEL_REVIEW]]
    for doc in [base,cap,abel,rev]:
        for p,h in doc['sources'].items():assert sha(p)==h,p
    assert rev['sources'][ABEL]==sha(ABEL)
    assert cap['status']=='PASS_EXACT_134_WALL_DOTTED_ARC_ESCAPE_CUTTER'
    S=base['support'];J=[[F(x)for x in row]for row in base['support_inverse']]
    assert S==cap['strict_escape']['support'] and len(S)==18
    assert all(x>=0 for row in J for x in row)
    lo,hi=map(F,base['parameter']['isolating_interval']);assert lo*lo>280
    Cmin=lo*lo/(lo*lo+2);assert 0<Cmin<1
    support_ceiling=(1+1/Cmin)/2
    assert support_ceiling<F(281,280)
    allowed=[i for i in range(18)if J[i][i]<support_ceiling]
    assert [S[i]for i in allowed]==[1,10,12,15,17,31,61,87]
    assert all(J[i][i]>=support_ceiling for i in range(18)if i not in allowed)
    ez=list(map(F,cap['strict_escape']['frame_coefficients']))
    eb=list(map(F,cap['strict_escape']['dual_coefficients']))
    assert ez==[F(1),F(-1,100),F(0)]
    assert eb==[F(-1)if i==5 else F(-1,100000)for i in S]
    z=[-sum(J[i][k]*eb[k]for k in range(18))for i in range(18)]
    K=ez[0]*ez[0]-ez[1]*ez[1]-ez[2]*ez[2];assert K==F(9999,10000)
    rows=[]
    for j in allowed:
        u=max(4*J[j][j]-2,1/(hi*hi))
        leading=J[j][j]*K-z[j]*z[j]
        margin=u*leading-2*z[j]*z[j]
        assert u>0 and leading>0 and margin>0
        rows.append(dict(wall=S[j],inverse_diagonal=str(J[j][j]),escape_dual_pair_magnitude=str(z[j]),
                         residual_timelike_norm_lower=str(u),CS_linear_coefficient=str(leading),
                         CS_strict_margin=str(margin)))
    paths=[BASE,CAP,ABEL,ABEL_REVIEW,f'scripts/{TAG}.py']
    return dict(status='PASS_ALL_LABEL_SYMMETRIC_SINGLE_SUPPORT_ARC_COMMON_ESCAPE',
        sources={p:sha(p)for p in paths},support=S,centers=list(range(108,132)),
        side_contacts='finite Coxeter or parallel, common magnitude B in {0} union {2cos(pi/m):m>=3} union {2}',
        middle_contact='finite Coxeter or parallel, magnitude A in the same set',
        norm_identity='y^2=A^2/2-B^2/(2*C); Jjj*c^2=2-y^2',
        field_reduction='If c is finite or parallel then C=B^2/(A^2+2*Jjj*c^2-4) lies in Q^ab, contradiction; B=0 is separately impossible by the all-nondotted-chain argument.',
        monotonic_reduction='B<=A makes all26 chain pairings >=-A>=-2, hence all nondotted, impossible. Thus B>A and B>=1.',
        residual_norm_lower='u=-y^2=(B^2-A^2)/2+B^2/a^2 >1/a^2; c>2 also gives u>4Jjj-2.',
        inverse_diagonal_ceiling=str(support_ceiling),allowed_support_walls=[S[i]for i in allowed],
        remaining_exact_margins=rows,escape_dependency=CAP,
        conclusion='Every coherent cap in this all-label symmetric three-contact, one-support-contact family has strictly negative pairing with the same saved simple escape. Any finite mutually compatible collection from this family therefore also misses it.',
        scope='All real normalized norm2 caps retaining all133 old walls, right to17 of18 S61 roots, with equal finite/parallel angles at a consecutive chain triple\'s outer walls and finite/parallel middle angle. All24 centers, all18 support releases, all integer labels and arbitrary real coefficients are included. Dotted contacts elsewhere are unrestricted. Unequal outer contacts, two or more nonright support contacts, and changed retained walls remain outside.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['remaining_exact_margins']),flush=True)
