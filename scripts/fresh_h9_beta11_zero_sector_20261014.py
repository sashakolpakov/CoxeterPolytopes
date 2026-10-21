#!/usr/bin/env python3
"""Complete height-zero reflection sector of the beta11 H9 lattice."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_beta11_zero_sector_20261014'
BASIS='results/fresh_h9_h8_cap_height_20261014_host_basis.json'
POOL='results/fresh_h9_h8_cap_height_20261014_norm_pool_review.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def integer_shell(G,target):
    n=len(G);S=[[F(x)for x in r]for r in G];D=[];U=[]
    for i in range(n):
        d=S[i][i];assert d>0;D.append(d)
        U.append([F(0)if j<=i else S[i][j]/d for j in range(n)])
        for j in range(i+1,n):
            for k in range(i+1,n):S[j][k]-=S[j][i]*S[i][k]/d
    roots=[];visits=0;y=[0]*n
    def walk(i,budget):
        nonlocal visits
        visits+=1
        if i<0:
            if not budget:roots.append(y[:])
            return
        shift=sum(U[i][j]*y[j]for j in range(i+1,n));rho=budget/D[i]
        assert rho>=0
        radius=isqrt(rho.numerator//rho.denominator)+2
        center=(-shift).numerator//(-shift).denominator
        for z in range(center-radius,center+radius+2):
            used=D[i]*(z+shift)**2
            if used<=budget:y[i]=z;walk(i-1,budget-used)
    walk(n-1,F(target))
    for y in roots:assert sum(F(G[i][j])*y[i]*y[j]for i in range(n)for j in range(n))==target
    return sorted(roots),visits
def run():
    d,p=[json.loads((ROOT/f).read_text())for f in(BASIS,POOL)]
    for obj in(d,p):
        for f,h in obj['sources'].items():assert sha(f)==h,f
    G=A.dm(d['gram']);phi=A.PHI;beta=3+phi;delta=1+4*phi
    H=[A.E(2)]+[A.Z]*8+[A.ONE]
    assert A.mv(G,H)==[-delta]+[A.Z]*9
    assert A.pair(H,H,G)==-2*delta and (-2*delta).sign()<0
    M=[[G[i][j]for j in range(1,10)]for i in range(1,10)]
    assert all(not x.b and not x.a.b and x.a.a.denominator==1 for row in M for x in row)
    Mi=[[int(x.a.a)for x in row]for row in M]
    assert all(Mi[i][i]==2 for i in range(9))
    roots,visits=integer_shell(Mi,2)
    assert len(roots)==242 and sum(bool(r[-1])for r in roots)==2
    assert all((not r[-1])or not any(r[:-1])for r in roots)
    # All full zero-height mirrors are supported on the nine listed walls.
    simple=[]
    for i in range(1,9):
        r=[A.Z]*10;r[i]=A.E(-1);simple.append(r)
    simple.append([A.Z]*9+[A.ONE])
    for r in roots:
        coeff=[-x for x in r[:-1]]+[r[-1]]
        assert all(x>=0 for x in coeff)or all(x<=0 for x in coeff)
    assert all(A.pair(H,r,G)==0 and A.pair(r,r,G)==2 for r in simple)
    # For q2beta, x0 forces every remaining coordinate into betaO.
    # Dividing the norm by beta² gives2/beta, not an even integral norm.
    reciprocal=1/beta
    assert reciprocal.a.a-reciprocal.a.b==F(4,11) and 2*reciprocal.a.b==F(-1,11)
    assert not((reciprocal.a.a-reciprocal.a.b).denominator==1 and (2*reciprocal.a.b).denominator==1)
    paths=[BASIS,POOL,f'scripts/{TAG}.py','scripts/fresh_h9_mixed_d4_coupling_20261014_review.py']
    return A.enc(dict(status='PASS_COMPLETE_BETA11_ZERO_HEIGHT_SECTOR',sources={f:sha(f)for f in paths},
        basepoint=H,basepoint_norm=-2*delta,height_formula='(r,H)^2/(-H²*q)=delta*x²/(2q)',
        zero_rational_gram=Mi,complete_integer_norm2_roots=roots,integer_shell_visits=visits,
        norm2_zero_roots=242,norm2beta_zero_roots=0,simple_zero_walls=simple,sector='E8+A1',
        norm2_completeness='Write the rational evenPD zero lattice vector asu+phi*v. Norm2 givesu²+v²=2 and2(u,v)+v²=0. Since the rational lattice isevenPD, eitheru²2,v0 oru0,v²2; the latter violates thephi coefficient. Hencev0 andu isone ofthe242exhaustively enumeratedordinaryroots.',
        norm2beta_completeness='The complete reflection residuefilter atx0 forces ally_i andt divisiblebybeta. Dividingr bybeta produces anintegral vectorofnorm2/beta, contradicting evenintegrality because1/beta=(4−phi)/11 isnotintegral.',
        scope='Complete heightzero reflection arrangement for this lattice, with exactsimple sector and Vinbergbasepoint. Completeness inherits the independently reviewed two-norm alphabet and residuefilter. No positiveheight root or finitevolume claim.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['norm2_zero_roots'],d['integer_shell_visits'])
