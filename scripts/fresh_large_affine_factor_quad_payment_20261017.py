#!/usr/bin/env python3
"""Exact new scalar checks for large affine factors with arbitrary spectators.

The report supplies the all-rank extension classification. No ambient bound.
"""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_large_affine_factor_quad_payment_20261017'
INPUTS=[
 'results/fresh_simple_cusp_affine_tree_completion_20261017.md',
 'results/fresh_simple_cusp_affine_tree_short_ends_20261017.md',
 'results/fresh_simple_cusp_affine_a_quadrilateral_20261017.md',
 'results/fresh_cusp_factor_payment_framework_20261017.json',
 'results/fresh_cusp_spectator_role_payments_20261017.json',
 'results/fresh_simple_finite_quad_vinberg_transfer_20261017.md']

def tree_arms(lengths):
    edges=[];leaves=[];j=1
    for length in lengths:
        a=0
        for _ in range(length): edges.append((a,j));a=j;j+=1
        leaves.append(a)
    G=2*s.eye(j)
    for a,b in edges:G[a,b]=G[b,a]=-1
    return G,leaves

def exact_positive(x):
    x=s.simplify(x)
    assert x.is_positive is True,x
    return str(x)

def run():
    rt=s.sqrt(2);Q=s.Rational
    gnext=4/(1+rt)
    assert s.simplify(gnext-4*(rt-1))==0
    bounds={
      'two_main_marks_P4':1-2*gnext/4,
      'double_terminal_and_main_P6':1-(Q(5,2)+gnext)/6,
      'fork_leaf_and_main_R6':1-(4+gnext)/6,
      'two_terminals_R17':1-Q(8,17),
      'endpoint4_below_next':gnext-2*rt/(1+rt),
      'D2_double_join_below_old14over9':Q(14,9)-2/(1+rt/2)}
    scalar={k:exact_positive(v)for k,v in bounds.items()}
    E8,leaves=tree_arms([1,2,4]);d=leaves[1]
    assert E8.inv()[d,d]==4

    # Literal D3 outer block in the main-path P4/E8 exception.
    # Order L, center, Lprime; d meets Lprime by a simple bond.
    G=s.Matrix([[2,-1,0],[-1,2,-1],[0,-1,2]]);J=G.inv()
    old_attach=s.Matrix([0,1,0]);wd=s.Matrix([0,0,1])
    u,v=s.symbols('u v',nonnegative=True)
    A=Q(1,2);D=Q(5,4);alpha=Q(1,4)
    rows=[]
    for name,w,umax in [
       ('outer_leaf3_other',s.Matrix([1,0,0]),1/(rt+Q(1,2))),
       ('outer_leaf3_same',s.Matrix([0,0,1]),1/(rt+Q(1,2))),
       ('outer_center3',s.Matrix([0,1,0]),1/(1+rt)),
       ('outer_leaf4_other',s.Matrix([rt,0,0]),1/(1+rt/2))]:
        C=2-(w.T*J*w)[0];E=(w.T*J*wd)[0]
        B=0
        F=s.expand((1-B*u-A*v)**2-alpha*(C*u*u+D*v*v+2*E*u*v))
        # C,D,E positive: F strictly decreases on the positive rectangle.
        assert C>0 and D>0 and E>0
        minimum=s.simplify(F.subs({u:umax,v:Q(2,3)}))
        rows.append(dict(profile=name,C=str(C),E=str(E),height_offset=str((old_attach.T*J*w)[0]),
            umax=str(s.simplify(umax)),polynomial=str(F),minimum=exact_positive(minimum)))
    # New inner affine profiles. Their support is separated from d's old D3 block.
    for name,l,umax in [('inner_endpoint4',rt,1/(1+rt)),('inner_second3',2,1/(1+rt))]:
        B=l/4;C=l*l/4;E=0
        F=s.expand((1-B*u-A*v)**2-alpha*(C*u*u+D*v*v))
        vals=[s.simplify(F.subs({u:uu,v:vv}))for uu,vv in product((0,umax),(0,Q(2,3)))]
        for z in vals:exact_positive(z)
        rows.append(dict(profile=name,B=str(B),C=str(C),E='0',umax=str(s.simplify(umax)),
            polynomial=str(F),rectangle_corner_values=list(map(str,vals))))

    # The two direct forbidden small stars added in the outer-only audit.
    affine=[]
    for name,edges in [
      ('D4_center_leaf',[(0,1,1),(0,2,1),(0,3,1),(0,4,1)]),
      ('B4_doubled_leaf_arm',[(0,1,1),(0,2,1),(0,3,1),(3,4,rt)])]:
        H=2*s.eye(5)
        for a,b,t in edges:H[a,b]=H[b,a]=-t
        ns=H.nullspace();assert len(ns)==1
        k=ns[0]
        if k[0]<0:k=-k
        assert all(x>0 for x in k)
        for j in range(5):
            M=H.minor_submatrix(j,j)
            assert all(M[:z,:z].det()>0 for z in range(1,5))
        affine.append(dict(name=name,gram=[[str(x)for x in r]for r in H.tolist()],kernel=list(map(str,k))))

    # All-rank available-simple-arm identities. Near-end exceptions have <=5 nodes.
    m,r,P=s.symbols('m r P',integer=True)
    arm_cases=[]
    for shift in (0,1,2):
        left=r-1+shift;right=m-3-left
        assert s.expand(left+right-(m-3))==0
        arm_cases.append(dict(inner_shift=shift,left=str(left),right=str(right),
            sum='m-3',near_end_exceptions=([1,2]if shift==0 else [1]if shift==1 else [])))
    paths=[Path(__file__)]+[ROOT/p for p in INPUTS]
    return dict(status='PASS_LARGE_AFFINE_FACTOR_QUADRILATERAL_SCALAR_CHECKS',
       sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
       generic_margins=scalar,E8_allowed_d_inverse_diagonal='4',
       new_E8_profile_checks=rows,small_affine_witnesses=affine,
       available_simple_arm_identities=arm_cases,
       scope='Quadrilaterals with a classical affine corner factor of order>=18, assuming no connected affine subdiagram of order<=9. Arbitrary acute PSD spectators are allowed. The all-finite-corner transfer and any ambient dimension conclusion remain separate.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    data=json.loads(json.dumps(run()));out=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],len(data['new_E8_profile_checks']),'new exceptional profiles',flush=True)
