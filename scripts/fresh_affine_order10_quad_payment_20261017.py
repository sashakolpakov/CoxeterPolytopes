#!/usr/bin/env python3
"""New boundary-order ideal-corner calculations, not a label/rank cutoff proof."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sympy as s
import fresh_private_affine_tree_exclusion_20261017 as P

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_affine_order10_quad_payment_20261017'
INPUTS=[
 'results/fresh_large_affine_factor_quad_payment_20261017.md',
 'results/fresh_large_affine_factor_quad_payment_20261017.json',
 'results/fresh_simple_cusp_affine_tree_completion_20261017.md',
 'results/fresh_simple_cusp_affine_tree_short_ends_20261017.md',
 'results/fresh_simple_cusp_affine_a_quadrilateral_20261017.md',
 'results/fresh_cusp_spectator_role_payments_20261017.md',
 'results/fresh_cusp_factor_payment_framework_20261017.json',
 'scripts/fresh_private_affine_tree_exclusion_20261017.py',
 'results/fresh_affine_carriers_order10_20261017.json',
 'results/fresh_affine_carriers_order10_20261017.md',
]

def build():
    # Full rank accounting: P>=4; D end order>=2, B end order>=1.
    # The only exceptional finite extensions of the inner classical block
    # are D6->E7 and D7->E8. Any connected join has full order m>=10.
    contexts=[];double=[]
    for m in range(10,18):
        for Pdist in range(4,m):
            for lk in 'BD':
                for rk in 'BD':
                    for r in range(2 if lk=='D'else 1,m):
                        t=m-r-Pdist-1
                        if t<(2 if rk=='D'else 1):continue
                        left='E7'if lk=='D'and r+Pdist==6 else'E8'if lk=='D'and r+Pdist==7 else None
                        right='E7'if rk=='D'and t+Pdist==6 else'E8'if rk=='D'and t+Pdist==7 else None
                        row=dict(order=m,distance=Pdist,left=[lk,r],right=[rk,t],left_exception=left,right_exception=right)
                        contexts.append(row)
                        if left==right=='E8':double.append(row)
    assert double==[
        dict(order=10,distance=5,left=['D',2],right=['D',2],left_exception='E8',right_exception='E8'),
        dict(order=11,distance=4,left=['D',3],right=['D',3],left_exception='E8',right_exception='E8')]
    # Literal remaining D11 dual-E8 configuration.  Ordinary nodes0..6,
    # left leaves7,8, right leaves9,10; marks a=1,b=5.
    edges,_=P.graph('D',11);G=P.gram(11,edges)
    a,b,c,d=1,5,11,12;x,y=s.symbols('x y',nonnegative=True)
    H=s.zeros(13);H[:11,:11]=G;H[c,c]=H[d,d]=2
    for i,j,t in[(c,10,1),(d,8,1),(a,c,x),(b,d,y)]:H[i,j]=H[j,i]=-t
    W=[i for i in range(11)if i not in(a,b)];T=[a,b,c,d]
    S=s.simplify(H.extract(T,T)-H.extract(T,W)*H.extract(W,W).inv()*H.extract(W,T))
    expected=s.Matrix([[s.Rational(1,4),-s.Rational(1,4),-x,-s.Rational(1,2)],
        [-s.Rational(1,4),s.Rational(1,4),-s.Rational(1,2),-y],
        [-x,-s.Rational(1,2),s.Rational(5,4),0],
        [-s.Rational(1,2),-y,0,s.Rational(5,4)]])
    assert S==expected
    u,v=s.symbols('u v',nonnegative=True);q=s.Rational
    F=s.expand(S.det().subs({x:1/u-q(1,2),y:1/v-q(1,2)})*u*u*v*v)
    assert s.simplify(F-(1-u-v-(u*u+v*v)/16+u*v/2))==0
    z=1/(s.sqrt(2)+q(1,2));margin=s.simplify(F.subs({u:z,v:z}))
    assert margin==q(181,98)-62*s.sqrt(2)/49 and margin.is_positive is True
    # Both derivatives <= -1+z/2 <0 throughout the closed rectangle.
    assert s.simplify(1-z/2).is_positive is True
    witnesses=[]
    for name,root,opp,leaf,verts0,verts1 in[
        ('c',c,a,10,[9,10,6,5,4,3,2,1,c],[1,2,3,4,5,6,10,c]),
        ('d',d,b,8,[7,8,0,1,2,3,4,5,d],[8,0,1,2,3,4,5,d])]:
        for value,verts in[(0,verts0),(1,verts1)]:
            K=H.subs({x:value if root==c else 0,y:value if root==d else 0}).extract(verts,verts)
            ns=K.nullspace();assert len(ns)==1
            nu=ns[0]
            if nu[0]<0:nu=-nu
            assert all(t>0 for t in nu)
            L=s.diag(*nu)*K*s.diag(*nu)
            assert all(sum(L.row(i))==0 for i in range(L.rows))
            assert all(L[i,j]<=0 for i in range(L.rows)for j in range(L.rows)if i!=j)
            witnesses.append(dict(root=name,opposite_product_magnitude=value,vertices=verts,
                                  order=len(verts),type='affine_E8'if value==0 else'affine_A7',
                                  positive_kernel=list(map(str,nu))))
    # The new c/d nodes really occupy the inverse-diagonal4 E8 role.
    E8rows=[]
    for root,verts in[(c,[2,3,4,5,6,9,10,c]),(d,[0,1,2,3,4,7,8,d])]:
        K=H.extract(verts,verts);assert K.det()==1 and K.inv()[-1,-1]==4
        E8rows.append(dict(new_root=root,vertices=verts,inverse_new_diagonal='4'))
    bound=s.Rational(1,9);assert 1-s.Rational(8,9)==bound>0
    paths=[Path(__file__)]+[ROOT/p for p in INPUTS]
    return dict(status='PASS_AFFINE_ORDER10_QUADRILATERAL_BOUNDARY_INGREDIENTS',
       sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
       minimum_affine_order=10,forbidden_affine_order_through=9,
       boundary_main_context_count=len(contexts),boundary_main_contexts=contexts,
       double_E8_contexts=double,dual_E8=dict(gram=[[str(t)for t in row]for row in H.tolist()],
         white=W,marked=T,schur=[[str(t)for t in row]for row in S.tolist()],
         heights=['x+1/2','y+1/2'],reciprocal_upper=str(s.simplify(z)),
         determinant_polynomial=str(F),strict_minimum=str(margin),small_affine_witnesses=witnesses,
         E8_roles=E8rows),two_terminal_minimum_margin=str(bound),
       scope='New boundary-order role accounting and literal dual-E8 determinant; all-rank ideal-corner classification/transfer is in report. No ambient dimension conclusion.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=json.loads(json.dumps(build()));out=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['boundary_main_context_count'],'contexts',data['dual_E8']['strict_minimum'],flush=True)
