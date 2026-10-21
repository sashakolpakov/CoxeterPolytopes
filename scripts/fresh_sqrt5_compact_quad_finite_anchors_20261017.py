#!/usr/bin/env python3
"""Exact finite anchors for one opposite product of a compact arithmetic quad."""
import argparse,hashlib,json
from fractions import Fraction as F
from pathlib import Path
import sympy as s
import fresh_sqrt5_adjacent_discount_witness_20261016 as A

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_compact_quad_finite_anchors_20261017'
K=A.K;PHI=K(F(1,2),F(1,2))

def run():
    # Generic Schur determinant identity; no diagram producer is imported.
    a,b,c,d,p,q,r,t,u,v=s.symbols('a b c d p q r t u v')
    M=s.Matrix([[a,p,u,t],[p,b,q,v],[u,q,c,r],[t,v,r,d]])
    det=s.Poly(M.det(),u,v)
    assert s.Poly(det.as_expr(),v).coeff_monomial(v**2)==u**2-a*c
    # If both full opposite products <=-6, white elimination only decreases
    # them. The two opposite-root sums have diagonal <=-8 and offdiagonal
    # strictly between -8 and0, hence a negative definite restriction.
    assert 4+2*(-6)==-8
    x=s.symbols('x',nonnegative=True)
    assert s.expand(64-(-8+x)**2)==16*x-x**2
    finite_squares=[K(0),K(1),K(2),1+PHI,K(3),2+PHI]
    direct=[];square=[];discarded_direct=[];discarded_square=[]
    for i in range(-6,7):
        for j in range(-3,1):
            g=K(i)+j*PHI
            if(g+6).sign()<=0 or g.sign()>0 or(g.conj()+2).sign()<=0 or(2-g.conj()).sign()<=0:continue
            if(g+2).sign()>=0 and g not in(K(0),K(-1),-PHI):continue
            record=dict(phi_coefficients=[i,j],value=g.enc(),square=(g*g).enc())
            if not g.b:discarded_direct.append(record)
            else:direct.append(record)
    for i in range(-34,39):
        for j in range(-1,18):
            z=K(i)+j*PHI
            if z.sign()<0 or(36-z).sign()<=0 or z.conj().sign()<0 or(4-z.conj()).sign()<=0:continue
            if(z-4).sign()<=0 and z not in finite_squares:continue
            record=dict(phi_coefficients=[i,j],square=z.enc())
            if not z.b:discarded_square.append(record)
            else:square.append(record)
    assert len(direct)==8 and len(square)==59
    assert len(discarded_direct)==2 and len(discarded_square)==4
    # Every direct product's square is also a member of the larger domain.
    assert all(r['square']in[x['square']for x in square]for r in direct)
    paths=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_adjacent_discount_witness_20261016.py',
        'results/fresh_sqrt5_mixed_d4_d5_migration_20261017.json']
    return dict(status='PASS_EXACT_COMPACT_ARITHMETIC_QUADRILATERAL_FINITE_ANCHOR_DOMAIN',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        physical_bound='At least one full opposite norm-two product is strictly greater than -6.',
        direct_product_anchors=direct,squared_product_anchors=square,
        rational_direct_excluded=discarded_direct,rational_squares_excluded=discarded_square,
        direct_coefficient_box=dict(a=[-6,6],b=[-3,0],proof='For g=a+b*phi, -6<g<=0 and -2<g_conjugate<2 imply -4<b<1, so b=-3..0. Since1<phi<2, then -6<a<6; the displayed box safely contains every integer possibility.'),
        square_coefficient_box=dict(a=[-34,38],b=[-1,17],proof='For z=a+b*phi,0<=z<36 and0<=z_conjugate<4 imply -2<b<18, so b=-1..17. Since1<phi<2, then -34<a<38; the displayed box safely contains every integer possibility.'),
        normalized_finite_squared_products=[x.enc()for x in finite_squares],
        rational_exclusion_scope='For the applications with all other edge contacts0/-1, a rational direct framing makes the edge Gram rational. If only its square is rational and the product is not inK, any fixed3-edge cycle is forbidden; without such a cycle its Schur residual is rational, so physical and conjugate signs agree. The8/59 filtered domains use this ADE-edge premise. The larger10/63 domains before these exclusions do not need that premise.',
        identity_specialization=dict(leading_coefficient='coefficient_v^2(det S) = -det S[{0,2}]',
            null_lift='At a singular opposite2-block, positive diagonals and nonpositive offdiagonal give a strictly positive residual kernel y. White coefficients -W_inverse*B*y are nonnegative.',
            nonzero='The resulting null combination of wall normals has negative product with a common strict timelike center, so it is nonzero and future.',
            feasibility='Every wall product is nonpositive by acuteness and the nonnegative coefficients. Hence this is a feasible ideal ray, forbidden by compactness.',
            conclusion='An identically-zero specialization in the other framing parameter cannot occur in an actual compact acute quadrilateral. Do not infer this for a single edge or a noncompact face.'),
        scope='A complete finite choice for at least ONE framing product, followed by its nonzero quadratic rank equation for the other. This package does not itself solve a diagram, bound both products independently, classify all quad patterns, or prove a dimension bound. The8/59 filtered domain additionally assumes all fixed edge contacts are0/-1.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    data=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],len(data['direct_product_anchors']),len(data['squared_product_anchors']),flush=True)
