#!/usr/bin/env python3
"""Separate symbolic reconstruction of the coupled-H3 angle identities.

Uses SymPy, but imports neither the producer nor its arithmetic helpers.
Geometric completeness and invariant-lattice arguments are reviewed in prose.
"""
from pathlib import Path
import argparse, hashlib, json
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_coupled_h3_obstruction_20261014'
SQ=s.sqrt(5);C5=(3+SQ)/8
X,U,R=s.symbols('X U R')
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def simp(x):return s.cancel(s.simplify(x),extension=SQ)
def monic(x,z):return s.Poly(x,z,extension=SQ).monic().as_expr()
def normcheck(p,z,record):
 p=monic(p,z);n=s.Poly(s.expand(p*p.xreplace({SQ:-SQ})),z,domain=s.QQ).monic()
 assert n.degree()==record['degree']and n.is_irreducible
 assert list(map(str,n.all_coeffs()))==record['coefficients_descending']
 return n.degree()
def run():
 packet=json.loads((ROOT/'results'/f'{TAG}.json').read_text())
 for p,h in packet['source_sha256'].items():assert sha(p)==h
 assert packet['status']=='PASS_EXACT_COUPLED_H3_ALGEBRA'
 degrees=[]
 for row in packet['joined_cases']:
  cold=s.sympify(row['old_cosine_squared']);xx=s.sympify(row['X']);yy=s.sympify(row['Y'])
  assert simp(cold/(1-yy)-C5)==0
  assert simp(xx*yy/((1-xx)*(1-yy))-s.Rational(1,4))==0
  ratio=simp(xx/yy);norm=simp(ratio*ratio.xreplace({SQ:-SQ}))
  assert norm==s.sympify(row['field_norm_of_ratio'])and not s.sqrt(norm).is_Rational
  assert simp((1-xx)/(1+3*xx)-yy)==0
  assert s.expand(s.radsimp(xx)).coeff(SQ)!=0
  assert s.expand(s.radsimp(yy*yy)).coeff(SQ)!=0
 for row in packet['two_contact_cases']:
  dc,df,dl,gf,gl=row['old_dc_df_dl_gcf_gcl']
  A=s.sympify(row['center_contact_squared'])
  assert simp(gl*gl/(dl*(dc-A))-C5)==0
  # Reconstruct from the target label3 Gram determinant rather than
  # assume the displayed quadratic coefficients.
  P=monic(4*(gf-A*R)**2-(dc-A)*(df-A*R*R),R)
  assert simp(P-s.sympify(row['ratio_polynomial_over_Qsqrt5']))==0
  degrees.append(normcheck(P,R,row['contact_ratio']))
  B=monic(s.resultant(P,U-A*R*R,R),U)
  degrees.append(normcheck(B,U,row['leaf_contact_squared']))
 for row in packet['three_contact_cases']:
  dc,df,dl,gf,gl=row['old_dc_df_dl_gcf_gcl']
  c1,c2=(C5,s.Rational(1,4))if row['five_on_leaf_pair']else(s.Rational(1,4),C5)
  Z=s.Symbol('Z')
  # fl=-XR and cl=gl-gf R after the exact cf cancellation.
  first=X*X*Z-c1*(df-X)*(dl-X*Z)
  r2=simp(s.solve(first,Z)[0])
  assert simp(r2-s.sympify(row['R_squared']))==0
  second=s.expand((gl-gf*R)**2-c2*(dc-gf*gf/X)*(dl-X*R*R))
  constant=simp(second.subs(R*R,r2).subs(R,0))
  assert simp(second.subs(R*R,r2)-constant+2*gl*gf*R)==0
  P=monic(simp(constant**2-4*gl*gl*gf*gf*r2).as_numer_denom()[0],X)
  assert simp(P-s.sympify(row['quartic_over_Qsqrt5']))==0
  degrees.append(normcheck(P,X,row['fixed_leaf_contact_squared']))
  quantity=r2 if row['five_on_leaf_pair']else simp(X*X*r2)
  num,den=quantity.as_numer_denom();B=monic(s.resultant(P,U*den-num,X),U)
  degrees.append(normcheck(B,U,row['five_edge_contact_ratio_gate']))
 A=simp(2-1/(2*C5));T=1/A;B=simp(2*(2-T)/(2+3*T))
 assert simp(4*B*T-(2-B)*(2-T))==0
 S=simp(((2-A)*(2-B)/4-1-A*B)/2)
 assert simp(1+2*S+A*B-(2-A)*(2-B)/4)==0
 defect=simp(S*S-A*B)
 assert simp(defect-s.sympify(packet['D4_extension_defect']))==0 and defect!=0
 for row in packet['B_endpoint_bounds']:
  da=row['anchor_norm'];d=row['wall_norm'];q=s.sympify(row['normal_norm']);g=row['smallest_nonzero_pairing']
  assert simp(1/(2*(1-1/(q*da)))-C5)==0
  assert s.sign(simp(d-q*g*g))==-1
 assert sorted(degrees)==[4]*10+[8]*20
 sources=[f'scripts/{TAG}.py',f'results/{TAG}.json',f'scripts/{TAG}_review.py']
 return {'status':'PASS_SEPARATE_SYMBOLIC_H3_ANGLE_RECONSTRUCTION',
  'sources':{p:sha(p)for p in sources},'irreducible_polynomials':len(degrees),
  'quartic_norm_polynomials':degrees.count(4),'degree8_norm_polynomials':degrees.count(8),
  'scope':'Independent formula reconstruction using the same SymPy CAS, without producer imports. Geometric completeness and actual compact-face/lattice implications are separately reviewed in prose.'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/'results'/f'{TAG}_review.json'
 if args.verify:assert json.loads(p.read_text())==out
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
