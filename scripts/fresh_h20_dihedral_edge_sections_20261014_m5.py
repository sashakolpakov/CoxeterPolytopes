#!/usr/bin/env python3
"""Exact mapping of both I2(5)->H3 orientations to certified obstructions.

This is a source-bound corollary, not a new implementation of the coupled
H3 irreducibility census. Its normalizations and complete old rank-three
possibilities are computed separately with SymPy.
"""
from pathlib import Path
import argparse, hashlib, json
import sympy as S
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_dihedral_edge_sections_20261014_m5'
H3='results/fresh_h20_coupled_h3_obstruction_20261014.json'
Q5='results/fresh_h20_h4_face_closure_20261014.json'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 h3=json.loads((ROOT/H3).read_text());old=json.loads((ROOT/Q5).read_text())
 assert h3['status']=='PASS_EXACT_COUPLED_H3_ALGEBRA'
 for p,h in h3['source_sha256'].items():assert sha(p)==h,p
 for p,h in old['sources'].items():assert sha(p)==h,p
 q=(1+S.sqrt(5))/4;C=(3+S.sqrt(5))/8;x=S.Symbol('x');assert S.simplify(1/q-(2-1/(2*C)))==0
 cases=[]
 for d in[2,4]:
  for g in([-1,0]if d==2 else[-2,0]):
   G=S.Matrix([[2,-1,g],[-1,2,0],[g,0,d]])
   if not all(G[:i,:i].det()>0 for i in range(1,4)):continue
   poly=(3+2*q)*x*x-8*g*q*x+4*g*g*q*q-d*q*(2*q-1)
   projected_norm_a=2-1/q
   equation=4*(g-x/q)**2-projected_norm_a*(d-x*x/q)
   assert S.simplify(q*q*equation-poly)==0
   if g==0:
    Y=1/(2*q);X=(1-Y)/(1+3*Y)
    assert S.simplify(poly.subs(x*x,d*q*X))==0
    status='JOINED_A1_PLUS_A2_H3_THEOREM'
   else:
    A=1/q;ratio_poly=A*(3*A+2)*x*x-8*g*A*x+4*g*g-(2-A)*d
    assert S.simplify(ratio_poly*q*q-poly)==0
    shape=[2,d,2,g,-1]
    assert shape in[r['old_dc_df_dl_gcf_gcl']for r in h3['two_contact_cases']]
    status='CONNECTED_TWO_CONTACT_H3_THEOREM'
   cases.append({'old_new_wall_norm':d,'old_anchor_pairing':g,'old_other_pairing':0,
    'old_rank3_Gram':list(map(list,G.tolist())),'necessary_contact_equation':str(S.expand(poly)),
    'acute_branch_inequality':f'x > ({g})*q','mapped_obstruction':status})
 assert [(r['old_new_wall_norm'],r['old_anchor_pairing'])for r in cases]==[(2,-1),(2,0),(4,-2),(4,0)]
 # The other orientation is orthogonal to the modified anchor and has
 # label3 against the unchanged old b. Its norm equation forces old A3.
 other=[]
 for d,k in[(2,0),(2,-1),(4,0),(4,-2)]:
  for gb in([-1]if d==2 else[-2]):
   if S.simplify(d-q*k*k-2*gb*gb)==0:other.append([d,k,gb])
 assert other==[[2,0,-1]]
 # In the absence of H3 boundary extensions, the previously checked q5
 # table makes every finite neighbor of a modified facet right and original.
 allowed={'FORBIDDEN_COXETER_ANGLE','NO_FINITE_INTERSECTION','ORTHOGONAL_TO_UNMODIFIED_RATIONAL_ROOT'}
 verdicts=sorted({r['verdict']for r in old['modified_facet_adjacency_table']})
 assert set(verdicts)<=allowed,verdicts
 sources=[H3,Q5,'scripts/fresh_h20_coupled_h3_obstruction_20261014.py',
  'scripts/fresh_h20_h4_face_closure_20261014.py',f'scripts/{TAG}.py']
 def clean(v):
  if isinstance(v,S.Basic):return int(v)if v.is_Integer else str(v)
  if isinstance(v,list):return[clean(x)for x in v]
  if isinstance(v,dict):return{k:clean(x)for k,x in v.items()}
  return v
 return clean({'status':'PASS_ALL_REAL_ONE_CONTACT_A2_TO_I2_5_COROLLARY',
  'sources':{p:sha(p)for p in sources},'section_norm':str(q),
  'extension_at_modified_anchor_cases':cases,'extension_at_unchanged_end_only_case':other,
  'orthogonal_boundary_table_verdicts':verdicts,
  'scope':'All real transverse one-contact A2->I2(5) sections in the literal source. Relies on the frozen one-contact A3 and coupled-H3 obstruction proofs and their exact algebra certificates; separately checks the complete orientation/normalization reduction.'})
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();z=run();out=ROOT/f'results/{TAG}.json'
 if args.verify:assert z==json.loads(out.read_text())
 else:out.write_text(json.dumps(z,indent=2)+'\n')
 print(z['status'])
