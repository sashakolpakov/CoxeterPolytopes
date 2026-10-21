#!/usr/bin/env python3
"""Exact complete spherical extension pool for a one-contact A2->B2 face.

The geometric finite-reduction theorem is in the matching report. This
file independently checks the literal source pool and all eight contact
recipes, without importing the floating-point scout.
"""
from pathlib import Path
from collections import Counter
import argparse, hashlib, json
import sympy as S
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_b2_edge_sections_20261014'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def sign(x):
 v=S.simplify(x)
 if v==0:return 0
 if v.is_positive:return 1
 if v.is_negative:return -1
 raise AssertionError(v)
def run():
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input'];Q=-S.Matrix(raw['basis_gram']).applyfunc(S.Rational)
 W=S.Matrix(raw['wall_coordinates_as_columns_transposed']).applyfunc(S.Rational);G=W*Q*W.T
 assert W.shape==(210,22)
 pools={}
 for d1,d2 in[(2,2),(2,4),(4,4)]:
  pools[d1,d2]=sorted({G[i,j]for i in range(210)for j in range(i)if sorted([G[i,i],G[j,j]])==[d1,d2]})
 assert pools=={(2,2):[-1,0],(2,4):[-2,0],(4,4):[-8,-6,-4,-2]}
 A=S.Matrix([[2,-1],[-1,2]]);assert A.inv()[0,0]==S.Rational(2,3)
 H=A-S.diag(1,0);assert H.det()==1 and 1/(H[0,0]*H[1,1])==S.Rational(1,2)
 x=S.Symbol('x');rows=[];excluded=[]
 for d in[2,4]:
  pool=pools[2,d]
  for ga in pool:
   for gb in pool:
    old=S.Matrix([[2,-1,ga],[-1,2,gb],[ga,gb,d]])
    if old.det()<=0:
     excluded.append({'d':d,'ga':int(ga),'gb':int(gb),'old_determinant':int(old.det())});continue
    proposals=[]
    if gb==0:
     proposals.append(('ORTHOGONAL',ga))
     poly=5*x*x-8*ga*x+4*ga*ga-d
     for root in S.solve(poly,x):
      if sign(root-ga)>0:proposals.append(('B3_AT_MODIFIED_END',root))
    if gb<0 and d-ga*ga==2*gb*gb:
     proposals.append(('B3_AT_UNCHANGED_END',ga))
    for kind,t in proposals:
     if sign(d-t*t)<=0:continue
     contact=S.Matrix([1,0,t]);K=S.simplify(old-contact*contact.T)
     assert all(sign(K[:j,:j].det())>0 for j in range(1,4))
     assert all(sign(K[i,j])<=0 for i in range(3)for j in range(i))
     ca=S.simplify(K[0,2]**2/(K[0,0]*K[2,2]));cb=S.simplify(K[1,2]**2/(K[1,1]*K[2,2]))
     expected={'ORTHOGONAL':(0,0),'B3_AT_MODIFIED_END':(S.Rational(1,4),0),'B3_AT_UNCHANGED_END':(0,S.Rational(1,4))}[kind]
     assert (ca,cb)==expected
     target=S.simplify(t-(2*ga+gb)/3)
     rows.append({'d':d,'ga':int(ga),'gb':int(gb),'kind':kind,'normal_contact':str(t),
      'edge_plane_contact':str(target),'positive_projected_norm':str(S.simplify(d-t*t)),
      'projected_Gram':[[str(S.simplify(v))for v in K.row(i)]for i in range(3)],
      'positive_LDL_determinants':[str(S.simplify(K[:j,:j].det()))for j in range(1,4)]})
 assert len(rows)==8
 assert Counter(r['kind']for r in rows)=={'ORTHOGONAL':3,'B3_AT_MODIFIED_END':4,'B3_AT_UNCHANGED_END':1}
 return {'status':'PASS_COMPLETE_EIGHT_B2_FACE_CONTACT_RECIPES','sources':{p:sha(p)for p in[RAW,f'scripts/{TAG}.py']},
  'section_norm':'1','old_dual_weight_norm':'2/3','edge_plane_norm':'1/3',
  'source_pairing_pools':{str(k):list(map(int,v))for k,v in pools.items()},
  'recipes':rows,'nonpositive_old_Gram_cases':excluded,
  'scope':'Complete possible extra bounding facet contacts of an actual B2 face in a hypothetical Coxeter section. Every such facet lies at a finite face point and extends the induced B2 to B2+A1 or B3. The numerical edge-position census and global ridge verification are separate.'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();z=run();out=ROOT/f'results/{TAG}.json'
 if args.verify:assert z==json.loads(out.read_text())
 else:out.write_text(json.dumps(z,indent=2)+'\n')
 print(z['status']);print(len(z['recipes']))
