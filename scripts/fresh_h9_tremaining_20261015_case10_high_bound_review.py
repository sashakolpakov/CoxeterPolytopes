#!/usr/bin/env python3
"""Independent exact determinant and residual-plane review of the T10 high-label gate."""
import argparse,json
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_tremaining_20261015_case10_high_bound_review'
PRIMARY='results/fresh_h9_tremaining_20261015_case10_high_bound.json'
def run():
 p=json.loads((ROOT/PRIMARY).read_text());sha=lambda f:sha256((ROOT/f).read_bytes()).hexdigest()
 for f,h in p['sources'].items():assert sha(f)==h,f
 gp='results/fresh_h9_tree_cube_wedges_20261014.json';g=json.loads((ROOT/gp).read_text())['representatives'][10]['geometry'];V=[set(r['incident'])for r in g['vertices']]
 sp='results/fresh_h9_tremaining_20261015_search_case_10_batch_001.json';state=json.loads((ROOT/sp).read_text());labels={(a,b):m for a,b,m in state['solution']['pair_labels']}
 blocks=[[3,4],[6,7],[9,10],[12,13]];S=sum(blocks,[]);assert p['positive_support']==S
 expected={tuple(x)for x in combinations(S,2)if list(x)not in blocks}|{(9,11),(10,11),(3,16),(4,16),(6,16),(7,16),(11,16)}
 assert len(expected)==31 and {tuple(r['pair'])for r in p['forced_right_guards']}==expected
 cos={2:s.Integer(0),3:s.Rational(1,2),4:s.sqrt(2)/2,6:s.sqrt(3)/2};records=[]
 for row in p['forced_right_guards']:
  ids=row['support'];assert set(ids)==V[row['actual_vertex']] and set(row['pair'])<=set(ids)
  M=s.eye(len(ids))
  for i,j in combinations(range(len(ids)),2):
   pair=tuple(sorted((ids[i],ids[j])));m=3 if pair==tuple(row['pair'])else labels.get(pair,2)
   if m==7:m=6
   M[i,j]=M[j,i]=-cos[m]
  minors=[]
  for k in range(1,len(ids)+1):
   z=s.simplify(M[:k,:k].det());sign=s.sign(z);assert sign in(-1,0,1);minors.append(str(z))
   if sign<=0:break
  assert s.sign(s.sympify(minors[-1]))<=0
  records.append(dict(pair=row['pair'],actual_vertex=row['actual_vertex'],leading_minors=minors))
 u=s.sqrt(2)/2;J=s.Matrix([[1,-u],[-u,1]]).inv();contact=s.Matrix([-s.Rational(1,2)]*2)
 assert s.simplify((contact.T*J*contact)[0]-(1+u))==0
 assert [labels[a,b]for a,b in blocks]==[4,4,7,4]
 assert all(labels[min(i,j),max(i,j)]==3 for i in(11,16)for j in(12,13))
 assert labels[9,16]==2 and labels[10,16]==3
 t,F=s.symbols('t F',nonnegative=True);d=2+s.sqrt(3)+t
 residual=s.Matrix([[-u-F,-1-u],[-1-u,-u-d]])
 det=s.expand(residual.det());base=s.sqrt(6)/2-1
 assert s.simplify(det-base-F*(u+d)-u*t)==0 and base.is_positive
 # All terms beyond the positive boundary margin are nonnegative.
 assert p['remaining_integer_high_labels']==[7,8,9,10,11]
 files=[PRIMARY,gp,sp,f'scripts/{TAG}.py']
 return dict(status='PASS_INDEPENDENT_T10_HIGH_LABEL_DETERMINANT_REVIEW',sources={f:sha(f)for f in files},forced_right_guards=records,
   residual_determinant=str(det),m12_positive_margin=str(base),remaining_high_labels=[7,8,9,10,11],
   scope='No producer arithmetic imports. Exact leading determinants recheck all31 actual-vertex guards at the relaxed highlabel6. Symmetric Z-matrix Perron monotonicity extends failure to larger contacts. The positive8root projection yields a negative residual plane for every m>=12. Only the saved mandatory family is covered; labels7..11 and full T10 remain open.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
