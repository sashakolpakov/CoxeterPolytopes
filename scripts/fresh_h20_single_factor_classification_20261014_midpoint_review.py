#!/usr/bin/env python3
"""Separate SymPy replay of the 28 delicate midpoint witnesses.

This is additional arithmetic corroboration of a source/theorem review and
fresh producer replay, not a second complete graph-cut census.
"""
from pathlib import Path
import argparse,hashlib,json,math
import sympy as S
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_single_factor_classification_20261014'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
WALK='results/fresh_h21_source_corner_walk_20261013.json'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 data=json.loads((ROOT/f'results/{TAG}.json').read_text())
 for p,h in data['sources'].items():assert sha(p)==h,p
 assert data['status']=='PASS_COMPLETE_ALL_RETAINED_SINGLE_FACTOR_CLASSIFICATION_AND_GLOBAL_EXCLUSION'
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input'];Q=-S.Matrix(raw['basis_gram']).applyfunc(S.Rational)
 W=S.Matrix(raw['wall_coordinates_as_columns_transposed']).applyfunc(S.Rational);C=W*Q
 h=sum((W.row(i).T for i in range(210)),S.zeros(22,1));g=math.gcd(*map(int,h));h=h/g
 assert (h.T*Q*h)[0]==-14
 walk=json.loads((ROOT/WALK).read_text())['corner_orbits']
 recipes=json.loads((ROOT/f'results/{TAG}_recipes.json').read_text());assert len(recipes)==28
 rows=[]
 for r in recipes:
  cusp=walk[r['cusp']];comp=next(x for x in cusp['affine_components']if x['walls']==r['factor'])
  marks=dict(zip(comp['walls'],map(S.Rational,comp['positive_kernel']),strict=True))
  D=sum((marks[i]*W.row(i).T for i in r['factor']),S.zeros(22,1))
  e=sum((marks[i]*W.row(i).T for i in r['active_part']),S.zeros(22,1))
  p=S.Matrix(r['midpoint']).applyfunc(S.Rational);values=C*p
  q=(p.T*Q*p)[0];future=(p.T*Q*h)[0]
  assert q==S.Rational(r['norm'])<0 and future<0
  assert max(values)==0 and [i for i,x in enumerate(values)if x==0]==sorted(r['ridge'])
  assert (p.T*Q*(e-D/2))[0]==0 and (p.T*Q*D)[0]<0
  i,j=r['ridge'];assert (W.row(i)*Q*W.row(j).T)[0]==0
  assert (W.row(i)*Q*e)[0]>0 and (W.row(j)*Q*e)[0]==(W.row(j)*Q*D)[0]<0
  rows.append({'cusp':r['cusp'],'active_part':r['active_part'],'ridge':r['ridge'],
   'norm':str(q),'future_pairing':str(future),'null_pairing':str((p.T*Q*D)[0]),
   'strict_other_inequalities':208})
 assert len(rows)==28 and data['counts']['two_half_interval_witnesses']==14
 # Check the determinant interpolation algebra independently.
 q,u,b,di,dj,t=S.symbols('q u b di dj t')
 ni=di-u*u/q;nj=dj-b*b*(1+t)**2/q;p=u*b*(1+t)/q
 assert S.factor(ni*nj-p*p-(ni*dj-b*b*di*(1+t)**2/q))==0
 sources=[RAW,WALK,f'scripts/{TAG}.py',f'results/{TAG}.json',f'results/{TAG}_recipes.json',f'scripts/{TAG}_midpoint_review.py']
 return {'status':'PASS_SOURCE_THEOREM_REVIEW_AND_SEPARATE_28_MIDPOINT_REPLAY',
  'sources':{p:sha(p)for p in sources},'producer_counts':data['counts'],'independent_midpoints':rows,
  'reviewed_premises':['two affine blocks from the all-retained bounded Euclidean link',
   'rank-one crossing support and exhaustive star-cut coloring',
   'positive affine marks force the normal modulo the null line',
   'rational squared-cosine gate includes every integer finite label and parallelism',
   'full and half-interval interpolations establish actual obtuse ridges'],
  'scope':'Source/proof review plus a fresh producer saved-result replay. Separate SymPy arithmetic checks every28 delicate midpoint and the determinant formula; no second complete287-cut enumeration or1209-face reconstruction.'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();z=run();out=ROOT/f'results/{TAG}_midpoint_review.json'
 if args.verify:assert z==json.loads(out.read_text())
 else:out.write_text(json.dumps(z,indent=2)+'\n')
 print(z['status'])
