#!/usr/bin/env python3
"""Numerical whole-chamber gap proposal for the symmetric mate replacement."""
import json
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_contact_closure_20261014'
def run():
 D=json.loads((ROOT/f'results/{TAG}.json').read_text());C=json.loads((ROOT/f'results/{TAG}_backtracking_scout.json').read_text());ff=lambda x:float(F(x['a']))+2**.5*float(F(x['b']));Q=np.array([[ff(x)for x in row]for row in D['form']]);N=np.array([[ff(x)for x in row]for row in D['normals']]);H=np.array([ff(x)for x in D['strict_center']]);E=np.array([ff(x)for x in D['strict_escape']['point']]);removed=list(range(61,76,2));keep=[i for i in range(76)if i not in removed];caps=np.array([C['candidates'][i]['normal']for i in[1,3,5,7]]);other=caps.copy();other[:,19]*=-1;caps=np.vstack([caps,other]);N=np.vstack([N[keep],caps]);norms=np.einsum('ij,jk,ik->i',N,Q,N);G=N@Q@N.T/np.sqrt(norms[:,None]*norms[None,:]);bad=[]
 for i in range(68,76):
  for j in range(i):
   x=G[i,j]
   if abs(x)<1e-7 or abs(x+.5)<1e-7 or x<=-1+1e-7:continue
   bad.append([i,j,float(x)])
 assert not bad
 h=H/np.sqrt(-H@Q@H);_,_,V=np.linalg.svd((h@Q)[None,:]);X=V[1:].T;L=np.linalg.cholesky(X.T@Q@X);X=X@np.linalg.inv(L.T);rng=np.random.default_rng(20261014);esc=[]
 for _ in range(64):
  direction=rng.normal(size=20);direction/=np.linalg.norm(direction)
  for it in range(4):
   lp=linprog(-direction,A_ub=N@Q@X,b_ub=-N@Q@h,bounds=[(None,None)]*20,method='highs')
   if not lp.success:break
   direction=lp.x/np.linalg.norm(lp.x)
  if lp.success:
   point=h+X@lp.x;q=float(point@Q@point)
   if q>1e-7:esc.append(dict(point=point.tolist(),q=q,max_pair=float(np.max(N@Q@point))))
 out=dict(status='NUMERICAL_76_WALL_SYMMETRIC_REPLACEMENT_GAP_PROPOSALS',removed_prior_indices=removed,retained_prior_indices=keep,candidate_indices=[1,3,5,7],append_height19_reflections=True,escapes=sorted(esc,key=lambda e:-e['q'])[:5],scope='Numerical proposals only; exact acceptance belongs to separate replacement.py.')
 print(out['status'],'escapes',[r['q']for r in out['escapes']]);return out
if __name__=='__main__':
 d=run();(ROOT/f'results/{TAG}_replacement_scout.json').write_text(json.dumps(d,indent=2)+'\n')
