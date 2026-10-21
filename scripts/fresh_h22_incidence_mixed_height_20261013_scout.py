#!/usr/bin/env python3
"""Bounded numerical volume scout for a new exact-formula mixed-height core.

The integer orbit/pair data are checked exactly. LP and geometry below are
proposals only; no finite-volume or certified escape conclusion is made here.
"""
import json, hashlib
from pathlib import Path
from collections import deque, Counter
import numpy as np
from scipy.optimize import linprog

ROOT=Path(__file__).resolve().parents[1]
PRE='fresh_h22_incidence_mixed_height_20261013'
SOURCE='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
BANK='data/enumerations/fresh_h22_incidence_general_binary_20261013_bank.json'
SYM='data/enumerations/fresh_h20_global_interior_sections_20261013.json'
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def setup():
 raw=read(SOURCE)['coxeter_input'];Q=-np.array(raw['basis_gram'],dtype=np.int64);W=np.array(raw['wall_coordinates_as_columns_transposed'],dtype=np.int64)
 rec=read(BANK)['candidates'];gens=[np.array(g['matrix'],dtype=np.int64) for g in read(SYM)['ambient_symmetry_generators']]
 for g in gens:assert np.array_equal(g.T@Q@g,Q)
 seeds=[tuple(r['twice_old_projection']) for r in rec if (r['A'],r['weight']) in ((3,9),(4,16))]
 seen=set(seeds);todo=deque(seeds)
 while todo:
  a=np.array(todo.popleft(),dtype=np.int64)
  for g in gens:
   b=tuple(g@a)
   if b not in seen:seen.add(b);todo.append(b)
 D2=np.array(sorted(seen),dtype=np.int64)
 X2=np.array([r['twice_old_projection'] for r in rec if (r['A'],r['weight']) in ((3,7),(4,14))],dtype=np.int64)
 assert D2.shape==(322,22) and X2.shape==(720,22)
 assert np.all(D2%2==0);D=D2//2
 assert max(abs(D).max(),abs(X2).max())<100000
 assert 22*int(abs(Q@D.T).max())*int(abs(D).max())<2**63
 assert np.all(np.einsum('ij,ij->i',D,D@Q)==0)
 assert np.all(np.einsum('ij,ij->i',X2,X2@Q)==-4)
 DX=D@Q@X2.T;assert set(DX.flat)=={-8,-6,-4}
 DD=D@Q@D.T;assert set(DD.flat)=={0,-2,-3,-4}
 XX4=X2@Q@X2.T
 # Equal-sign height-sqrt3 pair = XX4/4+3; norm two.
 pp=XX4+12;ok=(pp==0)|(pp==-4)|(pp<=-8)
 np.fill_diagonal(ok,False);assert set(ok.sum(axis=1))=={239}
 wn=np.einsum('ij,ij->i',W,W@Q)
 return Q,W,D,X2,ok,wn,XX4
def main():
 Q,W,D,X2,ok,wn,XX4=setup();a=(np.sqrt(6)+np.sqrt(3))/2
 core=read('results/fresh_h22_incidence_core_20261013.json');h=np.array(core['center'],float);h/=np.sqrt(14)
 Qn=np.zeros((23,23));Qn[:22,:22]=Q;Qn[22,22]=1
 old=np.column_stack([W,np.zeros(210)])
 null=np.vstack([np.column_stack([a*D,np.full(322,sgn*np.sqrt(2))]) for sgn in (1,-1)])
 vert=np.vstack([np.column_stack([X2/2,np.full(720,sgn*np.sqrt(3))]) for sgn in (1,-1)])
 center=np.r_[h,0];eq=(center@Qn)[None,:]
 def scout(R,objectives,seed):
  rng=np.random.default_rng(seed);A=R@Qn;best=None
  for k in range(objectives):
   lp=linprog(rng.normal(size=23),A_ub=A,b_ub=np.zeros(len(R)),A_eq=eq,b_eq=[-1],bounds=[(None,None)]*23,method='highs')
   if not lp.success:return {'lp_status':int(lp.status),'message':lp.message,'trial':k}
   q=float(lp.x@Qn@lp.x)
   if best is None or q>best['square']:best={'trial':k,'square':q,'point':lp.x.tolist(),'max_residual':float((A@lp.x).max()),'active':np.flatnonzero(abs(A@lp.x)<1e-7).tolist()}
   if q>1e-6:return best
  return best
 relaxed=scout(np.vstack([old,null,vert]),200,322720)
 print('relaxed',relaxed,flush=True)
 # Large maximal cliques are a bounded geometric scout, never completeness.
 rng=np.random.default_rng(322239);cliques=[];seen=set()
 for k in range(120):
  cand=list(range(720));cl=[]
  while cand:
   if k<60:
    ds=ok[np.ix_(cand,cand)].sum(axis=1);score=ds+rng.random(len(cand))*8
    v=cand[int(np.argmax(score))]
   else:v=cand[int(rng.integers(len(cand)))]
   cl.append(v);cand=[j for j in cand if ok[v,j]]
  key=tuple(sorted(cl))
  if key not in seen:seen.add(key);cliques.append(list(key))
 out=[]
 for i,cl in enumerate(sorted(cliques,key=lambda x:(-len(x),x))[:12]):
  ids=cl+[j+720 for j in cl];r=scout(np.vstack([old,null,vert[ids]]),80,322000+i)
  out.append({'ids':cl,'both_height_signs':True,'walls':854+2*len(cl),'volume_scout':r})
  print('clique',i,len(cl),r['square'] if 'square'in r else r,flush=True)
 result={'status':'NUMERICAL_MIXED_SCALED_NULL_SQRT3_SCOUT','sources':{p:sha(p) for p in (SOURCE,BANK,SYM)},'null_scaling':'(sqrt(6)+sqrt(3))/2','null_projections':D.tolist(),'sqrt3_twice_projections':X2.tolist(),'null_old_pair_counts':dict(Counter(map(str,(D@Q@W.T).flat))),'sqrt3_graph_degree':239,'clique_sizes':dict(Counter(map(str,map(len,cliques)))),'relaxed_all720':relaxed,'selected_cliques':out,'scope':'Integer projection data and combinatorial clique gates are exact; all LP volume outputs are numerical proposals only.'}
 (ROOT/('data/enumerations/'+PRE+'_scout.json')).write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
