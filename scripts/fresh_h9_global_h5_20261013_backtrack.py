#!/usr/bin/env python3
"""Full H5 geometry of coherent changed cores from opposite branches."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0,str(Path(__file__).resolve().parent))
import fresh_h9_global_h5_20261013_scout as S
def run():
 R,G,c,names,d,p,C,NR,refs,sources,counts=S.load();bank=json.loads((S.P/'data/enumerations'/f'{S.TAG}_opposite_scout.json').read_text());rows=bank['candidates'];results=[]
 for i,row in enumerate(rows):
  drop={x['wall']for x in row['conflicts']};keep=[j for j,n in enumerate(names)if n not in drop];r=np.array(row['normal']);W=np.vstack([R[keep],r]);vs=S.S.S.geometry(W,G,p,C);qq=np.array([v['norm']for v in vs]);finite=[v for v in vs if v['norm']<-1e-7];ids=[j for j in range(len(W))if abs(W[j]@G@r)<2-1e-7];edges={}
  fixed=set(d['H5_face']['fixed_walls']);newnames=[names[j]for j in keep]+['opposite']
  for vi,v in enumerate(finite):
   incident=[j for j in v['active']if newnames[j]not in fixed]
   if len(incident)!=5:continue
   for missing in range(5):edges.setdefault(tuple(incident[:missing]+incident[missing+1:]),[]).append(vi)
  results.append({'index':i,'removed':sorted(drop),'walls':len(W),'projective_vertices':len(vs),'finite_vertices':len(finite),'spacelike_vertices':int(sum(qq>1e-7)),'maximum_norm':float(max(qq)),'finite_edges':sum(len(v)==2 for v in edges.values()),'open_edges':sum(len(v)==1 for v in edges.values()),'leading_point':vs[int(np.argmax(qq))]['point']})
 results.sort(key=lambda r:(r['maximum_norm'],r['open_edges']));print('best',results[:3],flush=True)
 # Strongest relaxed bank, including mutually conflicting opposite branches.
 ALL=np.vstack([R,NR,np.array([row['normal']for row in rows])]);A=ALL@G@C;b=-ALL@G@p;keep=np.linalg.norm(A,axis=1)>1e-7;assert min(b[keep])>0;aa=A[keep]/b[keep,None];checks=[];escape=None
 old=S.S.S.geometry(R,G,p,C);old.sort(key=lambda v:v['norm'],reverse=True)
 for i,v in enumerate(old):
  direction=C.T@G@np.array(v['point']);res=linprog(-direction,A_ub=aa,b_ub=np.ones(len(aa)),bounds=[(None,None)]*5,method='highs')
  assert res.success;z=res.x;q=float(z@z-1);checks.append(q)
  if q>1e-5:
   zz=.999*z;escape={'norm':float(zz@zz-1),'point':(p+C@zz).tolist(),'maximum_normalized_residual':float(max(aa@zz-1))};break
 out={'status':'NUMERICAL_GLOBAL_CHANGED_CORE_COMPARISON','coherent_cases':results,'combined_relaxation_candidates':len(NR)+len(rows),'combined_bank_escape':escape,'lp_maxima':checks,'scope':'Numerical geometry and finite-bank relaxation only. Removal sets preserve the original seed and fixed H5 normals but may destroy earlier compact faces. No exact candidate or exhaustion claim.'}
 (S.P/'data/enumerations'/f'{S.TAG}_backtrack.json').write_text(json.dumps(out,indent=2)+'\n');print('combined',escape,flush=True)
if __name__=='__main__':run()
