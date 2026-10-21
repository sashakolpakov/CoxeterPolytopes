#!/usr/bin/env python3
"""Seek an exact conic obstruction to incident walls cutting the fixed escape.

An LP proposes a support; all successful identities and coefficient signs are
then reconstructed in Q(sqrt(11)). An infeasible selected support is merely
a failed scout, not a classification.
"""
import argparse,json
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
from fresh_h22_coupled_cusp_20261013_verify import P,K,read,sha,RAW,LOCAL,ESCAPE,PRE,lift

PREFIX='fresh_h22_incident_support_20261013'
SCOUT='data/enumerations/'+PRE+'_other_scout.json'
NULLS='data/enumerations/fresh_h22_cusp_orbits_20261013.npz'
LOWER='results/fresh_h22_new_cusp_completion_20261013_lower.json'

def sources():
 raw=read(RAW)['coxeter_input'];W=raw['wall_coordinates_as_columns_transposed']
 q=-np.array(raw['basis_gram'],dtype=np.int64);h=np.array(W,dtype=np.int64).sum(0);h//=np.gcd.reduce(h)
 def arrays(path):
  with np.load(P/path)as z:return {k:z[k]for k in z.files}
 fr=arrays('data/enumerations/fresh_h22_new_cusp_completion_20261013_frames.npz')
 hi=arrays('data/enumerations/fresh_h22_new_cusp_completion_20261013_highest_roots.npz')
 nc=arrays('data/enumerations/fresh_h22_cusp_orbits_20261013.npz')
 rr=arrays('data/enumerations/fresh_h22_gap_face_orbit_20261013.npz')
 delta=[K(F(x))for x in read(LOCAL)['representatives'][0]['null_ray']]
 def normal(ref):
  name,i,s=ref
  if name=='delta':return delta
  if name=='old':return [K(x)for x in W[i]]+[K()]
  if name=='old_null':return [K(int(x))for x in nc['old_null_projections'][i]]+[K(0,F(s*int(nc['four_heights'][i]),44))]
  if name=='old_k6':return [K(F(int(x),6))for x in rr['six_old_projections'][i]]+[K(s)]
  if name=='center':return [K(int(x))for x in h]+[K(0,F(4*s,11))]
  if name=='highest_roots':
   a=hi['old_bases'][i];j=int(hi['cusp_indices'][i]);size=int(hi['component_sizes'][i]);norm=int(a@q@a)
   t=60 if size>=4 else 80 if size==3 else 20 if norm==2 else 40
   old=[F(int(x))+F(t*int(y),120)for x,y in zip(a,fr['six_cusp_old'][j],strict=True)]
  else:
   j=i
   if name=='first+':a=fr['first_tangent_numerator'][i];den=28;t=200;sign=1
   elif name=='first-':a=fr['first_tangent_numerator'][i];den=28;t=100;sign=-1
   elif name=='second+':a=fr['second_tangent_numerator'][i];den=5208;t=1;sign=1
   elif name=='second-':a=fr['second_tangent_numerator'][i];den=5208;t=1;sign=-1
   else:raise ValueError(name)
   old=[F(sign*int(x),den)+F(t*int(y),120)for x,y in zip(a,fr['six_cusp_old'][j],strict=True)]
  return [K(x)for x in old]+[K(F(s*t,10))]
 return normal

def solve_exact(columns,target):
 n=len(columns);M=[[columns[j][i]for j in range(n)]+[target[i]]for i in range(23)];row=0;pivots=[]
 for j in range(n):
  k=next((k for k in range(row,23)if M[k][j]!=K()),None)
  if k is None:continue
  M[row],M[k]=M[k],M[row];a=M[row][j];M[row]=[x/a for x in M[row]]
  for k in range(23):
   if k!=row and M[k][j]!=K():
    a=M[k][j];M[k]=[x-a*y for x,y in zip(M[k],M[row],strict=True)]
  pivots.append(j);row+=1
 assert len(pivots)==n
 assert all(all(x==K()for x in r)for r in M[row:])
 sol=[K()for _ in range(n)]
 for i,j in enumerate(pivots):sol[j]=M[i][-1]
 assert all(sum((sol[j]*columns[j][i]for j in range(n)),K())==target[i]for i in range(23))
 return sol

def run(orthogonal=False):
 normal=sources();d=read(SCOUT)
 refs=sorted({tuple(x)for r in d['rows']for x in r['solution']['active_caps']if x[0]not in('old_timelike','center')}|{('highest_roots',4,1)})
 refs=[('old',i,1)for i in range(210)]+refs+[('delta',0,1)]
 ns=[normal(ref)for ref in refs];ex=read(ESCAPE)
 target=[K(F(x,ex['denominator']))for x in ex['numerator'][:-1]]+[K(0,F(ex['numerator'][-1],11*ex['denominator']))]
 val=lambda x:float(x.a)+np.sqrt(11)*float(x.b)
 M=np.array([[val(x)for x in n]for n in ns]).T;y=np.array(list(map(val,target)))
 free_old=set(read(LOCAL)['representatives'][0]['incident_old_wall_indices'])if orthogonal else set()
 costs=np.r_[np.ones(len(ns)-1),0];bounds=[(0,None)]*(len(ns)-1)+[(None,None)]
 for i in free_old:costs[i]=0;bounds[i]=(None,None)
 scales=np.linalg.norm(M,axis=0);lp=linprog(costs,A_eq=M/scales,b_eq=y,bounds=bounds,method='highs')
 print('LP',lp.success,lp.message,flush=True)
 if not lp.success:
  return {'status':'NUMERICAL_SELECTED_CONIC_SUPPORT_INFEASIBLE','normals_tested':len(ns),'scope':'A failed selected-support LP, not a nonexistence certificate.'}
 ids=np.flatnonzero(abs(lp.x)>1e-8);columns=[ns[i]for i in ids];sol=solve_exact(columns,target)
 assert all(a.sign()>0 for i,a in zip(ids,sol,strict=True)if refs[i][0]!='delta'and i not in free_old)
 terms=[{'source':list(refs[i]),'coefficient':a.json()}for i,a in zip(ids,sol,strict=True)]
 print('exact terms',len(terms),terms,flush=True)
 return {'status':'PASS_EXACT_'+('ORTHOGONAL_INCIDENT'if orthogonal else 'ALL_INCIDENT')+'_ACUTE_NORMALS_MISS_FIXED_ESCAPE','sources':{p:sha(p)for p in(RAW,LOCAL,ESCAPE,SCOUT,NULLS,LOWER)},'free_incident_old_coefficients':sorted(free_old),'terms':terms,'scope':'Exact identity X=sum positive coefficients times retained current wall normals plus a real multiple of delta4 and unrestricted multiples of the stated free old walls. Consequently any new normal r perpendicular to delta4 and every stated free old wall, and nonpositive against the remaining support walls, has r.X<=0. This concerns walls incident at this one cusp, not arbitrary additions or other cusps. No center wall is used, so the conclusion is uniform in the old center label m>=3.'}

def verify_saved(path):
 d=read(path);assert d['status']=='PASS_EXACT_ORTHOGONAL_INCIDENT_ACUTE_NORMALS_MISS_FIXED_ESCAPE'
 assert all(sha(p)==v for p,v in d['sources'].items())
 normal=sources();loc=read(LOCAL)['representatives'][0]
 assert d['free_incident_old_coefficients']==loc['incident_old_wall_indices']
 assert {tuple(t['source'])for t in d['terms']if t['source'][0]not in('old','delta')}=={('old_null',95611,1),('old_null',95665,1)}
 ex=read(ESCAPE);target=[K(F(x,ex['denominator']))for x in ex['numerator'][:-1]]+[K(0,F(ex['numerator'][-1],11*ex['denominator']))]
 total=[K()for _ in range(23)];positive=0
 for term in d['terms']:
  ref=term['source'];a=K(*map(F,term['coefficient']))
  if ref[0]!='delta'and not(ref[0]=='old'and ref[1]in d['free_incident_old_coefficients']):
   assert a.sign()>0;positive+=1
  r=normal(ref);total=[x+a*y for x,y in zip(total,r,strict=True)]
 assert positive==5 and total==target
 return d

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');p.add_argument('--orthogonal',action='store_true');a=p.parse_args()
 path='results/'+PREFIX+('_orthogonal'if a.orthogonal else '')+'.json'
 if a.verify:out=verify_saved(path)
 else:
  out=run(a.orthogonal);(P/path).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
