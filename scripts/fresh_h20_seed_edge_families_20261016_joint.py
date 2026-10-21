#!/usr/bin/env python3
"""On-demand finite cover with selected-facet incidence obligations."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from fresh_h20_continuous_sample_joint_cover_20261016 import bits
from fresh_h20_quad_adjacent_h3_closure_20261016 import geometry158,face_data,poly,gate_mask
from fresh_h20_release144_cap_20261016_residual_scout import close
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_seed_edge_families_20261016_joint'
SEED='results/fresh_h20_open_edge_closure_20261016.json';BASE='results/fresh_h20_seed_edge_families_20261016.json';PAIRS='results/fresh_h20_seed_edge_families_20261016_pairs.json'
def members(m):
 while m:
  v=m&-m;m-=v;yield v.bit_length()-1
class Search:
 def __init__(self,R,M,olddegree):self.R=R;self.M=M;self.degree=olddegree;self.cache={}
 def row(self,i):
  if i not in self.cache:
   c=-(self.R@(self.M@self.R[i]))/2;ok=gate_mask(c[:,None]);ok[i]=False
   self.cache[i]=(bits(ok),bits(ok&(abs(c)<1-1e-7)))
  return self.cache[i]
 def solve(self,cov):
  point=[bits(v)for v in cov];cover=[bits(v)for v in cov.T];failed=set();nodes=0
  def dfs(allowed,need,selected):
   nonlocal nodes
   nodes+=1
   if nodes>200000:raise RuntimeError('node_guard')
   key=(allowed,need,selected)
   if key in failed:return None
   choices=None
   for v in members(need):
    opts=point[v]&allowed
    if not opts:failed.add(key);return None
    if choices is None or opts.bit_count()<choices.bit_count():choices=opts
   for v in members(selected):
    finite=self.row(v)[1];deficit=3-int(self.degree[v])-(finite&selected).bit_count()
    if deficit>0:
     opts=finite&allowed
     if opts.bit_count()<deficit:failed.add(key);return None
     if choices is None or opts.bit_count()<choices.bit_count():choices=opts
   if choices is None:return list(members(selected))
   order=sorted(members(choices),key=lambda j:(-(cover[j]&need).bit_count(),-int(self.degree[j]),j))
   for j in order:
    ans=dfs(allowed&self.row(j)[0],need&~cover[j],selected|(1<<j))
    if ans is not None:return ans
   failed.add(key);return None
  try:ans=dfs((1<<len(self.R))-1,(1<<len(point))-1,0);complete=True
  except RuntimeError:ans=None;complete=False
  return ans,dict(nodes=nodes,failed=len(failed),complete=complete,cached_rows=len(self.cache),node_guard=200000)
def run():
 S=json.loads((ROOT/SEED).read_text());D=json.loads((ROOT/BASE).read_text());P=json.loads((ROOT/PAIRS).read_text());G=geometry158();M=G['metric']
 pool=[dict(origin='three_side',**r)for r in S['candidates']]+[dict(origin='endpoint_sample',**r)for r in D['candidates']]+[dict(origin='solved_partner',**r)for r in P['candidates']]
 R=np.array([r['root']for r in pool]);AR=-R[:,1:]/R[:,0,None]
 seed=np.array(S['seed']['roots']);X=np.r_[G['X'],seed];h=X@M@X[151];Xbar=X-h[:,None]*X[151]/2;q=np.einsum('ij,jk,ik->i',Xbar,M,Xbar);oi=np.where(q>1e-6)[0]
 olddegree=(2*q[oi][None,:]-(R@M@Xbar[oi].T)**2>1e-6).sum(axis=1)
 assert max(abs(R@M@X[151]))<1e-4
 search=Search(R,M,olddegree);base,V,A,b,ids=face_data(G);sa=-seed[:,1:]/seed[:,0,None];A=np.r_[A,sa@V];b=np.r_[b,1-sa@base];ids+=list(range(158,163))
 Z=[np.array(z)for z in D['exterior_points']];center=np.mean([v['z']for v in S['seed']['face']['vertices']if v['radius2']<1-1e-6],axis=0)
 steps=[];final=None;reason='iteration_guard'
 for k in range(30):
  cov=np.array(Z)@AR.T>1+1e-7
  if not cov.any(axis=1).all():reason='uncovered_sample';steps.append(dict(iteration=k,uncovered=np.where(~cov.any(axis=1))[0].tolist()));break
  selected,meta=search.solve(cov)
  if selected is None:reason='infeasible_bounded_incidence_cover'if meta['complete']else'node_guard';steps.append(dict(iteration=k,search=meta));break
  co=-R[selected]@M@R[selected].T/2;np.fill_diagonal(co,0);assert gate_mask(co).all()
  finite=abs(co)<1-1e-7;np.fill_diagonal(finite,False);assert min(olddegree[selected]+finite.sum(1))>=3
  aa=AR[selected];face=poly(base,V,np.r_[A,aa@V],np.r_[b,1-aa@base],ids+list(range(163,163+len(selected))))
  row=dict(iteration=k,selected=selected,roots=[pool[i]for i in selected],search=meta,face=face);steps.append(row)
  print('STEP',k,'caps',len(selected),'max',face['maximum_squared_radius'],'nodes',meta['nodes'],flush=True)
  if not face['exterior']:reason='compact_candidate'if not face['near_sphere']else'near_sphere';final=row;break
  worst=max(face['vertices'],key=lambda v:v['radius2']);d=np.array(worst['z'])-center
  t=next(float(t)for t in np.roots([d@d,2*center@d,center@center-1])if np.isreal(t)and 0<float(t)<1)
  z=center+t*d;assert abs(z@z-1)<1e-7 and max(aa@z-1)<1e-6;Z.append(z)
 src=[f'scripts/{TAG}.py',SEED,BASE,PAIRS,'scripts/fresh_h20_continuous_sample_joint_cover_20261016.py','scripts/fresh_h20_quad_adjacent_h3_closure_20261016.py']
 return dict(status='PASS_BOUNDED_EDGE_FAMILY_INCIDENCE_COVER',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},candidate_count=len(pool),termination=reason,steps=steps,final=final,samples=[z.tolist()for z in Z],scope='Numerical finite bank; selected facets require3 potential intrinsic finite neighbors. Exact full-endpoint closure is checked only by each full-face census.200000node and30iteration guards explicit; no continuous or arbitrary-label exclusion.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert close(d,json.loads(p.read_text()))
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['termination'],len(d['steps']),flush=True)
