#!/usr/bin/env python3
"""Rank bounded endpoint proposals by actual-edge closure and compactness gates."""
import itertools,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scripts.fresh_h9_face_bridge_20261013_scout import setup,PREFIX,V,C

def run():
 Q,R,names,h,center,y,prior=setup();d=V.read(PREFIX+'_scout.json');M=R@Q@R.T
 edges=[e for e in prior['new_vertex_edge_frontier']['edges']if e['status']=='open_to_light_cone'];ids=[[names.index(n)for n in e['common_facets']]for e in edges]
 inv=[np.linalg.inv(M[np.ix_(ss,ss)])for ss in ids];T=prior['compact_H3_face']['common_facets'];out=[]
 for ci,c in enumerate(d['candidates']):
  if c['parallel_walls']:continue
  r=np.array(c['normal']);p=R@Q@r;closed=[i for i,(ss,mi)in enumerate(zip(ids,inv))if max(p[ss])<1e-6 and min(p[ss])>-2+1e-6 and p[ss]@mi@p[ss]<2-1e-7]
  prisms=[a for a in T if all(i in closed for i,e in enumerate(edges)if a not in e['common_facets'])]
  if c['target_cut']<=1e-8 and not prisms and len(closed)<4:continue
  ns=[names.index(n)for n in c['finite_profile']];aff=[]
  # Every PSD diagram containing r is contained in this finite-neighbor set.
  # The test is numerical ranking only; successful selections need exact replay.
  for k in range(1,len(ns)+1):
   for ss in itertools.combinations(ns,k):
    if any(M[a,b]<-2+1e-7 for a,b in itertools.combinations(ss,2)):continue
    A=np.block([[M[np.ix_(ss,ss)],p[list(ss),None]],[p[None,list(ss)],np.array([[2.]])]])
    ev=np.linalg.eigvalsh(A)
    if abs(ev[0])<1e-6:aff.append([names[a]for a in ss])
   if aff:break
  out.append({'candidate':ci,'closed_edges':closed,'closed_edge_count':len(closed),'compact_prism_released_walls':prisms,'target_cut':c['target_cut'],'affine_obstructions':aff})
 out.sort(key=lambda r:(not bool(r['affine_obstructions']),bool(r['compact_prism_released_walls']),r['target_cut']>0,r['closed_edge_count'],r['target_cut']),reverse=True)
 result={'status':'NUMERICAL_ENDPOINT_QUALITY_AND_AFFINE_GATE','source':PREFIX+'_scout.json','source_sha256':V.sha(PREFIX+'_scout.json'),'ranked':out}
 (V.ROOT/(PREFIX+'_rank.json')).write_text(json.dumps(result,indent=2)+'\n')
 print('RANKED',len(out),'AFFINE_FREE',sum(not r['affine_obstructions']for r in out),'PRISMS',sum(bool(r['compact_prism_released_walls'])and not r['affine_obstructions']for r in out))
 for r in out[:15]:print(r)

if __name__=='__main__':run()
