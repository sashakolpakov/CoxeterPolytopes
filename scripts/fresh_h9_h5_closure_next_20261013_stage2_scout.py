#!/usr/bin/env python3
"""New finite contacts to one gap-cutting endpoint with only new labels2 and3."""
import json,itertools,copy
from collections import defaultdict
import numpy as np
import fresh_h9_h5_closure_next_20261013_scout as B
BASE=B.TAG;TAG=BASE+'_stage2';original_load=B.load
def load():
 R,G,c,names,doc,p,C,J,RF,phys,target=original_load();bank=json.loads((B.P/'data/enumerations'/f'{BASE}_scout.json').read_text());selected=[]
 for k,i in enumerate((7308,)):
  rec=copy.deepcopy(bank['candidates'][i]);rec['name']=f'endpoint{k}';selected.append(rec);names.append(rec['name']);R=np.vstack([R,np.array(rec['normal'])])
 RF=np.c_[-R@G@p,R@G@C];vs=B.S.geometry(R,G,p,C);fixed=doc['H5_face']['fixed_walls'];edges=defaultdict(list);vertices=[]
 for v in vs:
  if v['norm']>=-1e-7:continue
  inc=[names[i]for i in v['active']if names[i]not in fixed];assert len(inc)==5,(inc,v)
  vertices.append(inc)
  for e in itertools.combinations(inc,4):edges[tuple(sorted(e))].append(len(vertices)-1)
 assert all(len(z)in(1,2)for z in edges.values());leading=max(vs,key=lambda v:v['norm']);target=np.array(leading['point']);doc=copy.deepcopy(doc);doc['H5_face']['edges']=[{'walls':list(e),'endpoints':inds}for e,inds in sorted(edges.items())]
 # The generic scout iterates this field only for allowed metric anchors.
 # Here anchors must be one of the newly retained three normals.
 doc['H5_face']['facets']=[r['name']for r in selected]
 state={'source':B.OLD,'source_sha256':B.S.B.V.sha(B.OLD),'selected':selected,'retained_facets':len(R),'finite_vertices':len(vertices),'finite_edges':sum(len(e)==2 for e in edges.values()),'open_edges':sum(len(e)==1 for e in edges.values()),'maximum_norm':leading['norm'],'target_point':target.tolist(),'target_support':[names[i]for i in leading['active']if names[i]not in fixed]}
 (B.P/'data/enumerations'/f'{BASE}_state2.json').write_text(json.dumps(state,indent=2)+'\n')
 return R,G,c,names,doc,p,C,J,RF,phys,target
if __name__=='__main__':
 B.load=load;B.TAG=TAG;B.run();path=B.P/'data/enumerations'/f'{TAG}_scout.json';out=json.loads(path.read_text());out['retained_facets']=124;out['required_metric_anchors']=['endpoint0'];out['sources']['data/enumerations/'+BASE+'_state2.json']=B.S.B.V.sha('data/enumerations/'+BASE+'_state2.json');out['sources']['scripts/'+BASE+'_scout.py']=B.S.B.V.sha('scripts/'+BASE+'_scout.py');out['scope']='Bounded nonright spherical endpoint profiles at the124-wall numerical state. Every new candidate must have a finite metric contact to endpoint0. Strict dotted other walls; numerical affine gate. Exact reconstruction remains required.';path.write_text(json.dumps(out,indent=2)+'\n')
