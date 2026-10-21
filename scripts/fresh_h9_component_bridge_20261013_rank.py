#!/usr/bin/env python3
"""Rank fresh trivalent caps by actual finite-edge contacts (numeric scout)."""
import json,itertools
from pathlib import Path
from collections import Counter
import numpy as np
import networkx as nx
P=Path(__file__).resolve().parents[1];PRE='data/enumerations/fresh_h9_component_bridge_20261013'
d=json.loads((P/(PRE+'_faces_scout.json')).read_text());s=json.loads((P/(PRE+'_trivalent_scout.json')).read_text())
G=np.array(d['gram']);R=np.array(d['normals']);PP=R@G@R.T;faces=[tuple(x)for x in d['spherical8_faces']];names=d['names'];verts=set();fends=[]
for face in faces:
 js=[];A=R[list(face)];M=A@G@A.T;Mi=np.linalg.inv(M)
 for j in range(len(R)):
  if j in face:continue
  p=PP[j,list(face)]
  if min(p)<-2+1e-6 or 2-p@Mi@p<=1e-6:continue
  js.append(j);verts.add(tuple(sorted((*face,j))))
 fends.append(js)
vlist=sorted(verts);vg=nx.Graph();vg.add_nodes_from(range(len(vlist)));vind={v:i for i,v in enumerate(vlist)}
for face,js in zip(faces,fends):
 if len(js)==2:vg.add_edge(*(vind[tuple(sorted((*face,j)))]for j in js))
comps=list(nx.connected_components(vg));vc={v:i for i,c in enumerate(comps)for v in c}
print('old finite vertices',len(verts),'components',[len(c)for c in comps],'faces endpointcounts',Counter(map(len,fends)),flush=True)
rr=[];phi=(1+np.sqrt(5))/2
for ci,row in enumerate(s['candidates']):
 if row['target_cuts'][1]<=1e-6 or row['conflicts']:continue
 r=np.array(row['normal']);pp=R@G@r;finite=np.flatnonzero(pp>-2+1e-6);vals=np.array([0.,-1.,-phi]);lab_ok=np.min(abs(pp[finite,None]-vals),axis=1).max()<1e-6
 newverts=[];components=set();closed=0
 for fi,(face,js)in enumerate(zip(faces,fends)):
  if not set(face)<=set(finite):continue
  p=pp[list(face)];M=PP[np.ix_(face,face)]
  if 2-p@np.linalg.solve(M,p)<=1e-6:continue
  newverts.append(fi)
  if len(js)==1:
   closed+=1;components.add(vc[vind[tuple(sorted((*face,js[0])))]] )
 rr.append({'candidate':ci,'cut':row['target_cuts'],'anchor':row['anchor'],'anchor_label':row['anchor_label'],'finite_neighbors':len(finite),'golden_finite_labels_only':bool(lab_ok),'new_vertices':newverts,'closed_old_edges':closed,'components_touched':sorted(components),'edge':[names[j]for j in row['edge']],'conflicts':row['conflicts']})
rr.sort(key=lambda r:(-len(r['components_touched']),-r['closed_old_edges'],not r['golden_finite_labels_only'],-r['cut'][1]))
out={'old_vertices':vlist,'old_components':[sorted(c)for c in comps],'face_existing_end_walls':fends,'candidates':rr,'scope':'Numeric geometric prioritization; no exact closure theorem.'}
(P/(PRE+'_rank_scout.json')).write_text(json.dumps(out,indent=2)+'\n')
for r in rr[:20]:print(r,flush=True)
