#!/usr/bin/env python3
"""Numerical ranking by inherited finite-component contacts."""
import json,itertools
from pathlib import Path
import numpy as np
import networkx as nx
ROOT=Path(__file__).resolve().parents[1];OLD='data/enumerations/fresh_h9_component_bridge_20261013';BASE='data/enumerations/fresh_h9_remote_frontier_20261013'
d=json.loads((ROOT/(OLD+'_faces_scout.json')).read_text());c=json.loads((ROOT/(OLD+'_verify.json')).read_text());s=json.loads((ROOT/(BASE+'_scout.json')).read_text())
G=np.array(d['gram']);R=np.vstack([d['normals'],[[sum(map(int,x))/2e80 for x in r]for r in c['normal_enclosures']]]);names=d['names']+[r['name']for r in c['recipes']];PP=R@G@R.T
faces=[tuple(z)for z in d['spherical8_faces']]+[tuple(names.index(n)for n in r['edge'])for r in c['new_vertex_frontier']];faces=list({tuple(sorted(x))for x in faces});verts=set();fends={}
for face in faces:
 Mi=np.linalg.inv(PP[np.ix_(face,face)]);ends=[]
 for j in range(104):
  if j in face:continue
  p=PP[j,list(face)]
  if p.min()>-2+1e-6 and 2-p@Mi@p>1e-6:ends.append(j);verts.add(tuple(sorted((*face,j))))
 fends[face]=ends
vg=nx.Graph();vg.add_nodes_from(verts)
for face,ends in fends.items():
 if len(ends)==2:vg.add_edge(*(tuple(sorted((*face,j)))for j in ends))
comps=list(nx.connected_components(vg));vc={v:i for i,vs in enumerate(comps)for v in vs};out=[]
for i,row in enumerate(s['compatible_proposals']):
 r=np.array(row['normal']);pp=R@G@r;ids=list(map(int,np.flatnonzero(pp>-2+1e-6)));nv=[];cc=set()
 for face in itertools.combinations(ids,8):
  M=PP[np.ix_(face,face)];eig=np.linalg.eigvalsh(M)
  if eig[0]<1e-7:continue
  p=pp[list(face)]
  if 2-p@np.linalg.solve(M,p)>1e-7:
   nv.append(list(face))
   for j in fends.get(face,[]):cc.add(vc[tuple(sorted((*face,j)))])
 out.append({'candidate':i,'new_vertices':nv,'components':sorted(cc),'target_cut':row['target_cut'],'anchor_label':row['anchor_label']})
out.sort(key=lambda z:(-len(z['components']),-len(z['new_vertices']),z['anchor_label']not in[2,3,5],-z['target_cut']))
(ROOT/(BASE+'_rank.json')).write_text(json.dumps({'old_components':[list(map(list,x))for x in comps],'ranking':out,'scope':'Numerical prioritization only; exact chamber/component certificates needed.'},indent=2)+'\n')
print('components',[len(x)for x in comps])
for r in out[:20]:print(r)
