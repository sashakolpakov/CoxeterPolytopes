#!/usr/bin/env python3
"""Numerical outgoing-edge diagnostic from all23 saved actual finite vertices."""
import argparse,hashlib,json,itertools
from pathlib import Path
import numpy as np
from collections import Counter
import fresh_h22_incidence_prism_20261013_verify as V
P=V.P;TAG='fresh_h22_pi7_edgewalk_20261013';OUT=P/'results'/f'{TAG}_scout.json'
def load():
 M,R,h=V.setup();cert=V.read(V.OUT);geo=V.read('results/'+V.PRE+'_geometry.json');recipes=V.read(V.REC)['recipes']
 for d in(cert,geo):
  for p,v in d['sources'].items():assert V.sha(p)==v
 cap=[]
 for r in recipes:
  u=list(map(V.F,r['u']));v=list(map(V.F,r['v']));q=V.F(r['norm']);qu=V.pair(M,u,u);qv=V.pair(M,v,v);rad=(q/qv,V.F(0),-q*q*qu/(4*qv));root=V.scale(V.sqrt(V.ki(rad)),r['radical_branch']);p=V.scale(V.T,-q/2);cap.append([V.add(V.scale(p,a),V.scale(root,b))for a,b in zip(u,v)])
 normals=np.vstack([np.array(R,dtype=float),np.array([[(a+b)/(2*V.S)for a,b in r]for r in cap])]);M=np.array(M,dtype=float);G=normals@M@normals.T;n=np.diag(G)
 return M,normals,np.array(h,dtype=float),G,n,geo,cert

def run():
 M,W,h,G,n,geo,cert=load();bound=np.sqrt(n[:,None]*n[None,:]);finite=G> -bound+1e-7;non_dotted=G>=-bound-1e-7;np.fill_diagonal(finite,False);np.fill_diagonal(non_dotted,False);rows=[];status=Counter()
 for vi,v in enumerate(geo['new_finite_vertices']):
  incident=v['base_walls']+[812+v['new_cap_index']];A=G[np.ix_(incident,incident)];lam=np.linalg.solve(A,W[incident]@M@h);start=h-lam@W[incident];assert start@M@start<0
  for released in incident:
   edge=[i for i in incident if i!=released];C=G[np.ix_(edge,edge)];CI=np.linalg.inv(C);candidates=np.where(finite[edge].all(axis=0))[0].tolist();ends=[];near=[]
   for j in candidates:
    pp=G[edge,j];schur=n[j]-pp@CI@pp
    if schur>1e-6:ends.append({'wall':j,'schur':float(schur)})
    elif schur>-1e-6:near.append({'wall':j,'schur':float(schur)})
   assert released in[j['wall']for j in ends];other=[j for j in ends if j['wall']!=released];assert len(other)<=1
   # Edge tangent perpendicular to start; choose the released wall's interior.
   r=W[released];proj=r-(G[released,edge]@CI)@W[edge];t=proj-(proj@M@start)/(start@M@start)*start
   if t@M@r>0:t=-t
   t=t*np.sqrt(-(start@M@start)/(t@M@t));null=start+t;assert abs(null@M@null)<1e-6*max(1.,abs(start@M@start))
   ps=W@M@null;ii=int(np.argmax(ps));pp=float(ps[ii]);active=np.where(abs(ps)<1e-6)[0].tolist()
   if other:kind='EXISTING_FINITE_ENDPOINT'
   elif pp>1e-6:kind='UNRESOLVED_NUMERICAL_ENDPOINT'
   elif near or len(active)>21:kind='POSSIBLE_IDEAL_ENDPOINT_NEEDS_LINK'
   else:kind='LIGHT_CONE_ESCAPE_CANDIDATE'
   row={'vertex':vi,'released':released,'edge':edge,'finite_pair_candidates':candidates,'finite_endpoints':ends,'near_zero_schur_candidates':near,'status':kind,'null_point':null.tolist(),'maximum_null_pairing':pp,'maximum_null_wall':ii,'near_null_active_walls':active}
   rows.append(row);status[kind]+=1
 print('status',dict(status),flush=True)
 out={'status':'NUMERICAL_EDGE_FRONTIER_SCOUT','sources':{p:V.sha(p)for p in(V.OUT,V.REC,'results/'+V.PRE+'_geometry.json','scripts/'+Path(__file__).name)},'facets':835,'known_actual_vertices':23,'directed_edges':len(rows),'counts':dict(status),'rows':rows,'scope':'Numerical diagnostic only. Every one of22 releases at each of23 saved vertices is screened against all835 existing walls. Finite endpoint proposals use Schur complements. Light-cone boundary signs distinguish escaping rays from possible ideal endpoints; cusp feasibility/link conclusions require exact checks.'}
 OUT.write_text(json.dumps(out,indent=2)+'\n');print(OUT,flush=True)
if __name__=='__main__':run()
