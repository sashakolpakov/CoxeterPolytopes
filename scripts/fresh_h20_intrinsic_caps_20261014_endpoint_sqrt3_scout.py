#!/usr/bin/env python3
"""Numerical coherent R4/(2sqrt3phi) endpoint replacement covering old and3/4 probes."""
import argparse,json,time,os
from collections import Counter
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix,vstack
from fresh_h20_intrinsic_caps_20261014 import ROOT,TAG,data,sha
from fresh_h20_intrinsic_caps_20261014_central_contact_scout import allowed

def run(seconds):
 t=time.time();phi=(1+5**.5)/2;D=data(phi);V=np.array([v['point']for v in D['vertices']]);source=f'results/{TAG}_annulus_endpoint_all_scout.json';pfile=f'results/{TAG}_endpoint.json';s=json.loads((ROOT/source).read_text());P=json.loads((ROOT/pfile).read_text())
 rows=[]
 for r in s['surviving_candidates']:
  if min(abs(r['radius']-4),abs(r['radius']-2*3**.5*phi))>1e-8:continue
  tr,pi,_,_=r['recipe'];typ=P['literal_edge_transports'][tr]['type'];pr=next((a for a in P['types'][typ]['profiles']if a['id']==pi),None)
  if pr is None or any(m not in[3,4,5,6]for w,m in pr['contacts']):continue
  rows.append(r)
 C=np.concatenate([V,np.array([r['center']for r in rows])]);rr=np.concatenate([np.ones(252),np.array([r['radius']for r in rows])]);n=len(C)
 probes=np.concatenate([V*.75,V]);dist=np.maximum(np.sum(C*C,axis=1)[:,None]+np.sum(probes*probes,axis=1)[None,:]-2*C@probes.T,0);cover=dist<rr[:,None]**2-1e-8
 dd=np.maximum(np.sum(C*C,axis=1)[:,None]+np.sum(C*C,axis=1)[None,:]-2*C@C.T,0);co=(dd-rr[:,None]**2-rr[None,:]**2)/(2*rr[:,None]*rr[None,:]);ok=allowed(co,list(np.cos(np.pi/np.arange(2,65)))+[1.]);np.fill_diagonal(ok,True);bad=np.argwhere(np.triu(~ok,1))
 print('PID',os.getpid(),'bank',n,'radiuscounts',dict(Counter(rr)),'conflicts',len(bad),'uncovered',int((~cover.any(axis=0)).sum()),flush=True)
 row,col=np.where(cover.T);A=coo_matrix((np.ones(len(row)),(row,col)),shape=(504,n)).tocsr();extra=coo_matrix((np.ones(2*len(bad)),(np.repeat(np.arange(len(bad)),2),bad.flatten())),shape=(len(bad),n)).tocsr();A=vstack([A,extra]).tocsr();lo=np.concatenate([np.ones(504),np.full(len(bad),-np.inf)]);hi=np.concatenate([np.full(504,np.inf),np.ones(len(bad))])
 objective=np.ones(n);objective[:252]=.001
 res=milp(objective,integrality=np.ones(n),bounds=Bounds(0,1),constraints=LinearConstraint(A,lo,hi),options={'time_limit':seconds,'mip_rel_gap':.01})
 selected=None;details={}
 if res.x is not None:
  selected=np.flatnonzero(res.x>.5).tolist();assert ok[np.ix_(selected,selected)].all()and cover[selected].any(axis=0).all()
  ss=np.array(selected);selP=(dist[ss]-rr[ss,None]**2).min(axis=0)
  details=dict(old_vertex_caps_retained=[i for i in selected if i<252],new_candidates=[rows[i-252]for i in selected if i>=252],selected_radius_counts=dict(Counter(rr[ss])),worst_probe_minimum_power=float(selP.max()))
 out=dict(status='NUMERICAL_COHERENT_ENDPOINT_REPLACEMENT_FOUND'if selected is not None else('NUMERICAL_R4_SQRT3_COVER_INFEASIBLE'if res.status==2 else'NUMERICAL_R4_SQRT3_COVER_UNRESOLVED'),sources={p:sha(p)for p in[source,pfile,f'scripts/{TAG}_endpoint_sqrt3_scout.py',f'scripts/{TAG}_quarter_face_scout.py']},fixed_interior_caps='central+252halfvertex+60quarter, with22seedwalls',bank_size=n,old_vertex_options=252,new_endpoint_options=len(rows),radius_counts={str(k):v for k,v in Counter(rr).items()},conflict_pairs=len(bad),probe_count=504,all_probes_covered_by_relaxation=bool(cover.any(axis=0).all()),solver_status=int(res.status),solver_message=res.message,selected=selected,selection=details,elapsed_seconds=time.time()-t,scope='Floating complete compatibility graph of this finiteR4/(2sqrt3phi) bank; original vertex caps optional. Covers only504 named probes, not a volume theorem. Any positive selection needs exact arithmetic, facets and endpoint checks.')
 path=ROOT/f'results/{TAG}_endpoint_sqrt3_scout.json';path.write_text(json.dumps(out,indent=2)+'\n');np.savez_compressed(ROOT/f'data/enumerations/{TAG}_endpoint_sqrt3_scout.npz',centers=C,radii=rr,cover=cover,compatible=ok)
 print(out['status'],res.message,'selected',len(selected)if selected is not None else None,'seconds',time.time()-t,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--seconds',type=float,default=60);a=p.parse_args();run(a.seconds)
