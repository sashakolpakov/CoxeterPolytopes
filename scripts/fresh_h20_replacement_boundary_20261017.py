#!/usr/bin/env python3
"""Numerical complete T-H3 boundary ledger from frozen exact171 coordinates.

Plane intersections outside the ball are separated from actual hyperbolic
vertices/edges. No exact arithmetic or ambient-volume conclusion is added.
"""
import argparse,hashlib,itertools,json,math
from collections import Counter
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog
import networkx as nx
import fresh_h20_opposite_prism_cap_20261017 as P
ROOT=P.ROOT;TAG='fresh_h20_replacement_boundary_20261017';INPUT='results/fresh_h20_replacement_hexagon_20261017.json'
TOL=1e-8

def geometry():
 G,Y,CT,J,q,V,C=P.geometry();D=P.P.load(INPUT);B=Y[[95,142,145,31]].copy()
 mid=lambda z:float((Fraction(z[0])+Fraction(z[1]))/2)
 for z in D['exact_definitions']:
  assert z['root']==len(Y)
  Y=np.r_[Y,[np.array([mid(x)for x in z['residual_coordinates']])@B]]
  CT=np.r_[CT,[[mid(x)for x in z['T_target_products']]]];q=np.r_[q,2.]
 errs=np.abs(np.einsum('ij,jk,ik->i',Y,G['m'],Y)+np.einsum('ij,jk,ik->i',CT,J,CT)-q);scaled=errs/(1+np.sum(Y*Y,axis=1)+np.sum(CT*CT,axis=1));assert max(scaled)<1e-7;G['norm_evaluation_diagnostic']=dict(maximum_absolute_cancellation=float(max(errs)),wall=int(np.argmax(errs)),maximum_scaled_residual=float(max(scaled)),note='Double-precision cancellation on large inherited physical coordinates; this diagnostic is not an exact norm certificate.')
 return G,Y,CT,J,q,D

def boundary(Y,m,ret,T):
 norm=np.einsum('ij,jk,ik->i',Y,m,Y);sides=[i for i in ret if norm[i]>TOL];negative=[i for i in ret if i not in T and norm[i]<-TOL];uncertain=[i for i in ret if i not in T and abs(norm[i])<=TOL]
 assert not uncertain,uncertain
 A=Y[sides,1:];b=Y[sides,0];sc=np.linalg.norm(A,axis=1)
 lp=linprog([0,0,0,-1],A_ub=np.c_[A/sc[:,None],np.ones(len(sides))],b_ub=b/sc,bounds=[(None,None)]*3+[(0,None)],method='highs');assert lp.success and lp.x[3]>1e-9
 bounds=[]
 for ax in range(3):
  for sign in(-1,1):
   obj=np.zeros(3);obj[ax]=sign;r=linprog(obj,A_ub=A,b_ub=b,bounds=[(None,None)]*3,method='highs');assert r.success
   bounds.append(dict(axis=ax,sign=sign,objective=float(r.fun)))
 verts=[];trials=Counter()
 for ids in itertools.combinations(sides,3):
  aa=Y[list(ids),1:];bb=Y[list(ids),0]
  if abs(np.linalg.det(aa))<1e-11:trials['singular']+=1;continue
  x=np.linalg.solve(aa,bb)
  if max((A@x-b)/(1+sc))>1e-8:trials['infeasible']+=1;continue
  trials['feasible']+=1
  if any(np.linalg.norm(x-z['point'])<1e-7 for z in verts):continue
  rr=float(x@x);active=[j for j in sides if abs(Y[j,1:]@x-Y[j,0])/(1+np.linalg.norm(Y[j]))<1e-7]
  verts.append(dict(point=x.tolist(),radius_squared=rr,active=active,kind='finite'if rr<1-1e-6 else'exterior'if rr>1+1e-6 else'near_sphere'))
 verts.sort(key=lambda z:(z['active'],z['point']));assert verts
 # Verify each reported plane vertex against every positive-side inequality.
 edges=[];open_edges=[];ideal=[];skeleton=nx.Graph();skeleton.add_nodes_from(i for i,v in enumerate(verts)if v['kind']=='exterior')
 for i,j in itertools.combinations(sides,2):
  U=Y[[i,j]];aa=U[:,1:];bb=U[:,0];dd=np.cross(*aa)
  if np.linalg.norm(dd)<1e-11:continue
  d=dd/np.linalg.norm(dd);x=aa.T@np.linalg.solve(aa@aa.T,bb);lo,hi=-math.inf,math.inf;failed=False
  for k in sides:
   a=float(Y[k,1:]@d);v=float(Y[k,0]-Y[k,1:]@x)
   if a>1e-10:hi=min(hi,v/a)
   elif a< -1e-10:lo=max(lo,v/a)
   elif v< -1e-8:failed=True;break
  if failed or hi-lo<1e-8:continue
  vv=[k for k,z in enumerate(verts)if i in z['active']and j in z['active']]
  assert len(vv)==2,(i,j,vv)
  mn=max(lo,min(hi,0.));edge_min=float((x+mn*d)@(x+mn*d))
  if all(verts[k]['kind']=='exterior'for k in vv)and edge_min>1+1e-7:skeleton.add_edge(*vv)
  gram=U@m@U.T;eig=np.linalg.eigvalsh(gram);rr=1-float(x@x)
  rec=dict(pair=[i,j],poly_vertex_indices=vv,poly_parameter_bounds=[lo,hi],minimum_squared_radius=edge_min,projected_gram=gram.tolist(),actual_hyperbolic_edge=False)
  if min(eig)>1e-9 and rr>1e-9:
   sphere=math.sqrt(rr);l,h=max(lo,-sphere),min(hi,sphere)
   if h-l>1e-7:
    ends=[]
    for t in(l,h):
     xx=x+t*d;kind='ideal'if abs(abs(t)-sphere)<1e-7 else'finite';act=[k for k in sides if abs(Y[k,1:]@xx-Y[k,0])/(1+np.linalg.norm(Y[k]))<1e-7]
     assert max((Y[ret,1:]@xx-Y[ret,0])/(1+np.linalg.norm(Y[ret],axis=1)))<1e-7
     e=dict(kind=kind,point=xx.tolist(),active=act)
     if kind=='ideal':
      matches=[k for k,z in enumerate(ideal)if np.linalg.norm(xx-z['point'])<1e-7]
      if matches:idx=matches[0];ideal[idx]['incident_edges'].append([i,j])
      else:idx=len(ideal);ideal.append(dict(point=xx.tolist(),active=act,incident_edges=[[i,j]]))
      e['ideal_index']=idx
     else:
      matches=[k for k,z in enumerate(verts)if np.linalg.norm(xx-z['point'])<1e-7];assert len(matches)==1;e['finite_vertex_index']=matches[0]
     ends.append(e)
    rec.update(actual_hyperbolic_edge=True,hyperbolic_endpoints=ends,finite_endpoint_count=sum(z['kind']=='finite'for z in ends))
    if rec['finite_endpoint_count']<2:open_edges.append(rec)
  edges.append(rec)
 faces=[]
 for i in sides:
  vs=[k for k,z in enumerate(verts)if i in z['active']];ee=[z for z in edges if i in z['pair']and z['actual_hyperbolic_edge']]
  faces.append(dict(wall=i,projected_norm=float(norm[i]),poly_vertex_indices=vs,poly_vertex_kinds=dict(Counter(verts[k]['kind']for k in vs)),actual_edges=[z['pair']for z in ee],open_edges=[z['pair']for z in ee if z['finite_endpoint_count']<2],compact=all(verts[k]['kind']=='finite'for k in vs)))
 fg=nx.Graph();fg.add_edges_from(z['pair']for z in open_edges)
 # Exterior-skeleton components are intentionally only a combinatorial
 # diagnostic, not a claimed complete ideal-patch topology theorem.
 comps=[]
 for comp in sorted(nx.connected_components(skeleton),key=lambda c:sorted(c)):
  inc=[z['pair']for z in open_edges if set(z['poly_vertex_indices'])&comp]
  comps.append(dict(exterior_vertex_indices=sorted(comp),incident_open_edges=inc,maximum_squared_radius=max(verts[k]['radius_squared']for k in comp)))
 return dict(actual_facet_ids=sides,negative_projected_ids=negative,uncertain_projected_ids=uncertain,chart_center=lp.x[:3].tolist(),bounded_coordinate_LPs=bounds,triple_trials=dict(trials),poly_vertices=verts,poly_vertex_counts=dict(Counter(z['kind']for z in verts)),poly_edges=edges,actual_edge_count=sum(z['actual_hyperbolic_edge']for z in edges),open_edges=open_edges,ideal_edge_endpoints=ideal,facets=faces,open_facet_incidence_components=[sorted(c)for c in nx.connected_components(fg)],exterior_skeleton_components=comps,maximum_squared_radius=max(z['radius_squared']for z in verts))

def run():
 G,Y,CT,J,q,D=geometry();T=G['T'];stages=[('original168',list(range(168))),('delete158_only',[i for i in range(168)if i!=158]),('replace158_with168',[i for i in range(169)if i!=158]),('also_delete31',[i for i in range(169)if i not in(31,158)]),('add_nonright169',[i for i in range(170)if i not in(31,158)]),('full171',D['actual_wall_labels'])]
 out=[]
 for name,ret in stages:
  row=boundary(Y,G['m'],ret,T);out.append(dict(stage=name,actual_ambient_wall_count=len(ret),retained_ids=ret,**row));print(name,row['poly_vertex_counts'],'actual edges',row['actual_edge_count'],'open',len(row['open_edges']),'components',row['open_facet_incidence_components'],flush=True)
 current=out[-1];hexface=next(z for z in current['facets']if z['wall']==169);assert hexface['compact']and len(hexface['poly_vertex_indices'])==6
 paths=[f'scripts/{TAG}.py',INPUT,f'scripts/{P.TAG}.py','scripts/fresh_h20_e8_partner_endpoints_20261017.py','results/fresh_h20_e8_cusp_completion_20261017.json','results/fresh_h20_e8_f4_cusps_20261017.json']
 return dict(status='NUMERICAL_COMPLETE_T_H3_BOUNDARY_LEDGER_FOR_EXACT171',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},support=T,geometry=dict(norm_evaluation_diagnostic=G['norm_evaluation_diagnostic'],residual_form=G['m'].tolist(),native_norms=q.tolist(),residual_normals=Y.tolist(),T_products=CT.tolist(),T_inverse=J.tolist()),stages=out,scope='Complete numerical plane/edge intersections in one fixed target-H3 Klein chart, from midpoint evaluation of frozen exact171 residual definitions. Actual hyperbolic edges are clipped by every retained halfspace and the unit ball; exterior chart vertices are not hyperbolic vertices. Connected components of the open facet-incidence graph and of the exterior polyhedral skeleton are separate diagnostics, not a proof of ideal-boundary topology. No new exact arithmetic, ambient vertex census, or finite-volume claim.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(out.read_text()))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
