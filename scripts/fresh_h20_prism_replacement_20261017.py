#!/usr/bin/env python3
"""Numerical explicit two-cap replacement of the obstructing old158 wall.

The normals have an exact core recipe. Old gates and the replacement whole
H3 cone are presently numerical; no exact ambient promotion is asserted.
"""
import argparse,hashlib,json,math
import numpy as np
import sympy as s
from scipy.optimize import linprog
from scipy.spatial import HalfspaceIntersection
import fresh_h20_opposite_prism_cap_20261017 as P
ROOT=P.ROOT;TAG='fresh_h20_prism_replacement_20261017'
def classify(x):
 if x>1+1e-7:return 'dotted'
 if abs(x)<1e-7:return 2
 if abs(x-1)<1e-7:return 'parallel'
 if x<0:return 'obtuse'
 if 0<x<1:
  lab=round(math.pi/math.acos(x))
  if lab>=3 and abs(x-math.cos(math.pi/lab))<1e-7:return int(lab)
 return 'angle_gap'
def poly(Y,indices):
 live=[i for i in indices if np.linalg.norm(Y[i,1:])>1e-9];A=Y[live,1:];b=Y[live,0];norm=np.linalg.norm(A,axis=1)
 lp=linprog([0,0,0,-1],A_ub=np.c_[A/norm[:,None],np.ones(len(live))],b_ub=b/norm,bounds=[(None,None)]*3+[(0,None)],method='highs')
 assert lp.success and lp.x[3]>1e-9
 raw=HalfspaceIntersection(np.c_[A,-b],lp.x[:3]).intersections;unique=[]
 for z in raw:
  if not any(np.linalg.norm(z-w)<1e-7 for w in unique):unique.append(z)
 v=np.array(unique);rr=np.sum(v*v,axis=1)
 out=[]
 for z,r in zip(v,rr):
  act=[j for j in live if abs(Y[j,1:]@z-Y[j,0])/(1+np.linalg.norm(Y[j]))<1e-7]
  out.append(dict(point=z.tolist(),radius_squared=float(r),active=act))
 return dict(raw_vertex_count=len(raw),vertices=out,finite=int(sum(rr<1-1e-6)),exterior=int(sum(rr>1+1e-6)),near_sphere=int(sum(abs(rr-1)<=1e-6)),maximum_squared_radius=float(max(rr)),center=lp.x[:3].tolist())
def run():
 G,Y,CT,J,q,V,C=P.geometry();m=G['m'];ids=C['core_ids'];SG=s.Matrix(C['core_gram']);ii=[i for i,j in enumerate(ids)if j!=158];B=s.eye(21)[:,ii];HI=SG.extract(ii,ii).inv();f=s.zeros(20,1);f[ii.index(ids.index(160))]=-1
 p=(B*HI*f).applyfunc(s.simplify);ee=s.eye(21)[:,ids.index(158)];u=(ee-B*HI*B.T*SG*ee).applyfunc(s.simplify);p2=s.simplify((p.T*SG*p)[0]);u2=s.simplify((u.T*SG*u)[0]);amp=s.sqrt(s.simplify((2-p2)/u2));recipes=[(p+sgn*amp*u).applyfunc(s.simplify)for sgn in(-1,1)]
 rs=np.array([np.array(z,float).ravel()@Y[ids]for z in recipes]);cts=np.array([np.array(z,float).ravel()@CT[ids]for z in recipes]);assert max(abs(cts).ravel())<1e-7
 products=rs@m@Y.T+cts@J@CT.T;co=-products/np.sqrt(2*q)[None,:];labels=[[classify(x)for x in row]for row in co]
 removed=sorted({j for row in labels for j,lab in enumerate(row)if lab in('obtuse','angle_gap')});retained=[j for j in range(168)if j not in removed]
 assert not set(G['T']+[159,160,161])&set(removed)
 allY=np.r_[Y,rs];face=poly(allY,[i for i in retained+[168,169]if i not in G['T']]);prism=poly(allY,[159,160,161,168,169])
 individual=[]
 for k,row in enumerate(labels):
  drop=[j for j,lab in enumerate(row)if lab in('obtuse','angle_gap')]
  keep=[j for j in range(168)if j not in drop]
  individual.append(dict(branch=k,removed_old_wall_ids=drop,branch_wall_count=len(keep)+1,face=poly(allY,[i for i in keep+[168+k]if i not in G['T']])))
 bad=[dict(wall=j,types=[labels[k][j]for k in range(2)],cosines=[float(co[k,j])for k in range(2)])for j in removed]
 paths=[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py','results/fresh_h20_e8_cusp_completion_20261017.json','results/fresh_h20_e8_exact_cluster_20261017.json','results/fresh_h20_e8_f4_cusps_20261017.json']
 return dict(status='NUMERICAL_TWO_CAP_PRISM_REPLACEMENT_BRANCH',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},old_base=168,
  exact_recipe=dict(particular=list(map(str,p)),axis=list(map(str,u)),amplitude=str(amp),signs=[-1,1],core_ids=ids),native_new_norms=[2,2],mutual_pair=str(s.simplify(2*p2-2)),
  residual4=rs.tolist(),individual_replacements=individual,all_old_labels=labels,removed_old_wall_ids=removed,removed_gate_details=bad,retained_old_wall_ids=retained,new_branch_wall_count=len(retained)+2,
  old_finite_corner_new_pairs=(rs@m@V.T).tolist(),standalone_prism=prism,replacement_face=face,
  scope='Explicit alternative to frozen168: remove every old wall with a numerically non-Coxeter contact to either exact-recipe cap, retain all remaining old walls, add both norm2 caps. Full targetH3 cone enumerated numerically. This is not yet an exact gate/face certificate or a claim of wholeH20 finite volume.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(out.read_text()))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],'removed',d['removed_old_wall_ids'],'branch',d['new_branch_wall_count'],'face',d['replacement_face']['finite'],d['replacement_face']['exterior'],flush=True)
