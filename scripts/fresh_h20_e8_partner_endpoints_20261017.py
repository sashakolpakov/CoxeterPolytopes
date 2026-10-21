#!/usr/bin/env python3
"""Complete finite endpoint diagrams; numerical conic intervals on frozen E8 recipe.

The diagrams are exact. The inherited E8 recipe and conic inequalities are
numerical and are not promoted to an exact coherent chamber by this script.
"""
import argparse, hashlib, itertools, json, math
from pathlib import Path
from fractions import Fraction
import networkx as nx
import numpy as np
import sympy as sp
from scipy.linalg import null_space
from scipy.spatial import HalfspaceIntersection
import fresh_h20_spherical_attachment_profiles_20261016_scout as S
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_e8_partner_endpoints_20261017'
RECIPE='results/fresh_h20_e8_triangle_cluster_20261016_recipe.json'
EDGES=[(159,160),(159,161),(160,161)]
def load(p):return json.loads((ROOT/p).read_text())
def graph(n,edges):
 g=nx.Graph();g.add_nodes_from(range(n));g.add_edges_from((i,j,{'label':m})for i,j,m in edges);return g
def models():
 yield 'A1',graph(1,[])
 for n in range(2,21):yield f'A{n}',graph(n,[(i,i+1,3)for i in range(n-1)])
 for n in range(3,21):yield f'B{n}',graph(n,[(i,i+1,4 if i==0 else 3)for i in range(n-1)])
 for n in range(4,21):yield f'D{n}',graph(n,[(0,2,3),(1,2,3)]+[(i,i+1,3)for i in range(2,n-1)])
 for n in [6,7,8]:
  ee=[];v=1
  for arm in [1,2,n-4]:
   prev=0
   for _ in range(arm):ee.append((prev,v,3));prev=v;v+=1
  yield f'E{n}',graph(n,ee)
 yield 'F4',graph(4,[(0,1,3),(1,2,4),(2,3,3)])
 for n in [3,4]:yield f'H{n}',graph(n,[(i,i+1,5 if i==0 else 3)for i in range(n-1)])
def support(edge,T):
 W=T+list(edge);g=nx.Graph();g.add_nodes_from(W)
 for path in [list(range(1,11)),[12,13,14,15],[61,87]]:
  g.add_edges_from((i,j,{'label':3})for i,j in zip(path,path[1:]))
 for i,j,m in [(1,159,3),(10,160,3),(159,161,4),(160,161,3)]:
  if i in W and j in W:g.add_edge(i,j,label=m)
 return W,g
def census(edge,T):
 W,old=support(edge,T);components=sorted((sorted(c)for c in nx.connected_components(old)),key=lambda c:(-len(c),c))
 G=sp.eye(19)*2
 for i,j,a in old.edges(data=True):G[W.index(i),W.index(j)]=G[W.index(j),W.index(i)]=-2*sp.cos(sp.pi/a['label'])
 J=G.inv();rows={};matches=0
 for name,g in models():
  for deleted in g:
   h=g.copy();h.remove_node(deleted)
   for bits in itertools.product([0,1],repeat=len(components)):
    ids=sum((c for c,b in zip(components,bits)if b),[])
    if len(ids)!=len(h):continue
    gm=nx.algorithms.isomorphism.GraphMatcher(h,old.subgraph(ids),edge_match=lambda a,b:a['label']==b['label'])
    for mapping in gm.isomorphisms_iter():
     contacts=tuple(sorted((mapping[j],g[deleted][j]['label'])for j in g.neighbors(deleted)));matches+=1
     if len(contacts)==1 and contacts[0][0]==151:continue
     row=rows.setdefault(contacts,dict(contacts=[list(x)for x in contacts],models=[]))
     row['models'].append(dict(type=name,deleted=deleted,mapping={str(k):v for k,v in mapping.items()}))
 out=[]
 for contacts,row in sorted(rows.items(),key=lambda x:(-len(x[0]),x[0])):
  d=sp.zeros(19,1)
  for i,m in contacts:d[W.index(i)]=2*sp.cos(sp.pi/m)
  cost=sp.simplify((d.T*J*d)[0]);res=sp.simplify(2-cost);assert res.is_positive
  row.update(id=len(out),cost=str(cost),residual_norm=str(res),numeric_residual_norm=float(res));out.append(row)
 return dict(edge=list(edge),support=W,components=components,gram=[[str(x)for x in r]for r in G.tolist()],inverse=[[str(x)for x in r]for r in J.tolist()],profiles=out,raw_matches=matches,parametric_profile=dict(contacts=[[151,'m']],integer_domain='m>=3',cost='2*cos(pi/m)^2',residual_norm='2*sin(pi/m)^2'))
def intersections(intervals,allowed):
 return [(max(a,c),min(b,d))for a,b in intervals for c,d in allowed if max(a,c)<min(b,d)]
def quadratic_intervals(a,b,c):
 scale=max(abs(a),abs(b),abs(c),1.);a/=scale;b/=scale;c/=scale;tol=1e-12
 if abs(a)<tol:
  if abs(b)<tol:return [(0.,math.inf)]if c<=tol else[]
  r=-c/b
  return [(max(0.,r),math.inf)]if b<0 else ([(0.,r)]if r>0 else[])
 disc=b*b-4*a*c
 if disc < -tol:return [(0.,math.inf)]if a<0 else[]
 root=math.sqrt(max(0.,disc));r0,r1=sorted(((-b-root)/(2*a),(-b+root)/(2*a)))
 if a>0:return [(max(0.,r0),r1)]if r1>max(0.,r0)else[]
 out=[]
 if r0>0:out.append((0.,r0))
 out.append((max(0.,r1),math.inf));return out
def geometry():
 D,G=S.geometry();R=load(RECIPE);Y=np.r_[G['Y'],[r['residual4']for r in R['steps']]]
 CT=np.c_[G['CS'],G['c151']];new=[]
 for rec in R['steps']:
  v=np.zeros(17)
  for j,m in rec['support_contacts']:v[G['T'].index(j)]=-2*np.cos(np.pi/m)
  new.append(v)
 CT=np.r_[CT,new];J=np.zeros((17,17));J[:16,:16]=G['J'];J[16,16]=.5
 q=np.r_[G['q'],[2.]*4]
 live=[i for i in range(162)if i not in G['T'] and np.linalg.norm(Y[i,1:])>1e-9]
 A=Y[live,1:];b=Y[live,0]
 verts=HalfspaceIntersection(np.c_[A,-b],np.array(R['face']['center'])).intersections
 rad=np.sum(verts*verts,axis=1);center=verts[rad<1-1e-6].mean(axis=0)
 assert center@center<1 and max(A@center-b)<-1e-8
 return G,Y,CT,J,q,np.r_[1.,center],live,verts
def poly(Y,live,center,r):
 A=np.r_[Y[live,1:],r[None,1:]];b=np.r_[Y[live,0],r[0]]
 if max(A@center[1:]-b)>=0:return None
 v=HalfspaceIntersection(np.c_[A,-b],center[1:]).intersections;rr=np.sum(v*v,axis=1)
 return dict(vertices=len(v),finite=int(sum(rr<1-1e-6)),exterior=int(sum(rr>1+1e-6)),near_sphere=int(sum(abs(rr-1)<=1e-6)),maximum_squared_radius=float(max(rr)))
def conic(cen,row,G,Y,CT,JT,q,H,live):
 edge=cen['edge'];m=G['m'];W=cen['support'];U=Y[edge];inv=np.linalg.inv(U@m@U.T)
 dt=np.zeros(17);de=np.zeros(2)
 for i,label in row['contacts']:
  if i in G['T']:dt[G['T'].index(i)]=2*np.cos(np.pi/label)
  else:de[edge.index(i)]=2*np.cos(np.pi/label)
 offset=-CT@JT@dt
 f=-de-offset[edge];p=(inv@f)@U
 h=2-dt@JT@dt-p@m@p;assert h>0 and abs(h-row['numeric_residual_norm'])<1e-7
 ht=np.array([1.,0,0,0]);ht-=((ht@m@U.T)@inv)@U;ht/=math.sqrt(-(ht@m@ht))
 es=null_space(np.r_[U@m,[ht@m]])[:,0];es/=math.sqrt(es@m@es)
 out=[]
 for sign in [-1,1]:
  va=math.sqrt(h)/2*(ht+sign*es);vb=math.sqrt(h)/2*(-ht+sign*es)
  a=va@m@Y.T;b=vb@m@Y.T;c=p@m@Y.T+offset
  intervals=[(0.,math.inf)];constraints=[]
  for j in range(162):
   if j in W:continue
   coeff=(a[j],c[j]+math.sqrt(2*q[j]),b[j]);intervals=intersections(intervals,quadratic_intervals(*coeff))
   constraints.append(dict(root=j,quadratic=list(map(float,coeff))))
  hc=(float(va@m@H),float(p@m@H+1e-9),float(vb@m@H))
  intervals=intersections(intervals,quadratic_intervals(*hc))
  samples=[]
  for lo,hi in intervals:
   if math.isinf(hi):ks=[max(1.,lo)*x for x in [1.001,1.1,2]]
   elif lo==0:ks=[hi*x for x in [.1,.5,.999]]
   else:ks=[lo+(hi-lo)*f for f in [.01,.5,.99]]
   for k in ks:
    r=p+k*va+vb/k;co=-(r@m@Y.T+offset)/np.sqrt(2*q)
    other=[j for j in range(162)if j not in W]
    if min(co[other])<1-1e-6:continue
    assert abs(r@m@r+dt@JT@dt-2)<1e-6
    samples.append(dict(k=k,residual4=r.tolist(),minimum_unprescribed_cosine=float(min(co[other])),face=poly(Y,live,H,r)))
  out.append(dict(branch=sign,particular=p.tolist(),positive_k_coefficient=va.tolist(),inverse_k_coefficient=vb.tolist(),T_pairings=(-dt).tolist(),residual_norm=h,constraints=constraints,center_quadratic=list(hc),all_unprescribed_dotted_intervals=[[a,None if math.isinf(b)else b]for a,b in intervals],samples=samples))
 return out
def edge_ledger(G,Y,CT,J,q,verts):
 m=G['m'];sides=G['sides']+list(range(158,162));rows=[]
 for i,j in itertools.combinations(sides,2):
  Q=Y[[i,j]]@m@Y[[i,j]].T
  if min(np.linalg.eigvalsh(Q))<1e-8:continue
  active=[]
  for k,v in enumerate(verts):
   scale=1+np.linalg.norm(Y[[i,j]],axis=1)
   if max(abs(Y[[i,j],1:]@v-Y[[i,j],0])/scale)<1e-7:active.append(k)
  if not active:continue
  finite=[k for k in active if verts[k]@verts[k]<1-1e-6]
  rows.append(dict(pair=[i,j],vertex_indices=active,finite_vertex_indices=finite,open=len(finite)<2))
 return rows
def run():
 G,Y,CT,J,q,H,live,verts=geometry();out=[]
 for edge in EDGES:
  c=census(edge,G['T']);families=[]
  for row in c['profiles']:
   z=conic(c,row,G,Y,CT,J,q,H,live)
   families.append(dict(profile=row['id'],branches=z))
  c['families']=families;out.append(c)
  print('EDGE',edge,'profiles',len(c['profiles']),'interval branches',sum(bool(b['all_unprescribed_dotted_intervals'])for r in families for b in r['branches']),flush=True)
 partner=np.array([[sp.Rational(12,11),-sp.Rational(1,11),-sp.sqrt(2)],[-sp.Rational(1,11),sp.Rational(12,11),-1],[-sp.sqrt(2),-1,2]],dtype=object)
 determinant=sp.factor(sp.Matrix(partner).det());assert determinant.is_negative
 ledger=edge_ledger(G,Y,CT,J,q,verts)
 src=[f'scripts/{TAG}.py',RECIPE,S.INPUT,f'scripts/{S.TAG}.py']
 return dict(status='EXACT_ENDPOINT_PROFILE_CENSUS_WITH_NUMERICAL_CONIC_INTERVALS',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},base_walls=162,base_is_numerical_recipe=True,old_face_vertices=verts.tolist(),edge_ledger=ledger,endpoint_censuses=out,partner_projected_gram=[[str(x)for x in row]for row in partner.tolist()],partner_projected_determinant=str(determinant),center=H.tolist(),scope='Finite diagram endpoint profiles are complete, with one explicit unbounded I2(m) family per edge. Conic formulas are analytic; interval endpoints and old ambient pairings are numerical. Only the all-unprescribed-dotted open parameter subfamilies are solved here. Isolated finite-angle coincidences and the unbounded I2(m) family are not exhausted. No exact new wall or finite-volume claim.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,load(str(out.relative_to(ROOT))))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
