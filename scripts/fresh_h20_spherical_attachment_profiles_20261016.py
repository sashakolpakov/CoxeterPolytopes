#!/usr/bin/env python3
"""Finite spherical attachment census and numerical full-ambient geometry.

No new wall/volume theorem is asserted by the numerical geometry cache.
All old roots use their native positive scales; the trial root has norm two.
"""
import argparse, hashlib, itertools, json
from pathlib import Path
from fractions import Fraction
import numpy as np
import networkx as nx
import mpmath as mp
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_spherical_attachment_profiles_20261016'
S16=[1,2,3,4,5,6,7,8,9,10,12,13,14,15,61,87]
PATHS=[[1,2,3,4,5,6,7,8,9,10],[12,13,14,15],[61,87],[151]]
T=S16+[151]
FILES=['results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json','results/fresh_h20_e7_bridge_face_20261014_exact.json','results/fresh_h20_bridge_edge_closure_20261014.json','results/fresh_h20_polar_pair_20261014.json','results/fresh_h20_coupled_roles_20261014_exact.json']
RAD=[1,2,79,158,11,22,869,1738,3,6,237,474,33,66,2607,5214]
def load(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def graph(n,edges):
 g=nx.Graph();g.add_nodes_from(range(n));g.add_edges_from((i,j,{'label':m})for i,j,m in edges);return g
def finite_models():
 yield 'A1',graph(1,[])
 for n in range(2,19):yield f'A{n}',graph(n,[(i,i+1,3)for i in range(n-1)])
 for n in range(3,19):yield f'B{n}',graph(n,[(i,i+1,4 if i==0 else 3)for i in range(n-1)])
 for n in range(4,19):yield f'D{n}',graph(n,[(0,2,3),(1,2,3)]+[(i,i+1,3)for i in range(2,n-1)])
 for n in [6,7,8]:
  edges=[];v=1
  for arm in [1,2,n-4]:
   prev=0
   for _ in range(arm):edges.append((prev,v,3));prev=v;v+=1
  yield f'E{n}',graph(n,edges)
 yield 'F4',graph(4,[(0,1,3),(1,2,4),(2,3,3)])
 for n in [3,4]:yield f'H{n}',graph(n,[(i,i+1,5 if i==0 else 3)for i in range(n-1)])
def census():
 old=nx.Graph();old.add_nodes_from(T)
 for path in PATHS:old.add_edges_from((i,j,{'label':3})for i,j in zip(path,path[1:]))
 gram=sp.eye(17)*2
 for i,j in old.edges:gram[T.index(i),T.index(j)]=gram[T.index(j),T.index(i)]=-1
 inv=gram.inv();rows={};matches=0
 for name,g in finite_models():
  for deleted in g:
   h=g.copy();h.remove_node(deleted)
   for bits in itertools.product([0,1],repeat=4):
    ids=sum((p for p,b in zip(PATHS,bits)if b),[])
    if len(ids)!=len(h):continue
    target=old.subgraph(ids)
    gm=nx.algorithms.isomorphism.GraphMatcher(h,target,edge_match=lambda a,b:a['label']==b['label'])
    for mapping in gm.isomorphisms_iter():
     contacts=tuple(sorted((mapping[j],g[deleted][j]['label'])for j in g.neighbors(deleted)))
     matches+=1
     # All one-A1 contacts are one unbounded I2(m) family below.
     if len(contacts)==1 and contacts[0][0]==151:continue
     key=contacts
     row=rows.setdefault(key,dict(contacts=[list(x)for x in contacts],models=[]))
     row['models'].append(dict(type=name,deleted=deleted,mapping={str(k):v for k,v in mapping.items()}))
 for contacts,row in rows.items():
  d=sp.zeros(17,1)
  for i,m in contacts:d[T.index(i)]=2*sp.cos(sp.pi/m)
  cost=sp.simplify((d.T*inv*d)[0]);res=sp.simplify(2-cost)
  assert res.is_positive,(contacts,cost)
  row.update(cost=str(cost),residual_norm=str(res),numeric_cost=float(cost),numeric_residual_norm=float(res))
 result=sorted(rows.values(),key=lambda r:(-len(r['contacts']),r['contacts']))
 for i,r in enumerate(result):r['id']=i
 e8=[r['id']for r in result if any(m['type']=='E8'for m in r['models'])]
 assert len(e8)==4 and all(result[i]['residual_norm']=='1/30'for i in e8)
 return dict(support=T,paths=PATHS,gram=[[str(x)for x in row]for row in gram.tolist()],inverse=[[str(x)for x in row]for row in inv.tolist()],profiles=result,E8_profiles=e8,raw_marked_deletions=matches,parametric_family=dict(contacts=[[151,'m']],integer_domain='m >= 3',cost='2*cos(pi/m)^2',residual_norm='2*sin(pi/m)^2'),proof='The connected finite component containing the new root contains whole old irreducible components. Deleting it leaves exactly a subset of the four old paths. Rank is at most18. A/B/D/E/F/H and I2(m) exhaust finite irreducibles; all labelled deletion isomorphisms are retained. Every root is normalized to norm2 for this census.')
def geometry():
 mp.mp.dps=75
 old,bridge,edge,polar,cur=[load(p)for p in FILES]
 radicals=[]
 def scalar(x,k):
  if k==0:return sum(mp.mpf(str(Fraction(c)))*mp.sqrt(d)for c,d in zip(x,RAD))
  return scalar(x[0],k-1)+scalar(x[1],k-1)*radicals[k-1]
 for x,k in [(old['field']['discriminant'],0),(bridge['field']['new_positive_sqrt_radical'],1),(edge['field']['new_positive_sqrt_radical'],2)]:radicals.append(mp.sqrt(scalar(x,k)))
 def v(raw,k):return mp.matrix([scalar(x,k)for x in raw])
 def rat(x):return mp.mpf(str(Fraction(x)))
 def midpoint(x):return (rat(x[0])+rat(x[1]))/2
 def boundsvec(raw):return mp.matrix([midpoint(x)for x in raw])
 Q=mp.matrix([[rat(x)for x in row]for row in old['form']])
 def dot(v,w):return (v.T*Q*w)[0]
 N=[v(row,1)for row in old['normals']]+[v(bridge['new_normal'],2)]+[v(row,3)for row in edge['new_normals']]+[mp.matrix([rat(x)for x in row])for row in polar['new_normals']]+[v(cur['cap_R']['normal'],1)]
 N=[n*mp.sqrt(2/dot(n,n))for n in N];assert len(N)==101
 hv=[v(cur['hex_basis'][key],1)/mp.sqrt(scalar(q,1))for key,q in zip(['C','U','V'],cur['hex_basis']['squared_lengths'])]
 for co in cur['hexagon']['frame_coefficients']:N.append(sum((w*scalar(c,0)for w,c in zip(hv,co)),mp.zeros(21,1)))
 base=load('results/fresh_h20_s61_multicap_20261014_exact.json');S=base['support'];J=mp.matrix([[rat(x)for x in row]for row in base['support_inverse']])
 frame=[v(base['frame']['X'],2)/mp.sqrt(scalar(base['frame']['A'],2)),v(base['frame']['r95'],2)/mp.sqrt(2),v(base['frame']['W'],2)/mp.sqrt(scalar(base['frame']['W_norm'],2))]
 for co in base['chain']['frame_coefficient_bounds']:N.append(sum((w*midpoint(c)for w,c in zip(frame,co)),mp.zeros(21,1)))
 omega=[sum((N[j]*J[i,k]for k,j in enumerate(S)),mp.zeros(21,1))for i in range(18)]
 arc=load('results/fresh_h20_dotted_arc_cap_20261016.json');cap=load('results/fresh_h20_complementary_cap_20261016.json')
 N.append(sum((w*midpoint(c)for w,c in zip(frame,arc['definition']['frame_bounds'])),mp.zeros(21,1))-omega[S.index(31)]*midpoint(arc['definition']['dual_coefficient_bounds']))
 N.append(sum((w*midpoint(c)for w,c in zip(frame,cap['definition']['frame_bounds'])),mp.zeros(21,1))-omega[S.index(5)])
 assert len(N)==135
 later=[('fresh_h20_affine_a17_completion_20261016','new_normal_bounds'),('fresh_h20_off_cusp_cap_20261016','physical_normal_bounds'),('fresh_h20_mixed_isolated_cap_20261016','physical_normal_bounds'),('fresh_h20_residual_h4_corner_20261016','new_normal_bounds'),('fresh_h20_h4_parallel_partner_20261016','physical_normal_bounds'),('fresh_h20_second_cusp_completion_20261016','new_normal_bounds'),('fresh_h20_second_cusp_pi8_cap_20261016','physical_normal_bounds'),('fresh_h20_pi8_triangle_face_20261016','physical_normal_bounds'),('fresh_h20_dominant142_cap_20261016','physical_normal_bounds')]
 sources=FILES+['results/fresh_h20_s61_multicap_20261014_exact.json','results/fresh_h20_dotted_arc_cap_20261016.json','results/fresh_h20_complementary_cap_20261016.json']
 for name,key in later:
  p='results/'+name+'.json';sources.append(p);raw=load(p)['definition'][key]
  N.extend([boundsvec(row)for row in raw]if key=='new_normal_bounds'else[boundsvec(raw)])
 p='results/fresh_h20_dominant_axis_cluster_20261016.json';sources.append(p)
 N.extend(boundsvec(row['physical_normal_bounds'])for row in load(p)['definition']['new_roots'])
 for name in ['fresh_h20_release144_cap_20261016','fresh_h20_compact_quad_cap_20261016']:
  p='results/'+name+'.json';sources.append(p);N.append(boundsvec(load(p)['definition']['physical_normal_bounds']))
 assert len(N)==158,len(N)
 NN=mp.matrix([[x for x in n]for n in N]);FG=NN*Q*NN.T
 GS=mp.matrix([[FG[i,j]for j in S16]for i in S16]);JI=GS**-1
 C=mp.matrix([[FG[i,j]for j in S16]for i in range(158)])
 F5=mp.matrix([[x for x in n]for n in frame+[N[17]/mp.sqrt(2),N[31]/mp.sqrt(2)]])
 MM=mp.diag([-1,1,1,1,1]);XX=NN*Q*F5.T*MM
 rebuilt=XX*MM*XX.T+C*JI*C.T
 error=max(abs(FG[i,j]-rebuilt[i,j])/(1+abs(FG[i,j]))for i in range(158)for j in range(158))
 assert error<mp.mpf('1e-30'),mp.nstr(error)
 X=np.array(XX.tolist(),float);M=np.diag([-1.,1.,1.,1.,1.])
 norms=np.array([float(FG[i,i])for i in range(158)]);cn=np.array(C.tolist(),float);jn=np.array(JI.tolist(),float)
 # A double-precision norm subtraction is unsafe for some very large old roots.
 # Native norms and fixed support pairings are always supplied separately.
 from fresh_h20_quad_adjacent_h3_closure_20261016 import geometry158
 G=geometry158();assert np.max(abs(norms-G['norms']))<2e-7
 H=frame[0]-frame[1]/2
 HX=np.array([(H.T*Q*f)[0]*(-1 if i==0 else 1)for i,f in enumerate(frame+[N[17]/mp.sqrt(2),N[31]/mp.sqrt(2)])],float)
 return dict(source_files=sources,physical_form=[[str(x)for x in row]for row in old['form']],support=S16,support_gram=np.array(GS.tolist(),float).tolist(),support_inverse=jn.tolist(),old_support_pairings=cn.tolist(),residual_normals=X.tolist(),residual_metric=M.tolist(),native_norms=norms.tolist(),weak_residual_center=HX.tolist(),maximum_relative_Gram_reconstruction_error=str(error),physical_normal_bounds_used='Frozen primary coordinate enclosures evaluated at75decimaldigits; original107 normals re-evaluated in the literal radical tower. Full Gram decomposition checked at75digits before conversion to scout doubles. Numerical scout only.')

def run():
 c=census();g=geometry();paths=[f'scripts/{TAG}.py']+g.pop('source_files')+['scripts/fresh_h20_quad_adjacent_h3_closure_20261016.py']
 return dict(status='COMPLETE_FINITE_GRAPH_PROFILE_CENSUS_WITH_NUMERICAL_AMBIENT_CACHE',sources={p:sha(p)for p in paths},census=c,geometry=g,scope='Complete labelled finite one-node graph profiles, with explicit unbounded I2(m) family. Numerical physical geometry cache only; no new wall or finite-volume assertion.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();d=run();path=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==load(str(path.relative_to(ROOT)))
 else:path.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['census']['profiles']),d['census']['E8_profiles'],d['geometry']['maximum_relative_Gram_reconstruction_error'],flush=True)
