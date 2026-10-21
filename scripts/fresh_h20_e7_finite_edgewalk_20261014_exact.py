#!/usr/bin/env python3
"""Two coherent nonright edge caps on the frozen88; exact geometric replay."""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
from collections import Counter
from hashlib import sha256
import fresh_h20_e7_residual_gate_20261014_exact_data as X
import fresh_h20_e7_causal_residuals_20261014_bridge_review as B
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e7_finite_edgewalk_20261014_exact'
BRIDGE='results/fresh_h20_e7_causal_residuals_20261014_bridge.json';K=X.K

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return[enc(y)for y in x]
 if isinstance(x,dict):return{str(k):enc(y)for k,y in x.items()}
 return x

def run():
 data=X.load();parent=data['parent'];base=data['base'];N=list(data['N']);H=data['H'];E=data['E'];dot=data['dot'];S=data['S'];face=parent['face_facet_walls'];sq2=(F(0),F(1))+(F(0),)*6
 dbridge=json.loads((ROOT/BRIDGE).read_text())
 for p,h in dbridge['sources'].items():assert sha(p)==h,p
 scales={i:sq2 if i in face[11:13]else K.scalar(1)for i in range(90)}
 def total(xs):
  v=K.Z
  for x in xs:v=K.add(v,x)
  return v
 def lin(vs,cs):return tuple(total(K.multiply(c,v[k])for c,v in zip(cs,vs))for k in range(21))
 def rlin(vs,cs):return lin(vs,list(map(K.scalar,cs)))
 recipes=[];new=[]
 for vi,contacts,boost in [(4,[0,14],F(7,4)),(10,[10,14],F(3,2))]:
  row=parent['face_finite_vertices'][vi];released=13;ids=[i for i in row['ambient_incident']if i!=released]
  scaled=[lin([N[i]],[scales[i]])for i in ids];G=[[K.rational(dot(a,b))for b in scaled]for a in scaled];I=K.inverse(G)
  bv=[K.neg(scales[i])if i in contacts else K.Z for i in ids]
  coeff=[total(K.scale(bv[j],I[i][j])for j in range(19))for i in range(19)];p=lin(scaled,coeff)
  assert all(dot(p,N[i])==K.scalar(-int(i in contacts))for i in ids)
  V=K.vector(row['ambient_point']);omega=data['W'][S.index(released)];b=dot(p,N[released]);assert b==K.scalar(F(2,3))
  residual=K.rational(K.sub(K.scalar(2),dot(p,p)));q=-K.rational(dot(V,V));degree=data['I'][S.index(released)][S.index(released)]
  assert residual>0 and boost*boost*q==F(4,9)*degree-residual
  r=lin([p,omega,V],[K.scalar(1),K.scalar(F(-2,3)),K.scalar(boost)])
  assert dot(r,r)==K.scalar(2)and dot(r,N[released])==K.Z
  new.append(r);recipes.append(dict(vertex=vi,released=released,contacts=[[i,3]for i in contacts],omitted_angle=2,support=ids,projection=p,projection_norm=dot(p,p),residual=residual,dual_coefficient=F(2,3),vertex_coefficient=boost,normal=r))
 assert dot(new[0],new[1])==K.Z
 # Bridge's positive radical is evaluated by exact squaring over an
 # independently encoded 16-dimensional multiquadratic base field.
 DZ=tuple(map(F,dbridge['field']['discriminant']));assert B.sign(DZ)>0;Z=B.Z;EZ=(Z,Z)
 def ee(x):return(tuple(x)+(F(0),)*(16-len(x)),Z)
 def ep(x):return tuple(tuple(map(F,c))for c in x)
 def ea(a,b):return(B.add(a[0],b[0]),B.add(a[1],b[1]))
 def es(a,q):return(B.scale(a[0],q),B.scale(a[1],q))
 def em(a,b):return(B.add(B.mul(a[0],b[0]),B.mul(DZ,B.mul(a[1],b[1]))),B.add(B.mul(a[0],b[1]),B.mul(a[1],b[0])))
 def sg(a):
  x,y=a;sx,sy=B.sign(x),B.sign(y)
  if not sy:return sx
  if not sx or sx==sy:return sy
  return sx*B.sign(B.sub(B.mul(x,x),B.mul(DZ,B.mul(y,y))))
 def et(xs):
  a=EZ
  for x in xs:a=ea(a,x)
  return a
 sp=[[(j,q)for j,q in enumerate(row)if q]for row in data['Q']];cache={};simple_cache={}
 def simple(v):
  if v not in simple_cache:simple_cache[v]=tuple(x[0][:8]for x in v)if all(x[1]==Z and not any(x[0][8:])for x in v)else None
  return simple_cache[v]
 def edot(v,w):
  a,b=simple(v),simple(w)
  if a is not None and b is not None:return ee(dot(a,b))
  if w not in cache:cache[w]=tuple(et(es(w[j],q)for j,q in row)for row in sp)
  return et(em(a,b)for a,b in zip(v,cache[w]))
 def ev(v):return tuple(ee(x)for x in v)
 EN=[ev(v)for v in N]+[tuple(map(ep,dbridge['normals'][-1]))]+[ev(v)for v in new]
 norms=list(map(F,base['norms']))+[F(2)]*3;EH=ev(H);counts=Counter();pairs=[]
 for i in [88,89]:
  for j in range(i):
   g=edot(EN[i],EN[j]);prod=norms[i]*norms[j];assert sg(g)<=0
   if sg(g)==0:kind='right'
   elif em(g,g)==ee(K.scalar(prod/4)):kind='pi3'
   elif em(g,g)==ee(K.scalar(prod/2)):kind='pi4'
   elif em(g,g)==ee(K.scalar(3*prod/4)):kind='pi6'
   elif em(g,g)==ee(K.scalar(prod)):kind='parallel'
   else:assert sg(ea(em(g,g),ee(K.scalar(-prod))))>0;kind='dotted'
   pairs.append(dict(walls=[i,j],pairing=g,relation=kind));counts[kind]+=1
 assert len(pairs)==177 and sg(edot(EH,EH))<0 and all(sg(edot(EH,n))<0 for n in EN)
 # All90 facets are actual by strict-center projection and acuteness.
 facet_points=[]
 for i,n in enumerate(EN):
  c=es(edot(EH,n),-1/norms[i]);p=tuple(ea(h,em(c,x))for h,x in zip(EH,n));assert sg(edot(p,p))<0 and sg(edot(p,EH))<0 and sg(edot(p,n))==0
  facet_points.append(p)
 old_vertices=[]
 for vi,row in enumerate(parent['face_finite_vertices']):
  V=K.vector(row['ambient_point']);assert all(K.sign(dot(V,r))<0 for r in new)
  assert sg(edot(ev(V),EN[87]))<0
  old_vertices.append(dict(vertex=vi,point=V,incident=row['ambient_incident'],new_pairings=[dot(V,r)for r in new]))
 preserved_cusps=[]
 for row in base['completed_cusp_links']:
  nu=K.vector(row['point']);gg=[dot(nu,r)for r in new];assert all(K.sign(g)<0 for g in gg)
  preserved_cusps.append(dict(point=nu,new_pairings=gg,quotient_rank=19))
 # Reconstruct every forward finite endpoint on the308 selected rays;
 # accept only strict nonincident old/other-new signs. No float gate.
 endpoints=[]
 for vi,row in enumerate(parent['face_finite_vertices']):
  V=K.vector(row['ambient_point']);active=row['ambient_incident']
  for sj,j in enumerate(S):
   omega=data['W'][sj];ids=[i for i in active if i!=j]
   for ni,r in enumerate(new):
    numer=dot(r,V);denom=dot(r,omega)
    if K.sign(denom)>=0:continue
    # Both pairings are rational for these two existing-field caps.
    t=K.rational(numer)/K.rational(denom)
    if t<=0:continue
    point=rlin([V,omega],[1,-t])
    if K.sign(dot(point,point))>=0:continue
    active_new=ids+[88+ni]
    signs=[sg(edot(ev(point),n))for n in EN]
    if any(s!=0 if i in active_new else s>=0 for i,s in enumerate(signs)):continue
    assert K.sign(dot(point,H))<0 and dot(point,N[j])==K.scalar(-t)
    assert all(K.sign(dot(V,r0))<0 for r0 in new)
    endpoints.append(dict(old_vertex=vi,released=j,new_wall=88+ni,edge_support=ids,parameter=t,point=point,norm=dot(point,point),incident=active_new))
 assert Counter(row['new_wall']for row in endpoints)=={88:19,89:12}
 assert len({tuple(sorted(row['incident']))for row in endpoints})==31
 # Complete finite-vertex census restricted to those incident to a new cap.
 # The bridge87 is dotted to both, so all relevant roots remain in K.
 phys={i:v for i,v in enumerate(N)};phys.update({88:new[0],89:new[1]});pgram={}
 def gg(i,j):
  key=(min(i,j),max(i,j))
  if key not in pgram:pgram[key]=dot(phys[key[0]],phys[key[1]])
  return pgram[key]
 neigh={k:[i for i in phys if i!=k and K.sign(K.sub(K.multiply(gg(k,i),gg(k,i)),K.scalar(norms[k]*norms[i])))<0]for k in [88,89]}
 universe=sorted(set(neigh[88]+neigh[89]+[88,89]));adj={i:set()for i in universe}
 for i,j in combinations(universe,2):
  if K.sign(K.sub(K.multiply(gg(i,j),gg(i,j)),K.scalar(norms[i]*norms[j])))<0:adj[i].add(j);adj[j].add(i)
 def cliques(candidates,need,prefix=()):
  if not need:yield prefix;return
  for k,i in enumerate(candidates):
   if len(candidates)-k<need:return
   yield from cliques([j for j in candidates[k+1:]if j in adj[i]],need-1,prefix+(i,))
 supports=set()
 for i in [88,89]:
  for c in cliques(sorted(neigh[i]),19):supports.add(tuple(sorted(c+(i,))))
 verts=[];rejected=[]
 for ids in sorted(supports):
  rows=[lin([phys[i]],[scales[i]])for i in ids];G=[[K.rational(dot(a,b))for b in rows]for a in rows]
  if not K.pd(G):rejected.append(ids);continue
  GI=K.inverse(G);hp=[dot(H,r)for r in rows];co=[total(K.scale(hp[j],GI[i][j])for j in range(20))for i in range(20)]
  point=lin([H]+rows,[K.scalar(1)]+[K.neg(x)for x in co]);assert K.sign(dot(point,point))<0 and K.sign(dot(point,H))<0
  assert all(sg(edot(ev(point),n))==0 if i in ids else sg(edot(ev(point),n))<0 for i,n in enumerate(EN))
  verts.append(dict(incident=ids,point=point,norm=dot(point,point)))
 assert len(supports)==3080 and len(rejected)==3040 and len(verts)==40
 assert sum(88 in v['incident'] and 89 in v['incident']for v in verts)==9
 active_sets={tuple(sorted(v['incident'])):i for i,v in enumerate(verts)}
 for row in endpoints:assert tuple(sorted(row['incident']))in active_sets
 # List all finite edges with a new endpoint among these new vertices and
 # the44 inherited ones. This is not a complete ambient edge census.
 allverts=old_vertices+verts;finite_edges=[]
 for i,j in combinations(range(len(allverts)),2):
  if j<44:continue
  common=sorted(set(allverts[i]['incident'])&set(allverts[j]['incident']))
  if len(common)==19:finite_edges.append(dict(vertices=[i,j],support=common))
 assert len(finite_edges)==145 and sum(e['vertices'][0]<44 for e in finite_edges)==31
 escape=tuple(map(ep,dbridge['strict_escape']['point']));assert sg(edot(escape,escape))>0 and sg(edot(escape,EH))<0
 epairs=[edot(escape,n)for n in EN];assert all(sg(g)<0 for g in epairs)
 paths=[X.BASE,X.PARENT,X.FACE,BRIDGE,f'scripts/{TAG}.py','scripts/fresh_h20_e7_residual_gate_20261014_exact_data.py','scripts/fresh_h20_e7_causal_residuals_20261014_bridge_review.py']
 return enc(dict(status='PASS_90_ACTUAL_WALLS_TWO_NONRIGHT_CAPS_AND31_ADJACENT_ENDPOINTS_STILL_INFINITE',sources={p:sha(p)for p in paths},field=dbridge['field'],form=data['Q'],normals=EN,norms=norms,recipes=recipes,new_pair_counts=dict(counts),new_pairs=pairs,strict_center=EH,facet_points=facet_points,preserved_old_vertices=old_vertices,preserved_cusps=preserved_cusps,selected_edge_endpoints=endpoints,new_cap_finite_neighbor_sets=neigh,new_cap_candidate_cliques=len(supports),new_cap_nonPD_cliques=rejected,new_cap_finite_vertices=verts,finite_edges_with_new_endpoint=finite_edges,strict_escape=dict(point=escape,norm=edot(escape,escape),all90_pairings=epairs),scope='Retains frozen88 and adds two mutually orthogonal norm2 caps. Complete finite vertices incident to these two new caps are enumerated through their entire finite-neighbor graphs.31 specified released-E7 edges have first finite endpoints certified by future timelike points and strict signs. Finite edges listed only among44 specified inherited vertices and the new-cap vertices. No complete old88 vertex/edge census or finite-volume claim.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status']);print(d['new_pair_counts'],'newvertices',len(d['new_cap_finite_vertices']),'finiteedges',len(d['finite_edges_with_new_endpoint']))
