#!/usr/bin/env python3
"""Polar and positive-plane roof truncations; exact H5 face frontier and escape."""
import argparse,json,sys
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_higher_face_20261013_truncation_verify as T
F=T.F;V=F.V;K=F.K;N=F.N;PREFIX=F.PREFIX+'_final'
POLAR=['seed2','seed6','nextH4','crosspair0','compact1','compact2','compact3','roof5','roof6']
RIGHT=['seed2','wall17','crosspair0','compact1','compact2','compact3','roof1','roof4'];ANCHOR='compact0'
def perp(x,E,g):
 M=[[V.pair(a,g,b)for b in E]for a in E];Mi=F.matinv(M);lam=V.mv(Mi,[V.pair(x,g,a)for a in E]);p=N.combine(lam,E);return[V.sub(a,b)for a,b in zip(x,p,strict=True)]
def polar(g,roots,names,center):
 E=[roots[names.index(n)]for n in POLAR];z=perp(center,E,g);q=V.pair(z,g,z);assert q[0]>0;r=[V.neg(V.mul(V.sqrt(V.div(V.rat(2),q)),x))for x in z]
 return r,{'right':POLAR,'complement_norm':V.enc(q),'norm':2}
def metric(g,roots,names,center,branch):
 E=[roots[names.index(n)]for n in RIGHT];a=perp(roots[names.index(ANCHOR)],E,g);qa=V.pair(a,g,a);assert qa[0]>0;c=perp(center,E,g);lam=V.div(V.pair(c,g,a),qa);z=[V.sub(x,V.mul(lam,y))for x,y in zip(c,a,strict=True)];qz=V.pair(z,g,z);assert qz[0]>0;rad=V.div(V.sub(V.rat(2),V.div(V.rat(1),qa)),qz);assert rad[0]>0
 t=V.sqrt(rad);t=V.neg(t)if branch<0 else t;r=[V.add(V.neg(V.div(x,qa)),V.mul(t,y))for x,y in zip(a,z,strict=True)]
 return r,{'right':RIGHT,'anchor':ANCHOR,'label':3,'anchor_projection_norm':V.enc(qa),'last_complement_norm':V.enc(qz),'radicand':V.enc(rad),'branch':branch,'norm':2}
def valid(r,ff,g,roots,names,h):
 if V.dot(r,list(map(V.golden,h)))[0]<=0:return False
 for n,a in zip(names,roots,strict=True):
  p=V.pair(r,g,a)
  if n in ff:
   z=V.neg(F.cos2(ff[n]))
   if not(p[0]<=z[1]and z[0]<=p[1]):return False
  elif p[1]>=-2*V.SCALE:return False
 return True

def build():
 prior,(G,I,h,g,roots,names,center,pairs)=T.build();assert prior==V.read(T.PREFIX+'_verify.json');checks=[];stats=[];newsets=[];recipes=[]
 for which,new in enumerate(['roofPolar','roofMetric']):
  if which==0:r,rec=polar(g,roots,names,center);ff={n:2 for n in POLAR}
  else:
   ff={n:2 for n in RIGHT};ff[ANCHOR]=3;branches=[(b,*metric(g,roots,names,center,b))for b in(-1,1)];good=[x for x in branches if valid(x[1],ff,g,roots,names,h)];assert len(good)==1,[(x[0],valid(x[1],ff,g,roots,names,h))for x in branches];_,r,rec=good[0]
  assert valid(r,ff,g,roots,names,h);pp=[]
  for n,a in zip(names,roots,strict=True):
   p=V.pair(r,g,a);pairs[F.key(new,n)]={'kind':'finite','label':ff[n]}if n in ff else{'kind':'dotted'};pp.append({'wall':n,'kind':pairs[F.key(new,n)],'pairing':V.enc(p)})
  roots.append(r);names.append(new);checks.append(pp);recipes.append(rec);st={'name':new,'subsets':512,'dotted_excluded':0,'positive_definite':0,'indefinite':0,'affine':0};neigh=list(ff)
  for mask in range(512):
   ss=[new]+[n for j,n in enumerate(neigh)if mask>>j&1]
   if any(pairs[F.key(a,b)]['kind']=='dotted'for a,b in combinations(ss,2)):st['dotted_excluded']+=1;continue
   typ,_=F.inertia(ss,pairs);st[typ]+=1
   if len(ss)==9 and typ=='positive_definite':newsets.append(ss)
  stats.append(st);print(new,st,flush=True)
 fixed=[n for n in F.COMMON if n!=F.FREE];facefacets=[];excluded=[]
 for n in names:
  if n in fixed:continue
  ss=fixed+[n]
  if any(pairs[F.key(a,b)]['kind']=='dotted'for a,b in combinations(ss,2)):excluded.append({'wall':n,'reason':'dotted'});continue
  typ,_=F.inertia(ss,pairs)
  if typ=='positive_definite':facefacets.append(n)
  else:excluded.append({'wall':n,'reason':'indefinite'})
 assert len(facefacets)==19
 assert all(pairs[F.key(n,a)]=={'kind':'finite','label':2}for n in facefacets for a in fixed)
 facevertices=[];census={'subsets':0,'dotted_excluded':0,'positive_definite':0,'indefinite':0,'affine':0}
 for subset in combinations(facefacets,5):
  census['subsets']+=1;ss=list(subset)
  if any(pairs[F.key(a,b)]['kind']=='dotted'for a,b in combinations(ss,2)):census['dotted_excluded']+=1;continue
  typ,_=F.inertia(ss,pairs);census[typ]+=1
  if typ=='positive_definite':facevertices.append(F.vertex(fixed+ss,roots,names,center,g,h,pairs))
 edges={}
 for i,v in enumerate(facevertices):
  for ee in combinations([n for n in v['incident']if n not in fixed],4):edges.setdefault(tuple(sorted(ee)),[]).append(i)
 assert all(len(v)in(1,2)for v in edges.values())
 # An additional finite/ideal endpoint would supply a ninth normal. The full
 # face census has found all finite endpoints. Affine-freeness rules out an
 # ideal endpoint supported by an additional facet. Singleton edge records
 # therefore run to the light cone without another retained facet.
 graph=[set()for _ in facevertices]
 for ends in edges.values():
  if len(ends)==2:graph[ends[0]].add(ends[1]);graph[ends[1]].add(ends[0])
 components=[];seen=set()
 for i in range(len(graph)):
  if i in seen:continue
  comp={i};todo=[i];seen.add(i)
  while todo:
   j=todo.pop()
   for k in graph[j]-seen:seen.add(k);comp.add(k);todo.append(k)
  components.append(sorted(comp))
 out={'status':'PASS_EXACT_LINKED_H5_ROOF_WITH_TRUNCATIONS','sources':{T.PREFIX+'_verify.json':V.sha(T.PREFIX+'_verify.json'),'scripts/'+Path(__file__).name:V.sha('scripts/'+Path(__file__).name)},'facets':len(roots),'new_recipes':recipes,'new_pairs':checks,'new_affine_census':stats,'H5_face':{'fixed_walls':fixed,'facets':facefacets,'excluded_walls':excluded,'vertex_census':census,'finite_vertices':facevertices,'finite_components':components,'finite_edges':sum(len(x)==2 for x in edges.values()),'open_edges':sum(len(x)==1 for x in edges.values()),'edges':[{'walls':list(k),'endpoints':v}for k,v in sorted(edges.items())]},'normals':[list(map(V.enc,r))for r in roots[-2:]],'scope':'All118 walls are acute-compatible actual facets with a common strict timelike center. The retained four orthogonal walls define an actual H5 Coxeter face. Its19 facets and finite-vertex/edge census are complete for this literal face; nonclosure of this face does not exclude further caps or replacements.'}
 escape_path=PREFIX+'_escape.json'
 if(V.ROOT/escape_path).exists():
  dd=V.read(escape_path);y=tuple(map(K.parse,dd['point']));q=K.pair(y,I,y);ft=K.pair(y,I,h);assert K.sign(q)>0 and K.sign(ft)<0;pp=[V.dot(r,list(map(V.golden,y)))for r in roots];assert all(p[0]>0 for p in pp)
  xx=[V.golden(K.neg(a))for a in K.mv(I,y)];EE=[roots[names.index(n)]for n in fixed];facepoint=perp(xx,EE,g);faceq=V.pair(facepoint,g,facepoint);assert faceq[0]>0 and V.dot(facepoint,list(map(V.golden,h)))[0]>0
  facechecks=[{'wall':n,'pairing':V.enc(V.pair(facepoint,g,r))}for n,r in zip(names,roots,strict=True)if n not in fixed];assert all(int(z['pairing'][1])<0 for z in facechecks)
  out['H5_face']['strict_spacelike_escape']={'construction':'Orthogonal projection of the rational physical escape onto the four fixed walls; their exact Gram is2I4.','point':list(map(V.enc,facepoint)),'norm':V.enc(faceq),'nonincident_inequalities':facechecks}
  out['sources'][escape_path]=V.sha(escape_path);out['strict_escape']={'point':dd['point'],'norm':K.enc(q),'future_center_pairing':K.enc(ft),'pairings':list(map(V.enc,pp))};out['scope']+=' The saved strict future spacelike escape proves this118-wall chamber has infinite volume.'
 return out,(G,I,h,g,roots,names,center,pairs)
def run(verify=False):
 out,_=build()
 if verify:assert out==V.read(PREFIX+'_verify.json')
 else:(V.ROOT/(PREFIX+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['facets'],'H5 vertices',len(out['H5_face']['finite_vertices']),'edges',out['H5_face']['finite_edges'],out['H5_face']['open_edges'])
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
