#!/usr/bin/env python3
"""Two coupled endpoint caps and an actual compact H2 face inside H5."""
import argparse,json,sys
from pathlib import Path
from itertools import combinations
from collections import Counter,defaultdict
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_h5_closure_20261013_verify as O
F=O.F;V=F.V;K=F.K;N=F.N;PRE='data/enumerations/fresh_h9_h5_closure_next_20261013';OLD=O.PRE+'_verify.json'
def old_state():
 G,I,h,g,roots,names,c,pairs,_,sources=O.old_state();doc=V.read(OLD)
 for p,s in doc['sources'].items():assert V.sha(p)==s,p
 for row in doc['new_caps']:
  spec=row['specification'];new=spec['name'];assert {r['wall']for r in row['new_pairs']}==set(names)
  for r in row['new_pairs']:pairs[F.key(new,r['wall'])]=r['relation']
  roots.append([O.frozen_interval(z)for z in row['normal']]);names.append(new)
 assert len(roots)==123
 return G,I,h,g,roots,names,c,pairs,doc,sources+[OLD]
def make(spec,g,roots,names,center,pairs):
 # Use inherited exact finite-angle identities in the projection Gram.
 E=[roots[names.index(n)]for n in spec['incident']]
 M=[]
 for i,(a,na) in enumerate(zip(E,spec['incident'],strict=True)):
  row=[]
  for j,(b,nb) in enumerate(zip(E,spec['incident'],strict=True)):
   if i==j:z=V.rat(2)
   else:
    rel=pairs[F.key(na,nb)]
    z=V.neg(F.cos2(rel['label']))if rel['kind']=='finite'else V.pair(a,g,b)
   row.append(z)
  M.append(row)
 MI=F.matinv(M);b=[V.neg(F.cos2(m))for m in spec['profile']];co=V.mv(MI,b);part=N.combine(co,E);qp=V.dot(b,co)
 cp=[V.pair(center,g,a)for a in E];cc=V.mv(MI,cp);proj=N.combine(cc,E);z=[V.sub(x,y)for x,y in zip(center,proj,strict=True)]
 # Orthogonal-projection identity avoids cancellation from squaring z.
 qz=V.sub(V.pair(center,g,center),V.dot(cp,cc));assert qz[0]>0
 rad=V.div(V.sub(V.rat(2),qp),qz);assert rad[0]>0 and spec['branch']in(-1,1)
 t=V.sqrt(rad);t=V.neg(t)if spec['branch']<0 else t;r=[V.add(x,V.mul(t,y))for x,y in zip(part,z,strict=True)]
 return r,{'part_norm':V.enc(qp),'complement_norm':V.enc(qz),'radicand':V.enc(rad)}
def run(verify=False):
 G,I,h,g,roots,names,c,pairs,old,sources=old_state();data=V.read(PRE+'_recipes.json');assert data['source']==OLD and data['source_sha256']==V.sha(OLD);fixed=old['H5_face']['fixed_walls'];new=[];affine=[];new_vertices=[]
 for spec in data['recipes']:
  name=spec['name'];assert len(set(spec['incident']))==9;ff=dict(zip(spec['incident'],spec['profile'],strict=True));assert all(ff[n]==2 for n in fixed);r,recipe=make(spec,g,roots,names,c,pairs);cp=V.dot(r,list(map(V.golden,h)));assert cp[0]>0;checks=[]
  for n,w in zip(names,roots,strict=True):
   p=V.pair(r,g,w)
   if n in ff:
    wanted=V.neg(F.cos2(ff[n]));assert p[0]<=wanted[1]and wanted[0]<=p[1];rel={'kind':'finite','label':ff[n]}
   else:assert p[1]<-2*V.SCALE,(name,n,V.enc(p));rel={'kind':'dotted'}
   pairs[F.key(name,n)]=rel;checks.append({'wall':n,'relation':rel,'pairing':V.enc(p)})
  roots.append(r);names.append(name);neighbors=list(ff);stat={'wall':name,'subsets':512,'dotted_excluded':0,'positive_definite':0,'indefinite':0,'affine':0}
  for mask in range(512):
   ss=[name]+[n for i,n in enumerate(neighbors)if mask>>i&1]
   if any(pairs[F.key(a,b)]['kind']=='dotted'for a,b in combinations(ss,2)):stat['dotted_excluded']+=1;continue
   typ,_=F.inertia(ss,pairs);stat[typ]+=1
   if typ=='positive_definite'and len(ss)==9:new_vertices.append(ss)
  affine.append(stat);new.append({'specification':spec,'projection':recipe,'normal':list(map(V.enc,r)),'center_covector_pairing':V.enc(cp),'new_pairs':checks});print(name,stat,flush=True)
 assert len(roots)==125;facets=old['H5_face']['facets']+[x['specification']['name']for x in new];assert len(facets)==26;census=Counter();vertices=[]
 for ss in combinations(facets,5):
  census['subsets']+=1
  if any(pairs[F.key(a,b)]['kind']=='dotted'for a,b in combinations(ss,2)):census['dotted_excluded']+=1;continue
  typ,_=F.inertia(ss,pairs);census[typ]+=1
  if typ=='positive_definite':vertices.append(F.vertex(fixed+list(ss),roots,names,c,g,h,pairs))
 edges=defaultdict(list)
 for i,v in enumerate(vertices):
  for ee in combinations([n for n in v['incident']if n not in fixed],4):edges[tuple(sorted(ee))].append(i)
 assert all(len(v)in(1,2)for v in edges.values());graph=[set()for _ in vertices]
 for e in edges.values():
  if len(e)==2:graph[e[0]].add(e[1]);graph[e[1]].add(e[0])
 seen={0};todo=[0]
 while todo:
  for j in graph[todo.pop()]-seen:seen.add(j);todo.append(j)
 assert len(seen)==len(vertices)
 # Complete face-facet and finite-vertex test for the proposed H2 polygon.
 common=data['compact_polygon_fixed_walls'];assert len(set(common))==7 and set(fixed)<set(common);assert F.inertia(common,pairs)[0]=='positive_definite';polygon_facets=[];excluded=[]
 for n in names:
  if n in common:continue
  ss=common+[n]
  if any(pairs[F.key(a,b)]['kind']=='dotted'for a,b in combinations(ss,2)):excluded.append({'wall':n,'reason':'dotted_pair_in_incidence_gram'});continue
  typ,piv=F.inertia(ss,pairs)
  if typ=='positive_definite':polygon_facets.append(n)
  else:assert typ=='indefinite';excluded.append({'wall':n,'reason':'indefinite_incidence_gram','pivots':piv})
 assert set(polygon_facets)=={'seed6','wall17','roofPolar','roof5','endpoint0','linked0'};pv=[];incidence=defaultdict(list)
 for a,b in combinations(polygon_facets,2):
  ss=common+[a,b]
  if any(pairs[F.key(x,y)]['kind']=='dotted'for x,y in combinations(ss,2)):continue
  typ,_=F.inertia(ss,pairs)
  if typ=='positive_definite':pv.append(F.vertex(ss,roots,names,c,g,h,pairs));incidence[a].append(len(pv)-1);incidence[b].append(len(pv)-1)
 assert len(pv)==6 and set(incidence)==set(polygon_facets)and all(len(v)==2 for v in incidence.values());pg=[set()for _ in pv]
 for e in incidence.values():pg[e[0]].add(e[1]);pg[e[1]].add(e[0])
 seen={0};todo=[0]
 while todo:
  for j in pg[todo.pop()]-seen:seen.add(j);todo.append(j)
 assert len(seen)==6
 # All six actual facet intervals have two finite endpoints, and form one
 # closed boundary cycle of a convex H2 intersection; hence this face is
 # a compact hexagon. No intrinsic Coxeter-angle assertion is inferred.
 oldy=tuple(map(K.parse,old['strict_escape']['point']));cut=V.dot(roots[123],list(map(V.golden,oldy)));assert cut[1]<0
 y=tuple(map(K.parse,data['strict_escape_proposal']));q=K.pair(y,I,y);future=K.pair(y,I,h);assert K.sign(q)>0 and K.sign(future)<0;ps=[V.dot(r,list(map(V.golden,y)))for r in roots];assert all(v[0]>0 for v in ps)
 x=[V.golden(K.neg(t))for t in K.mv(I,y)];EE=[roots[names.index(n)]for n in fixed];cor=[V.div(V.pair(x,g,e),V.rat(2))for e in EE];xx=[V.sub(a,b)for a,b in zip(x,N.combine(cor,EE),strict=True)];fq=V.pair(xx,g,xx);assert fq[0]>0 and V.dot(xx,list(map(V.golden,h)))[0]>0;fp=[{'wall':n,'pairing':V.enc(V.pair(xx,g,r))}for n,r in zip(names,roots,strict=True)if n not in fixed];assert all(int(p['pairing'][1])<0 for p in fp)
 out={'status':'PASS_EXACT_125_COHERENT_ENDPOINTS_CLOSE_COMPACT_H2_HEXAGON_BUT_H5_ESCAPES','sources':{p:V.sha(p)for p in sources+[PRE+'_recipes.json','scripts/'+Path(__file__).name,'scripts/'+Path(O.__file__).name,'scripts/'+Path(F.__file__).name]},'facets':125,'new_pair_count':sum(len(r['new_pairs'])for r in new),'new_norms':[2,2],'new_caps':new,'positive_affine_census':affine,'actual_new_vertices':[F.vertex(ss,roots,names,c,g,h,pairs)for ss in new_vertices],'old_strict_escape_first_cap_covector_pairing':V.enc(cut),'H5_face':{'fixed_walls':fixed,'facets':facets,'vertex_census':dict(census),'finite_vertices':vertices,'finite_components':[list(range(len(vertices)))],'finite_edges':sum(len(v)==2 for v in edges.values()),'open_edges':sum(len(v)==1 for v in edges.values()),'edges':[{'walls':list(e),'endpoints':v}for e,v in sorted(edges.items())],'strict_escape':{'point':list(map(V.enc,xx)),'norm':V.enc(fq),'nonincident':fp}},'compact_H2_hexagon':{'fixed_walls':common,'facets':polygon_facets,'excluded_wall_incidence_tests':excluded,'finite_vertices':pv,'sides':[{'wall':n,'endpoints':v}for n,v in incidence.items()],'scope':'Actual compact six-sided face. Intrinsic angles are not asserted to be Coxeter angles.'},'strict_escape':{'point':data['strict_escape_proposal'],'norm':K.enc(q),'future_center_pairing':K.enc(future),'all_wall_covector_pairings':list(map(V.enc,ps))},'scope':'Two norm-two caps retain123 walls, cut the previous gap, and close an actual compact H2 face by metric coupling. Complete finite-neighbor affine and literal H5 finite-vertex gates pass. The H5 and H9 chambers still have strict spacelike escapes; no finite-volume or all-completion claim.'}
 if verify:assert out==V.read(PRE+'_verify.json')
 else:(V.ROOT/(PRE+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],'H5',len(vertices),out['H5_face']['finite_edges'],out['H5_face']['open_edges'],flush=True);return out
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');run(a.parse_args().verify)
