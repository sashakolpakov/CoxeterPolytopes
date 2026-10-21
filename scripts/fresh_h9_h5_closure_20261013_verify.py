#!/usr/bin/env python3
"""Exact coherent H5 gap cuts and one coupled finite endpoint; no scout imports."""
import argparse,json,sys
from pathlib import Path
from itertools import combinations
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_higher_face_20261013_truncation_verify as T
F=T.F;V=F.V;K=F.K;N=F.N;PRE='data/enumerations/fresh_h9_h5_closure_20261013';OLD=F.PREFIX
BASECOS=F.cos2
FROZEN_SCALE=10**80
def frozen_interval(z):
 # Earlier certificates enclose exact, source-defined normals with denominator
 # 10**80. Preserve those same rational bounds during higher-precision reviews.
 assert V.SCALE%FROZEN_SCALE==0
 k=V.SCALE//FROZEN_SCALE
 return tuple(int(a)*k for a in z)
def cos2(m):
 if m==6:return V.sqrt(V.rat(3))
 if m==10:return V.sqrt(V.add(V.rat(2),V.golden(K.parse({'a':'0','b':'1'}))))
 return BASECOS(m)
F.cos2=cos2
def old_state():
 S=F.S;d=V.read(S.PREFIX+'_caps.json');G,I,h,g,roots,names,c,pairs,_=S.old_state();r,_=S.make(d['branch'],g,roots,names,c,pairs)
 for n in names:pairs[F.key(S.NEW,n)]={'kind':'finite','label':S.PROFILE[n]}if n in S.PROFILE else{'kind':'dotted'}
 roots=roots+[r];names=names+[S.NEW]
 sources=[S.PREFIX+'_caps.json',S.PREFIX+'_verify.json']
 for suffix,added in[('_verify.json',[f'roof{i}'for i in range(7)]),('_truncation_verify.json',['adaptive0','adaptive1','secondAdaptive0']),('_final_verify.json',['roofPolar','roofMetric'])]:
  path=OLD+suffix;doc=V.read(path);sources.append(path)
  for s,sha in doc['sources'].items():assert V.sha(s)==sha
  assert len(doc['normals'])==len(added)==len(doc['new_pairs'])
  for new,rr,checks in zip(added,doc['normals'],doc['new_pairs'],strict=True):
   assert {q['wall']for q in checks}==set(names)
   for q in checks:pairs[F.key(new,q['wall'])]=q['kind']
   roots.append([frozen_interval(z)for z in rr]);names.append(new)
 assert len(roots)==118 and doc['status']=='PASS_EXACT_LINKED_H5_ROOF_WITH_TRUNCATIONS'
 return G,I,h,g,roots,names,c,pairs,doc,sources
def make(spec,g,roots,names,center):
 E=[roots[names.index(n)]for n in spec['incident']];M=[[V.pair(a,g,b)for b in E]for a in E];MI=F.matinv(M);b=[V.neg(cos2(m))for m in spec['profile']];co=V.mv(MI,b);part=N.combine(co,E);qp=V.dot(b,co)
 proj=N.combine(V.mv(MI,[V.pair(center,g,a)for a in E]),E);z=[V.sub(x,y)for x,y in zip(center,proj,strict=True)];qz=V.pair(z,g,z);assert qz[0]>0;rad=V.div(V.sub(V.rat(2),qp),qz);assert rad[0]>0 and spec['branch']in(-1,1)
 t=V.sqrt(rad);t=V.neg(t)if spec['branch']<0 else t;r=[V.add(x,V.mul(t,y))for x,y in zip(part,z,strict=True)]
 # The inverse defines exact projection expressions: part.E=b, z.E=0,
 # part.z=0, and r²=part²+rad*z²=2. Strict interval pivot/radicand
 # bounds certify that all defining expressions exist at the chosen embedding.
 return r,{'part_norm':V.enc(qp),'complement_norm':V.enc(qz),'radicand':V.enc(rad)}
def run(verify=False):
 G,I,h,g,roots,names,center,pairs,prior,sources=old_state();data=V.read(PRE+'_recipes.json');fixed=prior['H5_face']['fixed_walls'];added=[];affine=[];newverts=[]
 for spec in data['recipes']:
  new=spec['name'];assert len(spec['incident'])==len(spec['profile'])==9 and len(set(spec['incident']))==9;r,recipe=make(spec,g,roots,names,center);ff=dict(zip(spec['incident'],spec['profile'],strict=True));assert all(ff[n]==2 for n in fixed);cp=V.dot(r,list(map(V.golden,h)));assert cp[0]>0;checks=[]
  for n,a in zip(names,roots,strict=True):
   p=V.pair(r,g,a)
   if n in ff:
    wanted=V.neg(cos2(ff[n]));assert p[0]<=wanted[1]and wanted[0]<=p[1];kind={'kind':'finite','label':ff[n]}
   else:assert p[1]<-2*V.SCALE,(new,n,p);kind={'kind':'dotted'}
   pairs[F.key(new,n)]=kind;checks.append({'wall':n,'relation':kind,'pairing':V.enc(p)})
  roots.append(r);names.append(new);stats={'wall':new,'subsets':512,'dotted_excluded':0,'positive_definite':0,'indefinite':0,'affine':0};neighbors=list(ff)
  for mask in range(512):
   ss=[new]+[n for j,n in enumerate(neighbors)if mask>>j&1]
   if any(pairs[F.key(a,b)]['kind']=='dotted'for a,b in combinations(ss,2)):stats['dotted_excluded']+=1;continue
   typ,_=F.inertia(ss,pairs);stats[typ]+=1
   if len(ss)==9 and typ=='positive_definite':newverts.append(ss)
  affine.append(stats);added.append({'specification':spec,'projection':recipe,'center_covector_pairing':V.enc(cp),'new_pairs':checks,'normal':list(map(V.enc,r))});print(new,stats,flush=True)
 assert len(roots)==123
 facefacets=prior['H5_face']['facets']+[r['name']for r in data['recipes']];assert len(facefacets)==24
 census={'subsets':0,'dotted_excluded':0,'positive_definite':0,'indefinite':0,'affine':0};vertices=[]
 for ss in combinations(facefacets,5):
  census['subsets']+=1
  if any(pairs[F.key(a,b)]['kind']=='dotted'for a,b in combinations(ss,2)):census['dotted_excluded']+=1;continue
  typ,_=F.inertia(ss,pairs);census[typ]+=1
  if typ=='positive_definite':vertices.append(F.vertex(fixed+list(ss),roots,names,center,g,h,pairs))
 edge={}
 for i,v in enumerate(vertices):
  for ee in combinations([n for n in v['incident']if n not in fixed],4):edge.setdefault(tuple(sorted(ee)),[]).append(i)
 assert all(len(e)in(1,2)for e in edge.values());graph=[set()for _ in vertices]
 for vs in edge.values():
  if len(vs)==2:graph[vs[0]].add(vs[1]);graph[vs[1]].add(vs[0])
 comps=[];seen=set()
 for i in range(len(vertices)):
  if i in seen:continue
  component={i};seen.add(i);todo=[i]
  while todo:
   for j in graph[todo.pop()]-seen:seen.add(j);component.add(j);todo.append(j)
  comps.append(sorted(component))
 original=tuple(map(K.parse,prior['strict_escape']['point']));firstcut=V.dot(roots[118],list(map(V.golden,original)));assert firstcut[1]<0
 y=tuple(map(K.parse,data['strict_escape_proposal']));q=K.pair(y,I,y);ft=K.pair(y,I,h);assert K.sign(q)>0 and K.sign(ft)<0;checks=[V.dot(r,list(map(V.golden,y)))for r in roots];assert all(p[0]>0 for p in checks)
 x=[V.golden(K.neg(z))for z in K.mv(I,y)];EE=[roots[names.index(n)]for n in fixed];cor=[V.div(V.pair(x,g,e),V.rat(2))for e in EE];xx=[V.sub(a,b)for a,b in zip(x,N.combine(cor,EE),strict=True)];fq=V.pair(xx,g,xx);assert fq[0]>0 and V.dot(xx,list(map(V.golden,h)))[0]>0;fchecks=[{'wall':n,'pairing':V.enc(V.pair(xx,g,r))}for n,r in zip(names,roots,strict=True)if n not in fixed];assert all(int(p['pairing'][1])<0 for p in fchecks)
 out={'status':'PASS_EXACT_123_COHERENT_H5_GAP_CUTS_AND_COUPLED_ENDPOINT_WITH_STRICT_ESCAPE','sources':{s:V.sha(s)for s in sources+[PRE+'_recipes.json','scripts/'+Path(__file__).name,'scripts/'+Path(F.__file__).name]},'facets':123,'new_norms':[2]*5,'new_caps':added,'new_pair_count':sum(len(r['new_pairs'])for r in added),'new_positive_affine_census':affine,'original_strict_escape_first_cap_covector_pairing':V.enc(firstcut),'actual_new_vertices':[F.vertex(ss,roots,names,center,g,h,pairs)for ss in newverts],'H5_face':{'fixed_walls':fixed,'facets':facefacets,'vertex_census':census,'finite_vertices':vertices,'finite_components':comps,'finite_edges':sum(len(v)==2 for v in edge.values()),'open_edges':sum(len(v)==1 for v in edge.values()),'edges':[{'walls':list(k),'endpoints':v}for k,v in sorted(edge.items())],'strict_spacelike_escape':{'point':list(map(V.enc,xx)),'norm':V.enc(fq),'nonincident':fchecks}},'strict_escape':{'point':data['strict_escape_proposal'],'norm':K.enc(q),'future_center_pairing':K.enc(ft),'all_wall_covector_pairings':list(map(V.enc,checks))},'scope':'Five exactly defined norm-two caps retain all118 walls. Two caps have prescribed finite contacts with closure0, and the last supplies an actual finite endpoint attached to the inherited H5 vertex component. Full finite-neighbor affine tests and literal H5 finite-vertex census are exact. Both H9 and the H5 face still have strict spacelike escapes; neither compactness nor an all-parameter completion obstruction is claimed.'}
 if verify:assert out==V.read(PRE+'_verify.json')
 else:(V.ROOT/(PRE+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],'vertices',len(vertices),'finite/open',out['H5_face']['finite_edges'],out['H5_face']['open_edges'],flush=True)
 return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
