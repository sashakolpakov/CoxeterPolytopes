#!/usr/bin/env python3
"""Exact full-face finite-contact truncations of the linked roof."""
import argparse,json,sys
from itertools import combinations
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_higher_face_20261013_verify as F
V=F.V;K=F.K;N=F.N;PREFIX=F.PREFIX+'_truncation'
SPECS=[{'name':'adaptive0','incident':['seed2','compact1','compact2','compact3','nextH4','crosspair0','roof0','roof1','roof5'],'profile':[2,2,2,2,2,3,2,3,4],'branch':-1},
 {'name':'adaptive1','incident':['seed2','compact1','compact2','compact3','wall17','crosspair0','roof0','roof1','adaptive0'],'profile':[2,2,2,2,5,2,5,3,2],'branch':-1},
 {'name':'secondAdaptive0','incident':['seed2','compact1','compact2','compact3','seed6','crosspair0','roof5','roof6','adaptive0'],'profile':[2,2,2,2,3,2,5,8,8],'branch':-1}]
def make(spec,g,roots,names,center):
 E=[roots[names.index(n)]for n in spec['incident']];M=[[V.pair(a,g,b)for b in E]for a in E];Mi=F.matinv(M);bs=[V.neg(F.cos2(m))for m in spec['profile']];co=V.mv(Mi,bs);part=N.combine(co,E);qp=V.dot(bs,co)
 lam=V.mv(Mi,[V.pair(center,g,a)for a in E]);proj=N.combine(lam,E);z=[V.sub(x,y)for x,y in zip(center,proj,strict=True)];qz=V.pair(z,g,z);assert qz[0]>0
 rad=V.div(V.sub(V.rat(2),qp),qz);assert rad[0]>0 and spec['branch']in(-1,1);t=V.sqrt(rad);t=V.neg(t)if spec['branch']<0 else t
 r=[V.add(x,V.mul(t,y))for x,y in zip(part,z,strict=True)]
 # Matrix inversion defines exact rational expressions. Nonzero pivot bounds
 # certify invertibility, so part.E=bs, z.E=0, part.z=0. Hence r²=2 exactly.
 return r,{'incident':spec['incident'],'profile':spec['profile'],'part_norm':V.enc(qp),'complement_norm':V.enc(qz),'radicand':V.enc(rad),'branch':spec['branch']}

def build():
 prior,(G,I,h,g,roots,names,center,pairs)=F.build();assert prior==V.read(F.PREFIX+'_verify.json');newchecks=[];recipes=[];stats=[];sets=[]
 for spec in SPECS:
  new=spec['name'];r,rec=make(spec,g,roots,names,center);ff=dict(zip(spec['incident'],spec['profile'],strict=True));assert V.dot(r,list(map(V.golden,h)))[0]>0;checks=[]
  for n,a in zip(names,roots,strict=True):
   p=V.pair(r,g,a)
   if n in ff:
    want=V.neg(F.cos2(ff[n]));assert p[0]<=want[1]and want[0]<=p[1];pairs[F.key(new,n)]={'kind':'finite','label':ff[n]}
   else:assert p[1]<-2*V.SCALE,(new,n,p);pairs[F.key(new,n)]={'kind':'dotted'}
   checks.append({'wall':n,'kind':pairs[F.key(new,n)],'pairing':V.enc(p)})
  roots.append(r);names.append(new);recipes.append(rec);newchecks.append(checks);st={'name':new,'subsets':512,'dotted_excluded':0,'positive_definite':0,'indefinite':0,'affine':0}
  for mask in range(512):
   ss=[new]+[n for j,n in enumerate(spec['incident'])if mask>>j&1]
   if any(pairs[F.key(a,b)]['kind']=='dotted'for a,b in combinations(ss,2)):st['dotted_excluded']+=1;continue
   typ,_=F.inertia(ss,pairs);st[typ]+=1
   if len(ss)==9 and typ=='positive_definite':sets.append(ss)
  stats.append(st);print(new,st,flush=True)
 vertices=[F.vertex(ss,roots,names,center,g,h,pairs)for ss in sets];survival=[F.vertex(v['incident'],roots,names,center,g,h,pairs)for v in prior['actual_new_vertices']]
 out={'status':'PASS_EXACT_FULL_FACE_ROOF_TRUNCATIONS','sources':{F.PREFIX+'_verify.json':V.sha(F.PREFIX+'_verify.json'),'scripts/'+Path(__file__).name:V.sha('scripts/'+Path(__file__).name)},'facets':len(roots),'new_norms':[2]*len(SPECS),'recipes':recipes,'new_pairs':newchecks,'affine_census':stats,'actual_new_vertices':vertices,'roof_vertex_survival':survival,'normals':[list(map(V.enc,r))for r in roots[-len(SPECS):]],'scope':'Coherent finite-contact truncations retain the complete seven-cap roof and its12 certified new vertices. No compactness inference from these local checks.'}
 return out,(G,I,h,g,roots,names,center,pairs)
def run(verify=False):
 out,_=build()
 if verify:assert out==V.read(PREFIX+'_verify.json')
 else:(V.ROOT/(PREFIX+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['facets'],len(out['actual_new_vertices']))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
