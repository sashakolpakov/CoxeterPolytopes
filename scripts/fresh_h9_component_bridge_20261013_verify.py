#!/usr/bin/env python3
"""Producer-free exact trivalent endpoint recipes and outward interval replay.

The finite equalities follow from orthogonal projection identities, not from
an interval containing zero. Every interval operation encloses its exact
real algebraic expression; all divisors and radicands have certified signs.
"""
import argparse,json,sys,copy,inspect
from functools import lru_cache
from itertools import combinations
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_prism_completion_20261013_verify as T
V=T.V;K=T.K;C=T.C
PREFIX='data/enumerations/fresh_h9_component_bridge_20261013'
NEW=['trivalentRemote','trivalentLeading']
SPECS=[
 {'edge':['seed1','seed3','seed4','seed5','seed7','wall22','walk2720','second33'],
  'labels':[2,3,2,3,2,2,2,3],'anchor':'wall54','anchor_label':3},
 {'edge':['seed0','seed1','seed2','seed3','seed6','seed7','coupled0','coupled1'],
  'labels':[2,3,2,3,2,2,2,3],'anchor':'wall17','anchor_label':5}]

def memoize_dependencies():
 # Prior checkers reconstruct the same nested predecessors many times.
 # Reuse only identical literal inputs, returning deep copies because the
 # older callers append to normal lists and pair dictionaries.
 for mod in list(sys.modules.values()):
  if not getattr(mod,'__name__','').startswith('scripts.fresh_h9_'):continue
  for name in ('check','reconstruct'):
   fn=getattr(mod,name,None)
   if not callable(fn)or getattr(fn,'_deep_memoized',False):continue
   if len(inspect.signature(fn).parameters)!=1:continue
   def wrap(fn):
    @lru_cache(None)
    def stored(scale,key):return fn(json.loads(key))
    def call(d):return copy.deepcopy(stored(V.SCALE,json.dumps(d,sort_keys=True)))
    call._deep_memoized=True
    return call
   setattr(mod,name,wrap(fn))
memoize_dependencies()

def combine(c,roots):
 return[V.total(V.mul(a,r[i])for a,r in zip(c,roots,strict=True))for i in range(10)]

def inertia_gate(ss,pairs):
 # Old H4 endpoints can contain label10. Certified nonzero LDL pivots
 # classify them without pretending their Gram lies in Q(phi).
 vals={2:V.rat(0),3:V.rat(-1),5:V.neg(V.PHI),10:V.neg(V.sqrt(V.add(V.rat(2),V.PHI)))}
 M=[[V.rat(2)if a==b else vals[pairs['|'.join(sorted((a,b)))]['label']]for b in ss]for a in ss]
 pivots=[]
 for j in range(len(ss)):
  p=M[j][j];assert p[0]>0 or p[1]<0,('unresolved LDL pivot',ss,j,p)
  pivots.append(V.enc(p))
  if p[1]<0:return 'indefinite',pivots
  for a in range(j+1,len(ss)):
   for b in range(j+1,len(ss)):M[a][b]=V.sub(M[a][b],V.div(V.mul(M[a][j],M[j][b]),p))
 return 'positive_definite',pivots

def old_state():
 d=V.read(T.PREFIX+'_caps.json');prior=T.check(d)
 assert prior==V.read(T.PREFIX+'_verify.json')
 G,I,h,g,old,names,caps,center,pairs,_=T.reconstruct(d)
 for row in prior['new_normals']:
  for n in names:pairs['|'.join(sorted((row['name'],n)))]=row['finite_profile'].get(n,{'kind':'dotted'})
 for row in prior['new_cross_pairs']:pairs['|'.join(sorted(row['walls']))]=row['relation']
 assert len(pairs)==5151
 return G,I,h,g,old+caps,names+T.NEW,center,pairs,prior

def make(spec,branch,g,old,names,center,pairs):
 ss=spec['edge'];assert len(ss)==len(set(ss))==8 and spec in SPECS
 M=C.gram(ss,pairs);assert T.P.T.psd(M)=='positive_definite'
 Mi=K.inverse(M);p=[K.Z if m==2 else K.neg(K.O)for m in spec['labels']]
 co=K.mv(Mi,p);qb0=K.dot(p,co);assert K.sign(K.sub(K.scale(K.O,2),qb0))>0
 A=[old[names.index(n)]for n in ss];mii=[[V.golden(x)for x in row]for row in Mi]
 def project(v):
  z=combine(V.mv(mii,[V.pair(v,g,a)for a in A]),A)
  return[V.sub(x,y)for x,y in zip(v,z,strict=True)]
 base=combine(list(map(V.golden,co)),A)
 a=old[names.index(spec['anchor'])];ap=[V.pair(a,g,r)for r in A]
 b=project(a);qb=V.sub(V.rat(2),V.dot(ap,V.mv(mii,ap)))
 assert qb[0]>0 or qb[1]<0
 finite=K.neg(K.O if spec['anchor_label']==3 else K.P)
 gamma=V.sub(V.golden(finite),V.dot(list(map(V.golden,co)),ap))
 lam=V.div(gamma,qb)
 particular=[V.add(x,V.mul(lam,y))for x,y in zip(base,b,strict=True)]
 qp=V.add(V.golden(qb0),V.div(V.mul(gamma,gamma),qb))
 z=project(center);mu=V.div(V.pair(z,g,b),qb)
 z=[V.sub(x,V.mul(mu,y))for x,y in zip(z,b,strict=True)]
 qz=V.pair(z,g,z);assert qz[0]>0 or qz[1]<0
 rad=V.div(V.sub(V.rat(2),qp),qz);assert rad[0]>0
 assert branch in(-1,1);root=V.sqrt(rad)
 if branch<0:root=V.neg(root)
 r=[V.add(x,V.mul(root,y))for x,y in zip(particular,z,strict=True)]
 # Exact identities: base has prescribed edge pairings; b,z lie in E^perp;
 # (b,z)=0; (r,a)=finite; r^2=qp+rad*qz=2. Rank(E)=8,
 # qb*qz != 0 prove all ten defining directions independent.
 rec={'edge_gram':[[K.enc(x)for x in row]for row in M],
      'base_coefficients':list(map(K.enc,co)),'base_norm':K.enc(qb0),
      'anchor_complement_norm':V.enc(qb),'particular_norm':V.enc(qp),
      'last_complement_norm':V.enc(qz),'radicand':V.enc(rad),'branch':branch}
 return r,rec

def check(d):
 assert d['source']==T.PREFIX+'_caps.json' and d['source_sha256']==V.sha(d['source'])
 assert d['prior_verify_sha256']==V.sha(T.PREFIX+'_verify.json')
 G,I,h,g,old,names,center,pairs,prior=old_state();hi=list(map(V.golden,h))
 assert d['specifications']==SPECS and len(d['branches'])==2
 caps=[];rows=[];counts={'right':0,'pi3':0,'pi5':0,'dotted':0}
 for spec,branch,new in zip(SPECS,d['branches'],NEW,strict=True):
  r,recipe=make(spec,branch,g,old,names,center,pairs);caps.append(r)
  profile=dict(zip(spec['edge'],spec['labels'],strict=True));profile[spec['anchor']]=spec['anchor_label']
  assert V.dot(r,hi)[0]>0
  dotted=[]
  for n,a in zip(names,old,strict=True):
   pp=V.pair(r,g,a)
   if n in profile:
    m=profile[n];want=V.rat(0)if m==2 else V.golden(K.neg(K.O if m==3 else K.P))
    assert pp[0]<=want[1]and want[0]<=pp[1] # Sanity only; identity proved above.
    counts[{2:'right',3:'pi3',5:'pi5'}[m]]+=1
    pairs['|'.join(sorted((new,n))) ]={'kind':'finite','label':m}
   else:
    assert pp[1]<-2*V.SCALE,(new,n,pp);counts['dotted']+=1
    dotted.append({'wall':n,'pairing':V.enc(pp)});pairs['|'.join(sorted((new,n))) ]={'kind':'dotted'}
  assert len(dotted)==93
  rows.append({'name':new,'norm':2,'recipe':recipe,'finite_profile':profile,
               'center_pairing':V.enc(V.dot(r,hi)),'dotted_checks':dotted})
 cross=V.pair(caps[0],g,caps[1]);assert cross[1]<-2*V.SCALE
 pairs['|'.join(sorted(NEW))]={'kind':'dotted'};counts['dotted']+=1
 allroots=old+caps;allnames=names+NEW;positive=[];stats={'subsets':0,'dotted_excluded':0,'positive_definite':0,'indefinite':0,'affine':0}
 for row in rows:
  finite=list(row['finite_profile'])
  for mask in range(1<<len(finite)):
   ss=[row['name']]+[n for i,n in enumerate(finite)if mask>>i&1];stats['subsets']+=1
   if any(pairs['|'.join(sorted((a,b)))]['kind']=='dotted'for a,b in combinations(ss,2)):
    stats['dotted_excluded']+=1;continue
   typ=T.P.T.psd(C.gram(ss,pairs));assert typ!='positive_semidefinite',ss
   stats[typ if typ=='positive_definite'else'indefinite']+=1
   if typ=='positive_definite':positive.append(ss)
 vertices=[]
 for ss in positive:
  if len(ss)!=9:continue
  M=C.gram(ss,pairs);Mi=K.inverse(M);rs=[allroots[allnames.index(n)]for n in ss]
  assert all(K.sign(x)>=0 for row in Mi for x in row)
  lam=V.mv([[V.golden(x)for x in row]for row in Mi],[V.pair(center,g,r)for r in rs]);assert all(x[1]<0 for x in lam)
  pr=combine(lam,rs);point=[V.sub(x,y)for x,y in zip(center,pr,strict=True)]
  q=V.pair(point,g,point);assert q[1]<0 and V.dot(point,hi)[0]>0
  pp=[]
  for n,r in zip(allnames,allroots,strict=True):
   if n in ss:continue
   x=V.pair(point,g,r);assert x[1]<0;pp.append({'wall':n,'pairing':V.enc(x)})
  vertices.append({'incident':ss,'gram':[[K.enc(x)for x in row]for row in M],
                   'point':list(map(V.enc,point)),'norm':V.enc(q),'nonincident':pp})
 # Every edge from a new vertex is checked against every retained wall.
 frontier=[]
 for ss in sorted({tuple(sorted(s))for v in vertices for s in combinations(v['incident'],8)}):
  ends=[];gates=[]
  for n in allnames:
   if n in ss or any(pairs['|'.join(sorted((n,a)))]['kind']=='dotted'for a in ss):continue
   typ,pivots=inertia_gate(list(ss)+[n],pairs);gates.append({'wall':n,'type':typ,'pivots':pivots})
   if typ=='positive_definite':ends.append(n)
  assert len(ends)in(1,2);frontier.append({'edge':list(ss),'end_walls':ends,'ninth_wall_gates':gates})
 target_source=T.PREFIX+'_orbit_obstruction_point.json'
 assert d['orbit_target_source_sha256']==V.sha(target_source)
 orbit_path=T.PREFIX+'_orbit_obstruction.json';orbit=V.read(orbit_path)
 assert d['prior_orbit_certificate_sha256']==V.sha(orbit_path)
 assert orbit['status']=='PASS_EXACT_FULL_SEED_ORBIT_COMMON_ESCAPE' and orbit['literal_orbit_entries']==13824
 assert orbit['source_sha256']==d['source_sha256']and orbit['point']==V.read(target_source)['point']
 assert all(int(row['smallest_lower_bound_pairing'][0])>0 for row in orbit['orbit_minima'])
 targets=[V.read(d['source'])['escape'],V.read(target_source)['point']];target_checks=[]
 for raw in targets:
  y=tuple(map(K.parse,raw));q=K.pair(y,I,y);f=K.pair(y,I,h)
  assert K.sign(q)>0 and K.sign(f)<0
  pp=[V.dot(r,list(map(V.golden,y)))for r in old];assert min(x[0]for x in pp)>0
  cuts=[V.dot(r,list(map(V.golden,y)))for r in caps];assert min(x[1]for x in cuts)<0
  target_checks.append({'norm':K.enc(q),'future_pairing':K.enc(f),'old_strict_pairings':list(map(V.enc,pp)),'new_pairings':list(map(V.enc,cuts))})
 assert prior['affine_census']['new_affine_circuits']==prior['affine_census']['inherited_affine_circuits']==0
 out={'status':'PASS_EXACT_TRIVALENT_ENDPOINT_PAIR','facets':104,'interval_scale':str(V.SCALE),
      'source_sha256':V.sha(PREFIX+'_caps.json'),'recipes':rows,'new_pair_counts':counts,
      'cross_pairing':V.enc(cross),'affine_census':stats,'inherited_affine_circuits':0,
      'actual_new_finite_vertices':vertices,'new_vertex_frontier':frontier,'target_checks':target_checks,
      'outside_previous_family_dependency':{'certificate':orbit_path,'sha256':V.sha(orbit_path),'literal_images':13824},
      'normal_enclosures':[list(map(V.enc,r))for r in caps],
      'scope':'Two exact trivalent endpoint caps preserve all102 actual facets. Exact finite equalities follow by projection; all other pairs are rigorously dotted. All PSD subsets involving either new wall are censused; the dotted new pair forbids a PSD subset containing both. This is not a compactness or finite-volume claim.'}
 if 'escape'in d:
  y=tuple(map(K.parse,d['escape']));q=K.pair(y,I,y);f=K.pair(y,I,h)
  assert K.sign(q)>0 and K.sign(f)<0
  pp=[V.dot(r,list(map(V.golden,y)))for r in allroots];assert min(x[0]for x in pp)>0
  out['remaining_escape']={'norm':K.enc(q),'future_center_pairing':K.enc(f),'strict_pairings':list(map(V.enc,pp))}
  out['scope']+=' A strict future spacelike escape proves infinite volume.'
 return out

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=V.read(PREFIX+'_caps.json');out=check(d)
 if a.verify:assert out==V.read(PREFIX+'_verify.json')
 else:(V.ROOT/(PREFIX+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['new_pair_counts'],len(out['actual_new_finite_vertices']))
