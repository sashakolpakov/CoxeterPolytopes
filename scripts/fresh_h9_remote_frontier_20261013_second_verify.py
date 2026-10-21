#!/usr/bin/env python3
"""Exact H3 endpoint cap retaining the preceding 105-wall chamber."""
import argparse,json,sys
from itertools import combinations
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_remote_frontier_20261013_verify as F
V=F.V;K=F.K;N=F.N;C=F.C
PREFIX=F.PREFIX+'_second';NEW='remoteH3'
EDGE=['seed2','seed3','seed5','seed6','seed7','wall22','bridge0','bridge1']
PROFILE={n:3 if n=='seed5' else 5 if n=='bridge1' else 2 for n in EDGE};PROFILE[F.ANCHOR]=3
N.memoize_dependencies()

def old_state():
 d=V.read(F.PREFIX+'_caps.json');cert=F.check(d);assert cert==V.read(F.PREFIX+'_verify.json')
 G,I,h,g,old,names,center,pairs,_=F.old_state();r,_=F.make(d['branch'],g,old,names,center,pairs)
 for n in names:pairs['|'.join(sorted((F.NEW,n))) ]={'kind':'finite','label':F.PROFILE[n]}if n in F.PROFILE else{'kind':'dotted'}
 return G,I,h,g,old+[r],names+[F.NEW],center,pairs,cert

def make(branch,g,old,names,center,pairs):
 A=[old[names.index(n)]for n in EDGE]
 assert C.gram(EDGE,pairs)==[[K.scale(K.O,2)if i==j else K.Z for j in range(8)]for i in range(8)]
 anchor=old[names.index(F.ANCHOR)];wall=old[names.index('wall22')]
 assert all(pairs['|'.join(sorted((F.ANCHOR,n)))]=={'kind':'finite','label':2}for n in EDGE if n!='wall22')
 d=V.pair(anchor,g,wall);assert d[1]<-2*V.SCALE
 b=[V.sub(x,V.div(V.mul(d,y),V.rat(2)))for x,y in zip(anchor,wall,strict=True)]
 qb=V.sub(V.rat(2),V.div(V.mul(d,d),V.rat(2)));assert qb[1]<0
 phi=V.golden(K.parse({'a':'0','b':'1'}))
 base=[V.div(V.neg(V.add(x,V.mul(phi,y))),V.rat(2))for x,y in zip(A[2],A[7],strict=True)]
 qbase=V.div(V.add(V.rat(2),phi),V.rat(2))
 u=[V.sub(x,V.div(y,qb))for x,y in zip(base,b,strict=True)];qu=V.add(qbase,V.div(V.rat(1),qb))
 lam=[V.div(V.pair(center,g,a),V.rat(2))for a in A];p=N.combine(lam,A);z=[V.sub(x,y)for x,y in zip(center,p,strict=True)]
 lm=V.div(V.pair(z,g,b),qb);z=[V.sub(x,V.mul(lm,y))for x,y in zip(z,b,strict=True)];qz=V.pair(z,g,z);assert qz[0]>0
 rad=V.div(V.sub(V.rat(2),qu),qz);assert rad[0]>0 and branch in(-1,1)
 root=V.sqrt(rad);root=V.neg(root)if branch<0 else root;r=[V.add(x,V.mul(root,y))for x,y in zip(u,z,strict=True)]
 # b is orthogonal to A, z is orthogonal to A and b, base.anchor=0.
 # These exact identities give r^2=qbase+1/qb+rad*qz=2;
 # (r,seed5)=-1,(r,bridge1)=-phi,(r,anchor)=-1.
 return r,{'dotted_anchor_pairing':V.enc(d),'anchor_complement_norm':V.enc(qb),'base_norm':V.enc(qbase),'particular_norm':V.enc(qu),'last_complement_norm':V.enc(qz),'radicand':V.enc(rad),'branch':branch}

def check(d):
 assert d['source']==F.PREFIX+'_caps.json'and d['source_sha256']==V.sha(d['source'])
 assert d['prior_verify_sha256']==V.sha(F.PREFIX+'_verify.json')
 G,I,h,g,old,names,center,pairs,prior=old_state();r,recipe=make(d['branch'],g,old,names,center,pairs);hi=list(map(V.golden,h))
 assert V.dot(r,hi)[0]>0;assert V.pair(r,g,r)[0]<=2*V.SCALE<=V.pair(r,g,r)[1]
 counts={'right':0,'pi3':0,'pi5':0,'dotted':0};dotted=[]
 for n,a in zip(names,old,strict=True):
  p=V.pair(r,g,a)
  if n in PROFILE:
   m=PROFILE[n];want=V.rat(0)if m==2 else V.rat(-1)if m==3 else V.neg(V.golden(K.parse({'a':'0','b':'1'})))
   assert p[0]<=want[1]and want[0]<=p[1]
   pairs['|'.join(sorted((NEW,n))) ]={'kind':'finite','label':m};counts[{2:'right',3:'pi3',5:'pi5'}[m]]+=1
  else:
   assert p[1]<-2*V.SCALE,(n,p);counts['dotted']+=1;dotted.append({'wall':n,'pairing':V.enc(p)})
   pairs['|'.join(sorted((NEW,n))) ]={'kind':'dotted'}
 assert counts=={'right':6,'pi3':2,'pi5':1,'dotted':96}
 positive=[];stats={'subsets':512,'dotted_excluded':0,'positive_definite':0,'indefinite':0,'affine':0};ff=list(PROFILE)
 for mask in range(512):
  ss=[NEW]+[n for i,n in enumerate(ff)if mask>>i&1]
  if any(pairs['|'.join(sorted((a,b)))]['kind']=='dotted'for a,b in combinations(ss,2)):stats['dotted_excluded']+=1;continue
  typ=N.T.P.T.psd(C.gram(ss,pairs));assert typ!='positive_semidefinite',ss
  if typ=='positive_definite':positive.append(ss);stats['positive_definite']+=1
  else:stats['indefinite']+=1
 assert prior['affine_census']['affine']==prior['inherited_affine_circuits']==0
 roots=old+[r];allnames=names+[NEW]
 vertices=[F.vertex(ss,roots,allnames,center,g,h,pairs)for ss in positive if len(ss)==9];assert len(vertices)==1 and set(vertices[0]['incident'])==set(EDGE+[NEW])
 # Preserve all four actual corners of the previously closed quadrilateral.
 oldvertices=[F.vertex(v['incident'],roots,allnames,center,g,h,pairs)for v in prior['actual_new_vertices']+prior['old_vertex_survival']]
 frontier=[]
 for ss in sorted({tuple(sorted(x))for v in vertices for x in combinations(v['incident'],8)}):
  ends=[];gates=[]
  for n in allnames:
   if n in ss or any(pairs['|'.join(sorted((n,a)))]['kind']=='dotted'for a in ss):continue
   typ,pivots=N.inertia_gate(list(ss)+[n],pairs);gates.append({'wall':n,'type':typ,'pivots':pivots})
   if typ=='positive_definite':ends.append(n)
  assert len(ends)in(1,2);frontier.append({'edge':list(ss),'end_walls':ends,'ninth_wall_gates':gates})
 target=tuple(map(K.parse,V.read(d['source'])['escape']));q=K.pair(target,I,target);future=K.pair(target,I,h);assert K.sign(q)>0 and K.sign(future)<0
 checks=[V.dot(a,list(map(V.golden,target)))for a in old];assert min(x[0]for x in checks)>0
 cut=V.dot(r,list(map(V.golden,target)));assert cut[1]<0
 out={'status':'PASS_EXACT_REMOTE_H3_ENDPOINT','facets':106,'source_sha256':V.sha(PREFIX+'_caps.json'),'recipe':recipe,'norm':2,'finite_profile':PROFILE,'normal_coordinate_enclosures':list(map(V.enc,r)),'new_pair_counts':counts,'dotted_checks':dotted,'strict_center_pairing':V.enc(V.dot(r,hi)),'affine_census':stats,'inherited_affine_circuits':0,'actual_new_vertices':vertices,'quadrilateral_vertex_survival':oldvertices,'new_vertex_frontier':{'finite_edges':sum(len(e['end_walls'])==2 for e in frontier),'open_rays':sum(len(e['end_walls'])==1 for e in frontier),'edges':frontier},'old_escape_checks':list(map(V.enc,checks)),'target_cut':V.enc(cut),'scope':'A new actual H3+6A1 finite endpoint retains all105 preceding walls and cuts their strict escape. The inherited compact quadrilateral survives. All new finite-neighbor subsets are affine-free; no finite-volume claim follows from that fact.'}
 if 'escape'in d:
  y=tuple(map(K.parse,d['escape']));q=K.pair(y,I,y);f=K.pair(y,I,h);assert K.sign(q)>0 and K.sign(f)<0
  pp=[V.dot(a,list(map(V.golden,y)))for a in roots];assert min(x[0]for x in pp)>0
  out['remaining_escape']={'norm':K.enc(q),'future_center_pairing':K.enc(f),'strict_pairings':list(map(V.enc,pp))};out['scope']+=' An exact strict future spacelike escape proves this106-wall chamber still has infinite volume.'
 return out

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=V.read(PREFIX+'_caps.json');out=check(d)
 if a.verify:assert out==V.read(PREFIX+'_verify.json')
 else:(V.ROOT/(PREFIX+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['new_pair_counts'],out['new_vertex_frontier'])
