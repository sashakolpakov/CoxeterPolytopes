#!/usr/bin/env python3
"""Exact fourth-side closure of a different inherited H2 face."""
import argparse,json,sys
from itertools import combinations
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_component_bridge_20261013_verify as N
V=N.V;K=N.K;C=N.C;PREFIX='data/enumerations/fresh_h9_remote_frontier_20261013'
EDGE=['seed1','seed2','seed5','seed6','seed7','wall22','bridge0','bridge1']
COMMON=[n for n in EDGE if n!='wall22'];ANCHOR='nextH4';NEW='remoteQuadrilateral'
PROFILE={n:3 if n in('seed1','bridge1')else 2 for n in EDGE};PROFILE[ANCHOR]=3
N.memoize_dependencies()

def old_state():
 d=V.read(N.PREFIX+'_caps.json');cert=N.check(d);assert cert==V.read(N.PREFIX+'_verify.json')
 G,I,h,g,old,names,center,pairs,_=N.old_state()
 caps=[N.make(s,b,g,old,names,center,pairs)[0]for s,b in zip(N.SPECS,d['branches'],strict=True)]
 for spec,new in zip(N.SPECS,N.NEW,strict=True):
  ff=dict(zip(spec['edge'],spec['labels'],strict=True));ff[spec['anchor']]=spec['anchor_label']
  for n in names:pairs['|'.join(sorted((new,n))) ]={'kind':'finite','label':ff[n]}if n in ff else{'kind':'dotted'}
 pairs['|'.join(sorted(N.NEW))]={'kind':'dotted'}
 return G,I,h,g,old+caps,names+N.NEW,center,pairs,cert

def make(branch,g,old,names,center,pairs):
 A=[old[names.index(n)]for n in EDGE];M=C.gram(EDGE,pairs)
 assert M==[[K.scale(K.O,2)if i==j else K.Z for j in range(8)]for i in range(8)]
 anchor=old[names.index(ANCHOR)];wall=old[names.index('wall22')]
 assert all(pairs['|'.join(sorted((ANCHOR,n)))]=={'kind':'finite','label':2}for n in COMMON)
 d=V.pair(anchor,g,wall);assert d[1]<-2*V.SCALE
 b=[V.sub(x,V.div(V.mul(d,y),V.rat(2)))for x,y in zip(anchor,wall,strict=True)]
 qb=V.sub(V.rat(2),V.div(V.mul(d,d),V.rat(2)));assert qb[1]<0
 base=[V.div(V.neg(V.add(x,y)),V.rat(2))for x,y in zip(A[0],A[7],strict=True)]
 # base^2=1, base.anchor=0, b perpendicular to all8 edge normals.
 u=[V.sub(x,V.div(y,qb))for x,y in zip(base,b,strict=True)];qu=V.add(V.rat(1),V.div(V.rat(1),qb))
 lambdas=[V.div(V.pair(center,g,a),V.rat(2))for a in A];p=N.combine(lambdas,A);z=[V.sub(x,y)for x,y in zip(center,p,strict=True)]
 lam=V.div(V.pair(z,g,b),qb);z=[V.sub(x,V.mul(lam,y))for x,y in zip(z,b,strict=True)];qz=V.pair(z,g,z);assert qz[0]>0
 rad=V.div(V.sub(V.rat(2),qu),qz);assert rad[0]>0;assert branch in(-1,1)
 root=V.sqrt(rad);root=V.neg(root)if branch<0 else root;r=[V.add(x,V.mul(root,y))for x,y in zip(u,z,strict=True)]
 # Exact orthogonal decomposition: r^2=qu+rad*qz=2;
 # edge pairings are PROFILE and r.anchor=-1. No interval-zero inference.
 return r,{'dotted_anchor_pairing':V.enc(d),'anchor_complement_norm':V.enc(qb),'particular_norm':V.enc(qu),'last_complement_norm':V.enc(qz),'radicand':V.enc(rad),'branch':branch}

def vertex(ss,old,names,center,g,h,pairs):
 M=C.gram(ss,pairs);assert N.T.P.T.psd(M)=='positive_definite';Mi=K.inverse(M);assert all(K.sign(x)>=0 for row in Mi for x in row)
 rs=[old[names.index(n)]for n in ss];lam=V.mv([[V.golden(x)for x in row]for row in Mi],[V.pair(center,g,r)for r in rs]);assert all(x[1]<0 for x in lam)
 p=N.combine(lam,rs);point=[V.sub(x,y)for x,y in zip(center,p,strict=True)];q=V.pair(point,g,point);assert q[1]<0 and V.dot(point,list(map(V.golden,h)))[0]>0
 checks=[]
 for n,r in zip(names,old,strict=True):
  if n in ss:continue
  pp=V.pair(point,g,r);assert pp[1]<0;checks.append({'wall':n,'pairing':V.enc(pp)})
 return{'incident':ss,'gram':[[K.enc(x)for x in row]for row in M],'point':list(map(V.enc,point)),'norm':V.enc(q),'nonincident_pairings':checks}

def check(d):
 assert d['source']==N.PREFIX+'_caps.json'and d['source_sha256']==V.sha(d['source'])
 assert d['prior_verify_sha256']==V.sha(N.PREFIX+'_verify.json')
 G,I,h,g,old,names,center,pairs,prior=old_state();r,recipe=make(d['branch'],g,old,names,center,pairs);hi=list(map(V.golden,h))
 assert V.dot(r,hi)[0]>0;assert V.pair(r,g,r)[0]<=2*V.SCALE<=V.pair(r,g,r)[1]
 dotted=[];counts={'right':0,'pi3':0,'dotted':0}
 for n,a in zip(names,old,strict=True):
  p=V.pair(r,g,a)
  if n in PROFILE:
   m=PROFILE[n];want=0 if m==2 else -V.SCALE;assert p[0]<=want<=p[1]
   pairs['|'.join(sorted((NEW,n))) ]={'kind':'finite','label':m};counts['right'if m==2 else'pi3']+=1
  else:
   assert p[1]<-2*V.SCALE,(n,p);counts['dotted']+=1;dotted.append({'wall':n,'pairing':V.enc(p)})
   pairs['|'.join(sorted((NEW,n))) ]={'kind':'dotted'}
 assert counts=={'right':6,'pi3':3,'dotted':95}
 positive=[];stats={'subsets':512,'dotted_excluded':0,'positive_definite':0,'indefinite':0,'affine':0}
 ff=list(PROFILE)
 for mask in range(512):
  ss=[NEW]+[n for i,n in enumerate(ff)if mask>>i&1]
  if any(pairs['|'.join(sorted((a,b)))]['kind']=='dotted'for a,b in combinations(ss,2)):stats['dotted_excluded']+=1;continue
  typ=N.T.P.T.psd(C.gram(ss,pairs));assert typ!='positive_semidefinite',ss
  if typ=='positive_definite':positive.append(ss);stats['positive_definite']+=1
  else:stats['indefinite']+=1
 assert prior['affine_census']['affine']==prior['inherited_affine_circuits']==0
 allroots=old+[r];allnames=names+[NEW]
 vertices=[vertex(ss,allroots,allnames,center,g,h,pairs)for ss in positive if len(ss)==9];assert len(vertices)==2
 # The old and new endpoints form an actual closed quadrilateral in the
 # common rank7 face. All four actual points obey every other wall strictly.
 oldsets=[COMMON+['seed3',n]for n in('wall22',ANCHOR)]
 oldvertices=[vertex(ss,allroots,allnames,center,g,h,pairs)for ss in oldsets]
 assert {frozenset(v['incident'])for v in vertices}=={frozenset(COMMON+[NEW,n])for n in('wall22',ANCHOR)}
 assert C.gram(COMMON,pairs)==[[K.scale(K.O,2)if i==j else K.Z for j in range(7)]for i in range(7)]
 for n in('seed3','wall22',ANCHOR):
  assert all(pairs['|'.join(sorted((n,a)))]=={'kind':'finite','label':2}for a in COMMON)
 assert sum(PROFILE[n]==3 for n in COMMON)==2 # New projected norm2-2/2=1.
 assert pairs['|'.join(sorted(('seed3','wall22')))]==pairs['|'.join(sorted(('seed3',ANCHOR)))]=={'kind':'finite','label':2}
 assert pairs['|'.join(sorted((NEW,'seed3')))]['kind']=='dotted'
 # Intrinsic consecutive pairings are 0,0,-1,0 and norms2,2,1,2:
 # angles pi/2,pi/2,pi/4,pi/2. Nonadjacent planes are disjoint.
 frontier=[]
 for ss in sorted({tuple(sorted(x))for v in oldvertices+vertices for x in combinations(v['incident'],8)}):
  ends=[];gates=[]
  for n in allnames:
   if n in ss or any(pairs['|'.join(sorted((n,a)))]['kind']=='dotted'for a in ss):continue
   typ,pivots=N.inertia_gate(list(ss)+[n],pairs);gates.append({'wall':n,'type':typ,'pivots':pivots})
   if typ=='positive_definite':ends.append(n)
  assert len(ends)in(1,2);frontier.append({'edge':list(ss),'end_walls':ends,'ninth_wall_gates':gates})
 cycle=[e for e in frontier if set(COMMON)<=set(e['edge'])and len(e['end_walls'])==2];assert len(cycle)==4
 target=tuple(map(K.parse,V.read(d['source'])['escape']));q=K.pair(target,I,target);f=K.pair(target,I,h);assert K.sign(q)>0 and K.sign(f)<0
 oldchecks=[V.dot(a,list(map(V.golden,target)))for a in old];assert min(x[0]for x in oldchecks)>0
 cut=V.dot(r,list(map(V.golden,target)));assert cut[1]<0
 out={'status':'PASS_EXACT_REMOTE_QUADRILATERAL_FACE','facets':105,'source_sha256':V.sha(PREFIX+'_caps.json'),
      'recipe':recipe,'norm':2,'finite_profile':PROFILE,'normal_coordinate_enclosures':list(map(V.enc,r)),
      'new_pair_counts':counts,'dotted_checks':dotted,'strict_center_pairing':V.enc(V.dot(r,hi)),
      'affine_census':stats,'inherited_affine_circuits':0,'actual_new_vertices':vertices,'old_vertex_survival':oldvertices,
      'compact_H2_face':{'common':COMMON,'sides_in_cycle':['seed3','wall22',NEW,ANCHOR],
       'projected_norms':[2,2,1,2],'intrinsic_angle_labels':[2,2,4,2],'cycle_edges':cycle},
      'vertex_component_frontier':{'vertices':4,'finite_edges':sum(len(e['end_walls'])==2 for e in frontier),'open_rays':sum(len(e['end_walls'])==1 for e in frontier),'edges':frontier},
      'old_escape_checks':list(map(V.enc,oldchecks)),'target_cut':V.enc(cut),
      'scope':'A genuinely different inherited2-vertex component becomes a compact Lambert quadrilateral face. All104 old walls, including16seed and the compactH4prism, remain actual facets. The entire new finite-neighbor PSD census is affine-free. No claim of whole-chamber finite volume follows from this face.'}
 if 'escape'in d:
  y=tuple(map(K.parse,d['escape']));q=K.pair(y,I,y);f=K.pair(y,I,h);assert K.sign(q)>0 and K.sign(f)<0
  pp=[V.dot(a,list(map(V.golden,y)))for a in allroots];assert min(x[0]for x in pp)>0
  out['remaining_escape']={'norm':K.enc(q),'future_center_pairing':K.enc(f),'strict_pairings':list(map(V.enc,pp))};out['scope']+=' A new exact strict future spacelike escape proves infinite volume.'
 return out

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=V.read(PREFIX+'_caps.json');out=check(d)
 if a.verify:assert out==V.read(PREFIX+'_verify.json')
 else:(V.ROOT/(PREFIX+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['new_pair_counts'],out['vertex_component_frontier']['open_rays'])
