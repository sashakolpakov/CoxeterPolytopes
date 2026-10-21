#!/usr/bin/env python3
"""Exact endpoint caps and compact H4 prism inside the retained 94-wall chamber."""
import argparse,json,sys
from pathlib import Path
from itertools import combinations
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_compact_next_20261013_verify as C
V=C.V;K=C.K;P=C.P
PREFIX='data/enumerations/fresh_h9_face_bridge_20261013'
NEW=['facePrism','faceGap']

def table(G):
 names,pairs=C.table(G)
 for i,n in enumerate(C.NEW):
  for old in names:pairs['|'.join(sorted((n,old))) ]={'kind':'finite','label':3 if i==0 and old==C.ANCHOR else 2}if old in C.S+[C.ANCHOR]else{'kind':'dotted'}
  for m in C.NEW[:i]:pairs['|'.join(sorted((n,m))) ]={'kind':'finite','label':2}
 return names+C.NEW,pairs

def lincomb(coeff,roots):
 return [V.total(V.mul(V.golden(a),r[j])for a,r in zip(coeff,roots,strict=True))for j in range(10)]

def reconstruct(d):
 assert d['source']==C.PREFIX+'_caps.json'and d['source_sha256']==V.sha(d['source'])
 prior=V.read(d['source']);cert=C.check(prior)
 G,I,h,g,old,names,caps,center,_=C.reconstruct(prior);old+=caps;names+=C.NEW
 expect,pairs=table(G);assert expect==names
 roots=[];records=[];assert len(d['caps'])==2
 for ci,row in enumerate(d['caps']):
  ss=row['source_walls'];assert len(ss)==len(set(ss))==9
  assert row['complement_wall']==['wall17','compact2'][ci]and row['branch']==-1
  assert {n for n,m in zip(ss,row['labels'],strict=True)if m!=2}==[{'compact0'},{'crosspair0','seed2','compact1'}][ci]
  M=C.gram(ss,pairs);Mi=K.inverse(M);vv=[K.Z if m==2 else K.neg(K.O)for m in row['labels']];assert set(row['labels'])<={2,3}
  rc=[old[names.index(n)]for n in ss];co=K.mv(Mi,vv);w=lincomb(co,rc);qw=K.dot(vv,co)
  an=row['complement_wall'];a=old[names.index(an)];ap=[C.gram([an,n],pairs)[0][1]for n in ss]
  az=K.mv(Mi,ap);azv=lincomb(az,rc);z=[V.sub(x,y)for x,y in zip(a,azv,strict=True)]
  qz=K.sub(K.scale(K.O,2),K.dot(ap,az));assert K.sign(qz)>0
  rad=K.div(K.sub(K.scale(K.O,2),qw),qz);assert K.sign(rad)>0
  assert K.add(qw,K.mul(rad,qz))==K.scale(K.O,2)
  assert row['branch']in(-1,1);lam=V.sqrt(V.golden(rad));lam=V.neg(lam)if row['branch']<0 else lam
  r=[V.add(x,V.mul(lam,y))for x,y in zip(w,z,strict=True)]
  roots.append(r);records.append({'particular_coefficients':list(map(K.enc,co)),'complement_projection_coefficients':list(map(K.enc,az)),'particular_norm':K.enc(qw),'complement_norm':K.enc(qz),'radicand':K.enc(rad),'branch':row['branch']})
 return G,I,h,g,old,names,roots,center,records,cert

def check(d):
 G,I,h,g,old,names,caps,center,params,prior=reconstruct(d);_,pairs=table(G);hi=list(map(V.golden,h));counts={'right':0,'pi3':0,'dotted':0};rows=[]
 for j,(r,row)in enumerate(zip(caps,d['caps'],strict=True)):
  cp=V.dot(r,hi);assert cp[0]>0;finite=dict(zip(row['source_walls'],row['labels'],strict=True));dot=[]
  for a,n in zip(old,names,strict=True):
   pp=V.pair(r,g,a)
   if n in finite:
    lab=finite[n];v=0 if lab==2 else -V.SCALE;assert pp[0]<=v<=pp[1];counts['right'if lab==2 else'pi3']+=1
    pairs['|'.join(sorted((NEW[j],n))) ]={'kind':'finite','label':lab}
   else:
    assert pp[1]<-2*V.SCALE,(NEW[j],n,pp);counts['dotted']+=1;dot.append({'old':n,'pairing':V.enc(pp)})
    pairs['|'.join(sorted((NEW[j],n))) ]={'kind':'dotted'}
  assert len(dot)==85
  rows.append({'name':NEW[j],'norm':2,'center_pairing_dual':V.enc(cp),'construction':params[j],'dotted_checks':dot})
 cross=V.pair(caps[0],g,caps[1]);assert cross[1]<-2*V.SCALE;counts['dotted']+=1;pairs['|'.join(sorted(NEW))]={'kind':'dotted'}
 allroots=old+caps;allnames=names+NEW;positive=[];indefinite=0
 for j,row in enumerate(d['caps']):
  ss=row['source_walls']
  for mask in range(1<<9):
   subset=[NEW[j]]+[n for i,n in enumerate(ss)if mask>>i&1];typ=P.T.psd(C.gram(subset,pairs));assert typ!='positive_semidefinite',subset
   if typ=='positive_definite':positive.append(subset)
   else:indefinite+=1
 vertices=[]
 for ss in positive:
  if len(ss)!=9:continue
  M=C.gram(ss,pairs);Mi=K.inverse(M);rs=[allroots[allnames.index(n)]for n in ss]
  lam=V.mv([[V.golden(x)for x in row]for row in Mi],[V.pair(center,g,r)for r in rs]);assert all(x[1]<0 for x in lam)
  p=[V.sub(center[i],V.total(V.mul(t,r[i])for t,r in zip(lam,rs,strict=True)))for i in range(10)]
  q=V.pair(p,g,p);assert q[1]<0 and V.dot(p,hi)[0]>0
  checks=[]
  for r,n in zip(allroots,allnames,strict=True):
   if n in ss:continue
   pp=V.pair(p,g,r);assert pp[1]<0;checks.append({'wall':n,'pairing':V.enc(pp)})
  vertices.append({'incident':ss,'gram':[[K.enc(x)for x in row]for row in M],'norm':V.enc(q),'point':list(map(V.enc,p)),'strict_nonincident_pairings':checks})
 # The first new facet closes all four old edges transverse to wall17.
 face_common=[n for n in prior['compact_H3_face']['common_facets']if n!='wall17']
 oldface=prior['actual_finite_vertices'];prism_vertices=[v for v in vertices if NEW[0]in v['incident']]
 assert len(prism_vertices)==4 and all(set(face_common)<=set(v['incident'])for v in prism_vertices)
 tetra=prior['compact_H3_face']['bounding_facets'];assert len(face_common)==5
 assert C.gram(face_common,pairs)==[[K.scale(K.O,2)if i==j else K.Z for j in range(5)]for i in range(5)]
 assert all(pairs['|'.join(sorted((a,b)))]=={'kind':'finite','label':2}for a in face_common for b in tetra+[NEW[0]])
 assert all(pairs['|'.join(sorted(('wall17',b)))]=={'kind':'finite','label':2}for b in tetra)
 # Old base projected orthogonally to the common five-face has positive square.
 fi=K.inverse(C.gram(face_common,pairs));fp=[C.gram(['wall17',n],pairs)[0][1]for n in face_common]
 baseq=K.sub(K.scale(K.O,2),K.pair(fp,fi,fp));assert baseq==K.scale(K.O,2)
 old_survival=[]
 for v in oldface:
  p=[tuple(map(int,x))for x in v['point_enclosures']];ps=[V.pair(p,g,r)for r in caps];assert max(x[1]for x in ps)<0
  old_survival.append({'incident':v['incident'],'new_pairings':list(map(V.enc,ps))})
 assert {frozenset(v['incident'])for v in prism_vertices}=={frozenset(face_common+[NEW[0]]+[n for n in tetra if n!=omit])for omit in tetra}
 # The side cone is the old tetrahedral cone. In its orthogonal base coordinate,
 # the old base is z<=0 and the new base has negative z coefficient. At every
 # side-cone generator the new inequality is strict at z=0. Hence the six-wall
 # cone is generated by exactly the four old and four new endpoint rays.
 # All eight are timelike and satisfy every remaining wall: a compact H4 prism.
 # Full rank8 edge-frontier replay, including the inherited four vertices.
 vp=oldface+vertices;edge_sets=sorted({tuple(sorted(ss))for v in vp for ss in combinations(v['incident'],8)})
 frontier=[]
 for ss in edge_sets:
  ends=[];dotted=bad=0
  for n in allnames:
   if n in ss:continue
   if any(pairs['|'.join(sorted((n,a)))]['kind']=='dotted'for a in ss):dotted+=1;continue
   typ=P.T.psd(C.gram(list(ss)+[n],pairs));assert typ!='positive_semidefinite'
   if typ=='positive_definite':ends.append(n)
   else:bad+=1
  assert len(ends)in(1,2)
  assert all(frozenset(ss+(n,))in{frozenset(v['incident'])for v in vp}for n in ends)
  frontier.append({'common_facets':list(ss),'finite_end_walls':ends,'status':'finite'if len(ends)==2 else'open_to_light_cone','excluded_by_dotted_pair':dotted,'excluded_by_indefinite_gram':bad})
 vv={frozenset(v['incident'])for v in vp};reached={next(iter(vv))}
 while True:
  before=len(reached)
  for e in frontier:
   if len(e['finite_end_walls'])!=2:continue
   ends={frozenset(e['common_facets']+[n])for n in e['finite_end_walls']}
   if reached&ends:reached.update(ends)
  if len(reached)==before:break
 assert reached==vv
 target=tuple(map(K.parse,d['target']));assert d['target']==V.read(d['source'])['escape']
 assert K.sign(K.pair(target,I,target))>0 and K.sign(K.pair(target,I,h))<0
 ti=list(map(V.golden,target));assert min(V.dot(r,ti)[0]for r in old)>0
 cut=[V.dot(r,ti)for r in caps];assert cut[1][1]<0
 out={'status':'PASS_EXACT_96_WALL_COMPACT_H4_PRISM_AND_TARGET_CUT','source_sha256':V.sha(PREFIX+'_caps.json'),'facets':96,'retained_facets':94,'new_normals':rows,'new_pair_counts':counts,'new_cross_pair':V.enc(cross),'normal_enclosures':[list(map(V.enc,r))for r in caps],
  'affine_census':{'all_subsets_with_one_new_facet':1024,'positive_definite':len(positive),'indefinite':indefinite,'new_affine_circuits':0,'inherited_affine_circuits':0,'two_new_facets_excluded_by_dotted_pair':True},
  'actual_finite_vertices':vertices,'prism_common_facets':face_common,'prism_common_type':'5A1','prism_old_base':'wall17','prism_new_base':NEW[0],'prism_projected_old_base_norm':K.enc(baseq),'prism_vertices':8,'prism_old_vertex_survival':old_survival,'target_cut_pairings':list(map(V.enc,cut)),
  'vertex_component_frontier':{'vertices':len(vp),'finite_edges':sum(e['status']=='finite'for e in frontier),'open_rays':sum(e['status']=='open_to_light_cone'for e in frontier),'edges':frontier},
  'scope':'Two exact finite-endpoint caps retain94walls and cut the saved strict gap. All new pairings and every new positive-semidefinite candidate subset are checked. The first cap closes four transverse edges to produce a compact H4 prism face; the second closes two further edges. No finite-volume claim.'}
 if 'escape'in d:
  yy=tuple(map(K.parse,d['escape']));q=K.pair(yy,I,yy);f=K.pair(yy,I,h);assert K.sign(q)>0 and K.sign(f)<0
  pp=[V.dot(r,list(map(V.golden,yy)))for r in allroots];assert min(x[0]for x in pp)>0
  out['remaining_escape']={'norm':K.enc(q),'future_center_pairing':K.enc(f),'strict_pairings':list(map(V.enc,pp))};out['scope']+=' An exact strict future spacelike escape proves infinite volume.'
 return out

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=V.read(PREFIX+'_caps.json');out=check(d)
 if a.verify:assert out==V.read(PREFIX+'_verify.json')
 else:(V.ROOT/(PREFIX+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['new_pair_counts'],out['vertex_component_frontier']['vertices'])
