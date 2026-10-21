#!/usr/bin/env python3
"""Exact symmetry-quotient corner walk of the210-wall rational H21 source."""
import argparse,json,hashlib,math,time,os,itertools
from pathlib import Path
from fractions import Fraction as F
from collections import deque,Counter
import numpy as np
from sympy import Matrix
from fresh_h22_incidence_core_20261013 import inertia
P=Path(__file__).resolve().parents[1];TAG='fresh_h21_source_corner_walk_20261013';OUT=P/'results'/f'{TAG}.json';GROUP=P/'data/enumerations'/f'{TAG}_group.npz'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json';SYMS='data/enumerations/fresh_h20_global_interior_sections_20261013.json';SAMPLE='data/enumerations/fresh_h22_borcherds_extension_20261013_scout.json';CUSPS='data/enumerations/fresh_h22_cusp_orbits_20261013.json'
def read(p):return json.loads((P/p).read_text())
def sha(p):return hashlib.sha256((P/p).read_bytes()).hexdigest()
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def mv(A,v):return[dot(r,v)for r in A]
def inverse(A):return[[F(x)for x in row]for row in Matrix(A).inv().tolist()]
def primitive(v):
 den=math.lcm(*(F(x).denominator for x in v));a=[int(F(x)*den)for x in v];gg=math.gcd(*a);assert gg;return tuple(x//gg for x in a)
def checkgroup(perms,parents,labels,gens):
 assert perms.shape==(241920,210) and np.array_equal(perms[0],np.arange(210))
 rows={r.tobytes()for r in perms};assert len(rows)==241920
 assert int(parents[0])==int(labels[0])==-1
 for i in range(1,len(perms)):
  p=int(parents[i]);g=int(labels[i]);assert 0<=p<i and 0<=g<6 and np.array_equal(perms[i],gens[g][perms[p]])
 for g in gens:assert all(g[r].tobytes()in rows for r in perms)
 return perms,parents,labels
def buildgroup(gens,manifest):
 if GROUP.exists():
  z=np.load(GROUP);assert str(z['source_manifest'])==manifest;perms=z['permutations'];parents=z['parents'];labels=z['generator'];assert len(perms)==241920;return checkgroup(perms,parents,labels,gens)
 ident=bytes(range(210));seen={ident};rows=[ident];parents=[-1];labels=[-1];cursor=0
 while cursor<len(rows):
  r=np.frombuffer(rows[cursor],dtype=np.uint8)
  for j,g in enumerate(gens):
   t=g[r].tobytes()
   if t not in seen:seen.add(t);rows.append(t);parents.append(cursor);labels.append(j)
  cursor+=1;assert len(rows)<=241920
 assert len(rows)==241920
 perms=np.frombuffer(b''.join(rows),dtype=np.uint8).reshape((-1,210));parents=np.array(parents,dtype=np.int32);labels=np.array(labels,dtype=np.int8)
 np.savez_compressed(GROUP,permutations=perms,parents=parents,generator=labels,source_manifest=np.array(manifest));print('GROUP',len(perms),flush=True);return checkgroup(perms,parents,labels,gens)

def run(max_vertices=200,max_seconds=1200,resume=False):
 started=time.monotonic();print('PID',os.getpid(),flush=True)
 raw=read(RAW)['coxeter_input'];Q=[[-int(x)for x in row]for row in raw['basis_gram']];W=[tuple(map(int,r))for r in raw['wall_coordinates_as_columns_transposed']];COV=[mv(Q,r)for r in W];G=[[dot(r,c)for c in COV]for r in W];h=primitive([sum(r[j]for r in W)for j in range(22)]);hc=mv(Q,h);assert dot(h,hc)==-14 and max(dot(r,hc)for r in W)<0
 syms=read(SYMS)['ambient_symmetry_generators'];gens=[]
 for rec in syms:
  M=rec['matrix'];perm=rec['wall_permutation'];assert sorted(perm)==list(range(210))
  assert mv(M,h)==list(h) and all(mv(M,r)==list(W[perm[i]])for i,r in enumerate(W))
  assert Matrix(M).T*Matrix(Q)*Matrix(M)==Matrix(Q);gens.append(np.array(perm,dtype=np.uint8))
 manifest=json.dumps({p:sha(p)for p in(RAW,SYMS)},sort_keys=True);perms,parents,labels=buildgroup(gens,manifest)
 records=[];known={};literal={};queue=deque();seed_ids=[];previous_elapsed=0;resume_history=[]
 if resume:
  old=read(str(OUT.relative_to(P)));assert old['status']=='CHECKPOINT_EXACT_INCOMPLETE_CORNER_ORBIT_WALK'
  for pp in (RAW,SYMS,SAMPLE,CUSPS):assert old['sources'][pp]==sha(pp)
  assert old['group_sha256']==hashlib.sha256(GROUP.read_bytes()).hexdigest()
  records=old['corner_orbits'];queue=deque(old['pending_orbit_ids']);seed_ids=old['seed_orbit_ids'];previous_elapsed=old['elapsed_seconds'];resume_history=old.get('resume_history',[])+[{'previous_producer_sha256':old['sources']['scripts/'+Path(__file__).name],'processed_orbits':old['processed_orbits'],'pending_orbits':len(queue)}]
  assert list(queue)==[r['id']for r in records if not r['processed']]
  for r in records:
   known[tuple(r['canonical_active_walls'])]=r['id'];literal[tuple(r['active_walls'])]=(r['id'],r['canonical_group_index'])
   for e in r['edge_orbits']:
    av=tuple(i for i,c in enumerate(COV)if dot(e['endpoint_primitive_ray'],c)==0);literal[av]=(e['endpoint_orbit'],e['endpoint_canonical_group_index'])
  print('RESUME',len(records),'pending',len(queue),flush=True)
 def canon(active):
  aa=np.array(active,dtype=np.uint8);images=np.sort(perms[:,aa],axis=1);idx=np.arange(len(perms))
  for j in range(len(active)):
   vals=images[idx,j];idx=idx[vals==vals.min()]
   if len(idx)==1:break
  cg=int(idx[0]);key=tuple(map(int,images[cg]));stab=np.flatnonzero(np.all(images==aa,axis=1));assert len(stab)and len(perms)%len(stab)==0;return key,cg,stab
 def add(v,parent=None):
  v=primitive(v);vv=dot(v,mv(Q,v));hp=dot(v,hc);assert vv<=0 and hp<0;pp=[dot(v,c)for c in COV];assert max(pp)==0;active=tuple(i for i,x in enumerate(pp)if x==0)
  if active in literal:return literal[active][0],literal[active][1]
  key,cg,stab=canon(active)
  if key in known:idx=known[key];literal[active]=idx,cg;return idx,cg
  idx=len(records);known[key]=idx;literal[active]=idx,cg;queue.append(idx)
  comps=[]
  if vv<0:
   assert len(active)==21 and inertia([[G[i][j]for j in active]for i in active])==[21,0,0]
  else:
   unseen=set(active)
   while unseen:
    todo=[unseen.pop()];cc=[]
    while todo:
     i=todo.pop();cc.append(i);near={j for j in unseen if G[i][j]!=0};unseen-=near;todo.extend(near)
    cc.sort();A=[[G[i][j]for j in cc]for i in cc];assert inertia(A)==[len(cc)-1,0,1];ker=Matrix(A).nullspace();assert len(ker)==1;weights=[F(x)for x in ker[0]]
    if weights[0]<0:weights=[-x for x in weights]
    assert min(weights)>0;s=[sum(c*W[i][j]for c,i in zip(weights,cc))for j in range(22)];j=next(j for j,x in enumerate(v)if x);mu=s[j]/v[j];assert mu>0 and all(x==mu*y for x,y in zip(s,v));comps.append({'walls':cc,'rank':len(cc)-1,'positive_kernel':list(map(str,weights)),'ray_multiple':str(mu)})
   assert sum(c['rank']for c in comps)==20
  records.append({'id':idx,'type':'finite'if vv<0 else'ideal','primitive_ray':list(v),'norm':vv,'center_pairing':hp,'active_walls':list(active),'canonical_active_walls':list(key),'canonical_group_index':cg,'stabilizer_group_indices':stab.tolist(),'orbit_size':len(perms)//len(stab),'affine_components':comps,'discovery':parent,'processed':False,'edge_orbits':[]})
  print('NEW',idx,records[-1]['type'],'active',len(active),'orbit',records[-1]['orbit_size'],'queue',len(queue),flush=True);return idx,cg
 # Seed the known fifteen cusp orbits and a genuine finite corner. This
 # accelerates exploration; it does not presume they exhaust the chamber.
 if not resume:
  for i,row in enumerate(read(CUSPS)['orbits']):seed_ids.append(add(row['representative'],{'known_cusp_sample':i})[0])
  finite=next(r for r in read(SAMPLE)['old_vertices']if F(r['norm'])<0);seed_ids.append(add(finite['vector'],{'known_finite_sample':True})[0])
 def save(status):
  out={'status':status,'sources':{p:sha(p)for p in(RAW,SYMS,SAMPLE,CUSPS,'scripts/'+Path(__file__).name)},'group_array':str(GROUP.relative_to(P)),'group_sha256':hashlib.sha256(GROUP.read_bytes()).hexdigest(),'group_order':len(perms),'dimension':21,'wall_count':210,'seed_orbit_ids':seed_ids,'corner_orbits':records,'pending_orbit_ids':list(queue),'counts':dict(Counter(r['type']for r in records)),'processed_orbits':sum(r['processed']for r in records),'edge_orbits_processed':sum(len(r['edge_orbits'])for r in records),'elapsed_seconds':previous_elapsed+time.monotonic()-started,'resume_history':resume_history,'scope':'Exact orbit-level corner walk. All listed endpoints, active sets and affine links are certified. Completeness is asserted only if the orbit queue is empty; otherwise pending vertices and their outgoing directions remain unexplored. The six generated symmetries are used, not an assumed full automorphism group.'}
  OUT.write_text(json.dumps(out,indent=2)+'\n')
 processed=0
 while queue and processed<max_vertices and time.monotonic()-started<max_seconds:
  idx=queue.popleft();r=records[idx];active=r['active_walls'];v=list(map(F,r['primitive_ray']));stab=perms[r['stabilizer_group_indices']]
  if r['type']=='finite':choices={(i,)for i in active};A=[[G[i][j]for j in active]for i in active];AI=inverse(A)
  else:choices={tuple(sorted(t))for t in itertools.product(*(c['walls']for c in r['affine_components']))}
  total=len(choices);cover=0
  while choices:
   deleted=min(choices);images=np.sort(stab[:,np.array(deleted,dtype=np.uint8)],axis=1);orbit={tuple(map(int,x))for x in images};assert orbit<=choices;choices-=orbit;cover+=len(orbit);edge=[i for i in active if i not in deleted];assert len(edge)==20
   if r['type']=='finite':
    k=active.index(deleted[0]);u=[sum(AI[k][j]*W[i][a]for j,i in enumerate(active))for a in range(22)];assert dot(v,mv(Q,u))==0
    up=[dot(u,c)for c in COV];vp=[dot(v,c)for c in COV];bounds=[(vp[j]/up[j],j)for j in range(210)if up[j]<0];assert bounds;t=min(t for t,j in bounds);assert t>0;end=[x-t*y for x,y in zip(v,u)];formula='v-t*u';parameter=t
   else:
    A=[[G[i][j]for j in edge]for i in edge];assert inertia(A)==[20,0,0];co=mv(inverse(A),[dot(W[i],hc)for i in edge]);p=[F(h[a])-sum(c*W[i][a]for c,i in zip(co,edge))for a in range(22)];assert dot(p,mv(Q,p))<0
    pp=[dot(p,c)for c in COV];vp=[dot(v,c)for c in COV];assert max(pp)==0
    bounds=[(-pp[j]/vp[j],j)for j in range(210)if vp[j]<0];assert bounds;t=max(t for t,j in bounds);end=[x+t*y for x,y in zip(p,v)];formula='p+t*delta';parameter=t
   assert all(dot(end,COV[i])==0 for i in edge)
   target,cg=add(end,{'parent_orbit':idx,'deleted_walls':list(deleted)});r['edge_orbits'].append({'deleted_representative':list(deleted),'edge_walls':edge,'stabilizer_orbit_size':len(orbit),'formula':formula,'parameter':str(parameter),'endpoint_primitive_ray':list(primitive(end)),'endpoint_orbit':target,'endpoint_canonical_group_index':cg})
  assert cover==total;r['all_incident_edge_count']=total;r['processed']=True;processed+=1;save('RUNNING_EXACT_CORNER_ORBIT_WALK');print('DONE',idx,'edge_orbits',len(r['edge_orbits']),'total',total,'vertices',len(records),'pending',len(queue),flush=True)
 save('PASS_COMPLETE_CORNER_ORBIT_CENSUS'if not queue else'CHECKPOINT_EXACT_INCOMPLETE_CORNER_ORBIT_WALK');print('FINAL','CLOSED'if not queue else'CHECKPOINT',len(records),'pending',len(queue),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--max-vertices',type=int,default=200);p.add_argument('--max-seconds',type=int,default=1200);p.add_argument('--resume',action='store_true');aa=p.parse_args();run(aa.max_vertices,aa.max_seconds,aa.resume)
