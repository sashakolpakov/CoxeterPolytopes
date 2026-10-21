#!/usr/bin/env python3
"""Producer-free exact replay of the rational H21 symmetry-quotient corner walk.

No endpoint inverse/projection implementation is imported.  Endpoint validity is
checked geometrically: two distinct feasible future causal vertices share
exactly the twenty positive-definite walls defining their edge.  Convexity then
excludes any intervening wall or another endpoint.
"""
import argparse, hashlib, itertools, json, math, time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parents[1];TAG='fresh_h21_source_corner_walk_20261013';SRC=P/'results'/f'{TAG}.json';OUT=P/'results'/f'{TAG}_independent.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def mv(A,v):return [dot(r,v)for r in A]
def pd(A):
 a=[list(map(F,r))for r in A]
 for k in range(len(a)):
  d=a[k][k];assert d>0
  for i in range(k+1,len(a)):
   for j in range(i,len(a)):
    a[i][j]-=a[i][k]*a[j][k]/d;a[j][i]=a[i][j]
 return True

def replay():
 t=time.monotonic();d=json.loads(SRC.read_text());assert d['status']in('CHECKPOINT_EXACT_INCOMPLETE_CORNER_ORBIT_WALK','PASS_COMPLETE_CORNER_ORBIT_CENSUS')
 for p,h in d['sources'].items():assert sha(P/p)==h,p
 gp=P/d['group_array'];assert sha(gp)==d['group_sha256'];z=np.load(gp);perms=z['permutations'];parents=z['parents'];labels=z['generator'];assert perms.shape==(241920,210)
 raw=json.loads((P/'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json').read_text())['coxeter_input'];Q=[[-int(x)for x in r]for r in raw['basis_gram']];W=[tuple(map(int,r))for r in raw['wall_coordinates_as_columns_transposed']];C=[mv(Q,r)for r in W];G=[[dot(r,c)for c in C]for r in W]
 h=[sum(r[i]for r in W)for i in range(22)];gg=math.gcd(*h);h=[x//gg for x in h];hc=mv(Q,h);assert dot(h,hc)==-14 and max(dot(r,hc)for r in W)<0
 syms=json.loads((P/'data/enumerations/fresh_h20_global_interior_sections_20261013.json').read_text())['ambient_symmetry_generators'];gens=[]
 for rec in syms:
  M=rec['matrix'];g=np.array(rec['wall_permutation'],dtype=np.uint8);assert sorted(g)==list(range(210));assert mv(M,h)==h
  assert all(mv(M,w)==list(W[int(g[i])])for i,w in enumerate(W))
  QM=[mv(Q,[M[i][j]for i in range(22)])for j in range(22)]
  assert all(dot([M[k][i]for k in range(22)],QM[j])==Q[i][j]for i in range(22)for j in range(22));gens.append(g)
 assert np.array_equal(perms[0],np.arange(210));assert int(parents[0])==int(labels[0])==-1
 rows={r.tobytes():i for i,r in enumerate(perms)};assert len(rows)==len(perms)
 for i in range(1,len(perms)):
  p=int(parents[i]);g=int(labels[i]);assert 0<=p<i and 0<=g<6;assert np.array_equal(perms[i],gens[g][perms[p]])
 for g in gens:
  assert all(g[r].tobytes()in rows for r in perms)
 print('GROUP PASS',len(perms),flush=True)
 records=d['corner_orbits'];keys=set();stabs=[];actives=[];rays=[];norms=[];pd_cache=set();ideal=finite=0
 def verify_vertex(v,active):
  assert len(v)==22 and math.gcd(*v)==1
  pp=[dot(v,c)for c in C];vv=dot(v,mv(Q,v));assert vv<=0 and dot(v,hc)<0 and max(pp)==0
  assert active==[i for i,x in enumerate(pp)if x==0];return vv
 for idx,r in enumerate(records):
  assert r['id']==idx;v=r['primitive_ray'];active=r['active_walls'];vv=verify_vertex(v,active);assert vv==r['norm'] and dot(v,hc)==r['center_pairing']
  aa=np.array(active,dtype=np.uint8);images=np.sort(perms[:,aa],axis=1);cg=r['canonical_group_index'];key=tuple(r['canonical_active_walls']);assert tuple(map(int,images[cg]))==key and key not in keys;keys.add(key)
  candidates=np.arange(len(perms))
  for j in range(len(active)):
   val=images[candidates,j];candidates=candidates[val==val.min()]
  assert tuple(map(int,images[int(candidates[0])]))==key
  stab=np.flatnonzero(np.all(images==aa,axis=1));assert stab.tolist()==r['stabilizer_group_indices'];assert len(perms)//len(stab)==r['orbit_size'];stabs.append(stab);actives.append(active);rays.append(v);norms.append(vv)
  if vv<0:
   assert r['type']=='finite'and len(active)==21 and not r['affine_components'];pd([[G[i][j]for j in active]for i in active]);finite+=1
  else:
   assert r['type']=='ideal';ideal+=1;components=r['affine_components'];assert sorted(i for c in components for i in c['walls'])==active and sum(c['rank']for c in components)==20
   for c in components:
    cc=c['walls'];weights=list(map(F,c['positive_kernel']));mu=F(c['ray_multiple']);assert len(cc)==c['rank']+1 and min(weights)>0 and mu>0
    assert all(sum(w*W[i][j]for w,i in zip(weights,cc))==mu*v[j]for j in range(22))
    assert all(sum(G[i][j]*w for j,w in zip(cc,weights))==0 for i in cc)
    pd([[G[i][j]for j in cc[:-1]]for i in cc[:-1]])
    reached={cc[0]}
    while True:
     new={j for i in reached for j in cc if G[i][j]!=0};more=reached|new
     if more==reached:break
     reached=more
    assert reached==set(cc)
   assert all(G[i][j]==0 for a,c in enumerate(components)for b,e in enumerate(components)if a!=b for i in c['walls']for j in e['walls'])
  if idx%100==0:print('CORNERS',idx,flush=True)
 edge_count=literal_count=0
 for idx,r in enumerate(records):
  if not r['processed']:assert not r['edge_orbits'];continue
  active=actives[idx];v=rays[idx];stab=perms[stabs[idx]]
  if r['type']=='finite':remaining={(i,)for i in active}
  else:remaining={tuple(sorted(x))for x in itertools.product(*(c['walls']for c in r['affine_components']))}
  total=len(remaining);assert total==r['all_incident_edge_count']
  for e in r['edge_orbits']:
   deleted=tuple(e['deleted_representative']);assert deleted in remaining
   orbit={tuple(map(int,row))for row in np.sort(stab[:,np.array(deleted,dtype=np.uint8)],axis=1)};assert orbit<=remaining and len(orbit)==e['stabilizer_orbit_size'];remaining-=orbit;literal_count+=len(orbit)
   edge=tuple(i for i in active if i not in deleted);assert list(edge)==e['edge_walls']and len(edge)==20
   if edge not in pd_cache:pd([[G[i][j]for j in edge]for i in edge]);pd_cache.add(edge)
   end=e['endpoint_primitive_ray'];pp=[dot(end,c)for c in C];aa=[i for i,x in enumerate(pp)if x==0];verify_vertex(end,aa);assert set(aa)&set(active)==set(edge)
   assert end!=v and dot(v,mv(Q,end))<0
   cg=e['endpoint_canonical_group_index'];target=e['endpoint_orbit'];assert sorted(int(perms[cg,j])for j in aa)==records[target]['canonical_active_walls']
   assert F(e['parameter'])>0 if e['formula']=='v-t*u'else e['formula']=='p+t*delta'
   edge_count+=1
  assert not remaining
 pending=[r['id']for r in records if not r['processed']];assert pending==d['pending_orbit_ids'];assert len(records)-len(pending)==d['processed_orbits']and edge_count==d['edge_orbits_processed'];assert d['counts']=={'ideal':ideal,'finite':finite}
 assert bool(pending)==(d['status']=='CHECKPOINT_EXACT_INCOMPLETE_CORNER_ORBIT_WALK')
 literal_finite=sum(r['orbit_size']for r in records if r['type']=='finite');literal_ideal=sum(r['orbit_size']for r in records if r['type']=='ideal');directed=sum(r['orbit_size']*r.get('all_incident_edge_count',0)for r in records)
 if not pending:assert directed%2==0
 out={'status':'PASS_INDEPENDENT_EXACT_CORNER_WALK_REPLAY','producer_sha256':sha(SRC),'checker_sha256':sha(Path(__file__)),'group_order':len(perms),'corner_orbits':len(records),'finite_corner_orbits':finite,'cusp_orbits':ideal,'processed_corner_orbits':d['processed_orbits'],'pending_corner_orbits':len(pending),'edge_orbits':edge_count,'literal_edge_directions_at_representatives':literal_count,'distinct_PD_edge_supports_checked':len(pd_cache),'queue_closed':not pending,'literal_finite_vertices':literal_finite,'literal_cusps':literal_ideal,'literal_edges':directed//2 if not pending else None,'elapsed_seconds':time.monotonic()-t,'scope':'Independent rational vertex/edge and finite permutation-group replay. Endpoint validity is proved from common PD20 support and distinct feasible future causal endpoints, without importing the producer inverse/projection code. This replay verifies no numerical evidence and does not turn an incomplete queue into an exhaustive census.'}
 return out
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=replay()
 if a.verify:
  old=json.loads(OUT.read_text());assert {k:v for k,v in old.items()if k!='elapsed_seconds'}=={k:v for k,v in d.items()if k!='elapsed_seconds'}
 else:OUT.write_text(json.dumps(d,indent=2)+'\n')
 print(json.dumps(d,indent=2))
