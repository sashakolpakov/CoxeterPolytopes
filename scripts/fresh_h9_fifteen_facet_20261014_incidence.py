#!/usr/bin/env python3
"""All codimension2..5 one-cut extensions of the770 saved Gale nerves.

Old induced obstructions survive unless the cut face is contained in their
support. Exact cut orbits are reduced by verified parent automorphisms.
"""
import argparse,json,os,time
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import pynauty
from fresh_h9_bugaenko_surgery_20261014_incidence import canonical,minimal
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_fifteen_facet_20261014'
SOURCE='results/fresh_h9_gale_frontier_20261014_geometry.json'
PAT='results/fresh_h9_newfacet_frontier_20261014_patterns.json';REL='results/fresh_h9_newfacet_frontier_20261014_relations.json'
BLOCK='results/fresh_h9_bugaenko_remaining_20261014_block_gate.json';LOCAL='results/fresh_h9_newfacet_frontier_20261014_local_independent.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def bits(m):return [i for i in range(m.bit_length())if m>>i&1]
def load(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d.get('sources',{}).items():assert sha(f)==h,f
 return d

def patterns():
 data=load(PAT);load(REL);load(LOCAL);load(BLOCK)
 out={(p['order'],tuple(p['minimal_nonfaces'])):p['pattern']for p in data['patterns']if p['verdict']=='EXACT_LOCAL_OBSTRUCTION'}
 out[7,(31,99)]='seven';return out
CACHE={}
def witnesses(mm,bad):
 for a,b in combinations(mm,2):
  if min(a.bit_count(),b.bit_count())<3:continue
  u=a|b;inside=[m for m in mm if m&u==m]
  if not a&b and set(inside)=={a,b}and max(a.bit_count(),b.bit_count())>=4:
   yield dict(kind='two_block',support=u,blocks=[a,b]);continue
  if min(a.bit_count(),b.bit_count())<4 or u.bit_count()>9 or any(m.bit_count()==2 for m in inside):continue
  vv=bits(u);red=tuple(sorted(sum(1<<k for k,j in enumerate(vv)if m>>j&1)for m in inside));key=(len(vv),red)
  if key not in CACHE:CACHE[key]=canonical(red,len(vv))
  can,perm=CACHE[key]
  if(len(vv),can)in bad:yield dict(kind='local',support=u,pattern=bad[len(vv),can],permutation=[vv[j]for j in perm])

def automorphisms(mm,vv):
 n=14;adj={i:[]for i in range(n+len(mm))}
 for j,m in enumerate(mm,n):
  for i in bits(m):adj[i].append(j);adj[j].append(i)
 graph=pynauty.Graph(n+len(mm),adjacency_dict=adj,vertex_coloring=[set(range(n)),set(range(n,n+len(mm)))])
 generators,a,b,orbs,k=pynauty.autgrp(graph);generators=[g[:n]for g in generators]
 for g in generators:
  assert sorted(g)==list(range(n))
  move=lambda m:sum(1<<g[i]for i in bits(m))
  assert set(map(move,mm))==set(mm)and set(map(move,vv))==set(vv)
 return generators

def cut(mm,vv,F):
 assert all(m&F!=m for m in mm)and 2<=F.bit_count()<=5
 M=minimal([F,*[m for m in mm if m&F!=F],*[(1<<14)|(m&~F)for m in mm]])
 V=set()
 for v in vv:
  if v&F==F:
   for i in bits(F):V.add((v^(1<<i))|(1<<14))
  else:V.add(v)
 assert all(v.bit_count()==9 for v in V)and max(map(int.bit_count,M))<=5
 return M,sorted(V)

def expand_case(row,source):
 old=source['cases'][row['parent']]
 M,V=cut(old['minimal_nonfaces'],old['facets'],row['cut_face'])
 order=row['canonical_permutation'];assert sorted(order)==list(range(15))
 inv={v:i for i,v in enumerate(order)}
 change=lambda m:sum(1<<inv[i]for i in bits(m))
 assert sorted(map(change,M))==row['canonical_minimal_nonfaces']
 return dict(row,minimal_nonfaces=M,vertices=V,canonical_vertices=sorted(map(change,V)))

def build():
 d=load(SOURCE);bad=patterns();parents=[];cases=[];unique={};recipes=[];stats=Counter();start=time.monotonic()
 for row in d['cases']:
  mm=row['minimal_nonfaces'];vv=row['facets'];common=16383;proof=[]
  for w in witnesses(mm,bad):
   new=common&w['support']
   if new!=common:proof.append(w);common=new
   if common.bit_count()<2:break
  allfaces=[sum(1<<i for i in f)for size in range(2,6)for f in combinations(range(14),size)if all(m&sum(1<<i for i in f)!=m for m in mm)]
  eligible={f for f in allfaces if f&common==f};stats['literal_source_faces']+=len(allfaces);stats['old_obstruction_excluded_cuts']+=len(allfaces)-len(eligible)
  stats['eligible_literal_cuts']+=len(eligible)
  gen=automorphisms(mm,vv)if eligible else[];orbits=[];left=set(eligible)
  while left:
   first=min(left);orbit={first};todo=[first]
   for f in todo:
    for g in gen:
     moved=sum(1<<g[i]for i in bits(f));assert moved in eligible
     if moved not in orbit:orbit.add(moved);todo.append(moved)
   assert orbit<=left;left-=orbit;orbits.append(sorted(orbit))
  parent=dict(case=row['case'],common_support=common,reducing_obstructions=proof,source_face_count=len(allfaces),eligible_cuts=sorted(eligible),verified_automorphism_generators=gen,cut_orbits=orbits)
  parents.append(parent);stats['parents_with_eligible_cuts'if eligible else'parents_all_cuts_closed']+=1
  for orb in orbits:
   F=orb[0];M,V=cut(mm,vv,F);can,order=canonical(M,15);inv={v:i for i,v in enumerate(order)}
   cv=sorted(sum(1<<inv[i]for i in bits(v))for v in V)
   if can not in unique:
    idx=len(cases);unique[can]=idx;new=next(witnesses(M,bad),None)
    cases.append(dict(case=idx,parent=row['case'],cut_face=F,canonical_minimal_nonfaces=can,canonical_permutation=order,new_obstruction=new))
    stats['new_local_obstruction'if new else'needs_label_search']+=1
   else:idx=unique[can]
   recipes.append(dict(parent=row['case'],cut_face=F,orbit_size=len(orb),case=idx))
  if row['case']%100==0:print('parent',row['case'],'types',len(cases),'open',stats['needs_label_search'],'seconds',round(time.monotonic()-start,2),flush=True)
 stats['operation_orbits']=len(recipes);stats['distinct_fifteen_facet_nerves']=len(cases)
 paths=[SOURCE,PAT,REL,BLOCK,LOCAL,f'scripts/{TAG}_incidence.py','scripts/fresh_h9_bugaenko_surgery_20261014_incidence.py']
 return dict(status='EXACT_ONE_CUT_FIFTEEN_FACET_FRONTIER',sources={p:sha(p)for p in paths},counts=dict(stats),parents=parents,recipes=recipes,cases=cases,
 scope='Every codimension2..5 face of each of770 frozen realizable source nerves is either covered by an unchanged old induced obstruction or belongs to one saved verified symmetry orbit. Every eligible representative is reconstructed by stellar subdivision and deduplicated. A new induced obstruction applies only after recomputing the actual fifteen-node nerve. No missing-triangle or label>=6 exclusion special to fourteen facets is used. Polytopality follows from a sufficiently small rational face truncation of the certified source; this is not all fifteen-facet types.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();print('PID',os.getpid(),flush=True);out=json.loads(json.dumps(build()));p=ROOT/f'results/{TAG}_incidence.json'
 if a.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
 print(out['status'],out['counts'],flush=True)
