#!/usr/bin/env python3
"""Independent Fraction reconstruction of the changed H20 rational geometry.

Rebuild the doubled49 H19 normals, lifted71 core and six rational polars.
Retained radical normals and all their cross pairs are covered by the separate
primary replay; they are not silently claimed as independent here.
"""
import argparse,hashlib,json
from collections import Counter,deque
from fractions import Fraction as F
from pathlib import Path
P=Path(__file__).resolve().parents[1]
SOURCE='data/enumerations/fresh_h20_changed_assembly_20261013_exact.json'
BASE='data/enumerations/fresh_dimension_bridge_20261013.json'
OLD='data/enumerations/fresh_h20_coupled_completion_20261013_second_exact.json'
OUT='results/fresh_h20_changed_assembly_20261013_root_review.json'
def read(p):return json.loads((P/p).read_text())
def sha(p):return hashlib.sha256((P/p).read_bytes()).hexdigest()
def decode(x):assert F(x[1])==0;return F(x[0])
def vec(x):return list(map(decode,x))
def dot(x,y):return sum((a*b for a,b in zip(x,y,strict=True)),F(0))
def mv(M,x):return [dot(row,x)for row in M]
def label(p,a,b):
 assert p<=0
 if p==0:return '2'
 z=p*p/(a*b)
 if z>1:return 'dotted'
 if z==1:return 'parallel'
 assert z in(F(1,4),F(1,2),F(3,4))
 return {F(1,4):'3',F(1,2):'4',F(3,4):'6'}[z]
def run():
 d=read(SOURCE);base=read(BASE);assert d['wall_count']==179 and len(d['retained_old_indices'])==102
 Q=list(map(vec,base['form']));original=list(map(vec,base['old_normals']));center=vec(d['strict_center'])
 pair=lambda x,y:dot(x,mv(Q,y));a=original[24];assert pair(a,a)==6
 assert all(label(pair(a,r),6,pair(r,r))in('2','6','parallel','dotted')for i,r in enumerate(original)if i!=24)
 doubled=[];seen=set()
 for i,r in enumerate(original):
  if i==24:continue
  rr=[x-pair(r,a)*y/3 for x,y in zip(r,a,strict=True)]
  for v in(r,rr):
   if tuple(v)not in seen:seen.add(tuple(v));doubled.append(v)
 assert doubled==list(map(vec,d['doubled_H19_normals']))and len(doubled)==49
 dn=[pair(r,r)for r in doubled];short=[i for i,q in enumerate(dn)if q==2];long=[i for i,q in enumerate(dn)if q==6]
 assert len(short)==27 and len(long)==22
 adj={i:[j for j in long if i!=j and pair(doubled[i],doubled[j])==-6]for i in long}
 signs={long[0]:1};todo=deque([long[0]])
 while todo:
  i=todo.popleft()
  for j in adj[i]:
   if j not in signs:signs[j]=-signs[i];todo.append(j)
   assert signs[j]==-signs[i]
 assert len(signs)==22 and signs==dict(d['long_height_signs'])
 core=[doubled[i]for i in short]
 for i in long:
  core.extend([doubled[i][:-1]+[F(2*signs[i])],doubled[i][:-1]+[F(-signs[i])]])
 assert core==list(map(vec,d['new_core_normals']))and len(core)==71
 norms=[pair(r,r)for r in core];assert all(q>0 for q in norms)
 counts=Counter(label(pair(r,core[j]),norms[i],norms[j])for i,r in enumerate(core)for j in range(i))
 assert dict(counts)==d['new_core_pair_counts']
 assert pair(center,center)==decode(d['strict_center_norm'])<0 and pair(center,a)==0
 assert all(pair(center,r)<0 for r in core)
 oldtarget=read(OLD)['strict_escape']['vector'];cov=mv(Q,core[59])
 cut=dot(cov,[F(x[0])for x in oldtarget]);cut_radical=dot(cov,[F(x[1])for x in oldtarget]);assert cut_radical==0
 assert cut==decode(d['previous176_escape_cut']['pairing'])>0
 roots=core[:];polar_rows=[]
 for entry in d['new_polars']:
  r=vec(entry['normal']);q=pair(r,r);assert q==decode(entry['norm'])>0 and pair(r,center)<0
  kinds=[label(pair(r,v),q,pair(v,v))for v in roots]
  for i,k in enumerate(kinds):
   archive_index=i if i<71 else 173+i-71
   assert entry['old_pair_checks'][archive_index]['kind']==k
  x=vec(entry['strict_target']['vector']);assert pair(x,x)>0 and pair(x,center)<0
  assert pair(r,x)==decode(entry['strict_target']['new_pairing'])>0
  assert all(pair(x,v)<0 for v in roots)
  roots.append(r);polar_rows.append({'norm':str(q),'rational_pair_counts':dict(sorted(Counter(kinds).items()))})
 escape=vec(d['strict_escape']['vector']);square=pair(escape,escape)
 assert square==decode(d['strict_escape']['norm'])>0 and pair(escape,center)<0
 maximum=max(pair(escape,r)for r in roots);assert maximum<0
 return {'status':'PASS_INDEPENDENT_FRACTION_CHANGED_CORE_AND_SIX_POLARS',
  'sources':{p:sha(p)for p in(SOURCE,BASE,OLD)},'doubled_H19_facets':49,'lifted_core_facets':71,
  'new_core_pair_counts':dict(sorted(counts.items())),'core_target_cut':str(cut),
  'six_polars':polar_rows,'strict_escape_square':str(square),'escape_rational77_maximum_pairing':str(maximum),
  'scope':'Independent Fraction reconstruction of all49 doubled source normals,71 lifted core normals, their2485 pairs and strict center, old-target cut, six rational polars against core/previous polars, all sequential target cuts, and final escape against these77 rational walls. The102 retained old normals and cross pairs remain explicit dependencies of the separately rerun full primary certificate.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
 if a.verify:assert out==read(OUT)
 else:(P/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
