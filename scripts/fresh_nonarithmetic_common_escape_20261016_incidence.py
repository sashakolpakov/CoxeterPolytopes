#!/usr/bin/env python3
"""Complete finite-edge incidence derived from the exact49-wall cone."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_common_escape_20261016_incidence'
OLD='results/fresh_nonarithmetic_joint_regions_20261016_cone.json';NEW='results/fresh_nonarithmetic_common_escape_20261016_cone.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 old,new=[json.loads((ROOT/p).read_text())for p in[OLD,NEW]]
 for d in[old,new]:
  for p,h in d['sources'].items():assert sha(p)==h,p
 vertices=new['complete_census']['finite_vertices'];oldset={tuple(x)for x in old['complete_census']['finite_vertices']}
 assert len(vertices)==178 and all(len(x)==9 for x in vertices)
 allrays=vertices+new['complete_census']['cusps']+new['complete_census']['hyperideal_rays']+new['complete_census']['opposite_sheet_rays']
 incidence=[0]*49
 for i,row in enumerate(allrays):
  for j in row:incidence[j]|=1<<i
 edges=[];types=Counter()
 for i,j in combinations(range(len(vertices)),2):
  common=sorted(set(vertices[i])&set(vertices[j]))
  if len(common)!=8:continue
  face=(1<<len(allrays))-1
  for k in common:face&=incidence[k]
  assert face==(1<<i)|(1<<j)
  endpoints=int(tuple(vertices[i])in oldset)+int(tuple(vertices[j])in oldset);types[endpoints]+=1
  edges.append(dict(vertices=[i,j],support=common,old_endpoint_count=endpoints))
 assert types[2]==395 and oldset<={tuple(x)for x in vertices}
 oldg=9*172-2*395;newg=9*(178-172)-types[1]-2*types[0]
 return dict(status='PASS_COMPLETE_49_WALL_FINITE_EDGE_INCIDENCE',sources={p:sha(p)for p in[f'scripts/{TAG}.py',OLD,NEW]},finite_vertices=178,finite_edges=len(edges),edges=edges,open_finite_germs=9*178-2*len(edges),closed_old_germs=types[1],remaining_old_germs=oldg-types[1],new_vertex_open_germs=newg,edge_types={str(k):v for k,v in sorted(types.items())},
  per_recent_facet={str(w):{k:sum(w in s for s in rows)for k,rows in new['complete_census'].items()}for w in range(33,49)},
  proof='All finite vertices have exactly nine active walls. Every finite pair sharing eight has a common-zero face containing exactly those two complete-cone rays, so it is a finite edge. A one-face with two finite endpoints must be detected this way. Nine germs per finite vertex minus two per finite edge is the complete open finite-germ count; each old-new edge closes one old germ.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved49incidence replay')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
