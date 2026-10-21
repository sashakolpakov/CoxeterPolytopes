#!/usr/bin/env python3
"""Exact partial forest weights after the immutable BFS/two-chain additions.

The frozen three-count LP/dual/DP engine is reused read-only. All literal face
graphs are reconstructed and deduplicated by the independent marked matcher.
No claim is made about the remaining finite or dotted attachment frontier.
"""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_a1_count_weights_20261015 as C
import fresh_sqrt5_sealed_subcore_closure_20261015_review as S

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_enlarged_a1_weights_20261015'
B=C.B;G=B.G
BFS='fresh_sqrt5_nonsingular_attachment_search_20261015_batch610'
CHAINS='fresh_sqrt5_two_chain_star_closure_20261015'

def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def inputs():
 base=read('results/fresh_sqrt5_quad_exhaustion_20261015.json')['labelled_union_supports']
 old=read('results/fresh_sqrt5_forest_sensitive_weights_20261015.json')
 first=read('results/fresh_sqrt5_singular_extension_frontier_20261015.json')
 closure=read('results/fresh_sqrt5_singular_attachment_closure_20261015.json')
 rows=[('old112',r['class_index'],base[r['union_index']])for r in old['sealed_quad_rows']]
 rows +=[('new46',r['index'],r['row'])for r in first['records']if not r['accepted']]
 rows +=[('closure11',r['index'],r['row'])for r in closure['records']if r['status']=='NEW_SEALED_PAID']
 rows +=[('zero_touch5',r['index'],r['row'])for r in first['records']if r['accepted']]
 assert len(rows)==174
 bfs=read('results/'+BFS+'.json');chain=read('results/'+CHAINS+'.json');export=read('results/'+CHAINS+'_export.json')
 assert len(bfs['records'])==20
 for r in bfs['records']:
  assert r['result']['status']=='UNPAID_SEALED_FRONTIER'and not r['result']['accepted']
  rows.append(('bfs20',r['index'],r['row']))
 unpaid=chain['new_sealed_unpaid_marked_classes'];assert unpaid==[50,55,62]
 for r in export['records']:
  if r.get('necessary_weight_obligation'):
   assert r['source_status']=='NEW_SEALED_UNPAID'and r['source_index']in unpaid
   rows.append(('two_chain3',r['source_index'],r['row']))
 assert len(rows)==197
 # Retain all newly sealed paid states too: their old bonus payment need not
 # survive an arbitrary reoptimization. Closed starts also keep zero-touch cases.
 checks=[('two_chain_'+r['source_kind'],r['source_index'],r['row'])for r in export['records']if not r.get('necessary_weight_obligation')]
 checks +=[('two_chain_zero_touch',r['index'],r['row'])for r in chain['records']if r['status']=='OPEN_WITH_CONTACTS']
 assert len(checks)==30 and chain['open_attachment_states']==0
 rows+=checks;assert len(rows)==227
 return rows,checks

def corners(row,orbitdata=None):
 es=dict(zip(combinations(range(row['nodes']),2),row['labels']));out=[]
 for a,b in product((0,2),(1,3)):
  cs=G.pieces(G.graph(row['white']+[a,b],es));component=next(c for c in cs if a in c)
  connected=b in component;base=G.base(component,a,b)if connected else F(0);key=None
  if connected and base==0:
   assert G.distances(component,a)[b]>=3;key=B.code(component,(a,b));ids=sorted(component)
   if orbitdata is not None and key not in orbitdata:
    orbitdata[key]=dict(type=G.family(component),rank=len(component),nodes=ids,labels=[component[u].get(v,2)for u,v in combinations(ids,2)],marked=[a,b])
  out.append(dict(marked=[a,b],core_forest=sorted(map(G.family,cs)),baseline=str(base),orbit=key))
 return out

def build():
 raw,checks=inputs();pool=defaultdict(list);faces=[];literal_maps=[];orbitdata={}
 for family,index,row in raw:
  M=S.matrix(row);inv=S.invariant(M);hit=None
  for j,N in pool[inv]:
   mp=S.isomorphism(M,N)
   if mp is not None:S.check_map(M,N,mp);hit=(j,mp);break
  if hit is None:
   j=len(faces);pool[inv].append((j,M));mp=list(range(row['nodes']))
   faces.append(dict(family=family,index=index,row=row,corners=corners(row,orbitdata),origins=[]))
  else:j,mp=hit
  origin=dict(family=family,index=index,to_representative=mp);faces[j]['origins'].append(origin)
  literal_maps.append(dict(**origin,representative=j))
 orbitkeys=sorted(orbitdata);variables=[(k,c)for k in orbitkeys for c in range(3)];index={v:i for i,v in enumerate(variables)}
 constraints=[];unique={}
 for fi,face in enumerate(faces):
  for outside in range(3):
   co=Counter();baseline=F(0)
   for c in face['corners']:
    baseline+=F(c['baseline'])
    if c['orbit']is not None:co[index[c['orbit'],min(2,c['core_forest'].count('A1')+outside)]]+=1
   key=(tuple(sorted(co.items())),baseline-1)
   unique.setdefault(key,[]).append(dict(face=fi,origins=face['origins'],external_A1_count_capped2=outside))
   constraints.append(dict(face=fi,external_A1_count_capped2=outside,baseline=str(baseline),coefficients=[[j,v]for j,v in sorted(co.items())]))
 rows=[];rhs=[];meta=[]
 for(co,b),witnesses in sorted(unique.items()):rows.append({j:-v for j,v in co});rhs.append(b);meta.append(dict(kind='face',witnesses=witnesses))
 components=[]
 for name,g in G.models():
  if len(g)>28:continue
  baseline=sum((G.base(g,a,b)for a,b in combinations(g,2)),F(0));coef=[Counter()for _ in range(3)];orbits=defaultdict(list)
  for a,b in combinations(g,2):
   key=B.code(g,(a,b));orbits[key].append([a,b])
   if key in orbitdata:
    assert G.base(g,a,b)==0
    for c in range(3):coef[c][index[key,c]]+=1
  components.append(dict(type=name,rank=len(g),baseline=baseline,coeffs=coef,affected_orbits=[dict(key=k,multiplicity=len(v),pairs=v)for k,v in sorted(orbits.items())if k in orbitdata]))
 assert len(components)==89
 assert{k for c in components for k in[c['affected_orbits'][i]['key']for i in range(len(c['affected_orbits']))]}==set(orbitkeys)
 return dict(faces=faces,faceconstraints=constraints,orbits=[dict(key=k,**orbitdata[k])for k in orbitkeys],variables=variables,components=components,rows=rows,rhs=rhs,meta=meta,zi=len(variables),literal_maps=literal_maps,extra_checks=checks)

def additional_checks(data,solution):
 x=list(map(F,solution['exact_primal']));index={v:i for i,v in enumerate(data['variables'])};out=[]
 for family,i,row in data['extra_checks']:
  cs=corners(row);sums=[]
  for outside in range(3):
   value=F(0)
   for c in cs:
    value+=F(c['baseline']);key=c['orbit'],min(2,c['core_forest'].count('A1')+outside)
    if key in index:value+=x[index[key]]
   sums.append(str(value))
  out.append(dict(family=family,index=i,corners=cs,payments=sums,paid=min(map(F,sums))>=1))
 return out

def run(saved=None):
 data=build();print('literal227 distinct',len(data['faces']),'orbits',len(data['orbits']),'vars',data['zi'],'constraints',len(data['rows']),flush=True)
 solution=C.solve(data)if saved is None else saved['solution'];proof=C.verify(data,solution);rank29=C.rank29_check(data,solution);checks=additional_checks(data,solution)
 paths=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_a1_count_weights_20261015.py','results/fresh_sqrt5_a1_count_weights_20261015.json','scripts/fresh_sqrt5_sealed_subcore_closure_20261015_review.py','results/'+BFS+'.json','scripts/'+BFS+'.py','results/'+CHAINS+'_export.json','results/'+CHAINS+'.json','scripts/'+CHAINS+'.py']
 for p in paths:
  if p.endswith('.json'):
   for q,h in read(p).get('sources',{}).items():assert sha(q)==h,(p,q)
 assert all(r['paid']for r in checks)
 return dict(status='PASS_EXACT_ENLARGED_PARTIAL_A1_COUNT_WEIGHT_OPTIMUM',sources={p:sha(p)for p in paths},literal_inputs=227,distinct_marked_face_graphs=len(data['faces']),literal_maps=data['literal_maps'],face_cases=3*len(data['faces']),distinct_constraints=len(data['rows']),orbit_count=len(data['orbits']),variables=data['zi'],finite_component_types=89,orbits=data['orbits'],variable_order=[dict(orbit=k,A1_count_capped2=c)for k,c in data['variables']],faces=data['faces'],face_constraints=data['faceconstraints'],solution=solution,exact_certificate=proof,target=77,meets_target=F(solution['optimum'])<=77,rank29_same_rule=rank29,additional_two_chain_checks=checks,all_additional_checks_paid=all(r['paid']for r in checks),
  scope='Partial known-face optimization: old174 plus20 immutable BFS unpaid seals and all31 two-chain sealed states plus two closed zero-touch starts, deduplicated as full marked graphs. Exact primal/dual and complete rank28 forest DP; same literal rule evaluated at rank29. Sealing proofs are inherited; this is not their independent replay. Remaining attachment/quad coverage is incomplete; no dimension theorem.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();p=ROOT/'results'/f'{TAG}.json';saved=json.loads(p.read_text())if a.verify else None;out=json.loads(json.dumps(run(saved)))
 if a.verify:assert out==saved
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['solution']['optimum'],'rank29',out['rank29_same_rule']['maximum'],'extra_paid',out['all_additional_checks_paid'])
