#!/usr/bin/env python3
"""Conservative rank28/29 weighting stress test for additional parameter families.

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
import fresh_sqrt5_enlarged_a1_weights_20261015 as E
import fresh_sqrt5_sealed_subcore_closure_20261015_review as S

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_parameter_weight_stress_20261015'
B=C.B;G=B.G
BFS='fresh_sqrt5_nonsingular_attachment_search_20261015_batch610'
CHAINS='fresh_sqrt5_two_chain_star_closure_20261015'

def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def inputs():
 rows,_=E.inputs();rows=[(fam,i,row,None)for fam,i,row in rows]
 families=read('results/fresh_sqrt5_subdivision_families_20261015.json')
 for r in families['rank_at_most29_patterns']:
  assert r['corner_rank']==r['row']['nodes']-2 and r['corner_rank']<=29
  rows.append(('subdivision',str(r['source_index'])+':'+str(r['subdivisions']),r['row'],None))
 assert len(families['rank_at_most29_patterns'])==191
 one=read('results/fresh_sqrt5_rooted_component_completion_20261015.json')
 assert len(one['rank29_instances'])==1715
 for r in one['rank29_instances']:
  assert r['row']['nodes']<=31 and r['rank29_orthogonal_spectator_rank']==31-r['row']['nodes']
  rows.append(('one_touch',r['index'],r['row'],None))
 dotted=read('results/fresh_sqrt5_dotted_rank29_exhaustion_20261015.json')
 assert dotted['status']=='COMPLETE_DOTTED_ARITHMETIC_ATTACHMENT_QUEUE_THROUGH_DIMENSION29'and dotted['processed']==dotted['marked_states']
 for r in dotted['records']:
  if r['index']in dotted['weight_obligations']:
   assert r['result']['singular'];rows.append(('dotted',r['index'],r['row'],r['ratio']))
 return rows

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
 raw=inputs();pool=defaultdict(list);faces=[];literal_maps=[];orbitdata={}
 for family,index,row,parameter in raw:
  M=S.matrix(row);inv=(S.invariant(M),tuple(parameter)if parameter is not None else None);hit=None
  for j,N in pool[inv]:
   mp=S.isomorphism(M,N)
   if mp is not None:S.check_map(M,N,mp);hit=(j,mp);break
  if hit is None:
   j=len(faces);pool[inv].append((j,M));mp=list(range(row['nodes']))
   faces.append(dict(family=family,index=index,row=row,dotted_ratio=parameter,corners=corners(row,orbitdata),origins=[]))
  else:j,mp=hit
  origin=dict(family=family,index=index,to_representative=mp);faces[j]['origins'].append(origin)
  literal_maps.append(dict(**origin,representative=j))
 orbitkeys=sorted(orbitdata);variables=[(k,c)for k in orbitkeys for c in range(3)];index={v:i for i,v in enumerate(variables)}
 constraints=[];unique={}
 for fi,face in enumerate(faces):
  for outside in range(min(2,29-(face['row']['nodes']-2))+1):
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
  if len(g)>29:continue
  baseline=sum((G.base(g,a,b)for a,b in combinations(g,2)),F(0));coef=[Counter()for _ in range(3)];orbits=defaultdict(list)
  for a,b in combinations(g,2):
   key=B.code(g,(a,b));orbits[key].append([a,b])
   if key in orbitdata:
    assert G.base(g,a,b)==0
    for c in range(3):coef[c][index[key,c]]+=1
  components.append(dict(type=name,rank=len(g),baseline=baseline,coeffs=coef,affected_orbits=[dict(key=k,multiplicity=len(v),pairs=v)for k,v in sorted(orbits.items())if k in orbitdata]))
 assert len(components)==92
 assert{k for c in components for k in[c['affected_orbits'][i]['key']for i in range(len(c['affected_orbits']))]}==set(orbitkeys)
 return dict(faces=faces,faceconstraints=constraints,orbits=[dict(key=k,**orbitdata[k])for k in orbitkeys],variables=variables,components=components,rows=rows,rhs=rhs,meta=meta,zi=len(variables),literal_maps=literal_maps,literal_count=len(raw))

def rank29_check(data,solution):
 fronts=C.maxforests(data,list(map(F,solution['exact_primal'])),29)
 maximum=max(r['maximum']for r in fronts)
 return dict(rank=29,target=84,maximum=str(maximum),meets_target=maximum<=84,
   fronts=[dict(context=r['context'],maximum=str(r['maximum']),maximizing_forest=[data['components'][i]['type']for i in r['forest']])for r in fronts],
   scope='Same literal rule evaluated on all92 allowed finite types, including possibly nonzero A29/B29/D29 bonuses; no independent rank29 optimization.')

def elementary_barrier(data):
 selected=[('old112',i)for i in(940,958,982,984,993,995,1107,1143,1145)]+[('dotted',107)]
 aggregate=Counter();out=[]
 for family,index in selected:
  fi=next(i for i,r in enumerate(data['faces'])if any(o['family']==family and o['index']==index for o in r['origins']))
  face=data['faces'][fi];constraint=next(r for r in data['faceconstraints']if r['face']==fi and r['external_A1_count_capped2']==0)
  assert F(constraint['baseline'])==0
  for j,c in constraint['coefficients']:
   assert data['variables'][j][1]==0;aggregate[j]+=c
  out.append(dict(family=family,index=index,face=fi,coefficients=constraint['coefficients'],core_corner_forests=[r['core_forest']for r in face['corners']],dotted_ratio=face['dotted_ratio']))
 e8=next(r for r in data['components']if r['type']=='E8');h4=next(r for r in data['components']if r['type']=='H4')
 assert len(out)==10 and len(aggregate)==11 and set(aggregate.values())=={1}
 assert all(e8['coeffs'][0][j]>=v for j,v in aggregate.items())
 assert e8['baseline']==14 and h4['baseline']==7 and 3*8+4==28
 return dict(face_constraints=out,aggregate_bonus_lower_bound=10,E8_baseline='14',E8_forced_total='24',forest=['E8','E8','E8','H4'],rank=28,forced_forest_total='79',
   scope='Exact obstruction to this fixed-baseline connected-pair-only weighting model on the stated local inputs. Dotted state107 is a necessary arithmetic local completion, not a claimed global compact realization. No dimension or polytope nonexistence claim.')

def run(saved=None):
 data=build();print('literals',data['literal_count'],'distinct',len(data['faces']),'orbits',len(data['orbits']),'vars',data['zi'],'constraints',len(data['rows']),flush=True)
 solution=C.solve(data)if saved is None else saved['solution'];proof=C.verify(data,solution);rank29=rank29_check(data,solution);barrier=elementary_barrier(data)
 assert F(solution['optimum'])==F(barrier['forced_forest_total'])
 paths=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_enlarged_a1_weights_20261015.py','results/fresh_sqrt5_enlarged_a1_weights_20261015.json','results/fresh_sqrt5_enlarged_a1_weights_20261015_review.json','scripts/fresh_sqrt5_a1_count_weights_20261015.py','scripts/fresh_sqrt5_sealed_subcore_closure_20261015_review.py','results/fresh_sqrt5_subdivision_families_20261015.json','results/fresh_sqrt5_rooted_component_completion_20261015.json','results/fresh_sqrt5_dotted_rank29_exhaustion_20261015.json']
 for p in paths:
  if p.endswith('.json'):
   for q,h in read(p).get('sources',{}).items():assert sha(q)==h,(p,q)
 return dict(status='PASS_EXACT_PARTIAL_PARAMETER_FAMILY_WEIGHT_STRESS',sources={p:sha(p)for p in paths},literal_inputs=data['literal_count'],distinct_marked_face_graphs=len(data['faces']),literal_maps=data['literal_maps'],face_cases=len(data['faceconstraints']),distinct_constraints=len(data['rows']),orbit_count=len(data['orbits']),variables=data['zi'],finite_component_types=92,orbits=data['orbits'],variable_order=[dict(orbit=k,A1_count_capped2=c)for k,c in data['variables']],faces=data['faces'],face_constraints=data['faceconstraints'],solution=solution,exact_certificate=proof,target=77,meets_target=F(solution['optimum'])<=77,rank29_same_rule=rank29,elementary_E8_barrier=barrier,
  scope='Conservative partial-face diagnostic: previous224 marked graphs,191 rank<=29 subdivision patterns,1715 one-touch component patterns and every singular obligation in the complete fixed-parameter dotted queue through dimension29. Dotted ratios retained in deduplication. Only external A1 counts compatible with total rank<=29 imposed. Subdivision and one-touch families have no arbitrary-extension sealing theorem. No complete finite-family coverage or dimension theorem.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();p=ROOT/'results'/f'{TAG}.json';saved=json.loads(p.read_text())if a.verify else None;out=json.loads(json.dumps(run(saved)))
 if a.verify:assert out==saved
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['solution']['optimum'],'rank29',out['rank29_same_rule']['maximum'])
