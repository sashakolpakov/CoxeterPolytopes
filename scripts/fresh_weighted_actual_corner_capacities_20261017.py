#!/usr/bin/env python3
"""Complete nonnegative orbit-weight separation in frozen compatibility graphs.

Local LPs propose certificates only during generation. Saved verification
enumerates the same literal cliques and checks exact charges and covering
mixtures without optimization. Unknown geometric pairs remain compatible.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json,time
import networkx as nx

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_weighted_actual_corner_capacities_20261017'
GEOMETRY='fresh_raw33_295_mixed_actual_corners_20261017'
MODEL='fresh_raw33_295_mixed_full_model_20261017'
SUPPORT=MODEL+'_support'
MAX_CLIQUES=1000000
ENUM_SECONDS=30

@lru_cache(None)
def read(tag):return json.loads((ROOT/f'results/{tag}.json').read_text())
@lru_cache(None)
def digest(name):
 h=sha256()
 with (ROOT/name).open('rb')as f:
  for block in iter(lambda:f.read(1<<20),b''):h.update(block)
 return h.hexdigest()
@lru_cache(None)
def forests(tag):return {r['vertex']:r for r in read(tag)['forest_records']}
@lru_cache(None)
def relation(tag,v):
 r=forests(tag)[v];bad={tuple(c[:2])for c in r['incompatibility_certificates']};direct=set(bad)
 assert len(bad)==len(r['incompatibility_certificates'])
 owned=[]
 for dom in r.get('inherited_pair_domains',[]):
  old=forests(dom['source'])[dom['vertex']]
  assert dom['vertex']==v and old['forest']==r['forest']
  mp=dict(dom['new_to_old_star_indices']);rev={j:i for i,j in mp.items()}
  assert len(mp)==len(rev)==len(dom['new_to_old_star_indices'])
  # A direct/new exclusion cannot override an inherited UNKNOWN outcome.
  assert all(not(i in mp and j in mp)for i,j in direct)
  for i,j in mp.items():
   for key in('contacts','occupied_pairs','squared_timelike_heights','native_omitted_product'):
    assert r['literal_stars'][i][key]==old['literal_stars'][j][key]
  add={tuple(sorted((rev[i],rev[j])))for i,j in relation(dom['source'],v)
       if i in rev and j in rev and not any(rev[i]in q and rev[j]in q for q in owned)}
  assert not(bad&add);bad|=add
  count=sum(not any(i in q and j in q for q in owned)for i,j in combinations(mp,2))
  assert count==dom['pairs_reused'] and len(add)==dom['exclusions_reused']
  owned.append(set(mp))
 n=len(r['literal_stars'])
 assert r['tested_literal_pairs']==n*(n-1)//2
 assert len(bad)==r['tested_literal_pairs']-r['unknown_literal_pairs']
 return bad

def solve_exact(rows,rhs,n):
 """Solve a consistent full-column-rank rational linear system."""
 a=[[F(x)for x in r]+[F(b)]for r,b in zip(rows,rhs)];piv=[];k=0
 for j in range(n):
  pos=next((i for i in range(k,len(a))if a[i][j]),None)
  if pos is None:continue
  a[k],a[pos]=a[pos],a[k];z=a[k][j];a[k]=[x/z for x in a[k]]
  for i in range(len(a)):
   if i!=k and a[i][j]:
    z=a[i][j];a[i]=[x-z*y for x,y in zip(a[i],a[k])]
  piv.append(j);k+=1
 assert len(piv)==n
 assert all(any(r[:n])or not r[n]for r in a)
 out=[F()for _ in range(n)]
 for i,j in enumerate(piv):out[j]=a[i][n]
 return out

def exact_lp(vectors,demand):
 import numpy as np
 from scipy.optimize import linprog
 A=np.asarray(vectors,dtype=float);scale=max(demand)
 b=np.asarray([float(x/scale)for x in demand])
 res=linprog(-b,A_ub=A,b_ub=np.ones(len(A)),bounds=(0,None),method='highs',
             options={'time_limit':20,'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
 assert res.success,(res.status,res.message)
 ids=[i for i,x in enumerate(res.x)if x>1e-8]
 tight=[i for i,x in enumerate(A@res.x)if abs(x-1)<1e-7]
 w=[F()for _ in demand]
 if ids:
  vals=solve_exact([[vectors[j][i]for i in ids]for j in tight],[1]*len(tight),len(ids))
  for i,x in zip(ids,vals):w[i]=x
 assert all(x>=0 for x in w)
 assert all(sum(x*y for x,y in zip(v,w))<=1 for v in vectors)
 # Complementary marginals give a basic covering mixture. Its support is
 # recovered against exact demands, independently of the floating objective.
 lam=-res.ineqlin.marginals
 js=[j for j,x in enumerate(lam)if x>1e-12]
 covered=A.T@lam
 eq=[i for i,x in enumerate(covered)if abs(x-b[i])<1e-8*max(1e-4,abs(b[i]))]
 vals=solve_exact([[vectors[j][i]for j in js]for i in eq],[demand[i]for i in eq],len(js))
 mixture=[[j,str(x)]for j,x in zip(js,vals)if x]
 assert all(x>=0 for x in vals)
 cost=sum(vals,F());load=sum((x*y for x,y in zip(w,demand)),F())
 assert cost==load
 assert all(sum(vectors[j][i]*x for j,x in zip(js,vals))>=demand[i]for i in range(len(demand)))
 return dict(weights=list(map(str,w)),mixture=mixture,optimal_mass=str(cost))

def check_certificate(vectors,demand,c):
 w=list(map(F,c['weights']));mix=[(j,F(x))for j,x in c['mixture']]
 assert len(w)==len(demand) and all(x>=0 for x in w)
 assert len({j for j,x in mix})==len(mix)
 assert all(0<=j<len(vectors)and x>0 for j,x in mix)
 assert all(sum(x*y for x,y in zip(v,w))<=1 for v in vectors)
 assert all(sum(vectors[j][i]*x for j,x in mix)>=demand[i]for i in range(len(demand)))
 cost=sum((x for j,x in mix),F());load=sum((x*y for x,y in zip(w,demand)),F())
 assert cost==load==F(c['optimal_mass'])
 return cost

def main(verify=False,resume=False):
 start=time.monotonic();geometry=read(GEOMETRY);support=read(SUPPORT);model=read(MODEL)
 cert=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
 assert support['face_multipliers']==cert['face_multipliers']
 assert support['vertex_multipliers']==cert['vertex_multipliers']
 alpha={i:F(x)for i,x in cert['face_multipliers']};mu={v:F(x)for v,x in cert['vertex_multipliers']}
 restricted=support['exact_restrictions'][0]['certificate']
 assert support['exact_restrictions'][0]['name']=='latest_243_actual_corner_domain_Q'
 ralpha={i:F(x)for i,x in restricted['face_multipliers']};rmu={v:F(x)for v,x in restricted['vertex_multipliers']}
 known={r['face']for r in geometry['selected_faces']}
 assert len(known)==243 and len(known&set(alpha))==184
 assert all(ralpha.get(i,F())==alpha.get(i,F())for i in known)
 source_tags={GEOMETRY,MODEL,SUPPORT}
 def collect(tag):
  source_tags.add(tag)
  for row in read(tag)['forest_records']:
   for dom in row.get('inherited_pair_domains',[]):
    if dom['source']not in source_tags:collect(dom['source'])
 collect(GEOMETRY)
 for tag in sorted(source_tags):
  doc=read(tag)
  for name,h in doc.get('sources',doc.get('prepared',{}).get('sources',{})).items():assert digest(name)==h,name
 saved=read(TAG)if verify or resume else None
 savedrows={r['vertex']:r for r in saved['forest_records']}if saved else {}
 rows=[];columns=[];opened=[];counts=Counter();newmu=dict(mu);newrmu=dict(rmu)
 work=[]
 for row in geometry['forest_records']:
  groups=[g for g in row['corner_orbits']if any(alpha.get(i,F())for i,m in g['positive_face_coefficients'])]
  if groups:work.append((row,groups))
 assert len(work)==486
 for step,(row,groups)in enumerate(work,1):
  v=row['vertex'];bad=relation(GEOMETRY,v)
  group_ids=[g['index']for g in groups]
  nodes=sorted(i for g in groups for i in g['literal_stars']);owner={i:j for j,g in enumerate(groups)for i in g['literal_stars']}
  assert len(nodes)==len(owner)
  graph=nx.Graph();graph.add_nodes_from(nodes)
  graph.add_edges_from((i,j)for i,j in combinations(nodes,2)if(i,j)not in bad)
  demand=[sum((F(m)*alpha.get(i,F())for i,m in g['positive_face_coefficients']),F())for g in groups]
  assert all(x>0 for x in demand)
  patterns={};ncliques=0;t0=time.monotonic();guard=None
  for clique in nx.find_cliques(graph):
   ncliques+=1;clique=tuple(sorted(clique));vec=[0]*len(groups)
   for i in clique:vec[owner[i]]+=1
   vec=tuple(vec)
   if vec not in patterns or clique<patterns[vec]:patterns[vec]=clique
   if ncliques>MAX_CLIQUES or time.monotonic()-t0>ENUM_SECONDS:
    guard='literal_clique_enumeration_guard';break
  if guard:
   opened.append(dict(vertex=v,reason=guard,enumerated_cliques=ncliques));relation.cache_clear();continue
  # Coordinatewise dominated vectors cannot maximize a nonnegative charge.
  vv=sorted(patterns)
  maximal=[a for a in vv if not any(a!=b and all(x<=y for x,y in zip(a,b))for b in vv)]
  representatives=[patterns[a]for a in maximal]
  try:
   if v in savedrows:
    previous=savedrows[v]
    assert previous['orbit_indices']==group_ids and previous['demands']==list(map(str,demand))
    assert previous['maximal_count_vectors']==list(map(list,maximal))
    assert previous['literal_representatives']==list(map(list,representatives))
    c=previous['certificate']
   else:
    assert not verify
    c=exact_lp(maximal,demand)
  except Exception as e:
   if verify:raise
   opened.append(dict(vertex=v,reason='exact_local_certificate_recovery',detail=str(e),vectors=len(maximal)))
   relation.cache_clear();continue
  value=check_certificate(maximal,demand,c)
  for vector,clique in zip(maximal,representatives):
   assert all(tuple(sorted(p))not in bad for p in combinations(clique,2))
   assert tuple(Counter(owner[i]for i in clique).get(j,0)for j in range(len(groups)))==vector
  w=list(map(F,c['weights']));coef=Counter()
  for weight,g in zip(w,groups):
   for i,m in g['positive_face_coefficients']:coef[i]+=weight*m
  terms=sorted((i,x)for i,x in coef.items()if x)
  load=sum((alpha.get(i,F())*x for i,x in terms),F());assert load==value
  full_bad=value>mu.get(v,F());rest_bad=value>rmu.get(v,F())
  if full_bad or rest_bad:
   columns.append(dict(index=len(columns),vertex=v,vertex_cost='1',
      positive_face_coefficients=[[i,str(x)]for i,x in terms],
      orbit_indices=group_ids,orbit_weights=c['weights'],load=str(load),
      violates_full=full_bad,violates_restricted=rest_bad))
  newmu[v]=max(newmu.get(v,F()),value);newrmu[v]=max(newrmu.get(v,F()),value)
  record=dict(vertex=v,forest=row['forest'],orbit_indices=group_ids,
    orbit_coefficients=[g['positive_face_coefficients']for g in groups],demands=list(map(str,demand)),
    literal_nodes=nodes,compatibility_edges=graph.number_of_edges(),maximal_literal_cliques=ncliques,
    distinct_count_vectors=len(patterns),maximal_count_vectors=maximal,literal_representatives=representatives,
    certificate=c,baseline_full_mu=str(mu.get(v,F())),baseline_restricted_mu=str(rmu.get(v,F())))
  rows.append(record)
  counts.update(forests=1,orbits=len(groups),literal_stars=len(nodes),maximal_literal_cliques=ncliques,
    distinct_count_vectors=len(patterns),maximal_count_vectors=len(maximal),mixture_terms=len(c['mixture']),
    full_violations=int(full_bad),restricted_violations=int(rest_bad))
  relation.cache_clear()
  if step%25==0 or step==len(work):
   print('FORESTS',step,'/',len(work),'complete',len(rows),'open',len(opened),'vectors',counts['maximal_count_vectors'],
         'full/restricted raises',counts['full_violations'],counts['restricted_violations'],flush=True)
 repairs=[]
 for label,base,oldm,newm in(('full',cert,mu,newmu),('represented_Q_restriction',restricted,rmu,newrmu)):
  assert all(newm.get(v,F())>=x for v,x in oldm.items())
  objective=F(base['objective_numerator']);normalizer=sum(newm.values(),F())
  repairs.append(dict(domain=label,baseline_bound=base['bound'],objective_numerator=str(objective),
    normalizer=str(normalizer),bound=str(objective/normalizer),
    vertex_updates=[[v,str(oldm.get(v,F())),str(x)]for v,x in sorted(newm.items())if x>oldm.get(v,F())],
    vertex_multipliers=[[v,str(x)]for v,x in sorted(newm.items())if x]))
 sources={f'results/{tag}.json':digest(f'results/{tag}.json')for tag in sorted(source_tags)}
 sources[f'scripts/{TAG}.py']=digest(f'scripts/{TAG}.py')
 for tag in(GEOMETRY,MODEL,SUPPORT):
  for folder,ext in(('scripts','py'),('results','md')):sources[f'{folder}/{tag}.{ext}']=digest(f'{folder}/{tag}.{ext}')
 result=dict(status='PASS_COMPLETE_WEIGHTED_ACTUAL_CORNER_CONE'if not opened else'PARTIAL_WEIGHTED_ACTUAL_CORNER_CONE_WITH_OPEN_FORESTS',
  sources=sources,counts=dict(counts),face_multipliers=cert['face_multipliers'],
  restricted_face_multipliers=restricted['face_multipliers'],forest_records=rows,open_forests=opened,columns=columns,repairs=repairs,
  graph_domain='All complete saved corner orbits with positive current face mass; all unexcluded literal pairs remain compatible. First-owned inherited relations, including unknowns, are retained exactly.',
  theorem='An actual simultaneous set is a literal compatibility clique. Every clique extends to a maximal clique; nonnegative orbit charges therefore have exact graph capacity given by the saved nondominated orbit-count vectors. An exact charge and equal-mass nonnegative covering mixture prove the full orbit-weight cone optimum at these fixed demands.',
  preservation='All current face masses stay fixed. Only vertex allowances increase, so all439061 inherited finite-model columns, including96176 zero-budget and52 positive-cost signed columns, remain valid. No old-alpha implicit packing/joint cover is reused.',
  scope='A complete nonnegative orbit-weight cone for the stated conservative actual-corner graph domain, when no forests are open. Graph-compatible mixtures need not be geometrically realizable. No new field/provenance coverage, global optimizer, target12 feasibility or dimension theorem.')
 result=json.loads(json.dumps(result,sort_keys=True))
 out=ROOT/f'results/{TAG}.json'
 if verify:assert result==saved
 else:out.write_text(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
 print(result['status'],dict(counts),'bounds',[(r['domain'],float(F(r['bound'])))for r in repairs],
       'seconds',round(time.monotonic()-start,2),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');p.add_argument('--resume',action='store_true')
 a=p.parse_args();assert not(a.verify and a.resume);main(a.verify,a.resume)
