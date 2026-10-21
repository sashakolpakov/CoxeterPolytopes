#!/usr/bin/env python3
"""Exact partial rank28 weights with zero/one/at-least-two A1 components.
The face inventory and fixed baseline come from the frozen Boolean diagnostic.
The LP engine is copied read-only in provenance; all certificates are exact.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
import fresh_sqrt5_rank28_weight_optimization_20261015 as B
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_a1_count_weights_20261015';BOUND=28

def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def build():
 old=B.build();orbitkeys=[r['key']for r in old['orbits']];variables=[(k,c)for k in orbitkeys for c in range(3)];index={v:i for i,v in enumerate(variables)};zi=len(variables)
 constraints=[];unique={}
 for fi,f in enumerate(old['faces']):
  for outside in range(3):
   co=Counter();baseline=F(0)
   for corner in f['corners']:
    baseline+=F(corner['baseline'])
    if corner['orbit']is not None:
     context=min(2,corner['core_forest'].count('A1')+outside);co[index[corner['orbit'],context]]+=1
   key=(tuple(sorted(co.items())),baseline-1);unique.setdefault(key,[]).append(dict(face=fi,family=f['family'],index=f['index'],external_A1_count_capped2=outside))
   constraints.append(dict(family=f['family'],index=f['index'],external_A1_count_capped2=outside,baseline=str(baseline),coefficients=[[j,v]for j,v in sorted(co.items())]))
 rows=[];rhs=[];meta=[]
 for(co,b),w in sorted(unique.items()):rows.append({j:-v for j,v in co});rhs.append(b);meta.append(dict(kind='face',witnesses=w))
 components=[]
 for c in old['components']:
  coef=[Counter()for _ in range(3)]
  for o in c['affected_orbits']:
   for context in range(3):coef[context][index[o['key'],context]]+=o['multiplicity']
  components.append(dict(type=c['type'],rank=c['rank'],baseline=c['baseline'],coeffs=coef,affected_orbits=c['affected_orbits']))
 return dict(faces=old['faces'],faceconstraints=constraints,orbits=old['orbits'],variables=variables,components=components,rows=rows,rhs=rhs,meta=meta,zi=zi)

def maxforests(data,x,bound=BOUND):
 out=[]
 for context in range(3):
  totals=[c['baseline']+sum(v*x[j]for j,v in c['coeffs'][context].items())for c in data['components']]
  dp={(0,0):(F(0),())}
  for n in range(1,bound+1):
   for i,c in enumerate(data['components']):
    if c['rank']>n:continue
    increment=int(c['type']=='A1')
    for count in range(3):
     old=dp.get((n-c['rank'],count))
     if old is None:continue
     key=(n,min(2,count+increment));v=(old[0]+totals[i],tuple(sorted(old[1]+(i,))))
     if key not in dp or v>dp[key]:dp[key]=v
  value,forest=dp[bound,context]
  out.append(dict(context=context,maximum=value,forest=forest,component_totals=totals,
                  rank_values=[None if(n,context)not in dp else dp[n,context][0]for n in range(bound+1)]))
 return out

def rank29_check(data,solution):
 # The same literal rule, without reoptimization. Rank29 adds exactly A29,
 # B29 and D29 to the complete finite-component list. Their bonuses are zero.
 enlarged=dict(data);enlarged['components']=list(data['components']);added=[]
 orbitkeys={r['key']for r in data['orbits']}
 for name,g in B.G.models():
  if len(g)!=29:continue
  assert all(B.code(g,(a,b))not in orbitkeys for a,b in combinations(g,2))
  baseline=sum((B.G.base(g,a,b)for a,b in combinations(g,2)),F(0))
  c=dict(type=name,rank=29,baseline=baseline,coeffs=[Counter()for _ in range(3)],affected_orbits=[])
  enlarged['components'].append(c);added.append(dict(type=name,rank=29,baseline=str(baseline)))
 assert sorted(c['type']for c in added)==['A29','B29','D29']
 fronts=maxforests(enlarged,list(map(F,solution['exact_primal'])),29)
 maximum=max(r['maximum']for r in fronts)
 return dict(rank=29,target=84,maximum=str(maximum),meets_target=maximum<=84,
  finite_component_types=len(enlarged['components']),added_baseline_only_components=added,
  fronts=[dict(context=r['context'],maximum=str(r['maximum']),maximizing_forest=[enlarged['components'][i]['type']for i in r['forest']])for r in fronts],
  scope='Exact evaluation of the rank28-optimal rule on every allowed rank29 forest; no rank29 optimization or complete face-coverage claim.')

def forestrow(data,flag,forest):
 assert sum(data['components'][i]['rank']for i in forest)==BOUND
 assert min(2,sum(data['components'][i]['type']=='A1'for i in forest))==flag
 row=Counter({data['zi']:-1});b=F(0)
 for i in forest:
  c=data['components'][i];row.update(c['coeffs'][flag]);b-=c['baseline']
 return dict(row),b

def solve(data):
 rows=[dict(r)for r in data['rows']];rhs=list(data['rhs']);meta=list(data['meta']);cuts=[];cutkeys=set();zi=data['zi'];iterations=[]
 def add(flag,forest):
  key=(flag,tuple(forest));assert key not in cutkeys;cutkeys.add(key);row,b=forestrow(data,flag,forest);rows.append(row);rhs.append(b);m=dict(kind='forest',context=flag,components=list(forest));meta.append(m);cuts.append(m)
 for rec in maxforests(data,[F(0)]*zi):add(rec['context'],rec['forest'])
 for iteration in range(2000):
  rr=[];cc=[];vv=[]
  for i,r in enumerate(rows):
   for j,v in r.items():rr.append(i);cc.append(j);vv.append(float(v))
  mat=coo_matrix((vv,(rr,cc)),shape=(len(rows),zi+1)).tocsr();objective=np.zeros(zi+1);objective[-1]=1
  result=linprog(objective,A_ub=mat,b_ub=np.array([float(b)for b in rhs]),bounds=(0,None),method='highs');assert result.success,result.message
  x=[F(float(v)).limit_denominator(10**7)for v in result.x];assert min(x)>=0
  slack=[b-sum(v*x[j]for j,v in r.items())for r,b in zip(rows,rhs)];assert min(slack)>=0,('rational primal',min(slack))
  fronts=maxforests(data,x);iterations.append(dict(iteration=iteration,objective=str(x[zi]),maxima=[str(r['maximum'])for r in fronts],forests=[[data['components'][i]['type']for i in r['forest']]for r in fronts]))
  print('LP',iteration,'z',x[zi],'max',[r['maximum']for r in fronts],flush=True)
  bad=[r for r in fronts if r['maximum']>x[zi]]
  if bad:
   for r in bad:add(r['context'],r['forest'])
   continue
  y=[F(float(v)).limit_denominator(10**7)for v in result.ineqlin.marginals];assert all(v<=0 for v in y)
  dual=[F(0)]*(zi+1)
  for v,r in zip(y,rows):
   for j,w in r.items():dual[j]+=v*w
  assert all(v<=(1 if j==zi else 0)for j,v in enumerate(dual))
  assert sum(v*b for v,b in zip(y,rhs))==x[zi]
  return dict(optimum=str(x[zi]),exact_primal=list(map(str,x)),forest_constraints=cuts,nonzero_dual=[dict(row=i,multiplier=str(v),constraint=meta[i])for i,v in enumerate(y)if v],iterations=iterations)
 raise AssertionError('Unfinished cutting-plane frontier: iteration guard reached')

def verify(data,solution):
 rows=[dict(r)for r in data['rows']];rhs=list(data['rhs']);meta=list(data['meta']);seen=set()
 for c in solution['forest_constraints']:
  flag=c['context'];forest=tuple(c['components']);assert(flag,forest)not in seen;seen.add((flag,forest));row,b=forestrow(data,flag,forest);rows.append(row);rhs.append(b);meta.append(c)
 x=list(map(F,solution['exact_primal']));zi=data['zi'];assert len(x)==zi+1 and min(x)>=0
 slack=[b-sum(v*x[j]for j,v in r.items())for r,b in zip(rows,rhs)];assert min(slack)>=0
 y=[F(0)]*len(rows)
 for v in solution['nonzero_dual']:
  i=v['row'];assert y[i]==0 and v['constraint']==meta[i];y[i]=F(v['multiplier']);assert y[i]<0
 dual=[F(0)]*(zi+1)
 for v,r in zip(y,rows):
  for j,w in r.items():dual[j]+=v*w
 assert all(v<=(1 if j==zi else 0)for j,v in enumerate(dual))
 primal=x[zi];dv=sum(v*b for v,b in zip(y,rhs));assert primal==dv==F(solution['optimum'])
 fronts=maxforests(data,x);assert max(r['maximum']for r in fronts)==primal
 return dict(primal=str(primal),dual=str(dv),minimum_primal_slack=str(min(slack)),dual_constraints_checked=len(dual),constraints=len(rows),
   fronts=[dict(context=r['context'],maximum=str(r['maximum']),maximizing_forest=[data['components'][i]['type']for i in r['forest']],rank_maxima=[None if v is None else str(v)for v in r['rank_values']])for r in fronts],
   component_totals=[dict(type=c['type'],rank=c['rank'],baseline=str(c['baseline']),totals=[str(r['component_totals'][i])for r in fronts])for i,c in enumerate(data['components'])])

def run(saved=None):
 data=build();print('facegraphs',len(data['faces']),'variables',data['zi'],'distinctcases',len(data['rows']),flush=True)
 solution=solve(data)if saved is None else saved['solution'];proof=verify(data,solution);rank29=rank29_check(data,solution)
 paths=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_rank28_weight_optimization_20261015.py','results/fresh_sqrt5_rank28_weight_optimization_20261015.json','results/fresh_sqrt5_singular_attachment_closure_20261015_review.json','results/fresh_sqrt5_extended_forest_weights_20261015_review.json']
 for p in paths:
  if p.endswith('.json'):
   for q,h in read(p).get('sources',{}).items():assert sha(q)==h,(p,q)
 return dict(status='PASS_EXACT_PARTIAL_A1_COUNT_RANK28_WEIGHT_OPTIMUM',sources={p:sha(p)for p in paths},rank=28,A1_contexts=['zero','one','at_least_two'],face_graphs=174,face_cases=522,distinct_face_constraints=len(data['rows']),orbit_count=len(data['orbits']),variables=data['zi'],finite_component_types=len(data['components']),orbits=data['orbits'],variable_order=[dict(orbit=k,A1_count_capped2=c)for k,c in data['variables']],faces=data['faces'],face_constraints=data['faceconstraints'],solution=solution,exact_certificate=proof,target=77,meets_target=F(solution['optimum'])<=77,rank29_same_rule=rank29,
 rule='Fixed distance<=2/H4 baseline and independent nonnegative long-orbit bonuses for each capped entireA1 count0,1,2. Unlisted pair bonuses and disconnected pair weights are zero.',
 scope='Exact optimum for174 specified face graphs and all rank28 finite forests. The three externalA1-count cases cover arbitrary orthogonal spectators for reviewed sealed stars and the five closed starts. Remaining finite/dotted face coverage is incomplete; no dimension theorem.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();p=ROOT/'results'/f'{TAG}.json';saved=json.loads(p.read_text())if a.verify else None;out=json.loads(json.dumps(run(saved)))
 if a.verify:assert out==saved
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['solution']['optimum'],'target',out['target'],'meets',out['meets_target'])
