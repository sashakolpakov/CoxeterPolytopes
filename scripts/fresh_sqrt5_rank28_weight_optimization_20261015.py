#!/usr/bin/env python3
"""Exact A1-sensitive rank28 optimization for174 known face graphs.
SciPy proposes coefficients only. Fraction primal/dual checks and complete
finite-forest separation certify the optimum of this partial weighting system.
"""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
import fresh_sqrt5_forest_sensitive_weights_20261015_review as G
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_rank28_weight_optimization_20261015';BOUND=28

def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def code(g,marks):return json.dumps(G.code(g,marks),separators=(',',':'))
def edges(r):return dict(zip(combinations(range(r['nodes']),2),r['labels']))

def build():
 base=read('results/fresh_sqrt5_quad_exhaustion_20261015.json')['labelled_union_supports'];old=read('results/fresh_sqrt5_forest_sensitive_weights_20261015.json');first=read('results/fresh_sqrt5_singular_extension_frontier_20261015.json');closure=read('results/fresh_sqrt5_singular_attachment_closure_20261015.json')
 inputs=[('old112',r['class_index'],base[r['union_index']])for r in old['sealed_quad_rows']]
 inputs +=[('new46',r['index'],r['row'])for r in first['records']if not r['accepted']]
 inputs +=[('closure11',r['index'],r['row'])for r in closure['records']if r['status']=='NEW_SEALED_PAID']
 inputs +=[('zero_touch5',r['index'],r['row'])for r in first['records']if r['accepted']]
 assert Counter(fam for fam,_,_ in inputs)=={'old112':112,'new46':46,'closure11':11,'zero_touch5':5}
 rawfaces=[];orbitdata={}
 for fam,i,r in inputs:
  cornerdata=[];es=edges(r)
  for a,b in product((0,2),(1,3)):
   cs=G.pieces(G.graph(r['white']+[a,b],es));C=next(c for c in cs if a in c);has_a1=any(len(c)==1 for c in cs)
   connected=b in C;baseline=G.base(C,a,b)if connected else F(0);key=None
   if connected and baseline==0:
    assert G.distances(C,a)[b]>=3
    key=code(C,(a,b));ids=sorted(C)
    if key not in orbitdata:orbitdata[key]=dict(type=G.family(C),rank=len(C),nodes=ids,labels=[C[u].get(v,2)for u,v in combinations(ids,2)],marked=[a,b])
   cornerdata.append(dict(marked=[a,b],core_forest=sorted(map(G.family,cs)),core_has_a1=has_a1,baseline=str(baseline),orbit=key))
  rawfaces.append(dict(family=fam,index=i,corners=cornerdata))
 orbitkeys=sorted(orbitdata);variables=[(key,flag)for key in orbitkeys for flag in(0,1)];index={v:i for i,v in enumerate(variables)};zi=len(variables)
 faceconstraints=[];unique={}
 for fi,r in enumerate(rawfaces):
  for outside in(0,1):
   co=Counter();baseweight=F(0)
   for c in r['corners']:
    baseweight+=F(c['baseline'])
    if c['orbit']is not None:co[index[c['orbit'],int(c['core_has_a1']or outside)]]+=1
   key=(tuple(sorted(co.items())),baseweight-1)
   unique.setdefault(key,[]).append(dict(face=fi,family=r['family'],index=r['index'],external_A1=bool(outside)))
   faceconstraints.append(dict(family=r['family'],index=r['index'],external_A1=bool(outside),baseline=str(baseweight),coefficients=[[j,v]for j,v in sorted(co.items())]))
 rows=[];rhs=[];meta=[]
 for(co,b),witnesses in sorted(unique.items()):
  rows.append({j:-v for j,v in co});rhs.append(b);meta.append(dict(kind='face',witnesses=witnesses))
 components=[]
 for name,H in G.models():
  if len(H)>BOUND:continue
  baseweight=sum((G.base(H,a,b)for a,b in combinations(H,2)),F(0));coeffs=[Counter(),Counter()];orbits=defaultdict(list)
  for a,b in combinations(H,2):
   key=code(H,(a,b));orbits[key].append([a,b])
   if key in orbitdata:
    assert G.base(H,a,b)==0
    for flag in(0,1):coeffs[flag][index[key,flag]]+=1
  components.append(dict(type=name,rank=len(H),baseline=baseweight,coeffs=coeffs,
                         affected_orbits=[dict(key=k,multiplicity=len(v),pairs=v)for k,v in sorted(orbits.items())if k in orbitdata]))
 assert len(components)==89
 assert{key for c in components for key in[c['affected_orbits'][i]['key']for i in range(len(c['affected_orbits']))]}==set(orbitkeys)
 return dict(faces=rawfaces,faceconstraints=faceconstraints,orbits=[dict(key=k,**orbitdata[k])for k in orbitkeys],variables=variables,components=components,rows=rows,rhs=rhs,meta=meta,zi=zi)

def maxforests(data,x):
 out=[]
 for context in(0,1):
  values=[F(0)]+[None]*BOUND;forests=[()]+[None]*BOUND
  totals=[c['baseline']+sum(v*x[j]for j,v in c['coeffs'][context].items())for c in data['components']]
  for n in range(1,BOUND+1):
   options=[(values[n-c['rank']]+totals[i],tuple(sorted(forests[n-c['rank']]+(i,))))for i,c in enumerate(data['components'])if c['rank']<=n and values[n-c['rank']]is not None and(context or c['rank']!=1)]
   if options:values[n],forests[n]=max(options)
  # Context1 requires actual A1; its own pair contribution is zero.
  rank=BOUND if context==0 else BOUND-1
  a1=next(i for i,c in enumerate(data['components'])if c['type']=='A1')
  forest=forests[rank]if context==0 else tuple(sorted(forests[rank]+(a1,)))
  out.append(dict(context=context,maximum=values[rank],forest=forest,rank_values=values,component_totals=totals))
 return out

def forestrow(data,flag,forest):
 assert sum(data['components'][i]['rank']for i in forest)==BOUND
 assert any(data['components'][i]['type']=='A1'for i in forest)==bool(flag)
 row=Counter({data['zi']:-1});b=F(0)
 for i in forest:
  c=data['components'][i];row.update(c['coeffs'][flag]);b-=c['baseline']
 return dict(row),b

def elementary_barrier(data):
 parts=[]
 for name,classes in [('D6',(1018,543,545)),('D7',(535,563,843,967,1032,1052))]:
  comp=next(c for c in data['components']if c['type']==name);seen=set();terms=[]
  for ci in classes:
   face=next(f for f in data['faceconstraints']if f['family']=='old112'and f['index']==ci and f['external_A1'])
   assert F(face['baseline'])==0 and len(face['coefficients'])==1
   j,v=face['coefficients'][0];assert v==1 and j not in seen;seen.add(j)
   orbit,flag=data['variables'][j];assert flag==1
   multiplicity=comp['coeffs'][1][j];assert multiplicity>0
   literal=next(f for f in data['faces']if f['family']=='old112'and f['index']==ci)
   corner=next(c for c in literal['corners']if c['orbit']==orbit)
   assert name in corner['core_forest']
   terms.append(dict(class_index=ci,variable=j,minimum_bonus=1,multiplicity=multiplicity,forcing_corner_forest=corner['core_forest'],minimum_A1_spectators=corner['core_forest'].count('A1')))
  extra=sum(t['multiplicity']for t in terms);parts.append(dict(type=name,baseline=str(comp['baseline']),forced_extra=extra,forced_total=str(comp['baseline']+extra),terms=terms))
 assert [p['forced_total']for p in parts]==['15','21']
 assert 1+6+3*7==BOUND and F(parts[0]['forced_total'])+3*F(parts[1]['forced_total'])==78
 return dict(forest=['A1','D6','D7','D7','D7'],rank=BOUND,lower_bound='78',parts=parts,
  limitation='The D6 class1018 forcing corner has two entireA1 spectators, whereas this maximizing forest has only one. A boolean A1 flag identifies those contexts; a finer spectator-count rule is not excluded.')

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
 data=build();print('facegraphs',len(data['faces']),'orbits',len(data['orbits']),'variables',data['zi'],'distinctfaces',len(data['rows']),flush=True)
 solution=solve(data)if saved is None else saved['solution'];proof=verify(data,solution);barrier=elementary_barrier(data)
 assert F(solution['optimum'])==F(barrier['lower_bound'])
 paths=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py']
 for stem in('fresh_sqrt5_quad_exhaustion_20261015','fresh_sqrt5_forest_sensitive_weights_20261015','fresh_sqrt5_forest_sensitive_weights_20261015_review','fresh_sqrt5_singular_extension_frontier_20261015','fresh_sqrt5_extended_forest_weights_20261015_review','fresh_sqrt5_singular_attachment_closure_20261015','fresh_sqrt5_singular_attachment_closure_20261015_review','fresh_sqrt5_triangle_weights_20261015_review'):
  paths.append('results/'+stem+'.json')
 for p in paths:
  if p.endswith('.json'):
   for q,h in read(p).get('sources',{}).items():assert sha(q)==h,(p,q)
 return dict(status='PASS_EXACT_PARTIAL_A1_SENSITIVE_RANK28_WEIGHT_OPTIMUM',sources={p:sha(p)for p in paths},rank=28,face_graphs=174,face_cases=348,distinct_face_constraints=len(data['rows']),orbit_count=len(data['orbits']),variables=data['zi'],finite_component_types=len(data['components']),orbits=data['orbits'],variable_order=[dict(orbit=k,full_vertex_has_A1=bool(flag))for k,flag in data['variables']],faces=data['faces'],face_constraints=data['faceconstraints'],solution=solution,exact_certificate=proof,elementary_lower_bound=barrier,target=77,meets_target=F(solution['optimum'])<=77,
 rule='Distance<=2/H4 baseline retained. Independent nonnegative bonuses for each included long marked-pair orbit in each full-forest A1 context; every other pair bonus is zero. Disconnected pairs have weight zero.',
 scope='Exact optimum only for the174 specified face graphs with both external-A1 cases and all allowed rank28 finite forests. Complete remaining finite/dotted quadrilateral coverage is not established; no dimension theorem. Five nonsealed graphs are constrained in their zero-touch case, with their other attachments covered by the reviewed closure.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();p=ROOT/'results'/f'{TAG}.json';saved=json.loads(p.read_text())if a.verify else None;out=json.loads(json.dumps(run(saved)))
 if a.verify:assert out==saved
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['solution']['optimum'],'target',out['target'],'meets',out['meets_target'])
