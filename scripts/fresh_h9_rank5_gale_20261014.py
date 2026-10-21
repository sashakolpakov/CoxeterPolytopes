#!/usr/bin/env python3
"""Resumable bounded direct rank-five Gale walk with exact retained cases.

Annealing only selects proposals. Every recorded nerve and every gate are exact.
The walk may pass through large-missing-face states; those are never candidates.
"""
import argparse,json,math,os,random,time
from collections import Counter
import numpy as np
from fresh_h9_rank5_gale_20261014_core import ROOT,TAG,incidence,realize,simplex_link_factors,LocalGates,sha
from fresh_h9_fifteen_facet_20261014_engine import search
def tuples(x):return tuple(map(tuples,x))if isinstance(x,list)else x
def atomic(path,d):
 p=ROOT/path;tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(d,separators=(',',':'))+'\n');tmp.replace(p)
def run(trials,walks,steps,limit):
 gates=LocalGates();casefolder=ROOT/f'results/{TAG}_cases';casefolder.mkdir(exist_ok=True)
 files=[f'scripts/{TAG}.py',f'scripts/{TAG}_core.py','scripts/fresh_h9_fifteen_facet_20261014_engine.py',
  'scripts/fresh_h9_bugaenko_surgery_20261014_incidence.py','scripts/fresh_h9_bugaenko_surgery_20261014_incidence_independent.py',
  'scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/fresh_h9_bugaenko_closure_20261014.py']+gates.paths
 sources={p:sha(p)for p in files};outpath=f'results/{TAG}.json';statepath=f'results/{TAG}_state.json';start=time.monotonic()
 options=dict(random_seed=2026101406,random_trials=trials,walks=walks,steps_per_walk=steps,label_node_limit=limit)
 rows=[];seen=set();counts=Counter();rng=random.Random(options['random_seed']);pool=[];phase='random';index=0;wi=0;V=None;energy=None;survivor=None
 if(ROOT/statepath).exists():
  old=json.loads((ROOT/statepath).read_text());assert old['sources']==sources and old['options']==options
  rng.setstate(tuples(old['random_state']));counts=Counter(old['counts']);pool=old['pool'];phase=old['phase'];index=old['index'];wi=old['walk'];energy=old['energy'];survivor=old['survivor']
  V=np.array(old['current'],dtype=np.int64)if old['current']is not None else None
  for path in sorted(casefolder.glob('case_*.json')):
   row=json.loads(path.read_text());assert row['sources']==sources;rows.append(row);seen.add(tuple(row['geometry']['canonical_minimal_nonfaces']))
 def save(status):
  refs=[dict(case=r['case'],path=f'results/{TAG}_cases/case_{r["case"]:05d}.json',sha256=sha(f'results/{TAG}_cases/case_{r["case"]:05d}.json'),gate=r['gate']['gate'],verdict=r['gate'].get('search',{}).get('verdict'))for r in rows]
  data=dict(status=status,sources=sources,options=options,counts=dict(counts),distinct_exact_nerves=len(rows),first_full_category_survivor=survivor,
   phase=phase,phase_index=index,walk_index=wi,cases=refs,
   scope='Bounded direct15x5 integral Gale proposals; all saved5minors/circuit incidences are exact. Every saved nerve has minimal nonfaces<=5 and no simplex-boundary join factor in any vertex link, hence is not a standard one-face truncation. No exhaustive Gale or fifteen-facet claim. Only complete exact gates exclude a saved nerve; node-limited searches remain open.')
  atomic(outpath,data)
  atomic(statepath,dict(status=status,sources=sources,options=options,counts=dict(counts),random_state=rng.getstate(),pool=pool,phase=phase,index=index,walk=wi,current=V.tolist()if V is not None else None,energy=energy,survivor=survivor))
  return data
 def evaluate(matrix,origin):
  nonlocal survivor
  data,bad=incidence(matrix,require_small=False)
  if bad:counts[bad]+=1;return None
  mm=data['minimal_nonfaces'];large=sum(4**(m.bit_count()-5)-1 for m in mm if m.bit_count()>5)
  factors=simplex_link_factors(data['facets'])if not large else[]
  score=7*large+2*len(factors)+sum(.3 if m.bit_count()==4 else .1 if m.bit_count()==5 else -.1 if m.bit_count()==2 else 0 for m in mm)
  if large:counts['large_minimal_nonface']+=1;return(score,data)
  if factors:counts['simplex_link_factor']+=1;return(score,data)
  counts['admissible_novel_family_proposals']+=1
  row=realize(matrix,data);key=tuple(row['canonical_minimal_nonfaces'])
  if key in seen:return(score,data)
  seen.add(key);case=len(rows);row['case']=case;gate=gates.check(row['minimal_nonfaces'])
  if gate is None:
   print('LABEL_START',case,'origin',origin,'vertices',len(row['facets']),flush=True)
   result=search(row,limit);gate=dict(gate='complete_label_search',search=result)
   print('LABEL_END',case,result['verdict'],result['nodes'],flush=True)
   if result['verdict']=='FULL_CATEGORY_METRIC_REQUIRED':survivor=case
  counts[gate['gate']]+=1
  record=dict(case=case,sources=sources,origin=origin,geometry=row,gate=gate)
  atomic(f'results/{TAG}_cases/case_{case:05d}.json',record);rows.append(record)
  print('NEW',case,'origin',origin,'vertices',len(row['facets']),'missing',row['missing_size_counts'],'gate',gate['gate'],'seconds',round(time.monotonic()-start,1),flush=True)
  return(score,data)
 print('PID',os.getpid(),'resumed',phase,index,'cases',len(rows),'options',options,flush=True)
 if survivor is not None:return save('STOPPED_AT_FULL_CATEGORY_SURVIVOR')
 while phase=='random'and index<trials:
  trial=np.array([[rng.randint(-64,64)for j in range(5)]for i in range(14)],dtype=np.int64);trial=np.vstack([trial,-trial.sum(axis=0)])
  result=evaluate(trial,dict(phase='random',index=index));counts['random_proposals']+=1
  if result is not None:
   score,data=result;key=data['minimal_nonfaces']
   if not any(row['missing']==key for row in pool):
    pool.append(dict(score=score,V=trial.tolist(),missing=key));pool.sort(key=lambda r:(r['score'],r['missing']));pool=pool[:max(walks,8)]
  index+=1
  if index%250==0 or survivor is not None:save('RUNNING_RANDOM_INITIALIZATION');print('random',index,'cases',len(rows),'best',pool[0]['score']if pool else None,flush=True)
  if survivor is not None:return save('STOPPED_AT_FULL_CATEGORY_SURVIVOR')
 if phase=='random':assert pool;phase='walk';index=0;wi=0;V=None
 while wi<walks:
  if V is None:
   seed=pool[wi%len(pool)];V=np.array(seed['V'],dtype=np.int64);energy=seed['score'];index=0
  while index<steps:
   trial=V.copy();i,j=rng.sample(range(15),2);col=rng.randrange(5);delta=rng.choice([-1,1])*rng.choice([1,2,4,8,16,32,64,128])
   trial[i,col]+=delta;trial[j,col]-=delta;counts['walk_proposals']+=1
   if np.max(np.abs(trial))>2048:counts['coordinate_bound']+=1
   else:
    result=evaluate(trial,dict(phase='walk',walk=wi,index=index))
    if result is not None:
     score,data=result;temperature=2+6*(1-(index%2000)/2000)
     if score<=energy or rng.random()<math.exp((energy-score)/temperature):V=trial;energy=score;counts['accepted_walk_moves']+=1
   index+=1
   if index%250==0 or survivor is not None:
    save('RUNNING_BOUNDED_MUTATION_WALK');print('walk',wi,index,'cases',len(rows),'energy',round(energy,3),'seconds',round(time.monotonic()-start,1),flush=True)
   if survivor is not None:return save('STOPPED_AT_FULL_CATEGORY_SURVIVOR')
  wi+=1;index=0;V=None;energy=None;save('RUNNING_BOUNDED_MUTATION_WALK')
 phase='complete';return save('FINISHED_BOUNDED_DIRECT_RANK5_GALE_SEARCH')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--trials',type=int,default=3000);ap.add_argument('--walks',type=int,default=4);ap.add_argument('--steps',type=int,default=10000);ap.add_argument('--nodes',type=int,default=3000);a=ap.parse_args()
 d=run(a.trials,a.walks,a.steps,a.nodes);print(d['status'],d['distinct_exact_nerves'],d['counts'],flush=True)
