#!/usr/bin/env python3
"""Apply exact necessary Coxeter gates to the saved realizable Gale nerves."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from fresh_h9_bugaenko_surgery_20261014_incidence import canonical
from fresh_h9_newfacet_frontier_20261014_engine import search
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_gale_frontier_20261014'
GEOM=f'results/{TAG}_geometry.json'
PAT='results/fresh_h9_newfacet_frontier_20261014_patterns.json'
REL='results/fresh_h9_newfacet_frontier_20261014_relations.json'
LOCAL='results/fresh_h9_newfacet_frontier_20261014_local_independent.json'
BLOCK='results/fresh_h9_bugaenko_remaining_20261014_block_gate.json'
BLOCKAUDIT='results/fresh_h9_bugaenko_remaining_20261014_block_gate_review.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def checked(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d.get('sources',{}).items():assert sha(f)==h,f
 return d

def run():
 geom=checked(GEOM);patterns=checked(PAT);relation=checked(REL);audit=checked(LOCAL);checked(BLOCK);checked(BLOCKAUDIT)
 assert audit['obstruction_count']==44 and audit['six_node_relation_rows']==96
 bad={(p['order'],tuple(p['minimal_nonfaces'])):p['pattern']for p in patterns['patterns']if p['verdict']=='EXACT_LOCAL_OBSTRUCTION'}
 bad[7,(31,99)]='seven';cache={};stats=Counter();records=[]
 for row in geom['cases']:
  mm=row['minimal_nonfaces'];found=None
  for a,b in combinations(mm,2):
   if min(a.bit_count(),b.bit_count())<4:continue
   u=a|b;inside=[c for c in mm if c&u==c]
   if a&b==0 and set(inside)=={a,b}:
    found=dict(gate='two_block',blocks=[a,b]);break
   if u.bit_count()>9 or any(m.bit_count()==2 for m in inside):continue
   vv=[i for i in range(14)if u>>i&1]
   reduced=tuple(sorted(sum(1<<k for k,j in enumerate(vv)if m>>j&1)for m in inside));key=(len(vv),reduced)
   if key not in cache:cache[key]=canonical(reduced,len(vv))
   can,perm=cache[key]
   if(len(vv),can)in bad:
    found=dict(gate='local',pattern=bad[len(vv),can],permutation=[vv[j]for j in perm]);break
  if found is None:
   found=dict(gate='complete_label_search',search=search(row,5000))
   print('complete search',row['case'],found['search']['verdict'],found['search']['nodes'],flush=True)
  verdict=found.get('search',{}).get('verdict',found['gate']);stats[verdict]+=1
  records.append(dict(case=row['case'],**found))
 assert all(r.get('search',{}).get('verdict','EXACT_COMPLETE_OBSTRUCTION')=='EXACT_COMPLETE_OBSTRUCTION'for r in records)
 paths=[GEOM,PAT,REL,LOCAL,BLOCK,BLOCKAUDIT,f'scripts/{TAG}_labels.py',
  'scripts/fresh_h9_newfacet_frontier_20261014_engine.py','scripts/fresh_h9_bugaenko_closure_20261014.py',
  'scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/independent_h9_14_record30_a1_labels.py']
 return dict(status='PASS_ALL_SAMPLED_GALE_NERVES_EXCLUDED_COMPACT_H9',sources={p:sha(p)for p in paths},
  counts=dict(stats),node_limit=5000,complete_search_states=sum(r.get('search',{}).get('nodes',0)for r in records),
  cases=records,scope='Exact necessary local constraints, not numerical fitting. Inherited44 induced-nerve obstructions and the two-block negative-index theorem are applied only after literal induced-clutter checks. All remaining complete searches enforce shared Lanner labels, every actual spherical vertex, the96-row relation and early negative-index/modular-rank gates. Categories2..6 plus7 cover every finite integer label; none of the saved complete searches reaches a limit or leaves a necessary full-label survivor. Only these saved Gale matrices are excluded; this is not a whole-shadow, whole-Gale-family, or whole14-facet theorem.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}_labels.json'
 if a.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['counts'],out['complete_search_states'])
