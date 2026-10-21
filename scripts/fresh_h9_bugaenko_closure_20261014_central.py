#!/usr/bin/env python3
"""Track every eligible recipe of the central-singleton wedge carrier."""
import argparse,json,os,time
from collections import Counter
from pathlib import Path
import fresh_h9_bugaenko_surgery_20261014_incidence as I
from fresh_h9_bugaenko_closure_20261014_search import search
from fresh_h9_bugaenko_closure_20261014 import ROOT,OLD,sha
TAG='fresh_h9_bugaenko_closure_20261014_central'
def recipes(d):
 lookup={tuple(r['canonical_minimal_nonfaces']):r['case']for r in d['cases']}
 oldV=I.masks(json.loads((ROOT/I.SOURCE).read_text())['geometry']['vertices']);wb=1;carrier=I.BASE_BLOCKS[wb][0]
 blocks=[b+((11,)if i==wb else())for i,b in enumerate(I.BASE_BLOCKS)];M=I.masks([blocks[a]+blocks[b]for a,b in I.PAIRS]);V=[]
 for v in oldV:
  if v>>carrier&1:V.append(v|(1<<11))
  else:V.extend([v|(1<<carrier),v|(1<<11)])
 out=[];total=large=0
 for F,fc,H in I.first_faces(blocks,M):
  A,B=I.stellar(M,V,F,12)
  for S,sc,new in I.second_faces(blocks,F,H,A):
   total+=1
   if any(m.bit_count()>5 and m&S!=S for m in A):large+=1;continue
   C,W=I.stellar(A,B,S,13)
   if max(map(int.bit_count,C))>5:large+=1;continue
   canon,order=I.canonical(C,14);out.append(dict(first_face=list(I.tup(F)),second_face=list(I.tup(S)),case=lookup[canon]))
 assert total==15790 and len(out)==208 and large==15582 and len({r['case']for r in out})==128
 return out,large
def run(limit,verify):
 ip=f'results/{OLD}_incidence.json';d=json.loads((ROOT/ip).read_text());rr,large=recipes(d)
 files=[ip,f'results/{OLD}_labels_0_25.json']+[f'results/fresh_h9_bugaenko_closure_20261014_search_{s}.json'for s in['0_25','25_250','500_650','524_525']]
 known={}
 for p in files[1:]:
  data=json.loads((ROOT/p).read_text())
  for source,h in data['sources'].items():assert sha(source)==h,source
  for row in data['cases']:
   if row['verdict']in['EXACT_LANNER_VERTEX_OBSTRUCTION','EXACT_COMPLETE_OBSTRUCTION']:known[row['case']]=p
 ids=sorted({r['case']for r in rr}-set(known));assert len(ids)==16
 files+=['scripts/fresh_h9_bugaenko_closure_20261014_search.py','scripts/fresh_h9_bugaenko_closure_20261014.py',f'scripts/{OLD}_incidence.py',f'scripts/{TAG}.py']
 manifest={p:sha(p)for p in files};out=dict(status='RUNNING',sources=manifest,node_limit=limit,range=[min(ids),max(ids)+1],range_is_bounding_only=True,selected_case_ids=ids,eligible_recipes=rr,rejected_large_missing_face_recipes=large,distinct_carrier_types=128,prior_exact_exclusions=known,cases=[],scope='All central-singleton wedge recipes with two ordered codimension2..5 truncations. Search is required only on16 unresolved overlap types; range is merely a bounding interval. Only completed exact results exclude a type; numerical metrics are still required for a surviving full category.')
 partial=ROOT/f'results/{TAG}_partial.json';t=time.time()
 try:
  for i in ids:
   r=search(d['cases'][i],limit);out['cases'].append(r);print(i,r['verdict'],r['nodes'],'seconds',round(time.time()-t,2),flush=True)
   if not verify:partial.write_text(json.dumps(out,indent=2)+'\n')
 except KeyboardInterrupt:
  out['status']='INTERRUPTED_CHECKPOINT';partial.write_text(json.dumps(out,indent=2)+'\n');raise
 unresolved={r['case']for r in out['cases']if r['verdict']!='EXACT_COMPLETE_OBSTRUCTION'}
 out['status']='EXACT_COMPLETE_CENTRAL_CARRIER_OBSTRUCTION'if not unresolved else'BOUNDED_CENTRAL_CARRIER_CHECKPOINT'
 out['unresolved_types']=sorted(unresolved);out['verdict_counts']=dict(Counter(r['verdict']for r in out['cases']))
 out=json.loads(json.dumps(out))
 path=ROOT/f'results/{TAG}.json'
 if verify:assert out==json.loads(path.read_text())
 else:path.write_text(json.dumps(out,indent=2)+'\n')
 if partial.exists():partial.unlink()
 print(out['status'],out['verdict_counts'])
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--nodes',type=int,default=10000);ap.add_argument('--verify',action='store_true');a=ap.parse_args();print('PID',os.getpid(),flush=True);run(a.nodes,a.verify)
