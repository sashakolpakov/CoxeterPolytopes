#!/usr/bin/env python3
"""Necessary all-finite-label joins for the independently archived nerves."""
import argparse,contextlib,io,json,time,os
from collections import Counter
from hashlib import sha256
from pathlib import Path
from fresh_h9_highbase_construction_20261014_cell import join,codes

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_bugaenko_surgery_20261014'
SOURCE=f'results/{TAG}_incidence.json'
def tup(m):return tuple(i for i in range(m.bit_length())if m>>i&1)
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run(start,end,limit):
 t=time.time();src=json.loads((ROOT/SOURCE).read_text())
 for p,h in src['sources'].items():assert sha(p)==h,p
 end=min(end,len(src['cases']));rows=[];counts=Counter()
 print('PID',os.getpid(),'range',start,end,'state limit',limit,flush=True)
 for case in src['cases'][start:end]:
  M=list(map(tup,case['canonical_minimal_nonfaces']));V=list(map(tup,case['canonical_vertices']))
  with contextlib.redirect_stdout(io.StringIO()):out=join(M,V,limit)
  counts[out['verdict']]+=1;rows.append(dict(case=case['case'],vertices=len(V),missing_orders=dict(Counter(map(len,M))),**out))
  if len(rows)%25==0 or end-start<=25:print('done',len(rows),'last',case['case'],dict(counts),'seconds',round(time.time()-t,2),flush=True)
 paths=[SOURCE,f'scripts/{TAG}_labels.py','scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/independent_h9_14_record30_a1_labels.py']
 return dict(status='EXACT_BUGAENKO_SURGERY_NECESSARY_LABEL_CHECKPOINT',sources={p:sha(p)for p in paths},range=[start,end],total_incidence_types=len(src['cases']),all_incidence_types_processed=start==0 and end==len(src['cases']),state_limit=limit,verdict_counts=dict(counts),cases=rows,scope='Complete finite label categories on every joined minimal nonface, with necessary spherical tests at every actual vertex. Zero before a state limit excludes that incidence. Capped or surviving joins are unresolved, not counterexamples. Labels>=7 share a safe category in rank3/finite tests; infinity is separate.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--start',type=int,default=0);p.add_argument('--end',type=int,default=1000000);p.add_argument('--max-states',type=int,default=5000);p.add_argument('--verify',action='store_true');a=p.parse_args();out=run(a.start,a.end,a.max_states);out=json.loads(json.dumps(out));suffix=''if out['all_incidence_types_processed']else f'_{a.start}_{out["range"][1]}';path=ROOT/f'results/{TAG}_labels{suffix}.json'
 if a.verify:assert out==json.loads(path.read_text())
 else:path.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['verdict_counts'],flush=True)
