#!/usr/bin/env python3
"""Short saved-record audit and disjoint final family ledger.

This checks all literal case hashes and independent-replay records. It does
not perform another expensive DFS; that replay is separately archived.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_fifteen_completion_20261014'
BASE='results/fresh_h9_fifteen_facet_20261014_checkpoint.json'
PRIMARY=f'results/{TAG}.json';AUDIT=f'results/{TAG}_independent.json'
HASHES={}
def fingerprint(p):
 s=(ROOT/p).stat();return(s.st_size,s.st_mtime_ns,s.st_ino)
def sha(p):
 before=fingerprint(p)
 if p in HASHES:
  mark,digest=HASHES[p];assert before==mark,('changed during verification',p);return digest
 digest=sha256((ROOT/p).read_bytes()).hexdigest();assert fingerprint(p)==before
 HASHES[p]=(before,digest);return digest
def load(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d.get('sources',{}).items():assert sha(f)==h,f
 return d
def run():
 base=load(BASE);primary=load(PRIMARY);audit=load(AUDIT)
 assert not primary['status'].startswith(('RUNNING','STOPPED_DISPATCH'))
 assert audit['status']=='PASS_SEPARATE_ARITHMETIC_ALL_SAVED_COMPLETION_CASES'
 records={r['case']:r for r in primary['records']};checks={r['case']:r for r in audit['records']}
 assert len(records)==len(primary['records'])and set(records)==set(checks)
 assert set(records)<=set(base['open_case_ids'])and len(base['open_case_ids'])==1334
 inherited_closed=sum(base[k]for k in('inherited_local_or_block_exclusions','new_local_pattern_exclusions','complete_global_label_exclusions'))
 assert inherited_closed==44345 and base['distinct_nerves']==45679
 assert base['remaining_without_complete_global_search']==1334 and not(set(base['open_case_ids'])&set(base['excluded_by_global_search']))
 counts=Counter();gate_counts=Counter();states=0;minor_codes=0;closed=[];limited=[];survivors=[]
 for i,ref in records.items():
  path=ref['path'];assert sha(path)==ref['sha256'];d=load(path);r=d['result'];check=checks[i]
  assert d['sources']==primary['sources']and d['node_limit']==primary['node_limit_per_case']
  assert r['case']==i and r['verdict']==ref['verdict']and r['nodes']==ref['nodes']
  assert check['input']==path and check['input_sha256']==ref['sha256']
  assert check['verdict']==r['verdict']and check['states']==r['nodes']
  ap=f'results/{TAG}_audit_cases/case_{i:05d}.json';a=load(ap);assert a['result']==check
  assert a['sources']=={f:h for f,h in audit['sources'].items()if f!=PRIMARY}
  counts[r['verdict']]+=1;states+=r['nodes']
  gate_counts.update(r['counters']);minor_codes+=r['distinct_modular_minor_codes']
  if r['verdict']=='EXACT_COMPLETE_OBSTRUCTION':closed.append(i);assert r['nodes']<=d['node_limit']
  elif r['verdict']=='NODE_LIMIT_UNRESOLVED':limited.append(i)
  else:assert r['verdict']=='FULL_CATEGORY_METRIC_REQUIRED';survivors.append(i)
 assert dict(counts)==primary['verdict_counts']and states==primary['total_states']==audit['states']
 unsearched=sorted(set(base['open_case_ids'])-set(records));assert unsearched==primary['remaining_unsearched']
 unresolved=sorted(unsearched+limited+survivors);total_closed=inherited_closed+len(closed)
 assert total_closed+len(unresolved)==45679
 files=[BASE,PRIMARY,AUDIT,f'scripts/{TAG}_certificate.py']
 sources={p:sha(p)for p in files}
 for p,(mark,digest)in HASHES.items():assert fingerprint(p)==mark,('changed during verification',p)
 return dict(status='PASS_ALL_45679_ONE_CUT_NERVES_EXCLUDED'if not unresolved else'PASS_COMPLETED_CASES_WITH_EXPLICIT_OPEN_FRONTIER',
  sources=sources,starting_open_cases=1334,completed_cases=len(records),new_exact_exclusions=len(closed),
  states=states,largest_case_states=max((r['nodes']for r in primary['records']),default=0),
  total_completed_global_searches=base['complete_global_label_exclusions']+len(closed),
  total_completed_global_states=base['global_states']+states,
  search_counters=dict(gate_counts),distinct_modular_minor_codes_summed=minor_codes,
  verdict_counts=dict(counts),total_distinct_nerves=45679,total_exact_exclusions=total_closed,
  unsearched_cases=unsearched,node_limited_cases=sorted(limited),full_category_survivors=sorted(survivors),remaining_open_case_ids=unresolved,
  literal_case_hashes_checked=len(records),independent_case_records_checked=len(checks),
  scope='Every completed primary case is bound by its literal file hash, and its separate-arithmetic replay record is matched with the final independently bound aggregate and source manifest. The earlier45,679-nerve incidence coverage and44,345 disjoint exclusions are inherited from the frozen ledger. No node-limited, unsearched or full-category-surviving case is counted as excluded. Even an empty frontier excludes this entire one-truncation family only, not all fifteen-facet polytopes.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}_certificate.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['new_exact_exclusions'],len(d['remaining_open_case_ids']),flush=True)
