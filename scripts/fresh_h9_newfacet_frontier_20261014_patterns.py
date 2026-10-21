#!/usr/bin/env python3
"""Reusable induced-nerve Lanner/finite-diagram obstructions, at most9 nodes.

This is a necessary local relaxation. A complete local failure excludes all
literal occurrences; a survivor or node limit makes no negative claim.
"""
import argparse,json,os,time
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from fresh_h9_highbase_construction_20261014_cell import codes,finite
from fresh_h9_bugaenko_surgery_20261014_incidence import canonical
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_newfacet_frontier_20261014_patterns'
INC='results/fresh_h9_bugaenko_remaining_20261014_newfacet.json'
FRONT='results/fresh_h9_bugaenko_remaining_20261014_checkpoint.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def bits(m):return [i for i in range(m.bit_length())if m>>i&1]
def maximal_faces(n,missing):
 faces={v for v in range(1<<n)if all(v&m!=m for m in missing)}
 return sorted(v for v in faces if all((v|(1<<i))not in faces for i in range(n)if not v>>i&1))
def local_search(n,missing,limit):
 pairs=list(combinations(range(n),2));pi={p:i for i,p in enumerate(pairs)}
 faces=maximal_faces(n,missing)
 vv=[(v.bit_count(),tuple(pi[p]for p in combinations(bits(v),2)))for v in faces]
 tables=[(tuple(pi[p]for p in combinations(bits(m),2)),codes(m.bit_count()))for m in missing]
 nodes=0;dead=set();counts=Counter();candidate=None
 def good(labels):return all(finite(k,tuple(labels[p]for p in pp))for k,pp in vv)
 def rec(left,labels,known):
  nonlocal nodes,candidate
  nodes+=1
  if nodes>limit:return None
  key=(left,tuple(labels),tuple(sorted(known)))
  if key in dead:counts['cache_hits']+=1;return False
  if not good(labels):counts['nonspherical_face']+=1;return False
  pending=[i for i in range(len(tables))if left>>i&1]
  if pending:
   best=None
   for i in pending:
    pp,pool=tables[i];options=[c for c in pool if all(p not in known or labels[p]==m for p,m in zip(pp,c))]
    if not options:counts['empty_lanner_table']+=1;return False
    if best is None or len(options)<len(best[1]):best=(i,options)
   i,options=best;pp=tables[i][0]
   for code in options:
    child=labels[:]
    for p,m in zip(pp,code):child[p]=m
    result=rec(left^(1<<i),child,known|set(pp))
    if result is not False:return result
  else:
   unknown=set(range(len(pairs)))-known
   if not unknown:candidate=labels[:];return True
   p=min(unknown)
   for m in range(2,8):
    child=labels[:];child[p]=m
    result=rec(0,child,known|{p})
    if result is not False:return result
  dead.add(key);return False
 result=rec((1<<len(tables))-1,[2]*len(pairs),set())
 out=dict(verdict='LOCAL_LABEL_SURVIVOR'if result is True else'NODE_LIMIT_UNRESOLVED'if result is None else'EXACT_LOCAL_OBSTRUCTION',nodes=nodes,counters=dict(counts),maximal_faces=faces)
 if candidate is not None:out['label_vector']=candidate
 return out

def inventory(data,ids):
 patterns=defaultdict(dict);cache={}
 for case in ids:
  row=data['cases'][case];mm=row['canonical_minimal_nonfaces']
  unions={a|b for a,b in combinations(mm,2)if a.bit_count()>=4 and b.bit_count()>=4 and (a|b).bit_count()<=9}
  for union in sorted(unions):
   vertices=bits(union);missing=[m for m in mm if m&union==m]
   if min(map(int.bit_count,missing))==2:continue
   reduced=tuple(sorted(sum(1<<k for k,j in enumerate(vertices)if m>>j&1)for m in missing));n=len(vertices)
   if(n,reduced)not in cache:cache[n,reduced]=canonical(reduced,n)
   canon,order=cache[n,reduced]
   patterns[n,canon].setdefault(case,[vertices[j]for j in order])
 return patterns

def run(trials,limit):
 data=json.loads((ROOT/INC).read_text());front=json.loads((ROOT/FRONT).read_text())
 for d in(data,front):
  for p,h in d['sources'].items():assert sha(p)==h,p
 ids=front['changed_order']['open_case_ids'];assert len(ids)==3244
 bank=inventory(data,ids);todo=set(bank);open_ids=set(ids);records=[];covered={};start=time.monotonic()
 for trial in range(trials):
  useful=[k for k in todo if set(bank[k])&open_ids]
  if not useful:break
  key=min(useful,key=lambda k:(-len(set(bank[k])&open_ids),k[0],-len(k[1]),k[1]));todo.remove(key)
  n,missing=key;out=local_search(n,missing,limit)
  record=dict(pattern=len(records),order=n,minimal_nonfaces=missing,applications=[[i,p]for i,p in sorted(bank[key].items())],**out)
  if out['verdict']=='EXACT_LOCAL_OBSTRUCTION':
   for case,perm in record['applications']:
    # Verify every locally required spherical subset is an actual face.
    row=data['cases'][case];union=sum(1<<i for i in perm)
    lifted={sum(1<<perm[i]for i in bits(m))for m in missing}
    assert lifted=={m for m in row['canonical_minimal_nonfaces']if m&union==m}
    for face in out['maximal_faces']:
     lift=sum(1<<perm[i]for i in bits(face))
     assert any(v&lift==lift for v in row['canonical_vertices'])
    covered.setdefault(case,len(records));open_ids.discard(case)
  records.append(record)
  print(trial,'n',n,out['verdict'],'nodes',out['nodes'],'remaining',len(open_ids),'seconds',round(time.monotonic()-start,2),flush=True)
 files=[INC,FRONT,f'scripts/{TAG}.py','scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/independent_h9_14_record30_a1_labels.py','scripts/fresh_h9_bugaenko_surgery_20261014_incidence.py']
 return dict(status='EXACT_INDUCED_NERVE_OBSTRUCTIONS_WITH_OPEN_FRONTIER',sources={p:sha(p)for p in files},
  original_frontier_size=len(ids),pattern_inventory_size=len(bank),trial_limit=trials,node_limit=limit,
  verdict_counts=dict(Counter(r['verdict']for r in records)),patterns=records,
  excluded_cases=[[i,p]for i,p in sorted(covered.items())],remaining_case_ids=sorted(open_ids),
  scope='At most9-node induced nerve configurations with no dotted pairs. Complete Lanner tables, shared labels and every actual maximal spherical face are tested. Unknown meeting labels are relaxed to2 only for necessary Perron-Frobenius pruning; complete assignments cover finite labels2..6 and category7 for all finite m>=7. A fully exhausted local failure excludes every literal catalogue occurrence, independently of ambient metric/dotted distances. Surviving local categories and node limits are not exclusions or global constructions.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--trials',type=int,default=160);ap.add_argument('--nodes',type=int,default=10000);ap.add_argument('--verify',action='store_true');a=ap.parse_args()
 print('PID',os.getpid(),flush=True);out=json.loads(json.dumps(run(a.trials,a.nodes)))
 path=ROOT/f'results/{TAG}.json'
 if a.verify:assert out==json.loads(path.read_text())
 else:path.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],len(out['excluded_cases']),len(out['remaining_case_ids']),flush=True)
