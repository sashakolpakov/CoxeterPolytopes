#!/usr/bin/env python3
"""Separate stdlib audit of every saved H9 vertex list and base dimension.

No producer or Coxeter-table implementation is imported. This checks the
incidence inputs to the exact label gates, not the label joins themselves.
"""
from pathlib import Path
from itertools import combinations
from collections import Counter
from functools import lru_cache
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_nonpolygon_construction_20261014'
FILES=[f'results/{TAG}{s}.json'for s in['','_edgecut','_fourbase']]
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def mask(s):return sum(1<<i for i in s)
@lru_cache(None)
def allnine(n):return[mask(s)for s in combinations(range(n),9)]
def run():
    rows=[];counts=Counter()
    for source in FILES:
        z=json.loads((ROOT/source).read_text())
        for p,h in z['sources'].items():assert sha(p)==h,p
        for r in z['cases']:
            n=sum(r['weights']);M=list(map(mask,r['minimal_nonfaces']));stored={mask(v)for v in r['actual_vertices']}
            assert len(stored)==len(r['actual_vertices'])
            # This is a direct census of all C(n,9) sets, separate from
            # the producers' wedge-vertex construction.
            expected={v for v in allnine(n)if not any(v&m==m for m in M)}
            assert stored==expected
            assert all(not any(v&m==m for v in stored)and all(any(v&(m^(1<<i))==m^(1<<i)for v in stored)for i in range(n)if m>>i&1)for m in M)
            assert set(Counter(mask(t)for v in r['actual_vertices']for t in combinations(v,8)).values())=={2}
            twins={}
            for i in range(n):twins.setdefault(tuple(j for j,m in enumerate(M)if m>>i&1),[]).append(i)
            assert sorted(map(tuple,twins.values()))==sorted(map(tuple,r['blocks']))
            base_dimension=9-n+len(twins);assert base_dimension==(4 if source.endswith('_fourbase.json')else 3)
            adj=[set()for _ in range(n)]
            for m in r['minimal_nonfaces']:
                for i in m:adj[i].update(set(m)-{i})
            seen={0};todo=[0]
            for i in todo:
                for j in adj[i]-seen:seen.add(j);todo.append(j)
            assert len(seen)==n
            assert r['verdict']=='EXACT_MANDATORY_LANNER_VERTEX_OBSTRUCTION'and not r['state_limit_reached']
            assert r['stages'][-1]['passing_states']==0 and not any(s['truncated']for s in r['stages'])
            rows.append({'source':source,'case':r['case'],'facets':n,'actual_vertices':len(stored),
                'all_nine_sets_tested':len(allnine(n)),'recovered_base_dimension':base_dimension,
                'minimal_nonface_hypergraph_connected':True})
            counts[source]+=1
    assert list(counts.values())==[82,30,3]
    files=FILES+[f'scripts/{TAG}_incidence_audit.py']
    return {'status':'PASS_SEPARATE_ALL_115_VERTEX_LIST_AND_NONPOLYGON_AUDIT','sources':{p:sha(p)for p in files},
        'case_counts':dict(counts),'nine_subsets_checked':sum(r['all_nine_sets_tested']for r in rows),'cases':rows,
        'scope':'Independent full nine-subset incidence check, minimal-nonface and ridge checks, recovered base dimension and nonproduct criterion. Exact label joins are separately replayed by the primary checkers; convex realizations use the primary explicit base constructions and ordinary wedge theorem.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();z=run();out=ROOT/f'results/{TAG}_incidence_audit.json'
    if a.verify:assert z==json.loads(out.read_text())
    else:out.write_text(json.dumps(z,indent=2)+'\n')
    print(z['status']);print(z['nine_subsets_checked'])
