#!/usr/bin/env python3
"""Exact payment screen of minimal B1 rows and their saved first extensions.

This does not assert further attachment closure or geometric realizability
of a failed symbolic/dotted row. Every remaining spectator context is checked.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_global_discount_20261016 as G

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_global_b1_first_extensions_20261016'
INPUT='results/fresh_sqrt5_adjacent_discount_frontier_20261016.json'
RULE='results/fresh_sqrt5_global_discount_20261016.json'

def run():
    G.install();M=G.M;X=G.X
    d=json.loads((ROOT/INPUT).read_text());rule=json.loads((ROOT/RULE).read_text())
    assert rule['certificate']['kind']=='exact_feasible_restricted_rule'
    weights=list(map(F,rule['certificate']['primal']))
    assert weights[0]==0 and 0<=weights[1]<F(1,2)
    rows=[]
    for kind in ['records','outgoing_roots']:
        for r in d[kind]:
            if r['analysis']['baseline']!=1:continue
            row=r['row']
            if row['nodes']>23:continue
            rows.append((kind,r))
    values={};checks=0;summary=[];failures=[];overall=None
    for i,(kind,r)in enumerate(rows):
        row=r['row'];n=row['nodes'];edges=[(a,b,m)for(a,b),m in zip(combinations(range(n),2),row['labels'])]
        compiled=[X.compile_corner(row['white']+[a,b],edges,[a,b],[None])for a,b in product((0,2),(1,3))]
        lowest=None;worst=None;failed=0;count=0
        for ext in X.external_contexts(23-n):
            p,tables=X.add_bonus(compiled,ext)
            for table in tables:p=M.plus(p,table[None])
            key=M.key(p)
            if key not in values:values[key]=M.value(p,weights)
            payment=values[key];checks+=1;count+=1
            if lowest is None or payment<lowest:lowest=payment;worst=ext
            if payment<1:failed+=1
        assert lowest is not None
        identity=dict(family=kind,index=r.get('class_index',r.get('parent_class')),
                      contacts=r.get('contacts'),nodes=n,
                      inherited_status=r.get('inherited_status'),parameter_type=r.get('parameter_type'))
        summary.append(dict(**identity,contexts=count,failed_contexts=failed,minimum=str(lowest),worst_external_counts=list(worst)))
        if failed:failures.append(dict(**summary[-1],row=row,analysis=r['analysis']))
        if overall is None or lowest<overall:overall=lowest
        if i%100==0:print('B1',i,len(rows),'contexts',checks,'profiles',len(values),'failed_rows',len(failures),'minimum',str(overall),flush=True)
    failures.sort(key=lambda r:(F(r['minimum']),r['nodes'],r['index']))
    src=[f'scripts/{TAG}.py',f'scripts/{G.TAG}.py',INPUT,RULE,
         f'scripts/{X.TAG}.py','scripts/fresh_sqrt5_adjacent_discount_frontier_20261016.py']
    return dict(status='PASS_EXACT_GLOBAL_B1_MINIMAL_AND_FIRST_EXTENSION_PAYMENT_SCREEN',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in src},
                discount=str(weights[1]),rows=len(rows),family_counts=dict(Counter(k for k,r in rows)),
                common_spectator_cases=checks,distinct_payment_profiles=len(values),minimum_payment=str(overall),
                failed_rows=len(failures),failures=failures,summary=summary,
                scope='Exact rational payment evaluation for all saved baseline1 minimal-union rows and saved first-white extensions, with every current global spectator context. Includes overcomplete symbolic-dotted and rejected minimum rows as diagnostics; failures are not asserted realizable. No second-or-later extension closure, no global dimension bound.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['rows'],d['failed_rows'],d['minimum_payment'],flush=True)
