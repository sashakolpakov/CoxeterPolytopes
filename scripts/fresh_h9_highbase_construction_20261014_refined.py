#!/usr/bin/env python3
"""Coupled constraint order for the one-Lanner-four 120-cell recipe.

This retains every state; it does not enlarge the prior 20,000-state cap.
"""
import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_h9_highbase_construction_20261014_cell as C

ROOT=C.ROOT;TAG='fresh_h9_highbase_construction_20261014_refined'
SOURCE=f'results/{C.TAG}.json'
OUT=ROOT/f'results/{TAG}.json'


def run():
    source=json.loads((ROOT/SOURCE).read_text())
    for p,h in source['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    row=source['cases'][0];assert row['mandatory_lanner4']==1
    M=list(map(tuple,row['minimal_nonfaces']));V=row['actual_vertices']
    missing_edges={m for m in M if len(m)==2}
    supports=[m for m in M if len(m)>=3]
    sets=list(map(frozenset,supports));remaining=set(range(len(supports)))
    used=set();states=[()];stages=[]
    while states and remaining:
        def priority(k):
            b=sets[k];u=used|b
            closures=sum(sets[j]<=u for j in remaining)
            new_contacts=sum(tuple(sorted((i,j))) not in missing_edges
                             for i in b-used for j in used-b)
            return (len(b&used),len(b),closures,new_contacts,-k)
        k=max(remaining,key=priority);B=supports[k];used.update(B)
        options=[tuple(zip(combinations(B,2),code)) for code in C.codes(len(B))]
        restrictions=sorted({tuple(i for i in v if i in used) for v in V},key=lambda a:(-len(a),a))
        maximal=[]
        for a in restrictions:
            if not any(set(a)<=set(c) for c in maximal):maximal.append(a)
        passing=set();examined=compatible=0;truncated=False
        for state in states:
            old=dict(state)
            for opt in options:
                examined+=1
                if any(p in old and old[p]!=m for p,m in opt):continue
                compatible+=1;new=old|dict(opt)
                if all(C.finite(len(a),tuple(new.get(p,2) for p in combinations(a,2))) for a in maximal):passing.add(tuple(sorted(new.items())))
                if len(passing)>20000:truncated=True;break
            if truncated:break
        stages.append(dict(missing_face=B,input_states=len(states),option_count=len(options),
                           extensions_examined=examined,compatible_extensions=compatible,
                           maximal_actual_vertex_restrictions=maximal,passing_states=len(passing),
                           state_hash=C.digest(sorted(passing)),truncated=truncated))
        states=sorted(passing);remaining.remove(k)
        print('join',B,'states',len(states),'used',len(used),'truncated',truncated,flush=True)
        if truncated:break
    paths=[SOURCE,f'scripts/{TAG}.py',f'scripts/{C.TAG}.py','scripts/independent_h9_14_record30_a1_labels.py']
    return json.loads(json.dumps(dict(status='EXACT_REORDERED_120CELL_NECESSARY_JOIN',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        source_case=0,max_states=20000,stages=stages,unprocessed_missing_faces=len(remaining),
        state_limit_reached=bool(stages[-1]['truncated']),
        terminal_states=states if states and not remaining else [],
        verdict='EXACT_LANNER_VERTEX_OBSTRUCTION' if not states else 'LABEL_SURVIVOR_OR_EXPLICIT_LIMIT',
        scope='The same one-Lanner-four 125-facet recipe, with coupled triangle constraints prioritized. '
              'Every state retained until the explicit same 20000 limit; zero without truncation is an exact necessary obstruction.')))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run()
    if a.verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['verdict'])
