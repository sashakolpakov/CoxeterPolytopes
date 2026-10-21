#!/usr/bin/env python3
"""A bounded constructive 16-facet multiwedge branch over a rational 9-gon.

These incidences are polytopal by repeated ordinary wedge construction.
The label join is exact and uses all classical Lanner4/5 placements; labels
outside those missing faces are relaxed to right angles in vertex tests.
No surviving label state alone is a hyperbolic metric realization.
"""
from pathlib import Path
from itertools import combinations, product
from collections import Counter
import argparse,hashlib,json
from independent_h9_14_record30_a1_labels import lanner_codes,finite_code
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/fresh_h9_nonproduct_incidence_20261014_multiwedge.json'

def canonical(w):
    return min(tuple(w[(s+j)%9]for j in range(9))for s in range(9)) == tuple(w) and tuple(w)<=min(tuple(w[(s-j)%9]for j in range(9))for s in range(9))

def weight_vectors():
    out=[]
    def compositions(total,length,prefix=()):
        if length==1:
            yield prefix+(total,);return
        for first in range(1,total-length+2):yield from compositions(total-first,length-1,prefix+(first,))
    all_compositions=list(compositions(16,9));assert len(all_compositions)==6435
    for w in all_compositions:
        if not canonical(w):continue
        if any(w[i]+w[j]>5 for i,j in combinations(range(9),2)if(j-i)%9 not in(1,8)):continue
        out.append(w)
    assert len(out)==51 and all(max(w)<=3 for w in out)
    return sorted(out,key=lambda w:(-w.count(3),w))

def geometry(w):
    blocks=[];n=0
    for k in w:blocks.append(tuple(range(n,n+k)));n+=k
    M=sorted(tuple(sorted(blocks[i]+blocks[j]))for i,j in combinations(range(9),2)if(j-i)%9 not in(1,8))
    V=set()
    for i in range(9):
        F={i,(i+1)%9};outside=[blocks[j]for j in range(9)if j not in F]
        for omit in product(*outside):V.add(tuple(j for j in range(n)if j not in omit))
    assert all(len(A)==9 for A in V)
    assert all(not set(A)<=set(v)for A in M for v in V)
    masks=[sum(1<<i for i in A)for A in M]
    expected={A for A in combinations(range(n),9)if not any((sum(1<<i for i in A)&m)==m for m in masks)}
    assert V==expected
    counts=Counter(B for W in V for B in combinations(W,8));assert set(counts.values())=={2}
    return blocks,M,sorted(V)

def shaobject(o):return hashlib.sha256(json.dumps(o,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def join(M,V,maxstates):
    supports=[A for A in M if len(A)in(4,5)];states=[()];remaining=set(range(len(supports)));used=set();stages=[]
    while states and remaining:
        k=max(remaining,key=lambda i:(len(set(supports[i])&used),len(supports[i]),-i));B=supports[k];used.update(B)
        pairs=list(combinations(B,2));options=[tuple(zip(pairs,c))for c in lanner_codes(len(B))]
        restrictions=sorted({tuple(i for i in v if i in used)for v in V},key=lambda x:(-len(x),x))
        # Keep only maximal restrictions; every principal submatrix of a PD
        # matrix is PD, so this is an exact implication reduction.
        maximal=[]
        for A in restrictions:
            if not any(set(A)<=set(C)for C in maximal):maximal.append(A)
        passing=set();compatible=0;examined=0;truncated=False
        for state in states:
            old=dict(state)
            for option in options:
                examined+=1
                if any(p in old and old[p]!=m for p,m in option):continue
                compatible+=1;new=old|dict(option)
                if all(finite_code(len(A),tuple(new.get(p,2)for p in combinations(A,2)))for A in maximal):passing.add(tuple(sorted(new.items())))
                if len(passing)>maxstates:truncated=True;break
            if truncated:break
        states=sorted(passing)
        stages.append({'missing_face':B,'option_count':len(options),'input_states':examined//len(options)if not truncated else None,'extensions_examined':examined,'compatible_extensions':compatible,'actual_vertex_restrictions':len(maximal),'passing_states':len(states),'state_hash':shaobject(states),'truncated':truncated})
        remaining.remove(k)
        print('  join',B,'states',len(states),'compatible',compatible,'truncated',truncated,flush=True)
        if truncated:break
    return {'stages':stages,'terminal_states':states if states and not remaining else [],'unprocessed_missing_faces':len(remaining),'state_limit_reached':bool(stages and stages[-1]['truncated']),
            'verdict':'EXACT_MANDATORY_LANNER_VERTEX_OBSTRUCTION'if not states else 'BOUNDED_LABEL_JOIN_SURVIVOR_OR_LIMIT'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=8);p.add_argument('--max-states',type=int,default=20000);p.add_argument('--verify',action='store_true');a=p.parse_args();allw=weight_vectors();rows=[]
    settings={'limit':a.limit,'max_states':a.max_states}
    saved=json.loads(OUT.read_text())if a.verify else None
    if saved is not None:
        assert saved['settings']==settings
    print('complete bounded weight list',len(allw),flush=True)
    for i,w in enumerate(allw[:a.limit]):
        blocks,M,V=geometry(w);print('case',i,w,'vertices',len(V),'missing',dict(Counter(map(len,M))),flush=True)
        r={'case':i,'weights':w,'blocks':blocks,'minimal_nonfaces':M,'actual_vertices':V,**join(M,V,a.max_states)};rows.append(r)
        files=[Path(__file__).resolve(),ROOT/'scripts/independent_h9_14_record30_a1_labels.py']
        data={'status':'EXACT_NONPRODUCT_MULTIWEDGE_LABEL_JOIN','settings':settings,'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in files},'positive_weight_compositions_checked':6435,'available_weight_vectors':allw,'cases':rows,
          'scope':'Selected sixteen-facet H9 multiwedges over a9gon; full Lanner4/5 placements plus exact actual-vertex restrictions. Zero states exclude the selected combinatorial type for all Coxeter labels, because omitted meeting labels are relaxed to2. Positive states require metric rank/signature and compactness verification.'}
        data=json.loads(json.dumps(data))
        if not a.verify:OUT.write_text(json.dumps(data,indent=2)+'\n')
    if a.verify:assert data==saved;print('PASS saved exact multiwedge replay')
if __name__=='__main__':main()
