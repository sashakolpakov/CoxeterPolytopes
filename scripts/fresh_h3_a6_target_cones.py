#!/usr/bin/env python3
"""Cross-angle checked greedy H3+A6 slice completions and exact escapes."""
import argparse
import json
import random
from pathlib import Path
import time
from scripts.fresh_h3_a6_target_search import *
from scripts.fresh_rank3_h9_cones import escape

def plus(*xs):return (sum(x[0] for x in xs),sum(x[1] for x in xs))
def scale(x,k):return (x[0]*k,x[1]*k)
def ints(x):return int(x.a),int(x.b)
def covector(v,factor):
    x,s,t,u=v
    return scale(mul((0,factor),x),-1),plus(mul((6,6),s),mul((2,4),t)),plus(mul((2,4),s),mul((2,2),t)),scale(u,6)
def prod(v,w):return plus(*(mul(x,y) for x,y in zip(v,w)))
def allowed(p,nn):
    if ri_sign(*p)>0:return False
    z=scale(mul(p,p),4)
    return ri_sign(z[0]-4*nn,z[1])>=0 or z in {(nn*a,nn*b) for a,b in ((0,0),(1,0),(2,0),(3,0),(1,1),(2,1))}

def run(name):
    start=time.monotonic();source='results/fresh_h3_a6_target_'+name+'_neighbor_20261013.json';d=read(source)
    base=read(d['source']);gram=matrix(base['gram']);old=matrix(base['all_walls']);inside=tuple(map(G.from_json,base['strict_interior']))
    rows=d['cuts'];roots=tuple(tuple(map(G.from_json,r['root'])) for r in rows)
    params=tuple(tuple(ints(G.from_json(c)) for c in r['parameters']) for r in rows)
    cov=tuple(covector(r,1 if name=='phi' else 2) for r in params);qs=tuple(r['norm'] for r in rows)
    conflicts=[set() for _ in roots];edges=[]
    for i,v in enumerate(params):
        for j in range(i):
            if not allowed(prod(v,cov[j]),qs[i]*qs[j]):conflicts[i].add(j);conflicts[j].add(i);edges.append((j,i))
    heights=tuple(G.from_json(r['height']) for r in rows)
    orders=[sorted(range(len(roots)),key=lambda i:heights[i]),sorted(range(len(roots)),key=lambda i:(len(conflicts[i]),heights[i]))]
    rng=random.Random(20261013)
    for _ in range(3):
        order=list(range(len(roots)));rng.shuffle(order);orders.append(order)
    cases=[];seen=set()
    for order in orders:
        selected=[];blocked=set()
        for i in order:
            if i not in blocked:selected.append(i);blocked.update(conflicts[i])
        ids=tuple(sorted(selected))
        if ids in seen:continue
        seen.add(ids);assert all(any(j in conflicts[i] for j in ids) for i in range(len(roots)) if i not in ids)
        walls=(*old,*(roots[i] for i in ids))
        assert all(dot(inside,mv(gram,r)).sign()<0 for r in walls)
        witness=escape(walls,gram,20261013+len(cases))
        cases.append({'selection':ids,'wall_count':len(walls),'hyperideal_witness':witness,
            'status':'EXACT_SPACELIKE_ESCAPE' if witness else 'VOLUME_UNDECIDED'})
        print(json.dumps({'alpha':name,'case':len(cases)-1,'walls':len(walls),'escape_norm':witness['norm'] if witness else None}),flush=True)
    # This inequality intersection need not be Coxeter; a common escape
    # certifies every compatible subset simultaneously.
    dictionary=(*old,*roots);common=escape(dictionary,gram,20261020)
    out={'source':source,'source_walls':len(old),'dictionary_size':len(roots),'conflict_edges':edges,'cases':cases,
        'strict_interior':base['strict_interior'],'common_dictionary_escape':common,'seconds':round(time.monotonic()-start,3),
        'scope':'Five greedy orders, deduplicated, give maximal compatible subsets of the full finite-box dictionary. A common dictionary escape, when present, covers every subset, but not other roots in the slice.'}
    output=ROOT/('results/fresh_h3_a6_target_'+name+'_cones_20261013.json');output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'alpha':name,'conflicts':len(edges),'cases':len(cases),'common_escape':common['norm'] if common else None,'seconds':out['seconds']}),flush=True)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--alpha',choices=('phi','2phi'),required=True);a=p.parse_args();run(a.alpha)
