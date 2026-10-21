#!/usr/bin/env python3
"""Compatible completions from free-ratio rank-four cuts and their images."""
from __future__ import annotations
import argparse
import json
import time
import numpy as np
from scripts.fresh_multibasepoint_h9_merge import (
    ROOT, G, read, matrix, inverse, mv, dot, enc, encm,
)
from scripts.fresh_rank3_h9_search import ri_sign
from scripts.fresh_rank3_h9_cones import escape

SEARCH='results/fresh_rank4_h9_free_search_20261013.json'
BASE='results/fresh_rank3_h9_cones_20261013.json'
OUTPUT=ROOT/'results/fresh_rank4_h9_completion_20261013.json'
DICTIONARY=ROOT/'results/fresh_rank4_h9_dictionary_20261013.json'


def pair_conflict(a,b,q):
    if ri_sign(a,b)>0: return True
    x,y=4*(a*a+b*b),4*(2*a*b+b*b)
    if ri_sign(x-4*q,y)>=0:return False
    return (x,y) not in {(0,0),(q,0),(2*q,0),(3*q,0),(q,q),(2*q,q)}


def dictionary():
    source=read(BASE); search=read(SEARCH);gram=matrix(source['gram'])
    old=matrix(source['cases'][0]['all_walls']);swap=matrix(source['standard_block_swap_involution'])
    assert {mv(swap,r) for r in old}==set(old)
    roots=[];provenance=[];lookup={}
    for i,row in enumerate(search['cuts']):
        r=tuple(map(G.from_json,row['root']))
        for which,v in [('direct',r),('swap',mv(swap,r))]:
            assert v[0].sign()>0 and all(x.is_integral for x in v)
            if v not in lookup:lookup[v]=len(roots);roots.append(v);provenance.append([])
            provenance[lookup[v]].append([i,which])
    normals=tuple(mv(gram,r) for r in roots)
    A=np.array([[int(x.a) for x in r] for r in roots],dtype=object)
    B=np.array([[int(x.b) for x in r] for r in roots],dtype=object)
    C=np.array([[int(x.a) for x in r] for r in normals],dtype=object)
    D=np.array([[int(x.b) for x in r] for r in normals],dtype=object)
    pa=A@C.T+B@D.T;pb=A@D.T+B@C.T+B@D.T
    norms=[int(pa[i,i]) for i in range(len(roots))]
    assert all(q in (1,2) and pb[i,i]==0 for i,q in enumerate(norms))
    edges=[];neighbors=[set() for _ in roots]
    for i in range(len(roots)):
        for j in range(i):
            if pair_conflict(int(pa[i,j]),int(pb[i,j]),norms[i]*norms[j]):
                edges.append([j,i]);neighbors[i].add(j);neighbors[j].add(i)
    old_normals=tuple(mv(gram,r) for r in old)
    e8inv=inverse(tuple(tuple(row[1:9]) for row in gram[1:9]))
    rho=tuple(sum(e8inv[k],G(0)) for k in range(8))+(G(1),)
    den=100
    while True:
        inside=(G(1),*(x/den for x in rho))
        if dot(inside,mv(gram,inside)).sign()<0 and all(dot(inside,n).sign()<0 for n in (*old_normals,*normals)):break
        den*=2
    data={'base_source':BASE,'base_case':0,'search_source':SEARCH,'gram':source['gram'],
          'roots':encm(roots),'root_sources':provenance,'norms':norms,'conflict_edges':edges,
          'strict_interior':enc(inside),'scope':'All saved335 cutters and block-swap images; old382 inequalities are retained.'}
    DICTIONARY.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    return data,gram,old,tuple(roots),neighbors


def load_dictionary():
    d=json.loads(DICTIONARY.read_text());base=read(d['base_source'])
    gram=matrix(d['gram']);old=matrix(base['cases'][d['base_case']]['all_walls']);roots=matrix(d['roots'])
    neighbors=[set() for _ in roots]
    for i,j in d['conflict_edges']:neighbors[i].add(j);neighbors[j].add(i)
    return d,gram,old,roots,neighbors


def produce(rebuild=False):
    start=time.monotonic()
    d,gram,old,roots,neighbors=dictionary() if rebuild or not DICTIONARY.exists() else load_dictionary()
    heights=[2*G(0,1)*r[0]*r[0]/q for r,q in zip(roots,d['norms'])]
    increasing=sorted(range(len(roots)),key=lambda i:heights[i]);rng=np.random.default_rng(20261013)
    orders=[('increasing_height',increasing),('decreasing_height',list(reversed(increasing))),
            ('minimum_degree',sorted(increasing,key=lambda i:len(neighbors[i]))),
            ('maximum_degree',sorted(increasing,key=lambda i:-len(neighbors[i])))]
    for k in range(8):orders.append((f'random_{k}',list(map(int,rng.permutation(len(roots))))))
    seen=set();cases=[]
    for name,order in orders:
        ids=[];allowed=set(range(len(roots)))
        for i in order:
            if i not in allowed:continue
            ids.append(i);allowed.discard(i);allowed-=neighbors[i]
        ids=tuple(sorted(ids))
        if ids in seen:continue
        seen.add(ids);walls=(*old,*(roots[i] for i in ids))
        witness=escape(walls,gram,20261013+len(cases))
        row={'selection_method':name,'selected_root_indices':list(ids),'wall_count':len(walls),
             'hyperideal_witness':witness,
             'status':'EXACT_SPACELIKE_ESCAPE' if witness else 'VOLUME_UNDECIDED_NO_LP_ESCAPE'}
        cases.append(row);print(json.dumps({'method':name,'added':len(ids),'walls':len(walls),'escape_norm':witness['norm'] if witness else None}),flush=True)
        out={'dictionary':str(DICTIONARY.relative_to(ROOT)),'strict_interior':d['strict_interior'],'cases':cases,
             'scope':'Greedy compatible subsets of the finite cutter dictionary; no full-volume claim without a separate proof.',
             'seconds':round(time.monotonic()-start,3)}
        OUTPUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'dictionary_roots':len(roots),'conflicts':len(d['conflict_edges']),'cases':len(cases),'seconds':out['seconds']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--rebuild',action='store_true');a=p.parse_args();produce(a.rebuild)
