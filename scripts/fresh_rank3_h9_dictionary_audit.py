#!/usr/bin/env python3
"""Independent integer-pair conflict graph and twelve-cone replay."""
import argparse
import json
from pathlib import Path
from scripts.fresh_h4_blocks_independent_verify import (
    Z,O,P,parse,enc,scale,sign,integral,mv,tr,mm,pair,matrix,eye,
)
from scripts.fresh_multibasepoint_h9_verify import ROOT,read
from scripts.fresh_multibasepoint_h9_targeted_verify import escape
from scripts.fresh_rank3_h9_free_audit import ip,isum,times,kscale
from scripts.fresh_rank3_h9_newplane_verify import integer_sign

def conflict(product,den):
    if integer_sign(*product)>0:return True
    v=kscale(times(product,product),4)
    if integer_sign(v[0]-4*den,v[1])>=0:return False
    return v not in tuple(kscale(a,den) for a in ((0,0),(1,0),(2,0),(3,0),(1,1),(2,1)))

def verify():
    d=read('results/fresh_rank4_h9_dictionary_20261013.json');comp=read('results/fresh_rank4_h9_completion_20261013.json')
    b=read(d['base_source']);old=matrix(b['cases'][0]['all_walls']);gram=matrix(d['gram']);U=matrix(b['standard_block_swap_involution'])
    assert gram==matrix(b['gram']) and mm(U,U)==eye(10) and mm(tr(U),mm(gram,U))==gram
    assert all(integral(x) for row in U for x in row) and set(mv(U,r) for r in old)==set(old)
    roots=matrix(d['roots']);assert len(roots)==len(set(roots))==670
    source=read(d['search_source']);direct=matrix([r['root'] for r in source['cuts']])
    expect=[];provenance=[];lookup={}
    for i,r in enumerate(direct):
        for label,t in (('direct',r),('swap',mv(U,r))):
            if t not in lookup:lookup[t]=len(expect);expect.append(t);provenance.append([])
            provenance[lookup[t]].append([i,label])
    assert tuple(expect)==roots and provenance==d['root_sources']
    rint=tuple(tuple(map(ip,r)) for r in roots);cov=tuple(tuple(map(ip,mv(gram,r))) for r in roots)
    norms=[]
    for r,c in zip(rint,cov):
        q=isum(times(x,y) for x,y in zip(r,c));assert q in ((1,0),(2,0));norms.append(q[0])
    assert norms==d['norms']
    edges=[];neighbors=[set() for _ in roots]
    for i,r in enumerate(rint):
        for j in range(i):
            product=isum(times(x,y) for x,y in zip(r,cov[j]))
            if conflict(product,norms[i]*norms[j]):edges.append([j,i]);neighbors[i].add(j);neighbors[j].add(i)
    assert edges==d['conflict_edges'] and len(edges)==340
    inside=tuple(map(parse,d['strict_interior']));assert sign(inside[0])>0 and sign(pair(inside,gram,inside))<0
    assert tuple(map(parse,comp['strict_interior']))==inside
    assert all(sign(pair(inside,gram,r))<0 for r in (*old,*roots))
    reports=[];seen=set()
    for c in comp['cases']:
        ids=tuple(c['selected_root_indices']);chosen=set(ids)
        assert len(chosen)==len(ids) and all(0<=i<len(roots) for i in ids)
        assert all(not (neighbors[i]&chosen) for i in ids)
        assert all(neighbors[i]&chosen for i in range(len(roots)) if i not in chosen)
        assert ids not in seen;seen.add(ids)
        walls=(*old,*(roots[i] for i in ids));assert len(walls)==c['wall_count']
        v=escape(c['hyperideal_witness'],walls,gram)
        reports.append({'selection':c['selection_method'],'walls':len(walls),'escape_norm':enc(pair(v,gram,v))})
    assert len(reports)==12
    return {'status':'INDEPENDENT_670_NODE_CONFLICT_GRAPH_AND_12_CONES_PASS',
        'dictionary_roots':670,'root_pair_checks':670*669//2,'conflict_edges':340,
        'selections_maximal_in_dictionary':True,'cones':reports,
        'scope':'Exact finite-dictionary graph and all12 saved maximal compatible selections. Original-wall compatibility of335 direct roots was independently checked separately; their670 images follow by the verified integral isometry.'}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);a=p.parse_args();result=verify();print(json.dumps(result,sort_keys=True))
    if a.output:a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
