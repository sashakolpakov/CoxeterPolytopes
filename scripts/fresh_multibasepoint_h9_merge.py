#!/usr/bin/env python3
"""Merge compatible E8/H4 H9 prefixes and their transported block swaps.

This adds exact walls from different basepoints without raising any search
cutoff. Every merged cone receives an exact surviving escape certificate.
"""
from __future__ import annotations
import json
from pathlib import Path
from coxeter_polytopes.exact.golden import Golden as G, ONE as O, ZERO as Z
from scripts.fresh_h4_blocks_h9_search_20261013 import inverse, mv

ROOT=Path(__file__).resolve().parents[1]
BUGAENKO_TO_STANDARD=(0,1,8,2,3,4,5,6,7,9)
BLOCK_SWAP=(0,5,6,7,8,1,2,3,4,9)

def matrix(rows):return tuple(tuple(G.from_json(x) for x in row) for row in rows)
def enc(v):return [x.to_json() for x in v]
def encm(a):return [enc(row) for row in a]
def dot(a,b):return sum((x*y for x,y in zip(a,b,strict=True)),Z)
def tr(a):return tuple(zip(*a))
def mm(a,b):return tuple(tuple(dot(row,col) for col in tr(b)) for row in a)
def permutation(p):return tuple(tuple(O if i==p[j] else Z for j in range(len(p))) for i in range(len(p)))
def read(name):return json.loads((ROOT/name).read_text())
def coxeter_pair_value(value):
    return value>=4 or value in {G(0),G(1),G(2),G(3),G(1,1),G(2,1)}

def run(name,old_path,old_escape_path=None):
    h4_path=f'results/fresh_h4_blocks_h9_{name}_h30_20261013.json'
    iso_path=f'results/fresh_h4_blocks_h9_{name}_isometry_20261013.json'
    old=read(old_path);h4=read(h4_path);iso=read(iso_path)
    p=permutation(BUGAENKO_TO_STANDARD);s=permutation(BLOCK_SWAP)
    t=mm(p,matrix(iso['transform']));g=mm(p,mm(matrix(iso['target_gram']),tr(p)))
    assert mm(tr(t),mm(g,t))==matrix(h4['gram'])
    if 'gram' in old:assert g==matrix(old['gram'])
    u=mm(t,mm(s,inverse(t)))
    assert mm(u,u)==permutation(tuple(range(10))) and mm(tr(u),mm(g,u))==g
    assert all(x.is_integral for row in (*t,*inverse(t),*u) for x in row)
    alpha=G.from_json(h4['alpha'])
    old_walls=matrix(old['all_walls']);source_walls=matrix(h4['all_walls'])
    translated=tuple(mv(t,r) for r in source_walls)
    swapped=tuple(mv(u,r) for r in old_walls)
    assert {mv(s,r) for r in source_walls}==set(source_walls)
    provenance={r:[{'source':'old_prefix','index':i}] for i,r in enumerate(old_walls)}
    direct=[];symmetric=[]
    for source,roots,added in (('h4_translated',translated,direct),('old_block_swap',swapped,symmetric)):
        for i,r in enumerate(roots):
            if r not in provenance:provenance[r]=[];added.append(r)
            provenance[r].append({'source':source,'index':i})
    height=lambda r:alpha*r[0]*r[0]/dot(r,mv(g,r))
    added=sorted((*direct,*symmetric),key=height)
    walls=(*old_walls,*added);normals=tuple(mv(g,r) for r in walls)
    assert {mv(u,r) for r in walls}==set(walls)
    assert all(dot(r,normals[j]).sign()<=0 for i,r in enumerate(walls) for j in range(i))
    inside=tuple(map(G.from_json,old['strict_interior']))
    alternate_inside=mv(t,tuple(map(G.from_json,h4['strict_interior'])))
    for point in (inside,alternate_inside,mv(u,inside)):
        assert dot(point,mv(g,point)).sign()<0 and point[0].sign()>0
        assert all(dot(point,n).sign()<0 for n in normals)
    for r,n in zip(walls,normals):
        q=dot(r,n)
        assert q.sign()>0 and q.sign(1)>0 and all(x.is_integral for x in r)
        assert all((2*x/q).is_integral for x in n)
    bound=G.from_json(old.get('height_bound',old.get('completed_height')))
    assert all((height(r)-bound).sign()>0 for r in added)
    old_escape=read(old_escape_path) if old_escape_path else old['hyperideal_witness']
    escape=tuple(map(G.from_json,old_escape['vector']))
    products=tuple(dot(escape,n) for n in normals);escape_norm=dot(escape,mv(g,escape))
    assert escape_norm.sign()>0 and all(x.sign()<=0 for x in products)
    assert escape_norm==G.from_json(old_escape['norm'])
    witness={**old_escape,'wall_pairings':enc(products),'wall_count':len(walls),
        'source':'Unchanged old-prefix escape; all added inequalities verified exactly.'}
    new_records=[{'root':enc(r),'norm':dot(r,mv(g,r)).to_json(),'old_basepoint_height':height(r).to_json(),
        'sources':provenance[r]} for r in added]
    payload={
        'field':'Q(sqrt(5))','dimension':9,'alpha':alpha.to_json(),'gram':encm(g),
        'old_source':old_path,'old_escape_source':old_escape_path,'h4_source':h4_path,'isometry_source':iso_path,
        'bugaenko_to_standard_permutation':list(BUGAENKO_TO_STANDARD),
        'h4_block_swap_permutation':list(BLOCK_SWAP),'h4_to_standard_transform':encm(t),
        'standard_block_swap_involution':encm(u),'reflection_alignment_word':[],
        'old_wall_count':len(old_walls),'h4_translated_wall_count':len(translated),
        'h4_direct_new_walls':len(direct),'block_swap_additional_walls':len(symmetric),
        'direct_merged_wall_count':len(old_walls)+len(direct),'wall_count':len(walls),
        'all_walls':encm(walls),'new_walls':new_records,
        'complete_old_height':bound.to_json(),'source_h4_complete_height':h4['completed_height'],
        'new_height_min':height(added[0]).to_json(),'new_height_max':height(added[-1]).to_json(),
        'strict_interior':enc(inside),'translated_h4_strict_interior':enc(alternate_inside),
        'swapped_old_strict_interior':enc(mv(u,inside)),
        'basepoints':[enc(tuple(O if i==0 else Z for i in range(10))),enc(tuple(t[i][0] for i in range(10))),enc(tuple(u[i][0] for i in range(10)))],
        'hyperideal_witness':witness,'escape_fixed_by_block_swap':mv(u,escape)==escape,
        'geometry':{'status':'EXACT_SPACELIKE_ESCAPE'},
        'scope':'Compatible acute wall union from three basepoints; complete through the old height bound only. No full-group nonreflectivity or new-polytope claim.'
    }
    output=ROOT/f'results/fresh_multibasepoint_h9_{name}_20261013.json'
    output.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:payload[k] for k in ('alpha','old_wall_count','h4_direct_new_walls','block_swap_additional_walls','wall_count','new_height_min','new_height_max','escape_fixed_by_block_swap')},sort_keys=True),flush=True)

if __name__=='__main__':
    run('2phi','results/golden_h9_odd_a1_full_norms_h300_20261013.json','results/golden_h9_odd_a1_full_norms_h300_20261013_escape.json')
    run('sqrt5','results/fresh_golden_neighbor_sqrt5_h100.json')
