#!/usr/bin/env python3
"""Fraction-only replay of the H9 multi-basepoint merge certificates."""
import argparse
import json
from pathlib import Path
from scripts.fresh_h4_blocks_independent_verify import (
    Z,O,P,parse,enc,add,sub,scale,mul,div,sign,integral,dot,mv,tr,mm,pair,
    matrix,eye,inverse,rref,e8_host,host,
)

ROOT=Path(__file__).resolve().parents[1]
def permutation(p):return tuple(tuple(O if i==p[j] else Z for j in range(len(p))) for i in range(len(p)))
def read(path):return json.loads((ROOT/path).read_text())
def valid_coxeter_value(value):
    return sign(sub(value,scale(O,4)))>=0 or value in {Z,O,scale(O,2),scale(O,3),add(O,P),add(scale(O,2),P)}

def verify(path):
    d=json.loads(Path(path).read_text());alpha=parse(d['alpha']);g=matrix(d['gram'])
    assert d['dimension']==9 and alpha in (scale(P,2),sub(scale(P,2),O))
    assert d['bugaenko_to_standard_permutation']==[0,1,8,2,3,4,5,6,7,9]
    assert d['h4_block_swap_permutation']==[0,5,6,7,8,1,2,3,4,9]
    p=permutation(d['bugaenko_to_standard_permutation']);s=permutation(d['h4_block_swap_permutation'])
    assert g==mm(p,mm(e8_host(alpha,True),tr(p)))
    old=read(d['old_source']);h4=read(d['h4_source']);iso=read(d['isometry_source'])
    t=matrix(d['h4_to_standard_transform']);u=matrix(d['standard_block_swap_involution'])
    assert t==mm(p,matrix(iso['transform']))
    assert mm(tr(t),mm(g,t))==host(alpha,2,True)==matrix(h4['gram'])
    assert all(integral(x) for row in (*t,*inverse(t),*u) for x in row)
    assert u==mm(t,mm(s,inverse(t))) and mm(u,u)==eye(10) and mm(tr(u),mm(g,u))==g
    old_walls=matrix(old['all_walls']);h4_walls=matrix(h4['all_walls'])
    translated=tuple(mv(t,r) for r in h4_walls);swapped=tuple(mv(u,r) for r in old_walls)
    assert {mv(s,r) for r in h4_walls}==set(h4_walls)
    direct=set(translated)-set(old_walls);extra=set(swapped)-set(old_walls)-set(translated)
    walls=matrix(d['all_walls']);normals=tuple(mv(g,r) for r in walls)
    assert len(walls)==len(set(walls))==d['wall_count']
    assert set(walls)==set(old_walls)|set(translated)|set(swapped)
    assert walls[:len(old_walls)]==old_walls
    assert d['old_wall_count']==len(old_walls) and d['h4_direct_new_walls']==len(direct)
    assert d['block_swap_additional_walls']==len(extra)
    assert d['direct_merged_wall_count']==len(old_walls)+len(direct)
    assert set(mv(u,r) for r in walls)==set(walls)
    norms=tuple(dot(r,n) for r,n in zip(walls,normals))
    for i,(r,n) in enumerate(zip(walls,normals)):
        q=dot(r,n);assert sign(q)>0 and all(map(integral,r))
        assert all(integral(div(scale(x,2),q)) for x in n)
        for j in range(i):
            v=dot(r,normals[j]);assert sign(v)<=0
            assert valid_coxeter_value(div(scale(mul(v,v),4),mul(q,norms[j])))
    inside=tuple(map(parse,d['strict_interior']))
    h4_inside=tuple(map(parse,d['translated_h4_strict_interior']))
    swap_inside=tuple(map(parse,d['swapped_old_strict_interior']))
    assert inside==tuple(map(parse,old['strict_interior']))
    assert h4_inside==mv(t,tuple(map(parse,h4['strict_interior'])))
    assert swap_inside==mv(u,inside)
    for v in (inside,h4_inside,swap_inside):
        assert sign(pair(v,g,v))<0 and sign(v[0])>0
        assert all(sign(dot(v,n))<0 for n in normals)
    assert d['reflection_alignment_word']==[]
    bound=parse(d['complete_old_height']);previous=bound
    assert bound==parse(old.get('height_bound',old.get('completed_height')))
    records=d['new_walls'];assert len(records)==len(direct)+len(extra)
    assert [tuple(map(parse,z['root'])) for z in records]==list(walls[len(old_walls):])
    root_sources={'old_prefix':old_walls,'h4_translated':translated,'old_block_swap':swapped}
    for record in records:
        r=tuple(map(parse,record['root']));q=pair(r,g,r);h=div(mul(alpha,mul(r[0],r[0])),q)
        assert q==parse(record['norm']) and h==parse(record['old_basepoint_height'])
        assert sign(sub(h,bound))>0 and sign(sub(h,previous))>=0;previous=h
        assert record['sources'] and all(root_sources[z['source']][z['index']]==r for z in record['sources'])
    assert parse(d['new_height_min'])==parse(records[0]['old_basepoint_height'])
    assert parse(d['new_height_max'])==parse(records[-1]['old_basepoint_height'])
    witness=d['hyperideal_witness'];v=tuple(map(parse,witness['vector']));q=pair(v,g,v)
    assert q==parse(witness['norm']) and sign(q)>0 and sign(v[0])>0
    values=tuple(dot(v,n) for n in normals)
    assert values==tuple(map(parse,witness['wall_pairings'])) and all(sign(x)<=0 for x in values)
    active=witness['active_facets'];assert len(set(active))==len(active)==9
    assert all(values[i]==Z for i in active) and len(rref([normals[i] for i in active])[1])==9
    old_witness=read(d['old_escape_source']) if d['old_escape_source'] else old['hyperideal_witness']
    assert v==tuple(map(parse,old_witness['vector']))
    assert d['escape_fixed_by_block_swap']==(mv(u,v)==v)
    result={'artifact':str(path),'walls':len(walls),'direct_new':len(direct),'symmetry_new':len(extra),
        'status':'EXACT_COMPATIBLE_MERGE_WITH_SURVIVING_ESCAPE_PASS','escape_norm':enc(q)}
    print(json.dumps(result,sort_keys=True));return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('artifacts',type=Path,nargs='+');p.add_argument('--output',type=Path)
    a=p.parse_args();results=[verify(path) for path in a.artifacts]
    if a.output:a.output.write_text(json.dumps(results,indent=2,sort_keys=True)+'\n')
