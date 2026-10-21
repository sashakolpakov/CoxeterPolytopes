#!/usr/bin/env python3
"""Independent Fraction audit of the rank-three cut and cone certificates."""
import argparse
from itertools import combinations
import json
from pathlib import Path
from scripts.fresh_h4_blocks_independent_verify import (
    Z,O,P,parse,enc,add,sub,scale,mul,div,sign,integral,total,dot,mv,tr,mm,pair,matrix,eye,e8_host,
)
from scripts.fresh_multibasepoint_h9_verify import ROOT,read,permutation,valid_coxeter_value
from scripts.fresh_multibasepoint_h9_targeted_verify import escape

def compatible(r,s,gram):
    v=pair(r,gram,s)
    return sign(v)<=0 and valid_coxeter_value(div(scale(mul(v,v),4),mul(pair(r,gram,r),pair(s,gram,s))))

def verify_obstruction(cones,old,gram,u,z):
    obstruction=cones['slice_obstruction'];assert obstruction['case']==0
    assert set(mv(u,r) for r in old)==set(old)
    v=tuple(map(parse,cones['cases'][0]['hyperideal_witness']['vector']))
    assert pair(v,gram,v)==parse(obstruction['escape_norm']) and sign(pair(v,gram,v))>0
    project=lambda w:(w[0],div(pair(z,gram,w),scale(O,14)),w[9])
    for record,target in zip(obstruction['identities'],(v,mv(u,v)),strict=True):
        indices=record['old_wall_indices'];coefficients=tuple(map(parse,record['nonnegative_coefficients']))
        assert len(indices)==len(coefficients)==3 and all(sign(c)>=0 for c in coefficients)
        projected=project(target);assert projected==tuple(map(parse,record['projection']))
        recombined=tuple(total(mul(c,project(old[i])[j]) for i,c in zip(indices,coefficients)) for j in range(3))
        assert projected==recombined
    return True

def verify(search_path,cones_path):
    search=json.loads(Path(search_path).read_text());cones=json.loads(Path(cones_path).read_text())
    base=read(search['source']);assert cones['base_source']==search['source']
    alpha=parse(search['alpha']);assert alpha==scale(P,2)
    p=permutation([0,1,8,2,3,4,5,6,7,9]);gram=mm(p,mm(e8_host(alpha,True),tr(p)))
    assert gram==matrix(base['gram'])==matrix(cones['gram'])
    old=matrix(base['all_walls']);inside=tuple(map(parse,base['strict_interior']))
    assert sign(pair(inside,gram,inside))<0 and sign(inside[0])>0
    target=escape(search['targeted_escape'],old,gram)
    assert search['targeted_escape']==base['hyperideal_witness']
    direction=search['e8_direction'];assert direction==[7,14,20,16,12,8,4,10]
    # 3*7-20=1 proves the spatial rank-one module, and hence this ternary
    # coordinate module, is saturated in the ambient integral lattice.
    assert 3*direction[0]-direction[2]==1
    z=(Z,*(scale(O,k) for k in direction),Z);scalar=(Z,*(Z for _ in range(8)),O)
    assert pair(z,gram,z)==scale(O,14) and pair(scalar,gram,scalar)==O and pair(z,gram,scalar)==Z
    assert matrix(search['rank3_gram'])==((scale(P,-2),Z,Z),(Z,scale(O,14),Z),(Z,Z,O))
    u=matrix(cones['standard_block_swap_involution'])
    assert u==matrix(base['standard_block_swap_involution']) and mm(u,u)==eye(10) and mm(tr(u),mm(gram,u))==gram
    assert all(integral(x) for row in u for x in row)
    normals=tuple(mv(gram,r) for r in old);norms=tuple(dot(r,n) for r,n in zip(old,normals))
    for i,(r,n,q) in enumerate(zip(old,normals,norms)):
        assert q in (O,scale(O,2)) and all(map(integral,r))
        assert sign(dot(inside,n))<0 and all(integral(div(scale(x,2),q)) for x in n)
        for j in range(i):
            v=dot(r,normals[j]);assert sign(v)<=0
            assert valid_coxeter_value(div(scale(mul(v,v),4),mul(q,norms[j])))
    roots=[]
    for cut in search['cuts']:
        r=tuple(map(parse,cut['root']));q=parse(cut['norm']);x=parse(cut['x']);t=parse(cut['t']);s=parse(cut['u'])
        assert q in (O,scale(O,2)) and all(map(integral,r))
        assert r==(x,*(scale(t,k) for k in direction),s)
        assert sign(x)>0 and sign(t)>0 and sign(s)>=0
        assert sign(sub(scale(O,search['t_bound']),t))>=0 and sign(sub(scale(O,search['u_bound']),s))>=0
        assert sub(add(scale(mul(t,t),14),mul(s,s)),mul(alpha,mul(x,x)))==q==pair(r,gram,r)
        assert all(integral(div(scale(y,2),q)) for y in mv(gram,r))
        products=tuple(dot(r,n) for n in normals)
        assert products==tuple(map(parse,cut['wall_pairings'])) and all(sign(v)<=0 for v in products)
        assert pair(r,gram,target)==parse(cut['escape_pairing']) and sign(pair(r,gram,target))>0
        assert sign(pair(r,gram,inside))<0
        assert div(mul(alpha,mul(x,x)),q)==parse(cut['old_basepoint_height'])
        witness=sub(mul(s,target[1]),scale(mul(t,target[9]),7))
        assert witness==parse(cut['outside_previous_plane_witness']) and witness!=Z
        for value,oldq in zip(products,norms):assert valid_coxeter_value(div(scale(mul(value,value),4),mul(q,oldq)))
        roots.append(r)
    assert len(roots)==len(set(roots))==7
    orbits=tuple((r,mv(u,r)) for r in roots)
    assert tuple(tuple(tuple(map(parse,r)) for r in orbit) for orbit in cones['orbits'])==orbits
    assert all(r!=s and compatible(r,s,gram) for r,s in orbits)
    adjacency=tuple(tuple(i!=j and all(compatible(r,s,gram) for r in orbits[i] for s in orbits[j])
        for j in range(len(orbits))) for i in range(len(orbits)))
    assert adjacency==tuple(tuple(row) for row in cones['orbit_compatibility'])
    maximal=[]
    for size in range(1,len(orbits)+1):
        for selection in combinations(range(len(orbits)),size):
            if not all(adjacency[i][j] for i,j in combinations(selection,2)):continue
            if any(all(adjacency[k][i] for i in selection) for k in range(len(orbits)) if k not in selection):continue
            maximal.append(selection)
    assert maximal==[tuple(s) for s in cones['maximal_compatible_selections']]
    assert len(maximal)==len(cones['cases'])==4
    cases=[]
    for case,selection in zip(cones['cases'],maximal):
        assert case['selection']==list(selection)
        walls=tuple(dict.fromkeys((*old,*(r for i in selection for r in orbits[i]))))
        assert walls==matrix(case['all_walls']) and len(walls)==case['wall_count']==382
        assert tuple(map(parse,case['strict_interior']))==inside
        assert set(mv(u,r) for r in walls)==set(walls)
        for i in range(len(old),len(walls)):
            r=walls[i];q=pair(r,gram,r)
            assert q in (O,scale(O,2)) and all(map(integral,r))
            assert sign(pair(inside,gram,r))<0 and all(integral(div(scale(v,2),q)) for v in mv(gram,r))
            assert all(compatible(r,w,gram) for w in walls[:i])
        v=escape(case['hyperideal_witness'],walls,gram)
        assert case['status']=='EXACT_SPACELIKE_ESCAPE'
        cases.append({'selection':list(selection),'walls':len(walls),'escape_norm':enc(pair(v,gram,v))})
    verify_obstruction(cones,old,gram,u,z)
    result={'status':'EXACT_RANK3_CUTS_AND_FOUR_COXETER_CONES_PASS','cuts_outside_previous_plane':len(roots),
        'all_pairs_coxeter_verified':True,'exact_fixed_slice_obstruction':True,'maximal_selection_census_complete_for_seven_orbits':True,'cases':cases,
        'scope':'Independent certificate replay; does not independently repeat the full parameter-box enumeration.'}
    print(json.dumps(result,sort_keys=True));return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('search',type=Path);p.add_argument('cones',type=Path);p.add_argument('--output',type=Path)
    a=p.parse_args();result=verify(a.search,a.cones)
    if a.output:a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
