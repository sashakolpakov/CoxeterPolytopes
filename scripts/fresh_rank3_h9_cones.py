#!/usr/bin/env python3
"""Construct every maximal compatible subset of the seven new cut orbits."""
from __future__ import annotations
from itertools import combinations
import json
import time
import numpy as np
from scipy.optimize import linprog
from coxeter_polytopes.verify.chamber import verify_hyperideal_extreme_ray
from scripts.fresh_multibasepoint_h9_merge import ROOT,read,matrix,mv,dot,enc,encm,coxeter_pair_value,G

def compatible(a,b,gram):
    product=dot(a,mv(gram,b))
    return product.sign()<=0 and coxeter_pair_value(4*product*product/(dot(a,mv(gram,a))*dot(b,mv(gram,b))))

def escape(walls,gram,seed):
    phi=(1+5**.5)/2;real=lambda x:float(x.a)+float(x.b)*phi
    normals=tuple(mv(gram,r) for r in walls);f=np.array([[real(x) for x in n] for n in normals])
    gf=np.array([[real(x) for x in row] for row in gram]);rng=np.random.default_rng(seed)
    for _ in range(400):
        lp=linprog(rng.normal(size=9),A_ub=f[:,1:],b_ub=-f[:,0],bounds=[(None,None)]*9,method='highs')
        if not lp.success:continue
        point=np.r_[1.,lp.x]
        if point@gf@point<=1e-9:continue
        active=tuple(map(int,np.flatnonzero(abs(f@point)<1e-7)))
        if len(active)!=9:continue
        try:w=verify_hyperideal_extreme_ray(tuple(walls),gram,active,9)
        except Exception:continue
        return {'vector':enc(w.ray),'norm':w.norm.to_json(),'wall_pairings':enc(w.wall_pairings),'active_facets':list(w.active_facets)}
    return None

def slice_obstruction(cases,old,gram,u):
    from fractions import Fraction as F
    z=(G(0),*(G(k) for k in (7,14,20,16,12,8,4,10)),G(0))
    project=lambda w:(w[0],dot(z,mv(gram,w))/14,w[9])
    v=tuple(map(G.from_json,cases[0]['hyperideal_witness']['vector']))
    assert set(mv(u,r) for r in old)==set(old)
    records=[]
    for name,target,indices,coefficients in (
        ('escape',v,(107,1,8),(G(-8,5),G(F(-16,3),F(10,3)),G(F(-11,2),F(7,2)))),
        ('swapped_escape',mv(u,v),(101,1,8),(G(9,F(-11,2)),G(0),G(F(-21,2),F(13,2))))):
        projected=project(target)
        recombined=tuple(sum((c*project(old[i])[j] for i,c in zip(indices,coefficients)),G(0)) for j in range(3))
        assert projected==recombined and all(c.sign()>=0 for c in coefficients)
        records.append({'target':name,'projection':enc(projected),'old_wall_indices':list(indices),'nonnegative_coefficients':enc(coefficients)})
    return {'status':'EXACT_FIXED_SLICE_OBSTRUCTION','case':0,'escape_norm':dot(v,mv(gram,v)).to_json(),
        'projection_coordinates':'(time, pairing_with_z/14, last_scalar)','identities':records,
        'scope':'Every old-compatible real wall normal in S or U(S) has nonpositive pairing with this spacelike escape. This excludes completion by cuts confined to these slices, not the ambient host.'}

def run():
    start=time.monotonic();search_path='results/fresh_rank3_h9_search_20261013.json';d=read(search_path);base=read(d['source'])
    gram=matrix(base['gram']);old=matrix(base['all_walls']);u=matrix(base['standard_block_swap_involution'])
    inside=tuple(map(G.from_json,base['strict_interior']));target=tuple(map(G.from_json,d['targeted_escape']['vector']))
    roots=tuple(tuple(map(G.from_json,c['root'])) for c in d['cuts'])
    orbits=tuple(tuple(dict.fromkeys((r,mv(u,r)))) for r in roots)
    assert all(len(orbit)==2 and compatible(*orbit,gram) for orbit in orbits)
    adjacency=tuple(tuple(i!=j and all(compatible(a,b,gram) for a in orbits[i] for b in orbits[j]) for j in range(len(orbits))) for i in range(len(orbits)))
    maximal=[]
    for size in range(1,len(orbits)+1):
        for selection in combinations(range(len(orbits)),size):
            if not all(adjacency[i][j] for i,j in combinations(selection,2)):continue
            if any(all(adjacency[k][j] for j in selection) for k in range(len(orbits)) if k not in selection):continue
            maximal.append(selection)
    cases=[]
    for case,selection in enumerate(maximal):
        walls=tuple(dict.fromkeys((*old,*(r for i in selection for r in orbits[i]))))
        normals=tuple(mv(gram,r) for r in walls)
        assert len(walls)==len(old)+2*len(selection)
        assert all(dot(inside,n).sign()<0 for n in normals)
        assert any(dot(target,n).sign()>0 for n in normals)
        for i in range(len(old),len(walls)):
            assert all(compatible(walls[i],walls[j],gram) for j in range(i))
        witness=escape(walls,gram,20261013+case)
        cases.append({'selection':list(selection),'wall_count':len(walls),'all_walls':encm(walls),
            'strict_interior':base['strict_interior'],'hyperideal_witness':witness,
            'status':'EXACT_SPACELIKE_ESCAPE' if witness else 'VOLUME_UNDECIDED_NO_NUMERICAL_ESCAPE_FOUND'})
        print(json.dumps({'case':case,'selection':selection,'walls':len(walls),'escape_norm':witness['norm'] if witness else None}),flush=True)
    payload={'search_source':search_path,'base_source':d['source'],'alpha':base['alpha'],'dimension':9,'gram':base['gram'],
        'standard_block_swap_involution':base['standard_block_swap_involution'],'orbits':[[enc(r) for r in orbit] for orbit in orbits],
        'orbit_compatibility':adjacency,'maximal_compatible_selections':[list(s) for s in maximal],
        'cases':cases,'slice_obstruction':slice_obstruction(cases,old,gram,u),'seconds':round(time.monotonic()-start,3),
        'scope':'Every maximal compatible subset of these seven cut orbits, not every root of the ternary slice or ambient lattice.'}
    (ROOT/'results/fresh_rank3_h9_cones_20261013.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')

if __name__=='__main__':run()
