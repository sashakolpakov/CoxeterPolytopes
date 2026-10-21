#!/usr/bin/env python3
"""Bounded adaptive cuts of exact escape rays of the merged odd H9 cone."""
from __future__ import annotations
import argparse
import json
import time
import numpy as np
from scipy.optimize import linprog
from coxeter_polytopes.verify.chamber import verify_hyperideal_extreme_ray
from scripts.fresh_multibasepoint_h9_merge import ROOT,read,matrix,mv,dot,enc,encm,G,Z,coxeter_pair_value
from scripts.fresh_multibasepoint_h9_escape_plane import search_plane,plane_parameters

def escape_candidates(walls,gram,seed,maximum=12):
    phi=(1+5**.5)/2;real=lambda x:float(x.a)+float(x.b)*phi
    normals=tuple(mv(gram,r) for r in walls)
    f=np.array([[real(x) for x in n] for n in normals]);gf=np.array([[real(x) for x in r] for r in gram])
    rng=np.random.default_rng(seed);supports=set();out=[]
    for _ in range(200):
        lp=linprog(rng.normal(size=9),A_ub=f[:,1:],b_ub=-f[:,0],bounds=[(None,None)]*9,method='highs')
        if not lp.success:continue
        v=np.r_[1.,lp.x]
        if v@gf@v<=1e-9:continue
        active=tuple(map(int,np.flatnonzero(abs(f@v)<1e-7)))
        if len(active)!=9 or active in supports:continue
        supports.add(active)
        try:w=verify_hyperideal_extreme_ray(tuple(walls),gram,active,9)
        except Exception:continue
        _,_,beta=plane_parameters(gram,w.ray)
        out.append((int(beta.norm()),{'vector':enc(w.ray),'norm':w.norm.to_json(),
            'wall_pairings':enc(w.wall_pairings),'active_facets':list(w.active_facets)}))
        if len(out)>=maximum:break
    return sorted(out,key=lambda x:x[0])

def run(rounds,limit,trials):
    start=time.monotonic();source='results/fresh_multibasepoint_h9_2phi_20261013.json';base=read(source)
    gram=matrix(base['gram']);u=matrix(base['standard_block_swap_involution']);walls=list(matrix(base['all_walls']))
    inside=tuple(map(G.from_json,base['strict_interior']));records=[];added=[]
    def add_orbit(r,origin):
        orbit=[]
        for candidate in (r,mv(u,r)):
            if candidate in walls or candidate in orbit:continue
            assert all(dot(candidate,mv(gram,w)).sign()<=0 for w in (*walls,*orbit))
            assert dot(candidate,mv(gram,inside)).sign()<0
            q=dot(candidate,mv(gram,candidate))
            assert all(coxeter_pair_value(4*dot(candidate,mv(gram,w))**2/(q*dot(w,mv(gram,w)))) for w in (*walls,*orbit))
            orbit.append(candidate)
        for r in orbit:
            walls.append(r);added.append({'root':enc(r),'origin':origin})
        return len(orbit)
    initial=read('results/fresh_multibasepoint_h9_2phi_escape_plane_20261013.json')
    for cut in initial['cuts']:add_orbit(tuple(map(G.from_json,cut['root'])),'initial_escape_plane')
    final_escape=None
    output=ROOT/'results/fresh_multibasepoint_h9_2phi_targeted_rounds_20261013.json'
    def save(status):
        out={'source':source,'alpha':base['alpha'],'dimension':9,'gram':base['gram'],
            'all_walls':encm(walls),'wall_count':len(walls),'strict_interior':base['strict_interior'],
            'standard_block_swap_involution':base['standard_block_swap_involution'],
            'added_walls':added,'rounds':records,'status':status,'hyperideal_witness':final_escape,
            'maximum_rounds':rounds,'maximum_plane_trials_per_round':trials,'positive_t_bound':limit,
            'seconds':round(time.monotonic()-start,3),
            'scope':'Nonchronological acute integral reflection cuts; bounded adaptive discovery, not a complete Vinberg enumeration or a global obstruction.'}
        output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    save('SEARCH_IN_PROGRESS')
    for iteration in range(rounds):
        candidates=escape_candidates(walls,gram,20261013+iteration)
        if not candidates:
            records.append({'round':iteration,'status':'NO_NUMERICAL_ESCAPE_DISCOVERED_UNDECIDED'})
            save('VOLUME_UNDECIDED_REQUIRES_EXACT_CONE_AUDIT');return
        final_escape=candidates[0][1]
        progress=False
        round_record={'round':iteration,'candidate_escape_count':len(candidates),'trials':[]}
        for beta_norm,witness in candidates[:trials]:
            if beta_norm>10**10:
                round_record['trials'].append({'spatial_norm_field_norm':beta_norm,'status':'SKIPPED_LARGE_SLICE_DISCRIMINANT'});continue
            d={**base,'all_walls':encm(walls),'wall_count':len(walls),'hyperideal_witness':witness}
            result=search_plane(d,'2phi',limit,str(output.relative_to(ROOT)))
            result['targeted_escape']=witness
            round_record['trials'].append(result)
            if result['cuts']:
                cut=result['cuts'][0];count=add_orbit(tuple(map(G.from_json,cut['root'])),f'round_{iteration}')
                round_record['added_walls']=count;progress=True;final_escape=None;break
        records.append(round_record)
        print(json.dumps({'round':iteration,'walls':len(walls),'cut_found':progress}),flush=True)
        save('SEARCH_IN_PROGRESS' if progress else 'EXACT_ESCAPE_REMAINS_AFTER_BOUNDED_SLICE_SEARCH')
        if not progress:return
    candidates=escape_candidates(walls,gram,20261113,maximum=1)
    final_escape=candidates[0][1] if candidates else None
    save('EXACT_ESCAPE_REMAINS_AFTER_BOUNDED_ADAPTIVE_CUTS' if final_escape else 'VOLUME_UNDECIDED_REQUIRES_EXACT_CONE_AUDIT')
    print(json.dumps({'walls':len(walls),'seconds':round(time.monotonic()-start,3),'escape_remains':final_escape is not None}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--rounds',type=int,default=3)
    p.add_argument('--limit',type=int,default=100000);p.add_argument('--trials',type=int,default=4);a=p.parse_args()
    run(a.rounds,a.limit,a.trials)
