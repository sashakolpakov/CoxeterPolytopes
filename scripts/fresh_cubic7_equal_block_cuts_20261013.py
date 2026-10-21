#!/usr/bin/env python3
"""Target the two-block escapes by exact equal-coordinate block reflections.

This is a new restricted arithmetic slice, not a full Vinberg height increase.
For alpha=2theta, equal block coordinates (1+theta)*a have energy2a^2.
Primitive integral reflections in this slice have norm2, so solve
sum(a_i^2)-theta*x^2=1 with a_i>=0 in the distinguished embedding.
"""
from __future__ import annotations
import argparse
from itertools import combinations_with_replacement, permutations
import json
from pathlib import Path
import time
from fresh_cubic7_dihedral_h9_20261013 import K,T,Z,O,D,alphabet,decode_matrix,enc,mv,dot,pair,ray_key,finite_angle_counts,geometry
from fresh_cubic7_target_alphabet_20261013 import thin_alphabet as alphabet
from fresh_cubic7_asymmetric_cuts_20261013 import compatible as coxeter_compatible

ROOT=Path(__file__).resolve().parents[1]
def raw(x):
    assert x.is_integral
    return int(x.a),int(x.b),int(x.c)
def plus(a,b):return tuple(x+y for x,y in zip(a,b))
def minus(a,b):return tuple(x-y for x,y in zip(a,b))

def run(height,slots,output):
    start=time.monotonic()
    path=ROOT/'results/fresh_cubic7_dihedral_h9_2theta_h30_20261013.json'
    saved=json.loads(path.read_text());gram=decode_matrix(saved['gram']);old=decode_matrix(saved['all_walls'])
    escape=tuple(map(K.from_json,saved['geometry']['hyperideal_witness']['vector']))
    inside=tuple(map(K.from_json,saved['strict_interior']))
    times=tuple(x for x in alphabet((K(height)/T,-O/T,-O/T)) if x.sign()>0)
    coords=tuple(a for a in alphabet((K(height+1),O,O)) if a.sign()>=0)
    squares=tuple(raw(a*a) for a in coords)
    pairs={}
    for i,j in combinations_with_replacement(range(len(coords)),2):
        pairs.setdefault(plus(squares[i],squares[j]),[]).append((i,j))
    targets=sorted((T*x*x,x,raw(1+T*x*x)) for x in times)
    print(json.dumps({'phase':'start','slots':slots,'height':height,'time_controls':len(times),'coordinates':len(coords),'pair_energies':len(pairs)}),flush=True)
    roots=[];scount=0;compatible=0;cutters=[];keys={ray_key(r) for r in old};cov=tuple(mv(gram,r) for r in old)
    norms=tuple(dot(r,n) for r,n in zip(old,cov))
    for h,x,target in targets:
        matches=[]
        if slots==2:
            matches=[(0,0,i,j) for i,j in pairs.get(target,())]
        else:
            for energy,first in pairs.items():
                for i,j in first:
                    for k,l in pairs.get(minus(target,energy),()):
                        if j<=k:matches.append((i,j,k,l))
        scount+=len(matches)
        for ids in matches:
            for perm in sorted(set(permutations(ids))):
                r=(x,*(a for i in perm for a in (coords[i]/D,coords[i]/D)),Z)
                assert pair(r,gram,r)==2
                if not coxeter_compatible(r,K(2),cov,norms):continue
                key=ray_key(r)
                if key in keys:continue
                keys.add(key);compatible+=1
                p=pair(r,gram,escape)
                row={'vector':enc(r),'height':h.to_json(),'old_escape_pairing':p.to_json(),'coefficients':enc(tuple(coords[i] for i in perm))}
                roots.append(row)
                if p.sign()>0:cutters.append(row)
        if matches:print(json.dumps({'height':str(h),'multisets':len(matches),'compatible_total':compatible,'cutters_total':len(cutters)}),flush=True)
    walls=list(old);selected=[];wcov=list(cov);wnorms=list(norms)
    for row in roots if cutters else ():
        r=tuple(map(K.from_json,row['vector']))
        if coxeter_compatible(r,K(2),wcov,wnorms):
            walls.append(r);selected.append(row);wcov.append(mv(gram,r));wnorms.append(K(2))
    new_inside=inside
    if any(pair(new_inside,gram,w).sign()>=0 for w in walls):
        # Every positive-time wall is strictly negative at the time axis;
        # shrink all spatial coordinates to preserve the zero-height chamber.
        while any(pair(new_inside,gram,w).sign()>=0 for w in walls):
            new_inside=(O,*(x/2 for x in new_inside[1:]))
    assert pair(new_inside,gram,new_inside).sign()<0
    out={'status':'RESTRICTED_SLICE_COMPLETE','base':str(path.relative_to(ROOT)),'height_bound':height,'nonzero_block_limit':slots,
         'time_controls':len(times),'coordinate_alphabet':enc(coords),'time_alphabet':enc(times),'multiset_solutions':scount,
         'compatible_new_roots':roots,'old_escape_cutters':len(cutters),'selected_new_roots':selected,
         'selection_skipped_because_common_old_escape':not bool(cutters),
         'wall_count':len(walls),'walls':[enc(r) for r in walls],'strict_interior':enc(new_inside),
         'finite_angle_counts':finite_angle_counts(walls,gram),'geometry':geometry(walls,gram,9,inside=new_inside),
         'seconds':round(time.monotonic()-start,3),
         'scope':'Complete equal-coordinate block slice to the stated height; arbitrary block coordinates and the scalar direction are not searched.'}
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:out[k] for k in ['status','multiset_solutions','old_escape_cutters','wall_count','seconds']}),flush=True)
    print(json.dumps({'new_root_count':len(roots),'selected':len(selected),'geometry':out['geometry']['status']}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--height',type=int,default=10000)
    p.add_argument('--slots',type=int,choices=(2,4),default=4);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.height,a.slots,a.output)
