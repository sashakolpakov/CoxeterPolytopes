#!/usr/bin/env python3
"""Bounded targeted time/A4-direction/A1 slices of a changed prime-five host.

The five spatial directions are chosen from the fundamental-weight escapes,
not by increasing the full nine-dimensional Vinberg cutoff. Every slice has
its scalar coordinate free. All original-prefix and selected cross angles
are checked exactly. This is not a complete H9 root search.
"""
import argparse
from itertools import groupby
import json
from pathlib import Path
import time
from fresh_cubic9_positive_20261013 import *

SOURCE=ROOT/'results/fresh_cubic9_positive_20261013_h9_theta_h30.json'
OUTPUT=ROOT/'results/fresh_cubic9_positive_20261013_target.json'
DIRECTIONS=tuple(tuple(5*AINV[i][j] for i in range(4)) for j in range(4))+((O,K(2),K(2),O),)

def run(height=1000):
    start=time.monotonic();data=json.loads(SOURCE.read_text());host=Host(T,2);g=host.G
    old=tuple(tuple(map(K.from_json,r)) for r in data['all_walls']);cov=tuple(mv(g,r) for r in old)
    oldnorms=tuple(dot(r,c) for r,c in zip(old,cov));candidate_rows=[];dictionary={};control_count=0
    norms=(O,K(2),K(5),K(10))
    for q in norms:
        for x in alphabet((height*q/T,-q/T,-q/T)):
            if x.sign()<=0 or not (2*T*x/q).is_integral:continue
            control_count+=1;energy=q+T*x*x
            scalars=tuple(z for z in alphabet((energy/2,)*3) if z.sign()>=0 and (4*z/q).is_integral)
            squares={2*z*z:z for z in scalars}
            for di,d in enumerate(DIRECTIONS):
                dn=pair(d,A4,d);dc=mv(A4,d)
                for t in alphabet((energy/dn,)*3):
                    if t.sign()<0 or not all((2*t*w/q).is_integral for w in dc):continue
                    z=squares.get(energy-dn*t*t)
                    if z is None:continue
                    r=(x,Z,Z,Z,Z,*(t*v for v in d),z)
                    assert all(v.is_integral for v in r) and pair(r,g,r)==q
                    assert all((2*v/q).is_integral for v in mv(g,r))
                    key=ray_key(r)
                    if key in dictionary:continue
                    allowed=compatible(r,q,cov,oldnorms)
                    row={'norm':q.to_json(),'height':(T*x*x/q).to_json(),'direction_index':di,
                        'time':x.to_json(),'multiple':t.to_json(),'scalar':z.to_json(),'vector':enc(r),
                        'compatible_with_original_prefix':allowed}
                    dictionary[key]=row;candidate_rows.append(row)
    ordered=sorted((K.from_json(r['height']),r) for r in candidate_rows if r['compatible_with_original_prefix'])
    walls=list(old);kept=[];rejected=[]
    for h,packet in groupby(ordered,key=lambda row:row[0]):
        for _,row in packet:
            r=tuple(map(K.from_json,row['vector']));q=K.from_json(row['norm'])
            cc=tuple(mv(g,w) for w in walls);nn=tuple(dot(w,c) for w,c in zip(walls,cc))
            if ray_key(r) in {ray_key(w) for w in walls}:continue
            if compatible(r,q,cc,nn):walls.append(r);kept.append(row)
            else:rejected.append(row)
    inside=tuple(map(K.from_json,data['strict_interior']))
    # The original center can be beyond a new cap; choose a sufficiently small
    # displacement from the same regular zero-height vertex if necessary.
    while any(pair(inside,g,r).sign()>=0 for r in walls):inside=(O,*(v/2 for v in inside[1:]))
    assert pair(inside,g,inside).sign()<0 and all(pair(inside,g,r).sign()<0 for r in walls)
    out={'status':'COMPLETE_BOUNDED_TARGETED_SLICE_SEARCH','source':str(SOURCE.relative_to(ROOT)),
        'height_bound':height,'time_norm_controls':control_count,'norm_pool':enc(norms),
        'directions':encode_matrix(DIRECTIONS),'all_slice_candidates':candidate_rows,
        'compatible_candidates':len(ordered),'selected_new_roots':kept,'cross_incompatible_roots':rejected,
        'all_walls':encode_matrix(walls),'wall_count':len(walls),'strict_interior':enc(inside),
        'finite_angle_counts':finite_angle_counts(walls,g),'geometry':geometry(walls,g,9,inside=inside),
        'seconds':round(time.monotonic()-start,3),
        'scope':'Exact complete enumeration in five time/A4-ray/scalar slices through the stated height, with full original and selected cross-angle gates. No complete H9 enumeration or full-host conclusion.'}
    OUTPUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:out[k] for k in ('status','time_norm_controls','compatible_candidates','wall_count','seconds')}));return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--height',type=int,default=1000);a=p.parse_args();run(a.height)
