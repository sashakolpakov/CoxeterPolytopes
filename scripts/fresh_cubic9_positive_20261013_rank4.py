#!/usr/bin/env python3
"""Complete bounded saturated two-weight A4 plus scalar slices."""
import argparse
from itertools import combinations,groupby
import json
import time
from fresh_cubic9_positive_20261013 import *
SOURCE=ROOT/'results/fresh_cubic9_positive_20261013_h9_theta_h30.json'
OUTPUT=ROOT/'results/fresh_cubic9_positive_20261013_rank4.json'

def run(height=1000):
    start=time.monotonic();data=json.loads(SOURCE.read_text());host=Host(T,2);g=host.G
    old=tuple(tuple(map(K.from_json,r)) for r in data['all_walls']);cov=tuple(mv(g,r) for r in old)
    oldnorms=tuple(dot(r,c) for r,c in zip(old,cov));candidate_rows=[];seen=set();controls=[]
    for q in (O,K(2),K(5),K(10)):
        for x in alphabet((height*q/T,-q/T,-q/T)):
            if x.sign()<=0 or not (2*T*x/q).is_integral:continue
            energy=q+T*x*x;weights=tuple(q*k/2 for k in alphabet((8*energy/(q*q),)*3)
                if k.sign()>=0 and (q*k/2).is_integral)
            count=0
            ss={2*z*z:z for z in alphabet((energy/2,)*3) if z.sign()>=0 and (4*z/q).is_integral}
            for i,j in combinations(range(4),2):
                for u in weights:
                    for v in weights:
                        ww=tuple(u if k==i else v if k==j else Z for k in range(4))
                        # A4 inverse integrality is one exact O/(5) congruence.
                        if not (((i+1)*u+(j+1)*v)/5).is_integral:continue
                        e=AINV[i][i]*u*u+2*AINV[i][j]*u*v+AINV[j][j]*v*v
                        z=ss.get(energy-e)
                        if z is None:continue
                        y=mv(AINV,ww);r=(x,Z,Z,Z,Z,*y,z)
                        assert all(t.is_integral for t in r) and pair(r,g,r)==q
                        assert all((2*t/q).is_integral for t in mv(g,r))
                        key=ray_key(r)
                        if key in seen:continue
                        seen.add(key);count+=1
                        candidate_rows.append({'norm':q.to_json(),'height':(T*x*x/q).to_json(),'weight_indices':[i,j],
                            'time':x.to_json(),'weights':enc((u,v)),'scalar':z.to_json(),'vector':enc(r),
                            'compatible_with_original_prefix':compatible(r,q,cov,oldnorms)})
            controls.append({'height':(T*x*x/q).to_json(),'norm':q.to_json(),'time':x.to_json(),'alphabet_count':len(weights),'new_candidates':count})
    ordered=sorted(((K.from_json(r['height']),r) for r in candidate_rows if r['compatible_with_original_prefix']),key=lambda row:row[0])
    walls=list(old);kept=[];rejected=[]
    for _,row in ordered:
        r=tuple(map(K.from_json,row['vector']));q=K.from_json(row['norm'])
        if ray_key(r) in {ray_key(w) for w in walls}:continue
        cc=tuple(mv(g,w) for w in walls);nn=tuple(dot(w,c) for w,c in zip(walls,cc))
        if compatible(r,q,cc,nn):walls.append(r);kept.append(row)
        else:rejected.append(row)
    inside=tuple(map(K.from_json,data['strict_interior']))
    while any(pair(inside,g,r).sign()>=0 for r in walls):inside=(O,*(v/2 for v in inside[1:]))
    assert pair(inside,g,inside).sign()<0 and all(pair(inside,g,r).sign()<0 for r in walls)
    out={'status':'COMPLETE_BOUNDED_TWO_WEIGHT_SCALAR_SEARCH','source':str(SOURCE.relative_to(ROOT)),
        'height_bound':height,'controls':controls,'norm_pool':enc((O,K(2),K(5),K(10))),
        'all_slice_candidates':candidate_rows,'compatible_candidates':len(ordered),
        'selected_new_roots':kept,'cross_incompatible_roots':rejected,
        'all_walls':encode_matrix(walls),'wall_count':len(walls),'strict_interior':enc(inside),
        'finite_angle_counts':finite_angle_counts(walls,g),'geometry':geometry(walls,g,9,inside=inside),
        'seconds':round(time.monotonic()-start,3),
        'scope':'Every integral root supported on time, at most two nonnegative A4 Cartan covectors, and nonnegative scalar, with binary coordinates zero, through the stated height. All six saturated slices are included. This is not a complete H9 enumeration.'}
    OUTPUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:out[k] for k in ('status','compatible_candidates','wall_count','seconds')}));return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--height',type=int,default=1000);a=p.parse_args();run(a.height)
