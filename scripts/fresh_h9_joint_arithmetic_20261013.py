#!/usr/bin/env python3
"""Complete same-bound A4+A1 face completion in the prime-five H9 host.

All four A4 Cartan covectors and the scalar are free. A rank-one cross term
in the inverse A4 Cartan matrix permits an exact residue-indexed pair join.
The root dictionary is conditioned only on the original height-30 prefix,
not on later targeted walls; chronological replacement remains possible.
"""
import argparse
from itertools import groupby
import json
from pathlib import Path
import time
from fresh_cubic9_positive_20261013 import *
SOURCE=ROOT/'results/fresh_cubic9_positive_20261013_h9_theta_h30.json'
PREVIOUS=ROOT/'results/fresh_cubic9_positive_20261013_rank4.json'
OUTPUT=ROOT/'results/fresh_h9_joint_arithmetic_20261013.json'

@lru_cache(None)
def candidates(x,q):
    energy=q+T*x*x
    weights=tuple(q*k/2 for k in alphabet((8*energy/(q*q),)*3)
        if k.sign()>=0 and (q*k/2).is_integral)
    limit=T*x
    ss={2*z*z:z for z in alphabet((energy/2,)*3)
        if z.sign()>=0 and (4*z/q).is_integral and (limit-2*z).sign()>=0}
    pairs=[];groups={}
    residue=lambda v:tuple(int(t)%5 for t in (v.a,v.b,v.c))
    for a in weights:
        if (limit-a).sign()<0:continue
        for b in weights:
            if (limit-a-b).sign()<0:continue
            # Minimum possible energy with two adjacent covectors fixed.
            lower=K(F(2,3))*(a*a+a*b+b*b)
            if not le_all(lower,energy):continue
            A=a+2*b;el=(4*a*a+6*a*b+6*b*b)/5
            row=(a,b,A,el,a+b);pairs.append(row);groups.setdefault(residue(A),[]).append(row)
    out=[];joined=0
    for a,b,A,el,s in pairs:
        for d,c,B,er,t in groups[residue(A)]:
            if (limit-s-t).sign()<0:continue
            joined+=1;e=el+er+2*A*B/5
            z=ss.get(energy-e)
            if z is None:continue
            w=(a,b,c,d);y=mv(AINV,w)
            assert all(v.is_integral for v in y)
            r=(x,Z,Z,Z,Z,*y,z)
            assert pair(r,Host(T,2).G,r)==q
            assert all((2*v/q).is_integral for v in mv(Host(T,2).G,r))
            out.append(r)
    assert len(set(out))==len(out)
    return tuple(out),{'weights':len(weights),'pair_vectors':len(pairs),'joined_pairs':joined,'scalar_vectors':len(ss)}

def run(height=1000):
    started=time.monotonic();base=json.loads(SOURCE.read_text());previous=json.loads(PREVIOUS.read_text())
    host=Host(T,2);g=host.G;old=tuple(tuple(map(K.from_json,r))for r in base['all_walls'])
    prior=tuple(tuple(map(K.from_json,r))for r in previous['all_walls'])
    oldcov=tuple(mv(g,r)for r in old);oldnorms=tuple(dot(r,c)for r,c in zip(old,oldcov))
    pcov=tuple(mv(g,r)for r in prior);pnorms=tuple(dot(r,c)for r,c in zip(prior,pcov))
    oldescape=tuple(map(K.from_json,previous['geometry']['hyperideal_witness']['vector']))
    assert all(v==Z for v in oldescape[1:5])
    assert pair(oldescape,g,oldescape).sign()>0 and all(dot(oldescape,c).sign()<=0 for c in pcov)
    controls=sorted((T*x*x/q,q,x)for q in (O,K(2),K(5),K(10))
        for x in alphabet((height*q/T,-q/T,-q/T)) if x.sign()>0 and (2*T*x/q).is_integral)
    dictionary=[];levels=[];walls=list(old);keys={ray_key(r)for r in old}
    print(json.dumps({'phase':'start','controls':len(controls),'height':height}),flush=True)
    for h,packet in groupby(controls,key=lambda row:row[0]):
        accepted=[];cov=tuple(mv(g,r)for r in walls);nn=tuple(dot(r,c)for r,c in zip(walls,cov))
        for _,q,x in packet:
            roots,counts=candidates(x,q);allkept=0;chosen=0
            for r in roots:
                if ray_key(r) in {ray_key(w)for w in old}:continue
                if not compatible(r,q,oldcov,oldnorms):continue
                allkept+=1
                row={'height':h.to_json(),'norm':q.to_json(),'time':x.to_json(),'vector':enc(r),
                    'compatible_with_previous22':compatible(r,q,pcov,pnorms),
                    'previous_escape_pairing':pair(r,g,oldescape).to_json()}
                dictionary.append(row)
                if ray_key(r) in keys or not compatible(r,q,cov,nn):continue
                keys.add(ray_key(r));accepted.append(r);chosen+=1
            levels.append({'height':h.to_json(),'norm':q.to_json(),'time':x.to_json(),
                'raw_candidates_after_simple_prefix_bounds':len(roots),'original13_compatible':allkept,'accepted':chosen,**counts})
        finite_angle_counts(accepted,g);walls.extend(accepted)
        print(json.dumps({'height':str(h),'added':len(accepted),'walls':len(walls)}),flush=True)
    assert len({ray_key(tuple(map(K.from_json,r['vector'])))for r in dictionary})==len(dictionary)
    inside=tuple(map(K.from_json,base['strict_interior']))
    while any(pair(inside,g,r).sign()>=0 for r in walls):inside=(O,*(v/2 for v in inside[1:]))
    assert pair(inside,g,inside).sign()<0 and all(pair(inside,g,r).sign()<0 for r in walls)
    result={'status':'COMPLETE_SAME_BOUND_ALL_A4_AND_SCALAR_COMPLETION','source':str(SOURCE.relative_to(ROOT)),
        'previous_source':str(PREVIOUS.relative_to(ROOT)),'height_bound':height,
        'norm_pool':enc((O,K(2),K(5),K(10))),'control_count':len(controls),'levels':levels,
        'compatible_dictionary':dictionary,'all_walls':encode_matrix(walls),'wall_count':len(walls),
        'strict_interior':enc(inside),'finite_angle_counts':finite_angle_counts(walls,g),
        'geometry':geometry(walls,g,9,inside=inside)}
    # Extract the actual orthogonal face. Removed B simple walls vanish and
    # positive B-only walls restrict to strictly negative time inequalities.
    face_host=Host(T,0);fg=face_host.G
    face_walls=tuple((r[0],*r[5:])for r in walls if all(t==Z for t in r[1:5]))
    face_inside=(inside[0],*inside[5:])
    assert pair(face_inside,fg,face_inside).sign()<0 and all(pair(face_inside,fg,r).sign()<0 for r in face_walls)
    result['forced_h5_face']={'gram':encode_matrix(fg),'all_walls':encode_matrix(face_walls),'wall_count':len(face_walls),
        'strict_interior':enc(face_inside),'finite_angle_counts':finite_angle_counts(face_walls,fg),
        'geometry':geometry(face_walls,fg,5,inside=face_inside)}
    result['seconds']=round(time.monotonic()-started,3)
    result['scope']='Complete height-bounded integral root census in the full time+A4+A1 complement, conditioned only on the original13-wall chamber. Both missing A4 covectors and scalar are free. This is not a complete H9 root census.'
    OUTPUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'phase':'done','walls':len(walls),'face_walls':len(face_walls),'dictionary':len(dictionary),'seconds':result['seconds']}),flush=True)
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--height',type=int,default=1000);a=p.parse_args();run(a.height)
