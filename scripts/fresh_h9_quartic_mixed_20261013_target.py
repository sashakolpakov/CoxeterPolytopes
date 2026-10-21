#!/usr/bin/env python3
"""New rank-three, two-block slice aimed at the 62-wall escape."""
import argparse
from itertools import product
import json
import time
from fresh_h9_quartic_mixed_20261013 import *

SOURCE=ROOT/'results/fresh_h9_quartic_mixed_20261013_h1000.json'


def orbit(r):
    out=set()
    for reverse_a,reverse_b,swap in product((False,True),repeat=3):
        a=r[1:5][::-1]if reverse_a else r[1:5]
        b=r[5:9][::-1]if reverse_b else r[5:9]
        if swap:a,b=b,a
        out.add((r[0],*a,*b,r[9]))
    return out


def run(height,output):
    started=time.monotonic();data=json.loads(SOURCE.read_text());g=Host(2).g
    old=[tuple(map(K.from_json,r))for r in data['all_walls']]
    cov=[mv(g,r)for r in old];norms=[dot(r,c)for r,c in zip(old,cov)]
    keys={ray_key(r)for r in old};target=tuple(map(K.from_json,data['geometry']['vector']))
    controls=sorted((ALPHA*x*x/q,q,x)for q in NORMS
                   for x in balanced_alphabet((height*q/ALPHA,-q/ALPHA,-q/ALPHA,-q/ALPHA))if x>0)
    dictionary={};levels=[]
    print(json.dumps({'slice':'time + first-block omega0 + second-block omega1, scalar zero',
                      'height':height,'controls':len(controls)}),flush=True)
    for index,(h,q,x)in enumerate(controls):
        energy=q+ALPHA*x*x
        us=balanced_alphabet((energy/2,)*4)
        vs={6*v*v:v for v in balanced_alphabet((energy/6,)*4)}
        roots=[];accepted=0;cuts=0
        for u in us:
            rem=energy-2*u*u
            if rem not in vs:continue
            v=vs[rem]
            if not u or not v:continue
            r=(x,*(row[0]*u for row in FI),*(row[1]*v for row in FI),Z)
            assert pair(r,g,r)==q and all(z.is_integral for z in r)
            if all((z/T).is_integral for z in r):continue
            assert all((2*z/q).is_integral for z in mv(g,r))
            roots.append(r)
            for rr in orbit(r):
                key=ray_key(rr)
                if key in keys or key in dictionary:continue
                if not coxeter(rr,q,cov,norms):continue
                accepted+=1;cuts+=pair(rr,g,target)>0
                dictionary[key]={'height':h.to_json(),'norm':q.to_json(),'vector':enc(rr),
                                 'old62_escape_pairing':pair(rr,g,target).to_json()}
        levels.append({'height':h.to_json(),'norm':q.to_json(),'time':x.to_json(),
                       'u_alphabet':len(us),'v_alphabet':len(vs),'roots':menc(roots),
                       'accepted_orbit_images':accepted,'cutting_images':cuts})
        if roots or index%50==0:
            print(json.dumps({'index':index,'height':h.real(),'roots':len(roots),'accepted':accepted,
                              'cuts':cuts,'elapsed':round(time.monotonic()-started,2)}),flush=True)
    rows=sorted(dictionary.values(),key=lambda row:(K.from_json(row['height']),tuple(map(K.from_json,row['vector']))))
    walls=list(old);wc=list(cov);wq=list(norms);chosen=[];rejected=[]
    for i,row in enumerate(rows):
        r=tuple(map(K.from_json,row['vector']));q=K.from_json(row['norm'])
        if coxeter(r,q,wc,wq):chosen.append(i);walls.append(r);wc.append(mv(g,r));wq.append(q)
        else:rejected.append(i)
    inside=tuple(map(K.from_json,data['strict_interior']))
    while any(dot(inside,c)>=0 for c in wc):inside=(inside[0],*(z/2 for z in inside[1:]))
    out={'status':'COMPLETE_TARGETED_TWO_BLOCK_RANK_THREE_CENSUS','source':str(SOURCE.relative_to(ROOT)),
         'height_bound':height,'control_count':len(controls),'norm_pool':enc(NORMS),'levels':levels,
         'slice_gram':menc(((-ALPHA,Z,Z),(Z,K(2),Z),(Z,Z,K(6)))),
         'compatible_dictionary':rows,'selected_indices':chosen,'rejected_indices':rejected,
         'all_walls':menc(walls),'wall_count':len(walls),'strict_interior':enc(inside),
         'geometry':geometry(walls,g,inside,9),'seconds':round(time.monotonic()-started,3),
         'scope':'Complete bounded slice (time, omega0 in one F4 block, omega1 in the other), with both spatial weights strictly positive and scalar zero, including eight diagram-reversal/block-swap placements. Boundary one-weight planes are not included. Retains every original62-compatible new root in this stated bounded set; no full-height or full-host conclusion.'}
    output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'dictionary':len(rows),'chosen':len(chosen),'walls':len(walls),
                      'geometry':out['geometry']['status'],'seconds':out['seconds']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--height',type=int,default=100000)
    p.add_argument('--output',type=Path,default=ROOT/'results/fresh_h9_quartic_mixed_20261013_target.json')
    a=p.parse_args();run(a.height,a.output)
