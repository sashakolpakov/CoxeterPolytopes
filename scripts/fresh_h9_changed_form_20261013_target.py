#!/usr/bin/env python3
"""Bounded constructive cuts outside the forbidden end-weight binary ray.

Allow both B coordinates, two or three A6 end covectors, and the scalar.
Every necessary primitive norm class is included, including 14beta and 28beta.
"""
import argparse
import json
import time
from itertools import product
from fresh_h9_changed_form_20261013 import *
SOURCE=ROOT/'results/fresh_h9_changed_form_20261013_h9_2theta_h100.json'
OUTPUT=ROOT/'results/fresh_h9_changed_form_20261013_target.json'

def run(height=1000,all_end_pairs=False,three_adjacent=False):
    started=time.monotonic();data=json.loads(SOURCE.read_text());host=Host(2*T,1);g=host.G
    old=tuple(tuple(map(K.from_json,r))for r in data['all_walls'])
    cov=tuple(mv(g,r)for r in old);nn=tuple(pair(r,g,r)for r in old)
    target=tuple(map(K.from_json,data['geometry']['hyperideal_witness']['vector']))
    assert all(target[i]==Z for i in(1,2,9))
    controls=sorted((2*T*x*x/q,q,x)for q in host.pool
        for x in alphabet((height*q/(2*T),-q/(2*T),-q/(2*T)))
        if x.sign()>0 and(4*T*x/q).is_integral)
    levels=[];dictionary={};supports=((3,4,5),)if three_adjacent else tuple((i,5)for i in range(5))if all_end_pairs else((4,5),)
    for h,q,x in controls:
        energy=q+2*T*x*x
        weights=tuple(q*k/2 for k in alphabet((8*energy/(q*q),)*3)
                      if k.sign()>=0 and(q*k/2).is_integral)
        bs=[]
        for i,u in enumerate(weights):
            for v in weights[i:]:
                y=((2*u+T*v)/BETA,(T*u+2*v)/BETA);e=dot(y,(u,v))
                if all(t.is_integral for t in y)and le_all(e,energy):bs.append((e,y))
        ss=tuple((2*z*z,z)for z in alphabet((energy/2,)*3)
                 if z.sign()>=0 and(4*z/q).is_integral)
        remaining={}
        for eb,yb in bs:
            for ez,z in ss:
                if le_all(eb+ez,energy):remaining.setdefault(eb+ez,[]).append((*yb,z))
        aa=0;raw=0;acute_reject=0;angle_reject=0;imprimitive=0;kept=0
        for support in supports:
          residue=lambda w:tuple(int(t)%7 for t in(w.a,w.b,w.c))
          last={}
          for w in weights:last.setdefault(residue((support[-1]+1)*w),[]).append(w)
          for prefix in product(weights,repeat=len(support)-1):
            need=residue(-sum(((i+1)*w for i,w in zip(support[:-1],prefix)),Z))
            for v in last.get(need,()):
                ws=[Z]*6
                for i,w in zip(support,(*prefix,v)):ws[i]=w
                ws=tuple(ws);ys=mv(AI,ws)
                if not all(t.is_integral for t in ys):continue
                r0=(x,Z,Z,*ys,Z)
                if pair(r0,g,target).sign()<=0:continue
                ea=dot(ys,ws)
                if not le_all(ea,energy):continue
                aa+=1
                for b1,b2,z in remaining.get(energy-ea,()):
                    r=(x,b1,b2,*ys,z);raw+=1
                    assert pair(r,g,r)==q and all((2*t/q).is_integral for t in mv(g,r))
                    if all((t/2).is_integral for t in r):imprimitive+=1;continue
                    if any(dot(r,c).sign()>0 for c in cov):acute_reject+=1;continue
                    if not compatible(r,q,cov,nn):angle_reject+=1;continue
                    kept+=1
                    for rr in(r,(r[0],r[1],r[2],*reversed(r[3:9]),r[9])):
                        assert compatible(rr,q,cov,nn)
                        dictionary[ray_key(rr)]={'height':h.to_json(),'norm':q.to_json(),'vector':enc(rr),
                            'old_escape_pairing':pair(rr,g,target).to_json()}
        levels.append({'height':h.to_json(),'norm':q.to_json(),'time':x.to_json(),'weights':len(weights),
            'block_vectors':len(bs),'scalar_vectors':len(ss),'target_cutting_a6_vectors':aa,
            'raw_target_cutters':raw,'imprimitive':imprimitive,'positive_pair_rejected':acute_reject,
            'noncoxeter_rejected':angle_reject,'admissible_target_cutters':kept})
        if raw or kept:print(json.dumps({'height':str(h),'norm':str(q),'raw':raw,'kept':kept,'dictionary':len(dictionary)}),flush=True)
    rows=sorted(dictionary.values(),key=lambda row:(K.from_json(row['height']),tuple(map(K.from_json,row['vector']))))
    walls=list(old);chosen=[];rejected=[]
    for i,row in enumerate(rows):
        r=tuple(map(K.from_json,row['vector']));q=K.from_json(row['norm'])
        cc=tuple(mv(g,w)for w in walls);ns=tuple(pair(w,g,w)for w in walls)
        if compatible(r,q,cc,ns):chosen.append(i);walls.append(r)
        else:rejected.append(i)
    inside=tuple(map(K.from_json,data['strict_interior']))
    while any(pair(inside,g,r).sign()>=0 for r in walls):inside=(O,*(v/2 for v in inside[1:]))
    out={'status':'COMPLETE_BOUNDED_END_COVECTOR_B_SCALAR_TARGET','source':str(SOURCE.relative_to(ROOT)),
        'height_bound':height,'norm_pool':enc(host.pool),'control_count':len(controls),
        'a6_covector_supports':[list(s)for s in supports],'block_coordinates':'both free','scalar_coordinate':'free',
        'levels':levels,'compatible_dictionary':rows,'selected_indices':chosen,'rejected_indices':rejected,
        'all_walls':encode_matrix(walls),'wall_count':len(walls),'strict_interior':enc(inside),
        'finite_angle_counts':finite_angle_counts(walls,g),'geometry':geometry(walls,g,9,inside=inside),
        'seconds':round(time.monotonic()-started,3),
        'scope':'Complete search through the stated height in the saturated A6 covector slices listed, with full B and scalar freedom, filtered for cuts of the saved 23-wall escape. Reversed mirrors are added and a compatible height-ordered selection is checked. No full H9 height census or global nonreflectivity claim.'}
    suffix='_three_adjacent'if three_adjacent else '_all_end_pairs'if all_end_pairs else ''
    output=OUTPUT.with_name(OUTPUT.stem+suffix+'.json')
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps({'dictionary':len(rows),'selected':len(chosen),'walls':len(walls),'geometry':out['geometry']['status'],'seconds':out['seconds']}),flush=True);return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--all-end-pairs',action='store_true');p.add_argument('--three-adjacent',action='store_true');a=p.parse_args();run(all_end_pairs=a.all_end_pairs,three_adjacent=a.three_adjacent)
