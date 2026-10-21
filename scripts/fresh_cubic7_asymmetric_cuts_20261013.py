#!/usr/bin/env python3
"""Exact cubic-seven cuts with one free dihedral block and one equal block."""
import argparse
from itertools import permutations
import json
from pathlib import Path
import time
from fresh_cubic7_dihedral_h9_20261013 import K,T,Z,O,C,D,BETA,ALLOWED,decode_matrix,enc,mv,dot,pair,ray_key,finite_angle_counts,geometry
from fresh_cubic7_target_alphabet_20261013 import thin_alphabet

ROOT=Path(__file__).resolve().parents[1]

def compatible(r,q,cov,norms):
    for n,s in zip(cov,norms):
        p=dot(r,n)
        if p.sign()>0:return False
        value=4*p*p/(q*s)
        if (value-4).sign()<0 and value not in ALLOWED:return False
    return True

def run(height,output):
    started=time.monotonic();source=ROOT/'results/fresh_cubic7_dihedral_h9_2theta_h30_20261013.json'
    saved=json.loads(source.read_text());g=decode_matrix(saved['gram']);old=decode_matrix(saved['all_walls'])
    cov=tuple(mv(g,r) for r in old);norms=tuple(dot(r,n) for r,n in zip(old,cov));v=tuple(map(K.from_json,saved['geometry']['hyperideal_witness']['vector']))
    center=tuple(map(K.from_json,saved['strict_interior']));alpha=2*T
    rows=[];keys={ray_key(r) for r in old};census=[]
    for q in (K(2),2*BETA):
        times=tuple(x for x in thin_alphabet((height*q/alpha,-q/alpha,-q/alpha))
                    if x.sign()>0 and (2*alpha*x/q).is_integral)
        aa=tuple(a for a in thin_alphabet((q*(height+1)/2,q/2,q/2))
                 if a.sign()>=0 and (2*D*a/q).is_integral)
        aa_squares={2*a*a:a for a in aa}
        ww=tuple(w for w in thin_alphabet((2*q*(height+1),2*q,2*q))
                 if w.sign()>=0 and (2*w/q).is_integral)
        blocks=[]
        for i,u in enumerate(ww):
            for w in ww[i+1:]:
                y=((2*u+C*w)/BETA,(C*u+2*w)/BETA)
                if not all(t.is_integral for t in y):continue
                energy=u*y[0]+w*y[1]
                if (q*(height+1)-energy).sign()<0 or any((q-energy).sign(e)<0 for e in (1,2)):continue
                blocks.append((energy,y,(u,w)))
        count=0;compatible_count=0;cutters=0
        for x in times:
            target=q+alpha*x*x
            for energy,y,weights in blocks:
                a=aa_squares.get(target-energy)
                if a is None:continue
                count+=1
                for free,equal in permutations(range(4),2):
                    r=[x]+[Z]*9;r[1+2*free:3+2*free]=y;r[1+2*equal:3+2*equal]=(a/D,a/D);r=tuple(r)
                    assert pair(r,g,r)==q and all((2*t/q).is_integral for t in mv(g,r))
                    if not compatible(r,q,cov,norms):continue
                    key=ray_key(r)
                    if key in keys:continue
                    keys.add(key);compatible_count+=1;p=pair(r,g,v);cutters+=p.sign()>0
                    rows.append({'vector':enc(r),'norm':q.to_json(),'height':(alpha*x*x/q).to_json(),
                                 'old_escape_pairing':p.to_json(),'free_block':free,'equal_block':equal})
        census.append({'norm':q.to_json(),'time_controls':enc(times),'equal_coordinates':enc(aa),
                       'free_weights':enc(ww),'free_blocks':len(blocks),'norm_solutions':count,
                       'compatible_new_roots':compatible_count,'old_escape_cutters':cutters})
        print(json.dumps({k:(str(q) if k=='norm' else value) for k,value in census[-1].items() if k not in ('time_controls','equal_coordinates','free_weights')}),flush=True)
    rows.sort(key=lambda r:K.from_json(r['height']))
    walls=list(old);wcov=list(cov);wnorms=list(norms);selected=[]
    for row in rows:
        r=tuple(map(K.from_json,row['vector']))
        q=K.from_json(row['norm'])
        if compatible(r,q,wcov,wnorms):
            walls.append(r);wcov.append(mv(g,r));wnorms.append(q);selected.append(row)
    while any(dot(center,n).sign()>=0 for n in wcov):center=(O,*(x/2 for x in center[1:]))
    assert pair(center,g,center).sign()<0
    out={'status':'BOUNDED_ASYMMETRIC_SLICE_COMPLETE','source':str(source.relative_to(ROOT)),
         'height_bound':height,'census':census,'compatible_roots':rows,'selected_roots':selected,
         'wall_count':len(walls),'walls':[enc(r) for r in walls],'strict_interior':enc(center),
         'finite_angle_counts':finite_angle_counts(walls,g),'geometry':geometry(walls,g,9,inside=center),
         'scope':'All primitive reflections through the stated height supported on one genuinely asymmetric block and one equal-coordinate block, scalar zero. Full lattice remains undecided.',
         'seconds':round(time.monotonic()-started,3)}
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'walls':len(walls),'compatible':len(rows),'selected':len(selected),'geometry':out['geometry']['status'],'seconds':out['seconds']}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--height',type=int,default=1000);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.height,a.output)
