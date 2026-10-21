#!/usr/bin/env python3
"""Target a missing initial F4 middle edge in the dyadic-twisted host."""
import argparse, hashlib, json, time
from collections import Counter
from itertools import product
from fresh_h9_new_host_20261013 import *
from fresh_h9_quartic_mixed_20261013_endpoint import square_root
from fresh_h9_quartic_edge_walk_20261013 import finite_vertex

SOURCE=ROOT/'results/fresh_h9_new_host_20261013_h9_h100.json'
FACE=ROOT/'results/fresh_h9_new_host_20261013_h5_h100.json'
OUT=ROOT/'results/fresh_h9_new_host_20261013_edge.json'


def run(bound):
    started=time.monotonic();data=json.loads(SOURCE.read_text());fd=json.loads(FACE.read_text())
    g=Host(2).g;fg=Host(1).g;old=[tuple(map(K.from_json,r))for r in data['all_walls']]
    fw=[tuple(map(K.from_json,r))for r in fd['all_walls']];inside=tuple(map(K.from_json,data['strict_interior']))
    edge_ids=[0,2,3,4];edge=[fw[i]for i in edge_ids];eg=tuple(tuple(pair(r,fg,s)for s in edge)for r in edge);ei=inverse(eg)
    direction=(Z,*(row[1]for row in FI),Z);dn=pair(direction,fg,direction);assert dn==6
    profiles=[];candidates=[];passed=[];counts=Counter();alphabets=[]
    for q in NORMS:
        alph=[sorted({v for ratio in ALLOWED if(v:=square_root(q*pair(e,fg,e)*ratio/4))is not None and v.is_integral})for e in edge]
        alphabets.append({'norm':q.to_json(),'pairing_alphabets':[enc(a)for a in alph]})
        for fixed in product(*alph):
            co=mv(ei,tuple(-z for z in fixed));r0=tuple(sum((c*e[k]for c,e in zip(co,edge)),Z)for k in range(6))
            residual=q-pair(r0,fg,r0)
            if any(residual.sign(e)<=0 for e in range(4)):continue
            caps=(bound*q/ALPHA,*(-residual/ALPHA for _ in range(3)));times=balanced_alphabet(caps)
            idx=len(profiles);profiles.append({'norm':q.to_json(),'fixed_magnitudes':enc(fixed),'fixed_component':enc(r0),
                'residual_norm':residual.to_json(),'time_caps':enc(caps),'time_coordinates':enc(times)})
            for x in times:
                if x<=0:continue
                counts['quadratics']+=1;square=(residual+ALPHA*x*x)/dn;t=square_root(square)
                if t is None:continue
                counts['square_discriminants']+=1
                for ts in(t,-t):
                    r=tuple(a+ts*b for a,b in zip(r0,direction));r=(x,*r[1:]);assert pair(r,fg,r)==q
                    if not all(z.is_integral for z in r):continue
                    if all((z/T).is_integral for z in r):continue
                    if not all((2*z/q).is_integral for z in mv(fg,r)):continue
                    r9=(r[0],*r[1:5],*(Z,)*4,r[5]);item={'profile':idx,'vector':enc(r9),'norm':q.to_json(),'height':(ALPHA*x*x/q).to_json()}
                    candidates.append(item)
                    if not compatible(r9,q,old,g):continue
                    step=ALPHA*x/pair(r,fg,direction)
                    if step<=0:continue
                    point=tuple(a+step*b for a,b in zip((O,*(Z,)*5),direction))
                    assert pair(point,fg,point)<0
                    if any(pair(point,fg,w)>0 for w in fw):continue
                    item.update(endpoint_parameter=step.to_json(),endpoint=enc(point),endpoint_norm=pair(point,fg,point).to_json())
                    passed.append(item)
        print({'q':str(q),'profiles':len(profiles),'quadratics':counts['quadratics'],'roots':len(candidates),'compatible':len(passed),'seconds':round(time.monotonic()-started,2)},flush=True)
    passed.sort(key=lambda r:K.from_json(r['endpoint_parameter']))
    out={'status':'COMPLETE_BOUNDED_FINITE_COUPLING_EDGE_SEARCH','source':str(SOURCE.relative_to(ROOT)),
         'face_source':str(FACE.relative_to(ROOT)),'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(SOURCE,FACE)},
         'height_bound':bound,'alpha':ALPHA.to_json(),'edge_indices':edge_ids,'direction':enc(direction),'direction_norm':dn.to_json(),
         'finite_pairing_alphabets':alphabets,'profiles':profiles,'integral_candidates':candidates,'compatible_candidates':passed,'counts':dict(counts),
         'seconds':round(time.monotonic()-started,3),
         'scope':'Complete finite Coxeter coupling profiles on this fixed initial H5 edge, all five necessary primitive norm classes and all integral time coordinates under the stated root-height cap. Residual nonsquare classes are included. The H9 compatibility gates include all 14 current walls. A positive candidate does not by itself certify the first full-group mirror, and a bounded failure is not nonreflectivity.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print({'status':out['status'],'compatible':len(passed),'seconds':out['seconds']},flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--height',type=int,default=10000);a=p.parse_args();run(a.height)
