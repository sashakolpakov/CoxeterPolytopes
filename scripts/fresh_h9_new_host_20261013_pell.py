#!/usr/bin/env python3
"""One exact relative-unit step at a new host's open initial edge."""
import hashlib,json,time
from fresh_h9_new_host_20261013_edge import *
from fresh_h9_quartic_mixed_20261013_target_independent import field_sqrt

SOURCE_EDGE=OUT
OUT=ROOT/'results/fresh_h9_new_host_20261013_pell.json'

def run():
    start=time.monotonic();profiles=json.loads(SOURCE_EDGE.read_text());data=json.loads(SOURCE.read_text());fd=json.loads(FACE.read_text())
    g,fg=Host(2).g,Host(1).g;walls=[tuple(map(K.from_json,r))for r in data['all_walls']]
    fw=[tuple(map(K.from_json,r))for r in fd['all_walls']];direction=tuple(map(K.from_json,profiles['direction']))
    y=K(-8561928901343,-4633682531136,14616126887520,7910190870768)
    s=K(2736975019280,1481240206632,-4672308615352,-2528635201308)
    D=3*ALPHA;assert y*y-D*s*s==1 and y>1
    base=tuple(-z/6 for z in direction);jb=(-SQRT2/2,*(Z,)*5)
    assert pair(base,fg,base)==K(1)/6 and pair(jb,fg,jb)==-D/6 and pair(base,fg,jb)==0
    candidates=[];passed=[];squares=0
    for idx,p in enumerate(profiles['profiles']):
        q=K.from_json(p['norm']);r0=tuple(map(K.from_json,p['fixed_component']));res=K.from_json(p['residual_norm'])
        factor=square_root(6*res)
        if factor is None:continue
        squares+=1
        for power in(-1,1):
            translated=tuple(y*z+power*s*t for z,t in zip(base,jb))
            for sign in(-1,1):
                r=tuple(a+sign*factor*b for a,b in zip(r0,translated));assert pair(r,fg,r)==q
                if r[0]<=0 or not all(z.is_integral for z in r):continue
                if all((z/T).is_integral for z in r)or not all((2*z/q).is_integral for z in mv(fg,r)):continue
                r9=(r[0],*r[1:5],*(Z,)*4,r[5]);row={'profile':idx,'power':power,'sign':sign,'vector':enc(r9),'norm':q.to_json(),'height':(ALPHA*r[0]*r[0]/q).to_json()};candidates.append(row)
                if not compatible(r9,q,walls,g):continue
                step=ALPHA*r[0]/pair(r,fg,direction)
                if step<=0:continue
                point=tuple(a+step*b for a,b in zip((O,*(Z,)*5),direction));assert pair(point,fg,point)==-ALPHA/(y*y)<0
                if any(pair(point,fg,w)>0 for w in fw):continue
                row.update(endpoint_parameter=step.to_json(),endpoint=enc(point),endpoint_norm=pair(point,fg,point).to_json());passed.append(row)
    passed.sort(key=lambda r:K.from_json(r['endpoint_parameter']))
    out={'status':'EXACT_RELATIVE_UNIT_PROPOSALS','source':str(SOURCE_EDGE.relative_to(ROOT)),
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(SOURCE_EDGE,SOURCE,FACE)},
        'discriminant':D.to_json(),'unit_y':y.to_json(),'unit_s':s.to_json(),'seed_norm':(K(1)/6).to_json(),
        'square_ratio_profiles':squares,'integral_candidates':candidates,'compatible_candidates':passed,
        'scope':'Only the two powers -1,+1 of the displayed exact relative norm-one unit, both signs, and the finite-coupling profiles whose residual/seed ratio is square. No first-root or full relative-norm-orbit claim.',
        'seconds':round(time.monotonic()-start,3)}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print({'status':out['status'],'profiles':squares,'integral':len(candidates),'compatible':len(passed),'seconds':out['seconds']})

if __name__=='__main__':run()
