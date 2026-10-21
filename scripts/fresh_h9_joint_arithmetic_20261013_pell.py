#!/usr/bin/env python3
"""Exact first relative-norm cut of the fully freed A4+A1 prefix."""
import json
from fresh_cubic9_positive_20261013 import *
SOURCE=ROOT/'results/fresh_h9_joint_arithmetic_20261013.json'
OUTPUT=ROOT/'results/fresh_h9_joint_arithmetic_20261013_pell.json'

def run():
    data=json.loads(SOURCE.read_text());h=Host(T,2);g=h.G
    walls=[tuple(map(K.from_json,r))for r in data['all_walls']]
    x=8715+30785*T+16380*T*T;t=2181+7705*T+4100*T*T
    rows=[]
    for direction in ((2,4,6,3),(3,6,4,2)):
        r=(x,Z,Z,Z,Z,*(t*j for j in direction),Z);q=pair(r,g,r)
        assert q==5 and all((2*v/q).is_integral for v in mv(g,r))
        cc=[mv(g,w)for w in walls];nn=[dot(w,c)for w,c in zip(walls,cc)]
        assert compatible(r,q,cc,nn);walls.append(r)
        rows.append({'vector':enc(r),'norm':q.to_json(),'height':(T*x*x/q).to_json()})
    inside=tuple(map(K.from_json,data['strict_interior']))
    while any(pair(inside,g,r).sign()>=0 for r in walls):inside=(O,*(v/2 for v in inside[1:]))
    oldface=data['forced_h5_face'];oldv=tuple(map(K.from_json,oldface['geometry']['hyperideal_witness']['vector']))
    embedded=(oldv[0],Z,Z,Z,Z,*oldv[1:])
    cuts=[pair(embedded,g,tuple(map(K.from_json,row['vector'])))for row in rows]
    fg=Host(T,0).G;fw=tuple((r[0],*r[5:])for r in walls if all(t==Z for t in r[1:5]));fi=(inside[0],*inside[5:])
    out={'status':'EXACT_RELATIVE_NORM_TWO_ROOT_ORBIT','source':str(SOURCE.relative_to(ROOT)),
        'relative_norm_equation':'Norm_{K(sqrt(30T))/K}(T*x+t*sqrt(30T))=-5T',
        'binary_root_equation':'30*t^2-T*x^2=5','time':x.to_json(),'multiple':t.to_json(),
        'new_roots':rows,'old_forced_face_escape_pairings':enc(cuts),
        'all_walls':encode_matrix(walls),'wall_count':len(walls),'strict_interior':enc(inside),
        'finite_angle_counts':finite_angle_counts(walls,g),'geometry':geometry(walls,g,9,inside=inside),
        'forced_h5_face':{'gram':encode_matrix(fg),'all_walls':encode_matrix(fw),'wall_count':len(fw),
            'strict_interior':enc(fi),'finite_angle_counts':finite_angle_counts(fw,fg),'geometry':geometry(fw,fg,5,inside=fi)},
        'scope':'Two explicit primitive integral norm-five mirrors discovered by a relative norm solver, verified against every prior wall. Their enormous height is not a claim of chronological completeness.'}
    OUTPUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps({'walls':len(walls),'face_walls':len(fw),'geometry':out['geometry']['status'],'face_geometry':out['forced_h5_face']['geometry']['status']}));return out
if __name__=='__main__':run()
