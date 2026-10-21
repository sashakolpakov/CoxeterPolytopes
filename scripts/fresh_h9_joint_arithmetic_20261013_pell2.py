#!/usr/bin/env python3
"""An exact second relative-norm orbit cutting the end-weight/scalar gap."""
import json
from fresh_cubic9_positive_20261013 import *
SOURCE=ROOT/'results/fresh_h9_joint_arithmetic_20261013_pell.json'
OUTPUT=ROOT/'results/fresh_h9_joint_arithmetic_20261013_pell2.json'

def run():
    data=json.loads(SOURCE.read_text());h=Host(T,2);g=h.G
    walls=[tuple(map(K.from_json,r))for r in data['all_walls']]
    a=K(134319945814822091128381470365291574421001885898464615808428009390,
        474429987812617568624914142911519830010257575353034711272821386375,
        252438923813102754415739776663114654445420554320254413455694269015)
    t=K(19203896737533355913638233894574132022592195883486543641894083483,
        67829870239102909280966695754218190598032998643604019834208156260,
        36091520109189350016832053214072272266374679841132385656908886855)
    x=a/T;k=10+5*T-5*T*T;delta=120-50*T
    assert 20+2*k*k==delta and delta*t*t-T*x*x==5
    oldv=tuple(map(K.from_json,data['geometry']['hyperideal_witness']['vector']))
    rows=[];cuts=[]
    for direction in ((4,3,2,1),(1,2,3,4)):
        r=(x,Z,Z,Z,Z,*(t*j for j in direction),k*t);q=pair(r,g,r)
        assert q==5 and all(s.is_integral for s in r)
        assert all((2*v/q).is_integral for v in mv(g,r))
        cc=[mv(g,w)for w in walls];nn=[dot(w,c)for w,c in zip(walls,cc)]
        assert compatible(r,q,cc,nn);walls.append(r)
        rows.append({'vector':enc(r),'norm':q.to_json(),'height':(T*x*x/q).to_json()})
        cuts.append(pair(oldv,g,r))
    assert cuts[0].sign()>0
    inside=tuple(map(K.from_json,data['strict_interior']))
    while any(pair(inside,g,r).sign()>=0 for r in walls):inside=(O,*(v/2 for v in inside[1:]))
    fg=Host(T,0).G;fw=tuple((r[0],*r[5:])for r in walls if all(s==Z for s in r[1:5]));fi=(inside[0],*inside[5:])
    out={'status':'EXACT_SECOND_RELATIVE_NORM_TWO_ROOT_ORBIT','source':str(SOURCE.relative_to(ROOT)),
        'relative_norm_equation':'Norm(T*x+t*sqrt(T*(120-50T)))=-5T',
        'binary_root_equation':'(120-50T)*t^2-T*x^2=5','time':x.to_json(),'multiple':t.to_json(),
        'scalar_ratio':k.to_json(),'spatial_norm':delta.to_json(),
        'new_roots':rows,'old_escape_pairings':enc(cuts),
        'all_walls':encode_matrix(walls),'wall_count':len(walls),'strict_interior':enc(inside),
        'finite_angle_counts':finite_angle_counts(walls,g),'geometry':geometry(walls,g,9,inside=inside),
        'forced_h5_face':{'gram':encode_matrix(fg),'all_walls':encode_matrix(fw),'wall_count':len(fw),
            'strict_interior':enc(fi),'finite_angle_counts':finite_angle_counts(fw,fg),'geometry':geometry(fw,fg,5,inside=fi)},
        'scope':'Two explicit primitive norm-five mirrors at a new spatial direction. No minimal-height or complete Pell-orbit claim. Every pair is checked for an allowed Coxeter angle.'}
    OUTPUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'walls':len(walls),'face_walls':len(fw),'geometry':out['geometry']['status'],'face_geometry':out['forced_h5_face']['geometry']['status']}),flush=True)
    return out
if __name__=='__main__':run()
