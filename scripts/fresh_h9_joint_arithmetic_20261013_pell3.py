#!/usr/bin/env python3
"""Third certified Pell orbit, after a finite spatial-direction perturbation."""
import json
from fresh_cubic9_positive_20261013 import *
SOURCE=ROOT/'results/fresh_h9_joint_arithmetic_20261013_pell2.json'
PROPOSAL=ROOT/'results/fresh_h9_joint_arithmetic_20261013_pell3_proposal.json'
OUTPUT=ROOT/'results/fresh_h9_joint_arithmetic_20261013_pell3.json'

def witness(v,walls,g):
    ps=tuple(pair(v,g,r)for r in walls);q=pair(v,g,v)
    assert q.sign()>0 and all(p.sign()<=0 for p in ps)
    return {'status':'EXACT_ESCAPE','hyperideal_witness':{'vector':enc(v),'norm':q.to_json(),
        'wall_pairings':enc(ps),'active_facets':[i for i,p in enumerate(ps)if p==Z]},'vertices':[],
        'note':'A small rational escape is checked against every wall exactly.'}

def run():
    data=json.loads(SOURCE.read_text());proposal=json.loads(PROPOSAL.read_text());g=Host(T,2).G
    walls=[tuple(map(K.from_json,r))for r in data['all_walls']]
    x=K.from_json(proposal['time']);t=K.from_json(proposal['multiple'])
    d=tuple(map(K.from_json,proposal['spatial_vector']));delta=K.from_json(proposal['spatial_norm'])
    assert delta*t*t-T*x*x==5
    oldv=tuple(map(K.from_json,data['forced_h5_face']['geometry']['hyperideal_witness']['vector']))
    oldv=(oldv[0],Z,Z,Z,Z,*oldv[1:]);rows=[];cuts=[]
    for direction in (d,(*reversed(d[:4]),d[4])):
        r=(x,Z,Z,Z,Z,*(t*j for j in direction));q=pair(r,g,r)
        assert q==5 and all(s.is_integral for s in r)
        assert all((2*v/q).is_integral for v in mv(g,r))
        cc=[mv(g,w)for w in walls];nn=[dot(w,c)for w,c in zip(walls,cc)]
        assert compatible(r,q,cc,nn);walls.append(r)
        rows.append({'vector':enc(r),'norm':q.to_json(),'height':(T*x*x/q).to_json()});cuts.append(pair(oldv,g,r))
    assert cuts[0].sign()>0
    inside=tuple(map(K.from_json,data['strict_interior']))
    while any(pair(inside,g,r).sign()>=0 for r in walls):inside=(O,*(v/2 for v in inside[1:]))
    fg=Host(T,0).G;fw=tuple((r[0],*r[5:])for r in walls if all(s==Z for s in r[1:5]));fi=(inside[0],*inside[5:])
    # Rationally move the surviving end-weight/scalar escape into its
    # feasible face. This avoids carrying a large algebraic vertex in the
    # final certificate; its exact norm is1903957/1000000-T>0.
    v=(O,Z,Z,Z,Z,K(F(499,500)),K(F(1497,2000)),K(F(499,1000)),K(F(499,2000)),K(F(287,500)))
    assert pair(v,g,v)==K(F(1903957,1000000))-T
    fv=(v[0],*v[5:])
    out={'status':'EXACT_THIRD_RELATIVE_NORM_TWO_ROOT_ORBIT','source':str(SOURCE.relative_to(ROOT)),
        'proposal_source':str(PROPOSAL.relative_to(ROOT)),'target_kind':'forced_h5_face',
        'binary_root_equation':'delta*t^2-T*x^2=5','time':x.to_json(),'multiple':t.to_json(),
        'spatial_norm':delta.to_json(),'spatial_direction':enc(d),'new_roots':rows,'old_escape_pairings':enc(cuts),
        'all_walls':encode_matrix(walls),'wall_count':len(walls),'strict_interior':enc(inside),
        'finite_angle_counts':finite_angle_counts(walls,g),'geometry':witness(v,walls,g),
        'forced_h5_face':{'gram':encode_matrix(fg),'all_walls':encode_matrix(fw),'wall_count':len(fw),
            'strict_interior':enc(fi),'finite_angle_counts':finite_angle_counts(fw,fg),'geometry':witness(fv,fw,fg)},
        'scope':'Two explicit primitive norm-five mirrors in a perturbation of the third target direction. No completeness or periodic-isometry claim. The retained exact escape proves these203 walls still have infinite volume.'}
    OUTPUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'walls':len(walls),'face_walls':len(fw),'geometry':out['geometry']['status'],'face_geometry':out['forced_h5_face']['geometry']['status']}),flush=True);return out
if __name__=='__main__':run()
