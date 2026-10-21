#!/usr/bin/env python3
"""Fraction-only replay of later targeted relative-norm orbits."""
import argparse
import json
import time
from fresh_h9_joint_arithmetic_20261013_independent import f,ROOT,Census,Z,O,T,sq,root_checks
from fresh_cubic9_positive_20261013_independent import AI

def scout_replay(c):
    path=ROOT/'results/fresh_h9_joint_arithmetic_20261013_pell_scout.json'
    data=json.loads(path.read_text());prior=json.loads((ROOT/data['source']).read_text())
    fc=Census({'alpha':f.encoded(T),'blocks':0})
    walls=tuple(f.vecdecode(r)for r in prior['forced_h5_face']['all_walls'])
    v=f.vecdecode(data['target']);assert v==f.vecdecode(prior['forced_h5_face']['geometry']['hyperideal_witness']['vector'])
    d=tuple(f.divide(f.scale(t,30),f.sub(sq(T),O))for t in v[1:])
    assert d==f.vecdecode(data['base_spatial_vector'])
    directions=[(*tuple(f.scale(AI[j][i],5)for j in range(4)),Z)for i in range(4)]
    directions.append((Z,Z,Z,Z,f.scale(O,5)))
    unit=f.sub(f.add(f.scale(O,2),T),sq(T));up=O;found={};tests=0
    for j in range(1,5):
        up=f.mul(up,unit)
        for k,u in enumerate(directions):
            for sign in (-1,1):
                tests+=1;e=tuple(f.add(a,f.scale(f.mul(up,b),sign))for a,b in zip(d,u))
                r=(Z,*e);delta=fc.norm(r);b=f.dot(r,fc.cov(v))
                if f.sign(b)<=0 or f.sign(f.sub(sq(b),f.mul(T,delta)))<=0:continue
                okay=True
                for w in walls:
                    p=f.dot(r,fc.cov(w))
                    if w[0]==Z:
                        if f.sign(p)>0:okay=False;break
                    elif f.sign(p)>0 and f.sign(f.sub(f.mul(f.mul(T,delta),sq(w[0])),sq(p)))<=0:
                        okay=False;break
                if okay:found[j,k,sign]=(e,delta)
    assert tests==data['tested']==40 and len(found)==data['passed']==36
    assert set(found)=={(r['unit_power'],r['direction'],r['sign'])for r in data['directions']}
    for row in data['directions']:
        e,delta=found[row['unit_power'],row['direction'],row['sign']]
        assert e==f.vecdecode(row['spatial_vector']) and delta==f.decode(row['spatial_norm'])
    return {'exact_directions_tested':tests,'geometrically_feasible':len(found),
            'scope':'Only exact continuous direction inequalities; external PARI Hilbert observations are not replayed or used.'}

def run(number=2):
    start=time.monotonic()
    path=ROOT/f'results/fresh_h9_joint_arithmetic_20261013_pell{number}.json'
    data=json.loads(path.read_text());prior=json.loads((ROOT/data['source']).read_text())
    c=Census({'alpha':f.encoded(T),'blocks':2});fc=Census({'alpha':f.encoded(T),'blocks':0})
    old=tuple(f.vecdecode(r)for r in prior['all_walls'])
    added=tuple(f.vecdecode(row['vector'])for row in data['new_roots'])
    walls=tuple(f.vecdecode(r)for r in data['all_walls'])
    assert walls==(*old,*added)
    x=f.decode(data['time']);t=f.decode(data['multiple']);delta=f.decode(data['spatial_norm'])
    assert f.sub(f.mul(delta,sq(t)),f.mul(T,sq(x)))==f.scale(O,5)
    assert all(r[0]==x and c.norm(r)==f.scale(O,5)for r in added)
    if number==2:
        k=f.sub(f.add(f.scale(O,10),f.scale(T,5)),f.scale(sq(T),5))
        d=tuple(f.scale(O,i)for i in (4,3,2,1))+(k,)
    else:d=f.vecdecode(data['spatial_direction'])
    assert added[0]==(x,Z,Z,Z,Z,*(f.mul(t,s)for s in d))
    assert added[1]==(x,Z,Z,Z,Z,*(f.mul(t,s)for s in (*reversed(d[:4]),d[4])))
    for i,row in enumerate(data['new_roots']):
        assert f.decode(row['height'])==f.divide(f.mul(T,sq(x)),f.scale(O,5))
    escape=f.vecdecode(prior['geometry']['hyperideal_witness']['vector'])
    if data.get('target_kind')=='forced_h5_face':
        v=f.vecdecode(prior['forced_h5_face']['geometry']['hyperideal_witness']['vector'])
        escape=(v[0],Z,Z,Z,Z,*v[1:])
    cuts=tuple(f.dot(r,c.cov(escape))for r in added)
    assert cuts==f.vecdecode(data['old_escape_pairings']) and any(f.sign(p)>0 for p in cuts)
    h9=root_checks(c,data);face=data['forced_h5_face']
    projected=tuple((r[0],*r[5:])for r in walls if all(s==Z for s in r[1:5]))
    assert projected==tuple(f.vecdecode(r)for r in face['all_walls'])
    h5=root_checks(fc,face)
    result={'status':'PASS_TARGETED_RELATIVE_NORM_ORBIT','source':str(path.relative_to(ROOT)),
        'relative_norm_identity':True,'all_roots_integral_reflective':True,
        'old_escape_pairings':f.vecencode(cuts),'h9':h9,'forced_h5_face':h5,
        'seconds':round(time.monotonic()-start,3),
        'scope':'Independent Fraction field replay of roots, all pair Coxeter angles, strict centers, and exact escape witnesses. No chronological or norm-solver completeness claim.'}
    if number==3:result['direction_scout']=scout_replay(c)
    output=path.with_name(path.stem+'_independent.json')
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'walls':len(walls),'seconds':result['seconds']}))
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--round',type=int,default=2);args=p.parse_args();run(args.round)
