#!/usr/bin/env python3
"""Independent full A4/scalar completion and relative-norm cut replay.

The producer's residue-indexed pair join is replaced by sequential Schur
energy bounds and a final mod-five Cartesian join. Standard-library Fraction
arithmetic and the separately audited trace ellipsoid enumerate field integers.
"""
from fractions import Fraction as F
from itertools import groupby
import json
from pathlib import Path
import time
from fresh_cubic9_positive_20261013_independent import f,ROOT,Census,Z,O,T,AI,A4,sq
from fresh_cubic9_ellipsoid_audit_20261013 import ellipsoid
SOURCE=ROOT/'results/fresh_h9_joint_arithmetic_20261013.json'
PELL=ROOT/'results/fresh_h9_joint_arithmetic_20261013_pell.json'
OUTPUT=ROOT/'results/fresh_h9_joint_arithmetic_20261013_independent.json'


def candidates(c,x,q):
    energy=f.add(q,f.mul(T,sq(x)));limit=f.mul(T,x)
    weights=tuple(f.mul(f.scale(q,F(1,2)),k) for k in ellipsoid(f.divide(f.scale(energy,8),sq(q)))
        if f.sign(k)>=0 and f.integral(f.mul(f.scale(q,F(1,2)),k)))
    ss={f.scale(sq(z),2):z for z in ellipsoid(f.scale(energy,F(1,2)))
        if f.sign(z)>=0 and f.integral(f.divide(f.scale(z,4),q)) and f.sign(f.sub(limit,f.scale(z,2)))>=0}
    residue=lambda x:tuple(int(t)%5 for t in x)
    last={}
    for w in weights:last.setdefault(residue(f.scale(w,4)),[]).append(w)
    less=lambda a,b:all(f.sign(f.sub(b,a),e)>=0 for e in range(3))
    out=set()
    for a in weights:
        for b in weights:
            ab=f.add(a,b)
            if f.sign(f.sub(limit,ab))<0:continue
            aa=sq(a);bb=sq(b);axb=f.mul(a,b)
            lower=f.scale(f.add(f.add(aa,axb),bb),F(2,3))
            if not less(lower,energy):continue
            for d in weights:
                abd=f.add(ab,d)
                if f.sign(f.sub(limit,abd))<0:continue
                dd=sq(d)
                # Inverse A3 is1/4[[3,2,1],[2,4,2],[1,2,3]].
                lower=f.scale(f.sum_field((f.scale(aa,3),f.scale(bb,4),f.scale(dd,3),
                    f.scale(axb,4),f.scale(f.mul(a,d),2),f.scale(f.mul(b,d),4))),F(1,4))
                if not less(lower,energy):continue
                need=residue(f.neg(f.sum_field((a,f.scale(b,2),f.scale(d,3)))))
                for e in last.get(need,()):
                    if f.sign(f.sub(limit,f.add(abd,e)))<0:continue
                    w=(a,b,d,e);y=f.mv(AI,w);assert all(map(f.integral,y))
                    z=ss.get(f.sub(energy,f.dot(y,w)))
                    if z is None:continue
                    r=(x,Z,Z,Z,Z,*y,z);assert c.norm(r)==q
                    assert all(f.integral(f.divide(f.scale(t,2),q))for t in c.cov(r));out.add(r)
    return out


def root_checks(c,data):
    walls=tuple(f.vecdecode(r)for r in data['all_walls'])
    for r in walls:
        q=c.norm(r);assert q in tuple(f.scale(O,k)for k in (1,2,5,10)) or q==f.scale(f.BETA,2)
        assert all(map(f.integral,r)) and all(f.integral(f.divide(f.scale(t,2),q))for t in c.cov(r))
    return c.check_geometry(data)

def run():
    start=time.monotonic();data=json.loads(SOURCE.read_text());base=json.loads((ROOT/data['source']).read_text())
    previous=json.loads((ROOT/data['previous_source']).read_text());c=Census(base)
    old=tuple(f.vecdecode(r)for r in base['all_walls']);prior=tuple(f.vecdecode(r)for r in previous['all_walls'])
    pescape=f.vecdecode(previous['geometry']['hyperideal_witness']['vector'])
    assert all(t==Z for t in pescape[1:5]) and f.sign(c.norm(pescape))>0
    assert all(f.sign(f.dot(pescape,c.cov(r)))<=0 for r in prior)
    pool=tuple(f.scale(O,k)for k in (1,2,5,10));controls=c.controls(data['height_bound'],pool)
    saved={(f.decode(row['height']),f.decode(row['norm']),f.decode(row['time'])):row for row in data['levels']}
    assert set(controls)==set(saved) and len(controls)==data['control_count']
    selected=list(old);keys={f.ray_key(r)for r in old};dictionary={};checks=[]
    for h,packet in groupby(controls,key=lambda row:row[0]):
        accepted=[]
        for _,q,x in packet:
            roots=candidates(c,x,q);row=saved[h,q,x];assert len(roots)==row['raw_candidates_after_simple_prefix_bounds']
            allkept=0;added=0
            for r in sorted(roots):
                if f.ray_key(r)in {f.ray_key(w)for w in old} or not c.compatible(r,old):continue
                allkept+=1;dictionary[f.ray_key(r)]=(q,r)
                if f.ray_key(r)in keys or not c.compatible(r,selected):continue
                keys.add(f.ray_key(r));accepted.append(r);added+=1
            assert allkept==row['original13_compatible'] and added==row['accepted']
            checks.append({'height':f.encoded(h),'norm':f.encoded(q),'raw_candidates':len(roots),'original13_compatible':allkept,'accepted':added})
        assert all(c.compatible(r,accepted[:i])for i,r in enumerate(accepted));selected.extend(accepted)
        print(json.dumps({'controls':len(checks),'walls':len(selected)}),flush=True)
    walls=tuple(f.vecdecode(r)for r in data['all_walls'])
    assert {f.ray_key(r)for r in selected}=={f.ray_key(r)for r in walls}
    recorded={f.ray_key(f.vecdecode(r['vector'])):r for r in data['compatible_dictionary']}
    assert set(recorded)==set(dictionary)
    for key,(q,r) in dictionary.items():
        row=recorded[key];assert f.decode(row['norm'])==q and f.vecdecode(row['vector'])==r
        assert row['compatible_with_previous22']==c.compatible(r,prior)
        assert f.decode(row['previous_escape_pairing'])==f.dot(r,c.cov(pescape))
    full_geometry=root_checks(c,data)
    face=data['forced_h5_face'];fc=Census({'alpha':base['alpha'],'blocks':0})
    facewalls=tuple((r[0],*r[5:])for r in walls if all(t==Z for t in r[1:5]))
    assert facewalls==tuple(f.vecdecode(r)for r in face['all_walls'])
    face_geometry=root_checks(fc,face)
    pell=json.loads(PELL.read_text());assert pell['source']==str(SOURCE.relative_to(ROOT))
    added=tuple(f.vecdecode(row['vector'])for row in pell['new_roots'])
    assert tuple(f.vecdecode(row)for row in pell['all_walls'])==(*walls,*added)
    x=f.decode(pell['time']);t=f.decode(pell['multiple'])
    assert f.sub(f.scale(sq(t),30),f.mul(T,sq(x)))==f.scale(O,5)
    oldv=f.vecdecode(face['geometry']['hyperideal_witness']['vector']);embedded=(oldv[0],Z,Z,Z,Z,*oldv[1:])
    cuts=tuple(f.dot(r,c.cov(embedded))for r in added)
    assert cuts==f.vecdecode(pell['old_forced_face_escape_pairings']) and any(f.sign(p)>0 for p in cuts)
    pell_geometry=root_checks(c,pell);pell_face_geometry=root_checks(fc,pell['forced_h5_face'])
    out={'status':'PASS_COMPLETE_FACE_ENUMERATION_AND_PELL_CUTS','source':str(SOURCE.relative_to(ROOT)),
        'arithmetic':'Independent Fraction polynomial arithmetic; trace ellipsoid alphabets; sequential Schur bounds and final Cartesian mod5 join.',
        'controls':checks,'complete_candidate_census':True,'original13_compatible_dictionary':len(dictionary),
        'old22_escape_cutters':sum(f.sign(f.decode(row['previous_escape_pairing']))>0 for row in data['compatible_dictionary']),
        'h9_geometry':full_geometry,'h5_face_geometry':face_geometry,
        'pell_source':str(PELL.relative_to(ROOT)),'pell_equation_verified':True,'pell_old_escape_pairings':f.vecencode(cuts),
        'pell_h9_geometry':pell_geometry,'pell_h5_face_geometry':pell_face_geometry,
        'seconds':round(time.monotonic()-start,3),
        'scope':'The full time+A4+A1 complement is completely enumerated through1000 under the original13-wall constraints. Two additional targeted norm-five mirrors are verified independently without a height-completeness claim.'}
    OUTPUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':out['status'],'dictionary':len(dictionary),'seconds':out['seconds']}));return out
if __name__=='__main__':run()
