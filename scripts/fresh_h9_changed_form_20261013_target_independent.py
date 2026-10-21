#!/usr/bin/env python3
"""Independent saturated two- or three-covector/B/scalar targeted census."""
import argparse
import json
import time
from itertools import product
from fresh_h9_changed_form_20261013_independent import f,Census,Z,O,T,F,ROOT,sq,AI,ellipsoid,Path

def verify(path):
    data=json.loads(path.read_text());base=json.loads((ROOT/data['source']).read_text());c=Census(base)
    old=tuple(f.vecdecode(r)for r in base['all_walls']);target=f.vecdecode(base['geometry']['hyperideal_witness']['vector'])
    tc=c.cov(target);assert all(target[i]==Z for i in(1,2,9))
    supports=data.get('a6_covector_supports',[data.get('a6_covector_support',[4,5])])
    assert all(len(s)in(2,3)and s==sorted(set(s))and 0<=s[0]and s[-1]==5 for s in supports)
    assert f.vecdecode(data['norm_pool'])==c.pool
    controls=c.controls(data['height_bound']);saved={(f.decode(r['height']),f.decode(r['norm']),f.decode(r['time'])):r for r in data['levels']}
    assert set(controls)==set(saved)and len(controls)==data['control_count'];dictionary={};checks=[]
    for h,q,x in controls:
        energy=f.add(q,f.mul(c.alpha,sq(x)));less=lambda e:all(f.sign(f.sub(energy,e),i)>=0 for i in range(3))
        weights=tuple(f.mul(f.scale(q,F(1,2)),k)for k in ellipsoid(f.divide(f.scale(energy,8),sq(q)))
            if f.sign(k)>=0 and f.integral(f.mul(f.scale(q,F(1,2)),k)))
        bs=[]
        for u,v in product(weights,repeat=2):
            if f.sign(f.sub(v,u))<0:continue
            y=(f.divide(f.add(f.scale(u,2),f.mul(T,v)),f.BETA),f.divide(f.add(f.mul(T,u),f.scale(v,2)),f.BETA))
            e=f.dot(y,(u,v))
            if all(map(f.integral,y))and less(e):bs.append((e,y))
        ss={f.scale(sq(z),2):z for z in ellipsoid(f.scale(energy,F(1,2)))
            if f.sign(z)>=0 and f.integral(f.divide(f.scale(z,4),q))}
        counts=dict(weights=len(weights),block_vectors=len(bs),scalar_vectors=len(ss),target_cutting_a6_vectors=0,
                    raw_target_cutters=0,imprimitive=0,positive_pair_rejected=0,noncoxeter_rejected=0,admissible_target_cutters=0)
        for support in supports:
            residue=lambda v:tuple(int(t)%7 for t in v)
            byresidue={}
            for v in weights:byresidue.setdefault(residue(f.scale(v,support[-1]+1)),[]).append(v)
            for prefix in product(weights,repeat=len(support)-1):
                need=residue(f.neg(f.sum_field(f.scale(w,i+1)for i,w in zip(support[:-1],prefix))))
                for v in byresidue.get(need,()):
                    w=[Z]*6
                    for i,s in zip(support,(*prefix,v)):w[i]=s
                    w=tuple(w);y=f.mv(AI,w);assert all(map(f.integral,y))
                    r0=(x,Z,Z,*y,Z)
                    if f.sign(f.dot(r0,tc))<=0:continue
                    ea=f.dot(y,w)
                    if not less(ea):continue
                    counts['target_cutting_a6_vectors']+=1
                    for eb,yb in bs:
                        z=ss.get(f.sub(f.sub(energy,ea),eb))
                        if z is None:continue
                        r=(x,*yb,*y,z);counts['raw_target_cutters']+=1
                        assert c.norm(r)==q and all(f.integral(f.divide(f.scale(t,2),q))for t in c.cov(r))
                        if all(f.integral(f.scale(t,F(1,2)))for t in r):counts['imprimitive']+=1;continue
                        if any(f.sign(f.dot(r,c.cov(w)))>0 for w in old):counts['positive_pair_rejected']+=1;continue
                        if not c.compatible(r,old):counts['noncoxeter_rejected']+=1;continue
                        counts['admissible_target_cutters']+=1
                        for rr in(r,(r[0],r[1],r[2],*reversed(r[3:9]),r[9])):
                            assert c.compatible(rr,old);dictionary[f.ray_key(rr)]=(h,q,rr)
        row=saved[h,q,x];assert all(row[k]==v for k,v in counts.items()),(h,q,counts,row)
        checks.append({'height':f.encoded(h),'norm':f.encoded(q),**counts})
    rows=data['compatible_dictionary'];recorded={f.ray_key(f.vecdecode(r['vector'])):r for r in rows}
    assert set(dictionary)==set(recorded)
    for key,(h,q,r)in dictionary.items():
        row=recorded[key];assert f.decode(row['height'])==h and f.decode(row['norm'])==q and f.vecdecode(row['vector'])==r
        assert f.decode(row['old_escape_pairing'])==f.dot(r,tc)
    walls=list(old);chosen=[];rejected=[];last=Z
    for i,row in enumerate(rows):
        h=f.decode(row['height']);assert f.sign(f.sub(h,last))>=0;last=h;r=f.vecdecode(row['vector'])
        if c.compatible(r,walls):chosen.append(i);walls.append(r)
        else:rejected.append(i)
    assert chosen==data['selected_indices']and rejected==data['rejected_indices']
    assert tuple(walls)==tuple(f.vecdecode(r)for r in data['all_walls'])
    return {'source':str(path.relative_to(ROOT)),'status':'PASS_COMPLETE_TARGETED_CENSUS','controls':checks,
        'supports':supports,'dictionary_size':len(dictionary),'geometry':c.check_geometry(data)}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('sources',nargs='+',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    start=time.monotonic();rows=[verify(path.resolve())for path in a.sources]
    out={'status':'PASS','arithmetic':'Independent Fraction arithmetic, trace-ellipsoid alphabets, residue-filtered A6 covectors, direct B/scalar energy matching.',
        'cases':rows,'seconds':round(time.monotonic()-start,3)}
    a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','cases':len(rows),'seconds':out['seconds']}))
