#!/usr/bin/env python3
"""Independent unpruned integer-box and Cartesian root census for cubic-nine."""
import argparse
from functools import lru_cache
from itertools import product,groupby
import json
from pathlib import Path
import fresh_cubic9_dihedral_independent_audit_20261013 as f

BOUNDS=((f.F(1,12),f.F(7,6),f.F(1,6)),(f.F(1,4),f.F(1,6),f.F(5,12)),(f.F(1,6),f.F(5,12),f.F(1,4)))

def ceil(x):return -((-x.numerator)//x.denominator)

@lru_cache(None)
def alphabet(caps):
    if any(f.sign(c,e)<0 for e,c in enumerate(caps)):return ()
    ms=[]
    for e,c in enumerate(caps):
        m=0 if c==f.ZERO else 1
        while f.sign(f.sub(f.scale(f.ONE,m*m),c),e)<0:m*=2
        ms.append(m)
    bounds=[ceil(sum(b*m for b,m in zip(row,ms))) for row in BOUNDS]
    out=[]
    for coefficients in product(*(range(-b,b+1) for b in bounds)):
        x=tuple(map(f.F,coefficients));square=f.mul(x,x)
        if all(f.sign(f.sub(c,square),e)>=0 for e,c in enumerate(caps)):out.append(x)
    return tuple(out)

@lru_cache(None)
def block_list(target,q):
    values=[w for w in alphabet((f.scale(target,2),)*3) if f.sign(w)>=0 and f.integral(f.divide(f.scale(w,2),q))]
    result=[]
    for u,v in product(values,repeat=2):
        if f.sign(f.sub(v,u))<0:continue
        y=(f.divide(f.add(f.scale(u,2),f.mul(f.C,v)),f.BETA),f.divide(f.add(f.mul(f.C,u),f.scale(v,2)),f.BETA))
        if not all(f.integral(t) for t in y):continue
        energy=f.dot(y,(u,v))
        if all(f.sign(f.sub(target,energy),e)>=0 for e in range(3)):result.append((energy,y))
    return tuple(result)

def candidates(d,x,q):
    alpha=f.decode(d['alpha']);target=f.add(q,f.mul(alpha,f.mul(x,x)))
    blocks=block_list(target,q)
    scalars=[z for z in alphabet((target,)*3) if f.sign(z)>=0 and f.integral(f.divide(f.scale(z,2),q))] if d['scalar'] else [f.ZERO]
    squares={f.mul(z,z):z for z in scalars};out=set()
    for entries in product(blocks,repeat=d['blocks']):
        remainder=f.sub(target,f.sum_field(e for e,y in entries))
        if remainder not in squares:continue
        z=squares[remainder]
        out.add((x,*(a for e,y in entries for a in y),*((z,) if d['scalar'] else ())))
    return out,len(blocks)

def verify(d):
    alpha=f.decode(d['alpha']);height=f.decode(d['height_bound']);norms=f.vecdecode(d['norm_pool'])
    g=f.make_gram(alpha,d['blocks'],d['scalar']);walls=tuple(f.vecdecode(w) for w in d['all_walls'])
    controls=set()
    for q in norms:
        caps=(f.divide(f.mul(height,q),alpha),f.neg(f.divide(q,alpha)),f.neg(f.divide(q,alpha)))
        for x in alphabet(caps):
            if f.sign(x)>0 and f.integral(f.divide(f.scale(f.mul(alpha,x),2),q)):
                controls.add((f.divide(f.mul(alpha,f.mul(x,x)),q),q,x))
    saved={(f.decode(r['height']),f.decode(r['norm']),f.decode(r['time'])) for r in d['levels']}
    assert controls==saved and len(saved)==len(d['levels'])
    full=f.sign(f.sub(height,f.scale(f.ONE,30)))<=0
    out={'controls':len(controls),'complete_time_norm_census':True,'full_spatial_census_repeated':full}
    if not full:return out
    zero=[]
    for q in norms:
        roots,_=candidates(d,f.ZERO,q);zero.extend((q,r) for r in roots)
    assert set(zero)=={(f.decode(row['norm']),f.vecdecode(row['vector'])) for row in d['zero_height_dominant_census']}
    selected=list(walls[:d['dimension']]);keys={f.ray_key(r) for r in selected};counts=[]
    for _,packet in groupby(d['levels'],key=lambda row:row['height']):
        accepted=[];cov=tuple(f.mv(g,r) for r in selected)
        for row in packet:
            q=f.decode(row['norm']);x=f.decode(row['time']);roots,bc=candidates(d,x,q)
            assert len(roots)==row['candidate_roots'] and bc==row['block_vectors']
            added=[]
            for r in roots:
                if any(f.sign(f.dot(r,n))>0 for n in cov):continue
                key=f.ray_key(r)
                if key in keys:continue
                keys.add(key);added.append(r)
            assert len(added)==row['accepted'];accepted.extend(added);counts.append(len(roots))
        selected.extend(accepted)
    assert set(selected)==set(walls)
    return {**out,'zero_dominant_count':len(zero),'candidate_counts':counts,'walls':len(walls)}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('inputs',type=Path,nargs='+');p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    f.binary_census();out=[]
    for path in a.inputs:
        row=verify(json.loads(path.read_text()));row['source']=str(path);out.append(row);print(json.dumps(row),flush=True)
    a.output.write_text(json.dumps({'status':'PASS','arithmetic':'Independent stdlib Fraction; unpruned coefficient boxes and Cartesian spatial products.','cases':out},indent=2)+'\n')
