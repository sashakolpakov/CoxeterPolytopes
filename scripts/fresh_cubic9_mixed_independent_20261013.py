#!/usr/bin/env python3
"""Independent Fraction audit with complete bounded root reenumeration.

No producer or project arithmetic imports. A4 roots use a mod-5 Cartesian
join rather than the producer's prefix-energy recursion.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache,cmp_to_key
from itertools import product,groupby
from math import isqrt
import json
from pathlib import Path
import time
import fresh_cubic9_dihedral_independent_audit_20261013 as f
ROOT=Path(__file__).resolve().parents[1]
Z=f.ZERO;O=f.ONE;T=f.THETA;B=f.BETA;C=f.C;D=f.D
sq=lambda x:f.mul(x,x)
POOL=tuple(f.scale(b,n) for b in (O,B) for n in (2,4,10,20))
ALLOWED={Z,O,f.TWO,f.scale(O,3),sq(T),f.add(f.TWO,T)}
A4=tuple(tuple(f.scale(O,2 if i==j else -1 if abs(i-j)==1 else 0) for j in range(4)) for i in range(4))
AI=tuple(tuple(f.scale(O,F(min(i,j)*(5-max(i,j)),5)) for j in range(1,5)) for i in range(1,5))
G=tuple(tuple(f.neg(f.scale(T,2)) if i==j==0 else f.B[(i-1)%2][(j-1)%2]
    if 1<=i<=4 and 1<=j<=4 and (i-1)//2==(j-1)//2 else A4[i-5][j-5]
    if 5<=i<=8 and 5<=j<=8 else f.TWO if i==j==9 else Z for j in range(10)) for i in range(10))
DUAL=((F(11,9),F(2,9),F(-4,9)),(F(2,9),F(2,9),F(-1,9)),(F(-4,9),F(-1,9),F(2,9)))
BOUNDS=((F(1,12),F(7,6),F(1,6)),(F(1,4),F(1,6),F(5,12)),(F(1,6),F(5,12),F(1,4)))
SORT=cmp_to_key(lambda a,b:f.sign(f.sub(a,b)))

def ceil(x):return -((-x.numerator)//x.denominator)
def all_le(a,b):return all(f.sign(f.sub(b,a),e)>=0 for e in range(3))
@lru_cache(None)
def alphabet(caps):
    if any(f.sign(c,e)<0 for e,c in enumerate(caps)):return ()
    ms=[]
    for e,c in enumerate(caps):
        m=0 if c==Z else 1
        while f.sign(f.sub(f.scale(O,m*m),c),e)<0:m*=2
        ms.append(m)
    bs=[ceil(sum(b*m for b,m in zip(row,ms))) for row in BOUNDS];out=[]
    for raw in product(*(range(-b,b+1) for b in bs)):
        x=tuple(map(F,raw));xx=sq(x)
        if all(f.sign(f.sub(c,xx),e)>=0 for e,c in enumerate(caps)):out.append(x)
    return tuple(sorted(out,key=SORT))

@lru_cache(None)
def cov(v):return f.mv(G,v)
def norm(v):return f.dot(v,cov(v))
def compatible(v,walls):
    q=norm(v)
    for r in walls:
        p=f.dot(v,cov(r))
        if f.sign(p)>0:return False
        c=f.divide(f.scale(sq(p),4),f.mul(q,norm(r)))
        if f.sign(f.sub(c,f.FOUR))<0 and c not in ALLOWED:return False
    return True

def initial():
    out=[]
    for k in range(2):
        r=[Z]*10;r[1+2*k]=f.neg(O);out.append(tuple(r))
        r=[Z]*10;r[1+2*k]=D;r[2+2*k]=f.neg(D);out.append(tuple(r))
    for i in range(5,10):
        r=[Z]*10;r[i]=f.neg(O);out.append(tuple(r))
    return out

@lru_cache(None)
def candidates(x,q):
    target=f.add(q,f.mul(f.scale(T,2),sq(x)))
    weights=tuple(f.mul(f.scale(q,F(1,2)),w) for w in
        alphabet((f.divide(f.scale(target,8),sq(q)),)*3) if f.sign(w)>=0)
    assert all(map(f.integral,weights))
    bs=[]
    for u,v in product(weights,repeat=2):
        if f.sign(f.sub(v,u))<0:continue
        y=(f.divide(f.add(f.scale(u,2),f.mul(C,v)),B),f.divide(f.add(f.mul(C,u),f.scale(v,2)),B))
        e=f.dot(y,(u,v))
        if all(map(f.integral,y)) and all_le(e,target):bs.append((e,y))
    # A4 inverse integral iff w1+2w2+3w3+4w4 == 0 mod5.
    residue=lambda w:tuple(int(t)%5 for t in w)
    last={}
    for w in weights:last.setdefault(residue(f.scale(w,4)),[]).append(w)
    aa=[]
    for a,b,c in product(weights,repeat=3):
        need=residue(f.neg(f.add(f.add(a,f.scale(b,2)),f.scale(c,3))))
        for d in last.get(need,()):
            ww=(a,b,c,d);y=f.mv(AI,ww);assert all(map(f.integral,y))
            e=f.dot(y,ww)
            if all_le(e,target):aa.append((e,y))
    ss=tuple((f.scale(sq(z),2),z) for z in alphabet((f.scale(target,F(1,2)),)*3)
             if f.sign(z)>=0 and f.integral(f.divide(f.scale(z,4),q)))
    tails={}
    for e,y in aa:
        for s,z in ss:tails.setdefault(f.add(e,s),[]).append((*y,z))
    out=set();imprimitive=0
    for first,second in product(bs,repeat=2):
        left=f.sub(f.sub(target,first[0]),second[0])
        for tail in tails.get(left,()):
            r=(x,*first[1],*second[1],*tail)
            assert norm(r)==q and all(f.integral(f.divide(f.scale(t,2),q)) for t in cov(r))
            if all(f.integral(f.scale(t,F(1,2))) for t in r):imprimitive+=1;continue
            out.add(r)
    return out,{'weights':len(weights),'block_vectors':len(bs),'a4_vectors':len(aa),'scalar_vectors':len(ss),'imprimitive':imprimitive}

def verify(source,output):
    start=time.monotonic();data=json.loads(source.read_text())
    assert data['completed'] and data['dimension']==9
    assert tuple(f.vecdecode(r) for r in data['gram'])==G
    assert f.vecdecode(data['safe_primitive_norm_pool'])==POOL
    assert [[f.trace(f.mul(d,x)) for x in (O,T,sq(T))] for d in DUAL]==[[1,0,0],[0,1,0],[0,0,1]]
    assert all(f.sign(f.sub(f.scale(O,b),d),e)>0 and f.sign(f.add(f.scale(O,b),d),e)>0
        for d,bs in zip(DUAL,BOUNDS) for e,b in enumerate(bs))
    exponent=f.scale(B,10)
    assert f.integral(f.divide(exponent,f.scale(T,2)))
    assert f.integral(f.scale(exponent,F(1,2)))
    assert all(f.integral(f.mul(exponent,x)) for row in AI for x in row)
    assert all(f.integral(f.divide(f.mul(exponent,x),B)) for row in ((f.TWO,C),(C,f.TWO)) for x in row)
    assert all(f.integral(f.divide(f.scale(exponent,2),q)) for q in POOL)
    dominant,binary=f.binary_census();expected=set()
    binary={k:v for k,v in binary.items() if k not in ('full_H9_primitive_oriented_zero_roots','complete_zero_census_argument')}
    for q,tail in dominant.items():
        for k in range(2):
            r=[Z]*10;r[1+2*k:3+2*k]=tail;expected.add((q,tuple(r)))
    # Highest A4 root and positive A1 root; all20 A4 roots follow by Weyl action.
    expected.add((f.TWO,(Z,Z,Z,Z,Z,O,O,O,O,Z)))
    expected.add((f.TWO,(Z,Z,Z,Z,Z,Z,Z,Z,Z,O)))
    actual={(f.decode(r['norm']),f.vecdecode(r['vector'])) for r in data['zero_height_dominant_census']}
    assert actual==expected
    for q in POOL:
        roots,_=candidates(Z,q);assert {(q,r) for r in roots}=={row for row in expected if row[0]==q}
    controls=[]
    for q in POOL:
        cap=f.divide(q,f.scale(T,2))
        for x in alphabet((f.scale(cap,data['height_bound']),f.neg(cap),f.neg(cap))):
            if f.sign(x)>0 and f.integral(f.divide(f.mul(f.scale(T,4),x),q)):
                controls.append((f.divide(f.mul(f.scale(T,2),sq(x)),q),q,x))
    controls.sort(key=lambda row:SORT(row[0]));recorded={(f.decode(r['height']),f.decode(r['norm']),f.decode(r['time'])):r for r in data['levels']}
    assert set(controls)==set(recorded) and len(controls)==data['control_count']
    walls=initial();assert tuple(f.vecdecode(r) for r in data['initial_walls'])==tuple(walls)
    seen={f.ray_key(r) for r in walls};checks=[]
    for height,packet in groupby(controls,key=lambda row:row[0]):
        added=[]
        for _,q,x in packet:
            roots,counts=candidates(x,q);row=recorded[(height,q,x)]
            assert len(roots)==row['candidates'] and all(row[k]==v for k,v in counts.items())
            keep=[];rejected=0
            for r in sorted(roots):
                if any(f.sign(f.dot(r,cov(s)))>0 for s in walls):continue
                if not compatible(r,walls):rejected+=1;continue
                if f.ray_key(r) in seen:continue
                seen.add(f.ray_key(r));keep.append(r)
            assert len(keep)==row['accepted'] and rejected==row['acute_noncoxeter_rejected']
            added.extend(keep);checks.append({'height':f.encoded(height),'norm':f.encoded(q),'candidates':len(roots),'accepted':len(keep)})
        assert all(compatible(r,added[:i]) for i,r in enumerate(added));walls.extend(added)
    saved=tuple(f.vecdecode(r) for r in data['all_walls'])
    assert {f.ray_key(r) for r in walls}=={f.ray_key(r) for r in saved} and len(walls)==data['wall_count']
    assert all(compatible(r,saved[:i]) for i,r in enumerate(saved))
    for r in saved:
        q=norm(r);assert q in POOL and all(f.sign(q,e)>0 for e in range(3))
        assert all(map(f.integral,r)) and all(f.integral(f.divide(f.scale(x,2),q)) for x in cov(r))
        assert not all(f.integral(f.scale(x,F(1,2))) for x in r)
    center=f.vecdecode(data['strict_interior']);v=f.vecdecode(data['geometry']['hyperideal_witness']['vector'])
    assert f.sign(center[0])>0 and f.sign(norm(center))<0 and all(f.sign(f.dot(center,cov(r)))<0 for r in saved)
    assert f.sign(v[0])>0 and f.sign(norm(v))>0 and norm(v)==f.decode(data['geometry']['hyperideal_witness']['norm'])
    assert all(f.sign(f.dot(v,cov(r)))<=0 for r in saved)
    eps=O
    while f.sign(norm(tuple(f.add(x,f.mul(eps,y)) for x,y in zip(v,center))))<=0:eps=f.scale(eps,F(1,2))
    strict=tuple(f.add(x,f.mul(eps,y)) for x,y in zip(v,center))
    assert all(f.sign(f.dot(strict,cov(r)))<0 for r in saved)
    det=f.neg(f.mul(f.scale(T,20),sq(B)));comparisons=[]
    for name in ('theta','2theta','theta_beta','2theta_beta'):
        src=ROOT/f'results/fresh_cubic9_dihedral_h9_{name}_h30_20261013.json';old=json.loads(src.read_text())
        prior=f.neg(f.mul(f.decode(old['alpha']),sq(sq(B))));ratio=f.divide(det,prior);nn=f.field_norm(ratio)
        assert nn>0 and (isqrt(nn.numerator)**2!=nn.numerator or isqrt(nn.denominator)**2!=nn.denominator)
        comparisons.append({'source':str(src.relative_to(ROOT)),'determinant_ratio':f.encoded(ratio),'ratio_norm':str(nn),'nonsquare':True})
    out={'status':'PASS','source':str(source.relative_to(ROOT)),'wall_count':len(saved),
        'all_wall_pairs_checked':len(saved)*(len(saved)-1)//2,'height_bound':data['height_bound'],
        'bounded_levels_independently_reenumerated':checks,'determinant':f.encoded(det),
        'dual_exponent':f.encoded(exponent),'safe_primitive_norm_pool':f.vecencode(POOL),
        'compared_host_determinants':comparisons,
        'zero_census':{'binary_independent_census':binary,'dominant_orbits':6,'full_oriented_roots':94,'chamber':'I2(18)^2 x A4 x A1'},
        'escape':{'vector':f.vecencode(v),'norm':f.encoded(norm(v)),
            'strict_spacelike_vector':f.vecencode(strict),'strict_spacelike_norm':f.encoded(norm(strict))},
        'scope':'Four compared conductor-nine B^4+1 spaces have different determinant squareclasses. The bounded root census and every wall angle were independently replayed. The saved prefix has infinite volume; the full host remains open.',
        'seconds':round(time.monotonic()-start,3)}
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:out[k] for k in ('status','wall_count','all_wall_pairs_checked','seconds')}));return out

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();verify(a.source.resolve(),a.output)
