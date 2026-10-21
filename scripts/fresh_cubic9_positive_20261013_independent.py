#!/usr/bin/env python3
"""Independent Fraction reenumeration and geometry of changed prime-five hosts."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import product,groupby
from math import isqrt
import json
from pathlib import Path
import time
import fresh_cubic9_dihedral_independent_audit_20261013 as f
from fresh_cubic9_mixed_independent_20261013 import alphabet,all_le,AI,A4,SORT
ROOT=Path(__file__).resolve().parents[1]
Z=f.ZERO;O=f.ONE;T=f.THETA;B=f.BETA;sq=lambda x:f.mul(x,x)
ALLOWED={Z,O,f.TWO,f.scale(O,3),sq(T),f.add(f.TWO,T)}
class Census:
    def __init__(self,data):
        self.data=data;self.alpha=f.decode(data['alpha']);self.blocks=data['blocks'];self.n=2*self.blocks+5
        assert self.alpha in (T,f.mul(T,B)) and self.blocks in (0,2)
        self.pool=tuple(f.scale(b,k) for b in ((O,B) if self.blocks or self.alpha==f.mul(T,B) else (O,)) for k in (1,2,5,10))
        a=2*self.blocks+1
        self.G=tuple(tuple(f.neg(self.alpha) if i==j==0 else f.B[(i-1)%2][(j-1)%2]
            if 1<=i<a and 1<=j<a and (i-1)//2==(j-1)//2 else A4[i-a][j-a]
            if a<=i<a+4 and a<=j<a+4 else f.TWO if i==j==self.n else Z
            for j in range(self.n+1)) for i in range(self.n+1))
    @lru_cache(None)
    def cov(self,r):return f.mv(self.G,r)
    def norm(self,r):return f.dot(r,self.cov(r))
    def compatible(self,r,walls):
        q=self.norm(r)
        for s in walls:
            p=f.dot(r,self.cov(s))
            if f.sign(p)>0:return False
            c=f.divide(f.scale(sq(p),4),f.mul(q,self.norm(s)))
            if f.sign(f.sub(c,f.FOUR))<0 and c not in ALLOWED:return False
        return True
    def initial(self):
        rows=[]
        for k in range(self.blocks):
            r=[Z]*(self.n+1);r[1+2*k]=f.neg(O);rows.append(tuple(r))
            r=[Z]*(self.n+1);r[1+2*k]=f.D;r[2+2*k]=f.neg(f.D);rows.append(tuple(r))
        for i in range(2*self.blocks+1,self.n+1):
            r=[Z]*(self.n+1);r[i]=f.neg(O);rows.append(tuple(r))
        return rows
    @lru_cache(None)
    def candidates(self,x,q):
        target=f.add(q,f.mul(self.alpha,sq(x)))
        weights=tuple(f.mul(f.scale(q,F(1,2)),w) for w in alphabet((f.divide(f.scale(target,8),sq(q)),)*3)
            if f.sign(w)>=0 and f.integral(f.mul(f.scale(q,F(1,2)),w)))
        bs=[]
        for u,v in product(weights,repeat=2):
            if f.sign(f.sub(v,u))<0:continue
            y=(f.divide(f.add(f.scale(u,2),f.mul(T,v)),B),f.divide(f.add(f.mul(T,u),f.scale(v,2)),B))
            e=f.dot(y,(u,v))
            if all(map(f.integral,y)) and all_le(e,target):bs.append((e,y))
        residue=lambda w:tuple(int(t)%5 for t in w)
        last={}
        for w in weights:last.setdefault(residue(f.scale(w,4)),[]).append(w)
        aa=[]
        for a,b,c in product(weights,repeat=3):
            need=residue(f.neg(f.add(f.add(a,f.scale(b,2)),f.scale(c,3))))
            for d in last.get(need,()):
                ww=(a,b,c,d);y=f.mv(AI,ww);assert all(map(f.integral,y));e=f.dot(y,ww)
                if all_le(e,target):aa.append((e,y))
        ss=tuple((f.scale(sq(z),2),z) for z in alphabet((f.scale(target,F(1,2)),)*3)
            if f.sign(z)>=0 and f.integral(f.divide(f.scale(z,4),q)))
        tails={}
        for e,y in aa:
            for s,z in ss:tails.setdefault(f.add(e,s),[]).append((*y,z))
        out=set()
        for entries in product(bs,repeat=self.blocks):
            energy=f.sum_field(row[0] for row in entries)
            for tail in tails.get(f.sub(target,energy),()):
                r=(x,*(v for row in entries for v in row[1]),*tail)
                assert self.norm(r)==q and all(f.integral(f.divide(f.scale(t,2),q)) for t in self.cov(r))
                out.add(r)
        return out,{'weights':len(weights),'block_vectors':len(bs),'a4_vectors':len(aa),'scalar_vectors':len(ss)}
    def controls(self,height,pool=None):
        rows=[]
        for q in self.pool if pool is None else pool:
            cap=f.divide(q,self.alpha)
            for x in alphabet((f.scale(cap,height),f.neg(cap),f.neg(cap))):
                if f.sign(x)>0 and f.integral(f.divide(f.scale(f.mul(self.alpha,x),2),q)):
                    rows.append((f.divide(f.mul(self.alpha,sq(x)),q),q,x))
        return sorted(rows,key=lambda row:SORT(row[0]))
    def check_geometry(self,data):
        walls=tuple(f.vecdecode(r) for r in data['all_walls'])
        assert all(self.compatible(r,walls[:i]) for i,r in enumerate(walls))
        center=f.vecdecode(data['strict_interior']);v=f.vecdecode(data['geometry']['hyperideal_witness']['vector'])
        assert f.sign(center[0])>0 and f.sign(self.norm(center))<0 and all(f.sign(f.dot(center,self.cov(r)))<0 for r in walls)
        assert f.sign(v[0])>0 and f.sign(self.norm(v))>0 and self.norm(v)==f.decode(data['geometry']['hyperideal_witness']['norm'])
        assert all(f.sign(f.dot(v,self.cov(r)))<=0 for r in walls)
        eps=O
        while f.sign(self.norm(tuple(f.add(x,f.mul(eps,y)) for x,y in zip(v,center))))<=0:eps=f.scale(eps,F(1,2))
        strict=tuple(f.add(x,f.mul(eps,y)) for x,y in zip(v,center))
        assert all(f.sign(f.dot(strict,self.cov(r)))<0 for r in walls)
        return {'walls':len(walls),'all_pair_checks':len(walls)*(len(walls)-1)//2,'escape_norm':f.encoded(self.norm(v)),
            'strict_spacelike_vector':f.vecencode(strict),'strict_spacelike_norm':f.encoded(self.norm(strict))}

def verify(source):
    data=json.loads(source.read_text());c=Census(data)
    assert data['completed'] and tuple(f.vecdecode(r) for r in data['gram'])==c.G
    assert f.vecdecode(data['safe_primitive_norm_pool'])==c.pool
    dominant,_=f.binary_census();expected=set()
    for q,tail in dominant.items():
        for k in range(c.blocks):
            r=[Z]*(c.n+1);r[1+2*k:3+2*k]=tail;expected.add((q,tuple(r)))
    expected.add((f.TWO,(Z,*(Z for _ in range(2*c.blocks)),O,O,O,O,Z)))
    expected.add((f.TWO,(Z,*(Z for _ in range(2*c.blocks)),Z,Z,Z,Z,O)))
    actual={(f.decode(r['norm']),f.vecdecode(r['vector'])) for r in data['zero_height_dominant_census']}
    assert actual==expected
    for q in c.pool:
        roots,_=c.candidates(Z,q);assert {(q,r)for r in roots}=={row for row in expected if row[0]==q}
    controls=c.controls(data['height_bound']);recorded={(f.decode(r['height']),f.decode(r['norm']),f.decode(r['time'])):r for r in data['levels']}
    assert set(controls)==set(recorded) and len(controls)==data['control_count']
    walls=c.initial();assert tuple(f.vecdecode(r) for r in data['initial_walls'])==tuple(walls)
    seen={f.ray_key(r) for r in walls};checks=[]
    for height,packet in groupby(controls,key=lambda row:row[0]):
        added=[]
        for _,q,x in packet:
            roots,counts=c.candidates(x,q);row=recorded[(height,q,x)]
            assert len(roots)==row['candidates'] and all(row[k]==v for k,v in counts.items())
            keep=[];rejected=0
            for r in sorted(roots):
                if any(f.sign(f.dot(r,c.cov(s)))>0 for s in walls):continue
                if not c.compatible(r,walls):rejected+=1;continue
                if f.ray_key(r) in seen:continue
                seen.add(f.ray_key(r));keep.append(r)
            assert len(keep)==row['accepted'] and rejected==row['acute_noncoxeter_rejected']
            added.extend(keep);checks.append({'height':f.encoded(height),'norm':f.encoded(q),'candidates':len(roots),'accepted':len(keep)})
        assert all(c.compatible(r,added[:i]) for i,r in enumerate(added));walls.extend(added)
    saved=tuple(f.vecdecode(r) for r in data['all_walls'])
    assert {f.ray_key(r) for r in walls}=={f.ray_key(r) for r in saved} and len(walls)==data['wall_count']
    for r in saved:
        q=c.norm(r);assert q in c.pool and all(map(f.integral,r)) and all(f.integral(f.divide(f.scale(t,2),q)) for t in c.cov(r))
    return {'source':str(source.relative_to(ROOT)),'status':'PASS','complete_independent_census':True,
        'dimension':c.n,'alpha':f.encoded(c.alpha),'zero_dominant_orbits':len(expected),'controls':checks,'geometry':c.check_geometry(data)}

def run(sources,output):
    start=time.monotonic();rows=[verify(p.resolve()) for p in sources]
    out={'status':'PASS','arithmetic':'Independent standard-library Fraction field arithmetic; unpruned trace-dual boxes and mod5 Cartesian A4 enumeration.',
        'cases':rows,'seconds':round(time.monotonic()-start,3)}
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','cases':len(rows),'seconds':out['seconds']}));return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('sources',type=Path,nargs='+');p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.sources,a.output)
