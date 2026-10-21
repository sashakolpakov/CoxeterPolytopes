#!/usr/bin/env python3
"""Independent Fraction census of the new prime-seven cubic-nine hosts."""
import argparse
from functools import lru_cache
from itertools import product,groupby
import json
from pathlib import Path
import time
from fresh_cubic9_positive_20261013_independent import f,Census as OldCensus,Z,O,T,F,ROOT,sq,ALLOWED
from fresh_cubic9_mixed_independent_20261013 import SORT
from fresh_cubic9_ellipsoid_audit_20261013 import ellipsoid
AI=tuple(tuple(f.scale(O,F(min(i,j)*(7-max(i,j)),7))for j in range(1,7))for i in range(1,7))
PREFIX=tuple(tuple(tuple(f.scale(O,F(min(i,j)*(k+1-max(i,j)),k+1))for j in range(1,k+1))for i in range(1,k+1))for k in range(1,7))
A6=tuple(tuple(f.scale(O,2 if i==j else -1 if abs(i-j)==1 else 0)for j in range(6))for i in range(6))

class Census(OldCensus):
    def __init__(self,data):
        self.alpha=f.decode(data['alpha']);self.blocks=data['blocks'];self.n=7+2*self.blocks
        assert self.alpha in(T,f.scale(T,2))and self.blocks in(0,1)
        self.pool=tuple(f.scale(u,k)for u in((O,f.BETA)if self.blocks else(O,))
                        for k in((1,2,7,14)if self.alpha==T else(2,4,14,28))
                        if u!=O or k%7)
        a=1+2*self.blocks
        self.G=tuple(tuple(f.neg(self.alpha)if i==j==0 else f.B[i-1][j-1]
            if self.blocks and 1<=i<=2 and 1<=j<=2 else A6[i-a][j-a]
            if a<=i<a+6 and a<=j<a+6 else f.TWO if i==j==self.n else Z
            for j in range(self.n+1))for i in range(self.n+1))
    def initial(self):
        rows=[]
        if self.blocks:
            r=[Z]*(self.n+1);r[1]=f.neg(O);rows.append(tuple(r));r=[Z]*(self.n+1);r[1:3]=(f.D,f.neg(f.D));rows.append(tuple(r))
        for i in range(1+2*self.blocks,self.n+1):
            r=[Z]*(self.n+1);r[i]=f.neg(O);rows.append(tuple(r))
        return rows
    @lru_cache(None)
    def candidates(self,x,q):
        target=f.add(q,f.mul(self.alpha,sq(x)));less=lambda x:all(f.sign(f.sub(target,x),e)>=0 for e in range(3))
        weights=tuple(f.mul(f.scale(q,F(1,2)),k)for k in ellipsoid(f.divide(f.scale(target,8),sq(q)))
            if f.sign(k)>=0 and f.integral(f.mul(f.scale(q,F(1,2)),k)))
        residue=lambda x:tuple(int(t)%7 for t in x)
        last={}
        for w in weights:last.setdefault(residue(f.scale(w,6)),[]).append(w)
        aa=[];nodes=[0]*6
        def extend(ws):
            k=len(ws)+1
            if k==6:
                need=residue(f.neg(f.sum_field(f.scale(w,i+1)for i,w in enumerate(ws))))
                options=last.get(need,())
            else:options=weights
            for w in options:
                ww=(*ws,w);e=f.dot(ww,f.mv(PREFIX[k-1],ww))
                if not less(e):continue
                nodes[k-1]+=1
                if k<6:extend(ww)
                else:
                    y=f.mv(AI,ww);assert all(map(f.integral,y));aa.append((e,y))
        extend(())
        ss=tuple((f.scale(sq(z),2),z)for z in ellipsoid(f.scale(target,F(1,2)))
                 if f.sign(z)>=0 and f.integral(f.divide(f.scale(z,4),q)))
        tails={}
        for e,y in aa:
            for s,z in ss:tails.setdefault(f.add(e,s),[]).append((*y,z))
        bs=[(Z,())]
        if self.blocks:
            bs=[]
            for u,v in product(weights,repeat=2):
                if f.sign(f.sub(v,u))<0:continue
                y=(f.divide(f.add(f.scale(u,2),f.mul(T,v)),f.BETA),f.divide(f.add(f.mul(T,u),f.scale(v,2)),f.BETA))
                e=f.dot(y,(u,v))
                if all(map(f.integral,y))and less(e):bs.append((e,y))
        roots=set();imprimitive=0
        for e,y in bs:
            for tail in tails.get(f.sub(target,e),()):
                r=(x,*y,*tail);assert self.norm(r)==q
                assert all(f.integral(f.divide(f.scale(v,2),q))for v in self.cov(r))
                if all(f.integral(f.scale(v,F(1,2)))for v in r):imprimitive+=1;continue
                roots.add(r)
        return roots,{'weights':len(weights),'a6_vectors':len(aa),'scalar_vectors':len(ss),'block_vectors':len(bs),'imprimitive':imprimitive,'schur_nodes':nodes}

def verify(source):
    data=json.loads(source.read_text());c=Census(data)
    assert data['completed']and tuple(f.vecdecode(row)for row in data['gram'])==c.G
    assert f.vecdecode(data['safe_primitive_norm_pool'])==c.pool
    zero=set()
    for q in c.pool:
        roots,_=c.candidates(Z,q);zero.update((q,r)for r in roots)
    assert zero=={(f.decode(row['norm']),f.vecdecode(row['vector']))for row in data['zero_height_dominant_census']}
    expected={(f.TWO,(Z,*(Z for _ in range(2*c.blocks)),*(O for _ in range(6)),Z)),
              (f.TWO,(Z,*(Z for _ in range(2*c.blocks+6)),O))}
    if c.blocks:
        dominant,_=f.binary_census()
        for q,r in dominant.items():expected.add((q,(Z,*r,*(Z for _ in range(7)))))
    assert zero==expected
    controls=c.controls(data['height_bound']);saved={(f.decode(r['height']),f.decode(r['norm']),f.decode(r['time'])):r for r in data['levels']}
    assert set(controls)==set(saved)and len(controls)==data['control_count']
    walls=c.initial();assert walls==[f.vecdecode(r)for r in data['initial_walls']]
    keys={f.ray_key(r)for r in walls};checks=[]
    for height,packet in groupby(controls,key=lambda row:row[0]):
        added=[]
        for _,q,x in packet:
            roots,counts=c.candidates(x,q);row=saved[height,q,x]
            assert len(roots)==row['candidates']
            assert all(row[k]==v for k,v in counts.items()if k!='schur_nodes')
            keep=[];rejected=0
            for r in sorted(roots):
                if any(f.sign(f.dot(r,c.cov(w)))>0 for w in walls):continue
                if not c.compatible(r,walls):rejected+=1;continue
                if f.ray_key(r)in keys:continue
                keys.add(f.ray_key(r));keep.append(r)
            assert len(keep)==row['accepted']and rejected==row['acute_noncoxeter_rejected']
            added.extend(keep);checks.append({'height':f.encoded(height),'norm':f.encoded(q),'candidates':len(roots),'accepted':len(keep),'schur_nodes':counts['schur_nodes']})
        assert all(c.compatible(r,added[:i])for i,r in enumerate(added));walls.extend(added)
    actual=tuple(f.vecdecode(r)for r in data['all_walls'])
    assert len(actual)==len(walls)==data['wall_count']and {f.ray_key(r)for r in actual}==keys
    for r in actual:
        q=c.norm(r);assert q in c.pool and all(map(f.integral,r))
        assert all(f.integral(f.divide(f.scale(v,2),q))for v in c.cov(r))
        assert all(f.sign(q,e)>0 for e in range(3))
    det=f.determinant(c.G);expected_det=f.scale(f.mul(c.alpha,f.BETA if c.blocks else O),-14)
    assert det==expected_det
    return {'status':'PASS_COMPLETE_ENUMERATION','source':str(source.relative_to(ROOT)),
        'zero_dominant_orbits':len(zero),'complete_zero_census':True,'complete_time_census':True,
        'complete_spatial_census':True,'controls':checks,'determinant':f.encoded(det),'geometry':c.check_geometry(data)}

def run(sources,output):
    start=time.monotonic();rows=[verify(p.resolve())for p in sources]
    result={'status':'PASS','arithmetic':'Independent Fraction field arithmetic; spatial trace-ellipsoid alphabets; sequential Schur bounds with final mod7 join; independent time boxes.',
        'cases':rows,'seconds':round(time.monotonic()-start,3)}
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','cases':len(rows),'seconds':result['seconds']}));return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('sources',type=Path,nargs='+');p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.sources,a.output)
