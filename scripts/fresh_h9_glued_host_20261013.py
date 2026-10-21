#!/usr/bin/env python3
"""New conductor-sixteen arithmetic hosts <-T+1>+F4^k+<1>."""
import argparse
from functools import lru_cache
from itertools import combinations,groupby,product
import json
from pathlib import Path
import time
from fresh_h9_glued_host_20261013_arithmetic import *

ROOT=Path(__file__).resolve().parents[1]

@lru_cache(None)
def block_vectors(energy):
    weights=alphabet((2*energy,)*4);found=[];nodes=0
    def rec(prefix):
        nonlocal nodes
        nodes+=1;k=len(prefix)
        if k:
            lower=pair(prefix,SCHUR[k-1],prefix)
            if not le_all(lower,energy):return
            # Remaining covectors are nonnegative in the primary embedding.
            if pair(prefix,tuple(row[:k]for row in FI[:k]),prefix)>energy:return
        if k==4:
            v=mv(FI,prefix);assert all(x.is_integral for x in v)
            found.append((lower,v));return
        for w in weights:rec((*prefix,w))
    rec(())
    return tuple(found),{'alphabet':len(weights),'nodes':nodes,'block_vectors':len(found)}

def compatible(r,q,walls,g):
    for s in walls:
        ps=pair(r,g,s)
        if ps>0:return False
        ratio=4*ps*ps/(q*pair(s,g,s))
        if ratio<4 and ratio not in ALLOWED:return False
    return True

def geometry(walls,g,inside,n):
    import numpy as np
    from scipy.optimize import linprog
    cov=tuple(mv(g,r)for r in walls);M=np.array([[x.real()for x in c]for c in cov]);supports=[]
    spatial_inverse=inverse(tuple(tuple(g[i][j]for j in range(1,n+1))for i in range(1,n+1)))
    for j in range(n):
        direction=(Z,*(row[j]for row in spatial_inverse))
        if any(dot(c,direction)>0 for c in cov):continue
        step=O
        while pair(tuple(x+step*y for x,y in zip(inside,direction)),g,tuple(x+step*y for x,y in zip(inside,direction)))<=0:step*=2
        low=Z
        while True:
            v=tuple(x+step*y for x,y in zip(inside,direction))
            if pair(v,g,inside)<0:break
            mid=(low+step)/2;candidate=tuple(x+mid*y for x,y in zip(inside,direction))
            if pair(candidate,g,candidate)>0:step=mid
            else:low=mid
        assert pair(v,g,v)>0 and all(dot(c,v)<0 for c in cov)
        return {'status':'EXACT_STRICT_SPACELIKE_ESCAPE','vector':enc(v),'norm':pair(v,g,v).to_json(),
            'wall_pairings':enc(tuple(dot(c,v)for c in cov)),'active_facets':[],
            'recession_direction':enc(direction),'strict_vector':enc(v),'strict_norm':pair(v,g,v).to_json(),
            'strict_center_pairing':pair(v,g,inside).to_json()}
    rng=np.random.default_rng(20261013)
    for _ in range(160):
        lp=linprog(rng.normal(size=n),A_ub=M[:,1:],b_ub=-M[:,0],bounds=(None,None),method='highs')
        if not lp.success:continue
        ids=tuple(map(int,np.flatnonzero(abs(M[:,1:]@lp.x+M[:,0])<1e-7)))
        if len(ids)==n:supports.append(ids)
    seen=set()
    for ids in supports:
        if ids in seen:continue
        seen.add(ids)
        try:v=(O,*mv(inverse(tuple(tuple(cov[i][j]for j in range(1,n+1))for i in ids)),tuple(-cov[i][0]for i in ids)))
        except(StopIteration,ZeroDivisionError):continue
        values=tuple(dot(c,v)for c in cov)
        if any(x>0 for x in values)or pair(v,g,v)<=0:continue
        eps=O
        while pair(tuple(x+eps*y for x,y in zip(v,inside)),g,tuple(x+eps*y for x,y in zip(v,inside)))<=0:eps/=2
        strict=tuple(x+eps*y for x,y in zip(v,inside))
        assert all(dot(c,strict)<0 for c in cov)
        return {'status':'EXACT_STRICT_SPACELIKE_ESCAPE','vector':enc(v),'norm':pair(v,g,v).to_json(),
                'wall_pairings':enc(values),'active_facets':list(ids),'strict_vector':enc(strict),
                'strict_norm':pair(strict,g,strict).to_json(),'strict_center_pairing':pair(strict,g,inside).to_json()}
    return {'status':'UNDECIDED','exact_vertex_proposals':len(seen)}

class Host:
    def __init__(self,blocks):
        assert blocks in(0,1,2)
        self.blocks=blocks;self.n=4*blocks+1;n=self.n
        self.g=tuple(tuple(-ALPHA if i==j==0 else O if i==j==n else F4[(i-1)%4][(j-1)%4]
            if 1<=i<=4*blocks and 1<=j<=4*blocks and(i-1)//4==(j-1)//4 else Z for j in range(n+1))for i in range(n+1))
    def candidates(self,x,q):
        E=q+ALPHA*x*x
        vectors,stats=block_vectors(E)if self.blocks else(((Z,()),),{})
        scalars={z*z:z for z in alphabet((E,)*4)};found=[];imprimitive=0
        for packets in product(vectors,repeat=self.blocks):
            rem=E-sum((e for e,y in packets),Z)
            if rem not in scalars:continue
            r=(x,*(v for e,y in packets for v in y),scalars[rem])
            assert pair(r,self.g,r)==q
            if all((v/T).is_integral for v in r):imprimitive+=1;continue
            assert all((2*v/q).is_integral for v in mv(self.g,r))
            found.append(r)
        return found,{**stats,'scalar_vectors':len(scalars),'imprimitive':imprimitive}
    def run(self,height,output):
        started=time.monotonic();g=self.g;zero=[]
        for q in NORMS:
            roots,stats=self.candidates(Z,q)
            zero.extend({'norm':q.to_json(),'vector':enc(r)}for r in roots)
            print(json.dumps({'phase':'zero','q':str(q),'roots':len(roots),**stats}),flush=True)
        assert len(zero)==2*self.blocks+1
        walls=[tuple(-O if i==j else Z for i in range(self.n+1))for j in range(1,self.n+1)]
        keys={ray_key(r)for r in walls}
        controls=sorted((ALPHA*x*x/q,q,x)for q in NORMS for x in alphabet((height*q/ALPHA,-q/ALPHA,-q/ALPHA,-q/ALPHA))if x>0)
        out={'status':'RUNNING','field':'Q(T), T^4-4T^2+2=0, T=2cos(pi/8)',
            'host':'<-T+1>+F4^k+<1>','blocks':self.blocks,'dimension':self.n,'gram':menc(g),
            'alpha':ALPHA.to_json(),'norm_pool':enc(NORMS),'height_bound':height,'control_count':len(controls),
            'zero_dominant_roots':zero,'initial_walls':menc(walls),'levels':[],'completed':False}
        print(json.dumps({'phase':'positive','controls':len(controls),'zero_orbits':len(zero),'dimension':self.n}),flush=True)
        for h,packet in groupby(controls,key=lambda row:row[0]):
            added=[]
            for _,q,x in packet:
                roots,stats=self.candidates(x,q);kept=0;noncoxeter=0
                for r in roots:
                    if any(pair(r,g,s)>0 for s in walls):continue
                    if not compatible(r,q,walls,g):noncoxeter+=1;continue
                    key=ray_key(r)
                    if key in keys:continue
                    keys.add(key);added.append(r);kept+=1
                out['levels'].append({'height':h.to_json(),'norm':q.to_json(),'time':x.to_json(),
                    'candidates':len(roots),'accepted':kept,'noncoxeter':noncoxeter,**stats})
            assert all(compatible(r,pair(r,g,r),added[:i],g)for i,r in enumerate(added))
            walls.extend(added)
            out.update(all_walls=menc(walls),wall_count=len(walls))
            output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
            print(json.dumps({'height':h.real(),'new':len(added),'walls':len(walls)}),flush=True)
        tail=tuple(v for _ in range(self.blocks)for v in mv(FI,(O,)*4));eps=O
        while True:
            inside=(O,*(eps*v for v in tail),eps)
            if pair(inside,g,inside)<0 and all(pair(inside,g,r)<0 for r in walls):break
            eps/=2
        angle_counts={}
        for i,j in combinations(range(len(walls)),2):
            a,b=walls[i],walls[j];v=4*pair(a,g,b)**2/(pair(a,g,a)*pair(b,g,b))
            assert pair(a,g,b)<=0 and(v>=4 or v in ALLOWED)
            if v<4:angle_counts[str(v)]=angle_counts.get(str(v),0)+1
        out.update(status='COMPLETE_BOUNDED_QUARTIC_HOST',completed=True,all_walls=menc(walls),wall_count=len(walls),
            strict_interior=enc(inside),finite_angle_counts=angle_counts,geometry=geometry(walls,g,inside,self.n),
            seconds=round(time.monotonic()-started,3),
            scope='Complete dominant primitive-root census through the stated height in this new degree-four host. Exact saved cone geometry is separate from the unresolved full reflection group.')
        output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
        print(json.dumps({'status':out['status'],'walls':len(walls),'geometry':out['geometry']['status'],'seconds':out['seconds']}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--blocks',type=int,choices=(0,1,2),default=1)
    p.add_argument('--height',type=int,default=30);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    Host(a.blocks).run(a.height,a.output)
