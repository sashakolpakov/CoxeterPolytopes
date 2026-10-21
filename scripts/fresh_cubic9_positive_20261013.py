#!/usr/bin/env python3
"""Changed-coefficient prime-five hosts, with forced H5-face calibrations.

Search <-alpha>+B^k+A4+A1, k=0 or2, alpha=T orT*beta.
Safe primitive norm classes are divisors of10 beta (beta absent for k=0,
alpha=T). Norms divisible by4 are imprimitive by the dyadic argument.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations_with_replacement,groupby,permutations
import json
from pathlib import Path
import time
from fresh_cubic9_mixed_search_20261013 import (
    K,T,Z,O,C,D,BETA,BLOCK,ALLOWED,inverse,mv,dot,pair,enc,encode_matrix,
    ray_key,finite_angle_counts,geometry,cubic9_inertia,le_all,A4,AINV,PREFIX_INV,compatible)
from fresh_cubic9_alphabet_20261013 import alphabet
ROOT=Path(__file__).resolve().parents[1]

@lru_cache(None)
def spatial_lists(target,q):
    weights=tuple(q*k/2 for k in alphabet((8*target/(q*q),)*3)
                  if k.sign()>=0 and (q*k/2).is_integral)
    bs=[]
    for i,u in enumerate(weights):
        for v in weights[i:]:
            y=((2*u+C*v)/BETA,(C*u+2*v)/BETA)
            if not all(t.is_integral for t in y):continue
            e=dot(y,(u,v))
            if le_all(e,target):bs.append((e,y))
    aa=[]
    def extend(ws):
        k=len(ws)
        for w in weights:
            ww=(*ws,w);low=dot(ww,mv(PREFIX_INV[k],ww))
            if not le_all(low,target):continue
            if k<3:extend(ww)
            else:
                y=mv(AINV,ww)
                if all(t.is_integral for t in y):aa.append((low,y))
    extend(())
    ss=tuple((2*z*z,z) for z in alphabet((target/2,)*3) if z.sign()>=0 and (4*z/q).is_integral)
    return tuple(bs),tuple(aa),ss,len(weights)

class Host:
    def __init__(self,alpha,blocks):
        assert alpha in (T,T*BETA) and blocks in (0,2)
        self.alpha=alpha;self.blocks=blocks;self.n=2*blocks+5
        self.norms=tuple(K(k)*b for b in ((O,BETA) if blocks or alpha==T*BETA else (O,)) for k in (1,2,5,10))
        start=2*blocks+1
        self.G=tuple(tuple(-alpha if i==j==0 else BLOCK[(i-1)%2][(j-1)%2]
            if 1<=i<start and 1<=j<start and (i-1)//2==(j-1)//2 else A4[i-start][j-start]
            if start<=i<start+4 and start<=j<start+4 else K(2) if i==j==self.n else Z
            for j in range(self.n+1)) for i in range(self.n+1))
    def initial(self):
        out=[]
        for k in range(self.blocks):
            r=[Z]*(self.n+1);r[1+2*k]=-O;out.append(tuple(r))
            r=[Z]*(self.n+1);r[1+2*k]=D;r[2+2*k]=-D;out.append(tuple(r))
        for i in range(2*self.blocks+1,self.n+1):
            r=[Z]*(self.n+1);r[i]=-O;out.append(tuple(r))
        return out
    def candidates(self,x,q):
        target=q+self.alpha*x*x;bs,aa,ss,nw=spatial_lists(target,q);tails={}
        for e,y in aa:
            for s,z in ss:
                total=e+s
                if le_all(total,target):tails.setdefault(total,[]).append((*y,z))
        out=[]
        for ids in combinations_with_replacement(range(len(bs)),self.blocks):
            leftover=target-sum((bs[i][0] for i in ids),Z)
            for tail in tails.get(leftover,()):
                for order in sorted(set(permutations(ids))):
                    r=(x,*(v for i in order for v in bs[i][1]),*tail)
                    assert pair(r,self.G,r)==q and all(t.is_integral for t in r)
                    assert all((2*t/q).is_integral for t in mv(self.G,r))
                    out.append(r)
        return tuple(out),{'weights':nw,'block_vectors':len(bs),'a4_vectors':len(aa),'scalar_vectors':len(ss)}
    def run(self,height,output):
        start=time.monotonic();zero=[];g=self.G
        assert cubic9_inertia(g).negative==1 and all(cubic9_inertia(g,embedding=e).positive==self.n+1 for e in (1,2))
        for q in self.norms:
            roots,counts=self.candidates(Z,q)
            zero.extend({'norm':q.to_json(),'vector':enc(r)} for r in roots)
            print(json.dumps({'zero_norm':str(q),'roots':len(roots),**counts}),flush=True)
        assert len(zero)==2*self.blocks+2
        walls=self.initial();keys={ray_key(r) for r in walls}
        controls=sorted((self.alpha*x*x/q,q,x) for q in self.norms
            for x in alphabet((height*q/self.alpha,-q/self.alpha,-q/self.alpha))
            if x.sign()>0 and (2*self.alpha*x/q).is_integral)
        out={'status':'RUNNING','field':'Q(2cos(pi/9))','host':'<-alpha>+B^'+str(self.blocks)+'+A4+A1',
            'dimension':self.n,'blocks':self.blocks,'gram':encode_matrix(g),'alpha':self.alpha.to_json(),
            'safe_primitive_norm_pool':enc(self.norms),'zero_height_dominant_census':zero,
            'initial_walls':encode_matrix(walls),'height_bound':height,'control_count':len(controls),
            'levels':[],'completed':False,
            'scope':'Changed coefficient and forced-face calibration. Complete bounded primitive root search; no full-host reflectivity or worldwide novelty claim.'}
        print(json.dumps({'phase':'start','controls':len(controls),'zero_orbits':len(zero),'dimension':self.n}),flush=True)
        for h,packet in groupby(controls,key=lambda row:row[0]):
            accepted=[];cov=tuple(mv(g,r) for r in walls);norms=tuple(dot(r,c) for r,c in zip(walls,cov))
            for _,q,x in packet:
                roots,counts=self.candidates(x,q);kept=0;rejected=0
                for r in roots:
                    if any(dot(r,c).sign()>0 for c in cov):continue
                    if not compatible(r,q,cov,norms):rejected+=1;continue
                    key=ray_key(r)
                    if key in keys:continue
                    keys.add(key);accepted.append(r);kept+=1
                out['levels'].append({'height':h.to_json(),'norm':q.to_json(),'time':x.to_json(),
                    'candidates':len(roots),'accepted':kept,'acute_noncoxeter_rejected':rejected,**counts})
            finite_angle_counts(accepted,g);walls.extend(accepted)
            out.update(all_walls=encode_matrix(walls),wall_count=len(walls),completed_height=h.to_json())
            output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
            print(json.dumps({'height':str(h),'added':len(accepted),'walls':len(walls)}),flush=True)
        base=mv(inverse(BLOCK),(O,K(2)));tail=mv(AINV,(O,O,O,O));eps=K(F(1,2))
        while True:
            inside=(O,*(eps*y for _ in range(self.blocks) for y in base),*(eps*y for y in tail),eps)
            if pair(inside,g,inside).sign()<0 and all(pair(inside,g,r).sign()<0 for r in walls):break
            eps/=2
        out.update(status='COMPLETE_BOUNDED_CHANGED_COEFFICIENT_SEARCH',completed=True,
            completed_height=K(height).to_json(),strict_interior=enc(inside),finite_angle_counts=finite_angle_counts(walls,g))
        out['geometry']=geometry(walls,g,self.n,inside=inside)
        out['seconds']=round(time.monotonic()-start,3)
        output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
        print(json.dumps({'phase':'done','walls':len(walls),'geometry':out['geometry']['status'],'seconds':out['seconds']}),flush=True)
        return out
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--alpha',choices=('theta','theta_beta'),default='theta')
    p.add_argument('--blocks',type=int,choices=(0,2),default=2);p.add_argument('--height',type=int,default=30)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    Host(T if a.alpha=='theta' else T*BETA,a.blocks).run(a.height,a.output)
