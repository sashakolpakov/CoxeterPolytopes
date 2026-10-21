#!/usr/bin/env python3
"""Prime-seven conductor-nine hosts <-alpha>+B^k+A6+A1, k=0,1.

The rank-three cross term of inverse A6 is separable, permitting a complete
mod-seven join. The H7 complement is searched before adding the B block.
"""
import argparse
from itertools import groupby
import json
import time
from fresh_cubic9_positive_20261013 import (
    K,T,Z,O,C,D,BETA,BLOCK,ALLOWED,inverse,mv,dot,pair,enc,encode_matrix,
    ray_key,finite_angle_counts,geometry,cubic9_inertia,le_all,compatible,
    alphabet,ROOT,F,lru_cache,Path)
A6=tuple(tuple(K(2 if i==j else -1 if abs(i-j)==1 else 0)for j in range(6))for i in range(6))
AI=inverse(A6)
SCHUR=inverse(tuple(tuple(A6[i][j]for j in range(3))for i in range(3)))
LEFT=tuple(tuple(AI[i][j]for j in range(3))for i in range(3))

class Host:
    def __init__(self,alpha=T,blocks=0):
        assert alpha in (T,2*T) and blocks in (0,1)
        self.alpha=alpha;self.blocks=blocks;self.n=7+2*blocks;a=1+2*blocks
        self.pool=tuple(K(k)*u for u in ((O,BETA)if blocks else(O,))
                        for k in ((1,2,7,14)if alpha==T else(2,4,14,28))
                        if u!=O or k%7)
        self.G=tuple(tuple(-alpha if i==j==0 else BLOCK[i-1][j-1]
            if blocks and 1<=i<=2 and 1<=j<=2 else A6[i-a][j-a]
            if a<=i<a+6 and a<=j<a+6 else K(2)if i==j==self.n else Z
            for j in range(self.n+1))for i in range(self.n+1))
    def initial(self):
        rows=[]
        if self.blocks:
            r=[Z]*(self.n+1);r[1]=-O;rows.append(tuple(r));r=[Z]*(self.n+1);r[1:3]=(D,-D);rows.append(tuple(r))
        for i in range(1+2*self.blocks,self.n+1):
            r=[Z]*(self.n+1);r[i]=-O;rows.append(tuple(r))
        return rows
    @lru_cache(None)
    def candidates(self,x,q):
        target=q+self.alpha*x*x
        weights=tuple(q*k/2 for k in alphabet((8*target/(q*q),)*3)
                      if k.sign()>=0 and(q*k/2).is_integral)
        residue=lambda x:tuple(int(t)%7 for t in(x.a,x.b,x.c))
        triples=[];groups={}
        for a in weights:
            for b in weights:
                if not le_all(K(F(2,3))*(a*a+a*b+b*b),target):continue
                for c in weights:
                    w=(a,b,c)
                    if not le_all(dot(w,mv(SCHUR,w)),target):continue
                    u=a+2*b+3*c;e=dot(w,mv(LEFT,w));row=(w,u,e)
                    triples.append(row);groups.setdefault(residue(u),[]).append(row)
        aa=[];joins=0
        for left,u,el in triples:
            for right,v,er in groups[residue(u)]:
                joins+=1;e=el+er+2*u*v/7
                if not le_all(e,target):continue
                w=(*left,*reversed(right));y=mv(AI,w);assert all(t.is_integral for t in y)
                aa.append((e,y))
        ss=tuple((2*z*z,z)for z in alphabet((target/2,)*3)
                 if z.sign()>=0 and(4*z/q).is_integral)
        tails={}
        for e,y in aa:
            for s,z in ss:
                if le_all(e+s,target):tails.setdefault(e+s,[]).append((*y,z))
        bs=[(Z,())]
        if self.blocks:
            bs=[]
            for i,u in enumerate(weights):
                for v in weights[i:]:
                    y=((2*u+T*v)/BETA,(T*u+2*v)/BETA);e=dot(y,(u,v))
                    if all(t.is_integral for t in y)and le_all(e,target):bs.append((e,y))
        out=[];imprimitive=0
        for e,y in bs:
            for tail in tails.get(target-e,()):
                r=(x,*y,*tail);assert pair(r,self.G,r)==q
                assert all((2*t/q).is_integral for t in mv(self.G,r))
                if all((t/2).is_integral for t in r):imprimitive+=1;continue
                out.append(r)
        assert len(set(out))==len(out)
        return tuple(out),{'weights':len(weights),'triple_vectors':len(triples),'joined_pairs':joins,
            'a6_vectors':len(aa),'scalar_vectors':len(ss),'block_vectors':len(bs),'imprimitive':imprimitive}
    def run(self,height,output):
        started=time.monotonic();g=self.G;zero=[]
        assert cubic9_inertia(g).negative==1 and all(cubic9_inertia(g,embedding=e).positive==self.n+1 for e in(1,2))
        for q in self.pool:
            roots,counts=self.candidates(Z,q);zero.extend({'norm':q.to_json(),'vector':enc(r)}for r in roots)
            print(json.dumps({'zero_norm':str(q),'roots':len(roots),**counts}),flush=True)
        assert len(zero)==2+2*self.blocks
        walls=self.initial();seen={ray_key(r)for r in walls}
        controls=sorted((self.alpha*x*x/q,q,x)for q in self.pool
            for x in alphabet((height*q/self.alpha,-q/self.alpha,-q/self.alpha))
            if x.sign()>0 and(2*self.alpha*x/q).is_integral)
        out={'status':'RUNNING','host':'<-alpha>+B^k+A6+A1','dimension':self.n,'blocks':self.blocks,
            'alpha':self.alpha.to_json(),'gram':encode_matrix(g),'safe_primitive_norm_pool':enc(self.pool),
            'height_bound':height,'control_count':len(controls),'zero_height_dominant_census':zero,
            'initial_walls':encode_matrix(walls),'levels':[],'completed':False}
        print(json.dumps({'phase':'start','controls':len(controls),'dimension':self.n}),flush=True)
        for h,packet in groupby(controls,key=lambda row:row[0]):
            added=[];cov=tuple(mv(g,r)for r in walls);nn=tuple(dot(r,c)for r,c in zip(walls,cov))
            for _,q,x in packet:
                roots,counts=self.candidates(x,q);kept=0;rejected=0
                for r in roots:
                    if any(dot(r,c).sign()>0 for c in cov):continue
                    if not compatible(r,q,cov,nn):rejected+=1;continue
                    if ray_key(r)in seen:continue
                    seen.add(ray_key(r));added.append(r);kept+=1
                out['levels'].append({'height':h.to_json(),'norm':q.to_json(),'time':x.to_json(),
                    'candidates':len(roots),'accepted':kept,'acute_noncoxeter_rejected':rejected,**counts})
            finite_angle_counts(added,g);walls.extend(added)
            print(json.dumps({'height':str(h),'added':len(added),'walls':len(walls)}),flush=True)
        tail=mv(AI,(O,)*6);bs=mv(inverse(BLOCK),(O,K(2)));eps=O
        while True:
            inside=(O,*(eps*y for y in bs if self.blocks),*(eps*y for y in tail),eps)
            if pair(inside,g,inside).sign()<0 and all(pair(inside,g,r).sign()<0 for r in walls):break
            eps/=2
        out.update(status='COMPLETE_BOUNDED_PRIME_SEVEN_HOST',completed=True,all_walls=encode_matrix(walls),
            wall_count=len(walls),strict_interior=enc(inside),finite_angle_counts=finite_angle_counts(walls,g),
            geometry=geometry(walls,g,self.n,inside=inside),seconds=round(time.monotonic()-started,3),
            scope='Complete primitive bounded root census on the stated host. Prime7 distinguishes this quadratic space from the compared conductor-nine B4 and prime5 forms; no worldwide novelty or full-host reflectivity claim.')
        output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps({'phase':'done','walls':len(walls),'geometry':out['geometry']['status'],'seconds':out['seconds']}),flush=True)
        return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--alpha',choices=('theta','2theta'),default='theta')
    p.add_argument('--blocks',type=int,choices=(0,1),default=0);p.add_argument('--height',type=int,default=30)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    Host(T if args.alpha=='theta'else 2*T,args.blocks).run(args.height,args.output)
