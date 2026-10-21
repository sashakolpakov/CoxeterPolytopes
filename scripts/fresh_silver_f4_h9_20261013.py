#!/usr/bin/env python3
"""Exact bounded Vinberg searches on <-a>+F4^k+(<1>).

This uses equal-norm F4 with middle pairing -sqrt(2), not the rational
long/short F4 root lattice. Candidate norm alphabets include ramified
prime classes; positive parallel duplicates are identified explicitly.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, groupby, product
import json
from pathlib import Path
import time

from coxeter_polytopes.exact.silver import Silver as S, SQRT2, silver_determinant, silver_inertia
from coxeter_polytopes.search.silver_path_host import bounded_silver_integers

Z=S(0); O=S(1); ETA=2+SQRT2; UNIT=1+SQRT2

def inverse(matrix):
    n=len(matrix); a=[list(row)+[S(int(i==j)) for j in range(n)] for i,row in enumerate(matrix)]
    for j in range(n):
        k=next(i for i in range(j,n) if a[i][j]); a[j],a[k]=a[k],a[j]
        q=a[j][j]; a[j]=[v/q for v in a[j]]
        for i in range(n):
            if i!=j:
                q=a[i][j]; a[i]=[x-q*y for x,y in zip(a[i],a[j])]
    return tuple(tuple(row[n:]) for row in a)

def mv(matrix,v): return tuple(sum((a*b for a,b in zip(row,v)),Z) for row in matrix)
def dot(x,y): return sum((a*b for a,b in zip(x,y)),Z)
def pair(x,matrix,y): return dot(x,mv(matrix,y))
def enc(v): return [x.to_json() for x in v]
def ray_key(v):
    first=next(x for x in v if x)
    return (first.sign(),tuple(x/first for x in v))

F4=tuple(tuple(S(2) if i==j else -SQRT2 if {i,j}=={1,2}
               else -O if abs(i-j)==1 else Z for j in range(4)) for i in range(4))
FI=inverse(F4)
SCHUR=tuple(inverse(tuple(row[:k] for row in F4[:k])) for k in range(1,5))
assert silver_determinant(F4)==1

@lru_cache(None)
def block_vectors(target,q):
    alphabet=tuple(x for x in bounded_silver_integers((2*target,2*target))
                   if x.sign()>=0 and (2*x/q).is_integral)
    found=[]; nodes=0
    def rec(prefix):
        nonlocal nodes
        nodes+=1; k=len(prefix)
        if k:
            lower=pair(prefix,SCHUR[k-1],prefix)
            if any((target-lower).sign(e)<0 for e in (0,1)): return
            if (target-pair(prefix,tuple(row[:k] for row in FI[:k]),prefix)).sign()<0:return
        if k==4:
            v=mv(FI,prefix); assert all(x.is_integral for x in v)
            found.append((lower,v)); return
        for x in alphabet: rec((*prefix,x))
    rec(())
    return tuple(found),nodes

def angle_gate(walls,gram):
    cov=[mv(gram,r) for r in walls]; norms=[dot(r,c) for r,c in zip(walls,cov)]
    allowed={Z,O,S(2),S(3),2+SQRT2}; counts={}
    for i,j in combinations(range(len(walls)),2):
        p=dot(walls[i],cov[j]); assert p.sign()<=0
        v=4*p*p/(norms[i]*norms[j])
        if (v-4).sign()<0:
            assert v in allowed, ('non-Coxeter acute angle',i,j,str(v))
            counts[str(v)]=counts.get(str(v),0)+1
    return counts

def geometry(walls,gram,inside,dimension,exhaust=False):
    import numpy as np
    from scipy.optimize import linprog
    cov=tuple(mv(gram,r) for r in walls)
    real=lambda x:float(x.a)+float(x.b)*2**.5
    M=np.array([[real(x) for x in r] for r in cov])
    supports=[]
    if exhaust: supports=combinations(range(len(walls)),dimension)
    else:
        rng=np.random.default_rng(20261013)
        for _ in range(160):
            lp=linprog(rng.normal(size=dimension),A_ub=M[:,1:],b_ub=-M[:,0],
                       bounds=[(None,None)]*dimension,method='highs')
            if not lp.success:continue
            ids=tuple(map(int,np.flatnonzero(abs(M[:,1:]@lp.x+M[:,0])<1e-7)))
            if len(ids)==dimension: supports.append(ids)
    seen=set(); finite=0; cusps=0; escape=None
    for ids in supports:
        try:
            B=tuple(tuple(cov[i][j] for j in range(1,dimension+1)) for i in ids)
            spatial=mv(inverse(B),tuple(-cov[i][0] for i in ids))
        except (StopIteration,ZeroDivisionError):continue
        v=(O,*spatial)
        if v in seen:continue
        values=tuple(dot(c,v) for c in cov)
        if any(x.sign()>0 for x in values):continue
        seen.add(v); norm=pair(v,gram,v)
        if norm.sign()>0:
            escape={'vector':enc(v),'norm':norm.to_json(),'active_facets':list(ids),
                    'wall_pairings':enc(values)}
            if not exhaust:break
        elif norm.sign()<0:finite+=1
        else:cusps+=1
    return {'status':'EXACT_ESCAPE' if escape else 'COMPLETE_VERTEX_CENSUS' if exhaust else 'UNDECIDED',
            'finite_vertices':finite if exhaust else None,'cusps':cusps if exhaust else None,
            'tested_distinct_feasible_vertices':len(seen),'hyperideal_witness':escape}

class Host:
    def __init__(self,alpha,blocks,scalar):
        self.alpha=alpha; self.blocks=blocks; self.scalar=scalar; self.n=4*blocks+int(scalar)
        self.gram=tuple(tuple(-alpha if i==j==0 else O if scalar and i==j==self.n
            else F4[(i-1)%4][(j-1)%4] if 1<=i<=4*blocks and 1<=j<=4*blocks and (i-1)//4==(j-1)//4
            else Z for j in range(self.n+1)) for i in range(self.n+1))
        assert silver_inertia(self.gram).positive==self.n and silver_inertia(self.gram).negative==1
        assert silver_inertia(self.gram,embedding=1).positive==self.n+1
        # Primitive q divides2det(B); 2 is ramified and ETA is its TP prime.
        # TP units are squares. These lists represent every divisor class.
        pools={UNIT:(O,ETA,S(2)), SQRT2:(O,ETA,S(2),2*ETA),
               2*UNIT:(O,ETA,S(2),2*ETA,S(4))}
        self.norms=pools[alpha]

    def candidates(self,x,q):
        target=q+self.alpha*x*x
        vectors,nodes=block_vectors(target,q)
        scalars=tuple(z for z in bounded_silver_integers((target,target))
                      if z.sign()>=0 and (2*z/q).is_integral) if self.scalar else (Z,)
        by_square={z*z:z for z in scalars}; out=[]
        for blocks in product(vectors,repeat=self.blocks):
            left=target-sum((r[0] for r in blocks),Z)
            if left not in by_square:continue
            root=(x,*(c for _,v in blocks for c in v),*((by_square[left],) if self.scalar else ()))
            assert pair(root,self.gram,root)==q
            out.append(root)
        return out,nodes

    def run(self,height,output):
        started=time.monotonic(); zero=[]
        for q in self.norms:
            roots,_=self.candidates(Z,q)
            zero.extend({'norm':q.to_json(),'vector':enc(r)} for r in roots)
        # Completeness is also checked by a full finite root-system audit.
        walls=[tuple(-O if i==j else Z for i in range(self.n+1)) for j in range(1,self.n+1)]
        keys={ray_key(r) for r in walls}
        controls=sorted((self.alpha*x*x/q,q,x) for q in self.norms
            for x in bounded_silver_integers((height*q/self.alpha,-q/self.alpha))
            if x.sign()>0 and (2*self.alpha*x/q).is_integral)
        out={'field':'Q(sqrt2)','host':'<-alpha>+F4^'+str(self.blocks)+('+<1>' if self.scalar else ''),
             'alpha':self.alpha.to_json(),'dimension':self.n,'blocks':self.blocks,'scalar':self.scalar,
             'gram':[enc(row) for row in self.gram],'norm_pool':enc(self.norms),'height':height.to_json(),
             'zero_height_dominant_census':zero,'levels':[],'completed':False,
             'scope':'Bounded arithmetic construction. No global reflection-group conclusion.'}
        print(json.dumps({'phase':'start','dimension':self.n,'controls':len(controls),'zero_roots':len(zero)}),flush=True)
        for h,packet in groupby(controls,key=lambda row:row[0]):
            accepted=[]
            for _,q,x in packet:
                roots,nodes=self.candidates(x,q); kept=0; duplicates=0
                for r in roots:
                    assert all((2*c/q).is_integral for c in mv(self.gram,r))
                    if not all(pair(r,self.gram,w).sign()<=0 for w in walls):continue
                    key=ray_key(r)
                    if key in keys: duplicates+=1; continue
                    keys.add(key); accepted.append(r); kept+=1
                out['levels'].append({'height':h.to_json(),'norm':q.to_json(),'control':x.to_json(),
                                      'nodes':nodes,'candidates':len(roots),'roots':kept,'parallel_duplicates':duplicates})
            assert all(pair(r,self.gram,s).sign()<=0 for i,r in enumerate(accepted) for s in accepted[:i])
            walls.extend(accepted)
            out.update(all_walls=[enc(r) for r in walls],wall_count=len(walls),completed_height=h.to_json())
            output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
            print(json.dumps({'height':str(h),'new':len(accepted),'walls':len(walls)}),flush=True)
        spatial=tuple(c for _ in range(self.blocks) for c in mv(FI,(O,)*4))
        eps=S(F(1,2))
        while True:
            inside=(O,*(eps*c for c in spatial),*((eps,) if self.scalar else ()))
            if pair(inside,self.gram,inside).sign()<0 and all(pair(inside,self.gram,w).sign()<0 for w in walls):break
            eps/=2
        out.update(completed=True,completed_height=height.to_json(),strict_interior=enc(inside),
                   finite_angle_counts=angle_gate(walls,self.gram),root_seconds=round(time.monotonic()-started,3))
        out['geometry']=geometry(walls,self.gram,inside,self.n,exhaust=self.n<=4 and len(walls)<25)
        out['seconds']=round(time.monotonic()-started,3)
        output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
        print(json.dumps({'phase':'done','walls':len(walls),'geometry':out['geometry']['status'],'seconds':out['seconds']}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--alpha',default='1,1')
    p.add_argument('--blocks',type=int,choices=(1,2),default=2); p.add_argument('--no-scalar',action='store_true')
    p.add_argument('--height',type=int,default=30); p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); a.output.parent.mkdir(parents=True,exist_ok=True)
    Host(S(*(int(x) for x in a.alpha.split(','))),a.blocks,not a.no_scalar).run(S(a.height),a.output)
