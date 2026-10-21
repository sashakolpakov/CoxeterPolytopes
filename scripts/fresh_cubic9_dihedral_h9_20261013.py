#!/usr/bin/env python3
"""Cubic-nine Vinberg search with genuine dihedral-18 zero-height walls.

Host <-a> + I2(9)^k + optional <1>. The integral coordinate swap in each
dihedral block refines the initial chamber from I2(9) to I2(18).
All enumeration bounds and acceptance decisions are exact.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, groupby, permutations
import json
from pathlib import Path
import time

from coxeter_polytopes.exact.cubic9 import Cubic9 as K, THETA9 as T, cubic9_inertia
from fresh_cubic9_alphabet_20261013 import alphabet

ROOT=Path(__file__).resolve().parents[1]
Z=K(0);O=K(1);C=T;D=2+T-T*T;BETA=4-C*C
BLOCK=((K(2),-C),(-C,K(2)))
ALLOWED={Z,O,K(2),K(3),C*C,2+C}
assert D*D==2-C and BETA.norm()==3 and D.norm() in (1,-1)


def inverse(matrix):
    n=len(matrix);a=[list(row)+[K(int(i==j)) for j in range(n)] for i,row in enumerate(matrix)]
    for j in range(n):
        k=next(i for i in range(j,n) if a[i][j]);a[j],a[k]=a[k],a[j]
        q=a[j][j];a[j]=[x/q for x in a[j]]
        for i in range(n):
            if i==j:continue
            q=a[i][j]
            if q:a[i]=[x-q*y for x,y in zip(a[i],a[j])]
    return tuple(tuple(row[n:]) for row in a)


def dot(a,b):return sum((x*y for x,y in zip(a,b,strict=True) if x and y),Z)
def mv(A,v):return tuple(dot(row,v) for row in A)
def pair(a,B,b):return dot(a,mv(B,b))
def enc(v):return [x.to_json() for x in v]
def encode_matrix(A):return [enc(row) for row in A]
def decode_matrix(A):return tuple(tuple(map(K.from_json,row)) for row in A)
def le_all(a,b):return all((b-a).sign(e)>=0 for e in range(3))
def ray_key(v):
    x=next(t for t in v if t)
    return x.sign(),tuple(y/x for y in v)


@lru_cache(None)
def block_vectors(target,q):
    weights=tuple(w for w in alphabet((2*target,)*3) if w.sign()>=0 and (2*w/q).is_integral)
    out=[]
    for i,u in enumerate(weights):
        for v in weights[i:]:
            # w1>=0, w2>=w1 is the I2(18) chamber, not merely I2(9).
            y=((2*u+C*v)/BETA,(C*u+2*v)/BETA)
            if not all(x.is_integral for x in y):continue
            energy=u*y[0]+v*y[1]
            if le_all(energy,target):out.append((energy,y,(u,v)))
    return tuple(out)


def finite_angle_counts(walls,gram):
    cov=tuple(mv(gram,r) for r in walls);norms=tuple(dot(r,n) for r,n in zip(walls,cov));counts={}
    for i,j in combinations(range(len(walls)),2):
        p=dot(walls[i],cov[j]);assert p.sign()<=0
        value=4*p*p/(norms[i]*norms[j])
        if (value-4).sign()<0:
            assert value in ALLOWED,(i,j,str(value))
            counts[str(value)]=counts.get(str(value),0)+1
    return counts


def geometry(walls,gram,dimension,exhaust=False,inside=None):
    import numpy as np
    from scipy.optimize import linprog
    real=lambda x:float(x.a)+float(x.b)*1.879385241571817+float(x.c)*3.532088886237956
    cov=tuple(mv(gram,r) for r in walls)
    A=np.array([[real(x) for x in row] for row in cov])
    norms=np.linalg.norm(A[:,1:],axis=1);M=A/norms[:,None]
    supports=[]
    if exhaust:supports=combinations(range(len(walls)),dimension)
    else:
        rng=np.random.default_rng(20261013)
        for _ in range(100):
            prop=linprog(rng.normal(size=dimension),A_ub=M[:,1:],b_ub=-M[:,0],bounds=(None,None),method='highs')
            if not prop.success:continue
            ids=tuple(map(int,np.flatnonzero(abs(M[:,1:]@prop.x+M[:,0])<1e-7)))
            if len(ids)==dimension:supports.append(ids)
    vertices=[];seen=set();escape=None
    for ids in supports:
        try:v=(O,*mv(inverse(tuple(tuple(cov[i][j] for j in range(1,dimension+1)) for i in ids)),tuple(-cov[i][0] for i in ids)))
        except (StopIteration,ZeroDivisionError):continue
        if v in seen:continue
        ps=tuple(dot(v,n) for n in cov)
        if any(p.sign()>0 for p in ps):continue
        seen.add(v);q=pair(v,gram,v)
        row={'vector':enc(v),'norm':q.to_json(),'active_facets':list(ids),'wall_pairings':enc(ps)}
        if q.sign()>0:
            escape=row
            if not exhaust:break
        vertices.append(row)
    recession=None
    if not escape and inside is not None:
        # Affine sections can be unbounded without any hyperideal vertex.
        # Scout simple rays, but certify every pairing and the resulting
        # spacelike point exactly; failure to find a ray proves nothing.
        directions=[tuple(O if j in ids else Z for j in range(dimension+1))
                    for ids in [*( (i,) for i in range(1,dimension+1)),
                                *combinations(range(1,dimension+1),2),
                                tuple(range(1,dimension+1))]]
        for d in directions:
            if any(dot(d,n).sign()>0 for n in cov):continue
            assert pair(d,gram,d).sign()>0
            scale=O
            while pair(tuple(x+scale*y for x,y in zip(inside,d)),gram,
                       tuple(x+scale*y for x,y in zip(inside,d))).sign()<=0:
                scale*=2
            v=tuple(x+scale*y for x,y in zip(inside,d));ps=tuple(dot(v,n) for n in cov)
            assert all(p.sign()<0 for p in ps)
            escape={'vector':enc(v),'norm':pair(v,gram,v).to_json(),
                    'active_facets':[],'wall_pairings':enc(ps)}
            recession={'vector':enc(d),'norm':pair(d,gram,d).to_json(),
                       'wall_pairings':enc(tuple(dot(d,n) for n in cov))}
            break
    result={'status':'EXACT_ESCAPE' if escape else 'EXHAUSTIVE_AFFINE_VERTEX_CENSUS' if exhaust else 'UNDECIDED',
            'hyperideal_witness':escape,'vertices':vertices if exhaust else [],
            'note':'A complete vertex census alone requires a bounded affine section before asserting compactness.'}
    if recession:result['exact_affine_recession_ray']=recession
    # A positive dependence among spanning spatial normals proves that the
    # affine section has no nonzero recession direction.  Combined with the
    # exhaustive, timelike vertex census, this proves compactness.
    if exhaust and len(walls)==dimension+1:
        try:
            dependence=(*mv(inverse(tuple(tuple(cov[i][j] for i in range(dimension))
                                          for j in range(1,dimension+1))),
                            tuple(-cov[-1][j] for j in range(1,dimension+1))),O)
        except (StopIteration,ZeroDivisionError):
            dependence=()
        if dependence and all(t.sign()>0 for t in dependence):
            assert all(sum((t*r[j] for t,r in zip(dependence,cov)),Z)==Z
                       for j in range(1,dimension+1))
            result['positive_affine_normal_dependence']=enc(dependence)
            if vertices and all(K.from_json(v['norm']).sign()<0 for v in vertices):
                result['status']='EXACT_COMPACT_TRIANGLE_CALIBRATION' if dimension==2 else 'EXACT_COMPACT_SIMPLEX_CALIBRATION'
                result['note']='Positive dependence proves affine boundedness; every affine vertex is exactly timelike.'
                labels={Z:2,O:3,K(2):4,K(3):6,C*C:9,2+C:18}
                result['coxeter_labels']=sorted(labels[4*pair(r,gram,s)**2/(pair(r,gram,r)*pair(s,gram,s))]
                                                for r,s in combinations(walls,2))
    return result


class Host:
    def __init__(self,alpha,blocks,scalar):
        self.alpha=alpha;self.blocks=blocks;self.scalar=scalar;self.n=2*blocks+int(scalar)
        self.gram=tuple(tuple(-alpha if i==j==0 else O if scalar and i==j==self.n
                  else BLOCK[(i-1)%2][(j-1)%2] if 1<=i<=2*blocks and 1<=j<=2*blocks and (i-1)//2==(j-1)//2
                  else Z for j in range(self.n+1)) for i in range(self.n+1))
        assert cubic9_inertia(self.gram).negative==1
        assert all(cubic9_inertia(self.gram,embedding=e).positive==self.n+1 for e in (1,2))
        self.norms=(O,K(2),BETA,2*BETA) if alpha in (T,T*BETA) else (O,K(2),K(4),BETA,2*BETA,4*BETA)
        assert alpha in (T,2*T,T*BETA,2*T*BETA)

    def candidates(self,x,q):
        target=q+self.alpha*x*x;blocks=block_vectors(target,q)
        scalars=tuple(z for z in alphabet((target,)*3) if z.sign()>=0 and (2*z/q).is_integral) if self.scalar else (Z,)
        scalar_squares={z*z:z for z in scalars};matches=[]
        if self.blocks<=2:
            for ids in combinations_with_replacement(range(len(blocks)),self.blocks):
                leftover=target-sum((blocks[i][0] for i in ids),Z)
                if leftover in scalar_squares:matches.append((ids,scalar_squares[leftover]))
        else:
            assert self.blocks==4
            pairs={}
            for i,j in combinations_with_replacement(range(len(blocks)),2):
                energy=blocks[i][0]+blocks[j][0]
                if le_all(energy,target):pairs.setdefault(energy,[]).append((i,j))
            for ss,z in scalar_squares.items():
                for energy,first in pairs.items():
                    second=pairs.get(target-ss-energy,())
                    for i,j in first:
                        for k,l in second:
                            if j<=k:matches.append(((i,j,k,l),z))
        roots=[]
        for ids,z in matches:
            for perm in sorted(set(permutations(ids))):
                r=(x,*(a for i in perm for a in blocks[i][1]),*((z,) if self.scalar else ()))
                assert pair(r,self.gram,r)==q and all((2*a/q).is_integral for a in mv(self.gram,r))
                roots.append(r)
        return roots,len(blocks)

    def initial(self):
        walls=[]
        for k in range(self.blocks):
            e=[Z]*(self.n+1);e[1+2*k]=-O;walls.append(tuple(e))
            e=[Z]*(self.n+1);e[1+2*k]=D;e[2+2*k]=-D;walls.append(tuple(e))
        if self.scalar:walls.append(tuple(-O if j==self.n else Z for j in range(self.n+1)))
        return walls

    def run(self,height,output):
        started=time.monotonic();zero=[]
        for q in self.norms:
            roots,count=self.candidates(Z,q)
            zero.extend({'norm':q.to_json(),'vector':enc(r)} for r in roots)
            print(json.dumps({'zero_norm':str(q),'dominant_roots':len(roots),'block_vectors':count}),flush=True)
        walls=self.initial();keys={ray_key(r) for r in walls}
        controls=sorted((self.alpha*x*x/q,q,x) for q in self.norms
                        for x in alphabet((height*q/self.alpha,-q/self.alpha,-q/self.alpha))
                        if x.sign()>0 and (2*self.alpha*x/q).is_integral)
        out={'field':'Q(2cos(pi/9))','alpha':self.alpha.to_json(),'blocks':self.blocks,'scalar':self.scalar,
             'dimension':self.n,'gram':encode_matrix(self.gram),'norm_pool':enc(self.norms),
             'zero_height_dominant_census':zero,'height_bound':height.to_json(),'control_count':len(controls),
             'levels':[],'completed':False,
             'scope':'Bounded full candidate norm-class Vinberg search; no full-host reflectivity conclusion.'}
        print(json.dumps({'phase':'start','dimension':self.n,'controls':len(controls),'zero_roots':len(zero)}),flush=True)
        for h,packet in groupby(controls,key=lambda r:r[0]):
            accepted=[];cov=tuple(mv(self.gram,w) for w in walls)
            for _,q,x in packet:
                roots,count=self.candidates(x,q);kept=0
                for r in roots:
                    if any(dot(r,n).sign()>0 for n in cov):continue
                    for old,n in zip(walls,cov):
                        value=4*dot(r,n)**2/(q*dot(old,n))
                        assert (value-4).sign()>=0 or value in ALLOWED,('acute_forbidden_angle',str(value))
                    key=ray_key(r)
                    if key in keys:continue
                    keys.add(key);accepted.append(r);kept+=1
                out['levels'].append({'height':h.to_json(),'norm':q.to_json(),'time':x.to_json(),
                                      'block_vectors':count,'candidate_roots':len(roots),'accepted':kept})
            assert all(pair(r,self.gram,s).sign()<=0 for i,r in enumerate(accepted) for s in accepted[:i])
            walls.extend(accepted);out.update(all_walls=encode_matrix(walls),wall_count=len(walls),completed_height=h.to_json())
            output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
            print(json.dumps({'height':str(h),'added':len(accepted),'walls':len(walls)}),flush=True)
        base=mv(inverse(BLOCK),(O,K(2)));eps=K(F(1,2))
        while True:
            inside=(O,*(eps*y for _ in range(self.blocks) for y in base),*((eps,) if self.scalar else ()))
            if pair(inside,self.gram,inside).sign()<0 and all(pair(inside,self.gram,w).sign()<0 for w in walls):break
            eps/=2
        out.update(completed=True,completed_height=height.to_json(),strict_interior=enc(inside),
                   finite_angle_counts=finite_angle_counts(walls,self.gram),root_seconds=round(time.monotonic()-started,3))
        out['geometry']=geometry(walls,self.gram,self.n,exhaust=self.n<=3 and len(walls)<=30,inside=inside)
        out['seconds']=round(time.monotonic()-started,3)
        output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
        print(json.dumps({'phase':'done','walls':len(walls),'geometry':out['geometry']['status'],'seconds':out['seconds']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--alpha',choices=('theta','2theta','theta_beta','2theta_beta'),default='theta')
    p.add_argument('--blocks',type=int,choices=(1,2,4),default=4);p.add_argument('--no-scalar',action='store_true')
    p.add_argument('--height',type=int,default=10);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();alpha={'theta':T,'2theta':2*T,'theta_beta':T*BETA,'2theta_beta':2*T*BETA}[a.alpha]
    Host(alpha,a.blocks,not a.no_scalar).run(K(a.height),a.output)
