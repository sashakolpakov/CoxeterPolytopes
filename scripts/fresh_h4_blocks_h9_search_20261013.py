#!/usr/bin/env python3
"""Exact Vinberg exploration with field-valued H4 positive blocks.

Default target: <-a> + H4 + H4 + <1>, in dimension nine over Q(phi).
All arithmetic gates are exact. A finite prefix is not a full-group verdict.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import groupby, product
import json
from pathlib import Path
import time

from coxeter_polytopes.exact.golden import Golden as G, ZERO as Z, ONE as O, PHI
from coxeter_polytopes.exact.matrix import determinant, inertia
from coxeter_polytopes.search.golden_path_host import bounded_golden_integers
from coxeter_polytopes.verify.chamber import diagnose_chamber_cone, verify_hyperideal_extreme_ray
from scripts.search_golden_path_h9_20261013 import complete_norm_pool


def inverse(matrix):
    n=len(matrix)
    a=[list(row)+[G(i==j) for j in range(n)] for i,row in enumerate(matrix)]
    for j in range(n):
        pivot=next(i for i in range(j,n) if a[i][j])
        a[j],a[pivot]=a[pivot],a[j]
        q=a[j][j];a[j]=[v/q for v in a[j]]
        for i in range(n):
            if i!=j:
                q=a[i][j];a[i]=[x-q*y for x,y in zip(a[i],a[j])]
    return tuple(tuple(row[n:]) for row in a)


def mv(matrix,v):
    return tuple(sum((a*b for a,b in zip(row,v)),Z) for row in matrix)


def pair(x,matrix,y):
    return sum((a*b for a,b in zip(x,mv(matrix,y))),Z)


H4=tuple(tuple(G(2) if i==j else -PHI if {i,j}=={0,1} else -O
              if abs(i-j)==1 else Z for j in range(4)) for i in range(4))
H4_INV=inverse(H4)
SCHUR=tuple(inverse(tuple(row[:k] for row in H4[:k])) for k in range(1,5))


@lru_cache(None)
def block_vectors(target,q):
    """Every primary-dominant H4 vector with energy <= target at both places."""
    alphabet=tuple(w for w in bounded_golden_integers((2*target,2*target))
                   if w.sign()>=0 and (2*w/q).is_integral)
    vectors=[];nodes=0
    def rec(prefix):
        nonlocal nodes
        nodes+=1;k=len(prefix)
        if k:
            minimum=pair(prefix,SCHUR[k-1],prefix)
            if any((target-minimum).sign(p)<0 for p in (0,1)):
                return
            # Primary inverse entries and all dominant weights are nonnegative.
            if (target-pair(prefix,tuple(row[:k] for row in H4_INV[:k]),prefix)).sign()<0:
                return
        if k==4:
            y=mv(H4_INV,prefix)
            assert all(c.is_integral for c in y)
            assert pair(y,H4,y)==minimum
            vectors.append((minimum,y));return
        for w in alphabet:
            rec((*prefix,w))
    rec(())
    return tuple(vectors),nodes


class Host:
    def __init__(self,alpha,blocks=2,scalar=True):
        self.alpha=alpha;self.blocks=blocks;self.scalar=scalar
        self.n=4*blocks+int(scalar)
        self.gram=tuple(tuple(-alpha if i==j==0 else
            O if scalar and i==j==self.n else
            H4[(i-1)%4][(j-1)%4] if 1<=i<=4*blocks and 1<=j<=4*blocks
                and (i-1)//4==(j-1)//4 else Z
            for j in range(self.n+1)) for i in range(self.n+1))
        assert inertia(self.gram).positive==self.n and inertia(self.gram).negative==1
        assert inertia(self.gram,embedding=1).positive==self.n+1
        assert determinant(H4)==5-3*PHI and (1/determinant(H4)).is_integral
        self.norms=complete_norm_pool(determinant(self.gram))
        # For a=2phi, spatial unimodularity plus q=-2phi*x² mod4
        # excludes primitive norm4, just as in the cubic and E8 odd hosts.
        if alpha==2*PHI:
            self.norms=tuple(q for q in self.norms if q!=4)

    def candidates(self,x,q):
        target=q+self.alpha*x*x
        block,nodes=block_vectors(target,q)
        if self.scalar:
            scalars=tuple(u for u in bounded_golden_integers((target,target))
                          if u.sign()>=0 and (2*u/q).is_integral)
        else:
            scalars=(Z,)
        scalar_buckets={u*u:u for u in scalars}
        output=[]
        for choices in product(block,repeat=self.blocks):
            remaining=target-sum((v[0] for v in choices),Z)
            if remaining not in scalar_buckets:
                continue
            vector=(x,*(c for _,y in choices for c in y),
                    *((scalar_buckets[remaining],) if self.scalar else ()))
            assert pair(vector,self.gram,vector)==q
            output.append(vector)
        return tuple(output),nodes,len(block)

    def geometry(self,walls,inside):
        if (len(walls)<=22 and self.n<=9) or self.n<=5:
            g=diagnose_chamber_cone(walls,self.gram,inside,self.n)
            summary={'finite_vertices':len(g.finite_vertices),'cusps':len(g.cusps),
                     'hyperideal_rays':len(g.hyperideal_rays),'opposite_sheet_rays':len(g.opposite_sheet_rays)}
            support=next((r for r in g.hyperideal_rays if len(r)==self.n),None)
            witness=None
            if support is not None:
                w=verify_hyperideal_extreme_ray(walls,self.gram,support,self.n)
                witness={'vector':enc(w.ray),'norm':w.norm.to_json(),'wall_pairings':enc(w.wall_pairings),
                         'active_facets':list(w.active_facets)}
            return summary,witness
        import numpy as np
        from scipy.optimize import linprog
        phi=(1+5**.5)/2
        real=lambda z:float(z.a)+float(z.b)*phi
        f=np.array([[real(z) for z in mv(self.gram,r)] for r in walls])
        rng=np.random.default_rng(20261013)
        for _ in range(100):
            lp=linprog(rng.normal(size=self.n),A_ub=f[:,1:],b_ub=-f[:,0],
                       bounds=[(None,None)]*self.n,method='highs')
            if not lp.success:continue
            active=np.flatnonzero(abs(f[:,1:]@lp.x+f[:,0])<1e-7)
            if len(active)!=self.n:continue
            try:w=verify_hyperideal_extreme_ray(walls,self.gram,tuple(map(int,active)),self.n)
            except Exception:continue
            return {'status':'EXACT_SPACELIKE_ESCAPE'}, {'vector':enc(w.ray),'norm':w.norm.to_json(),
                    'wall_pairings':enc(w.wall_pairings),'active_facets':list(w.active_facets)}
        return {'status':'UNDECIDED_NO_ESCAPE_FOUND'},None

    def run(self,height,output):
        start=time.monotonic()
        zero=[]
        for q in self.norms:
            roots,_,_=self.candidates(Z,q)
            zero.extend((q,r) for r in roots)
        assert len(zero)==self.blocks+int(self.scalar)
        assert sum(q==2 for q,_ in zero)==self.blocks
        assert sum(q==1 for q,_ in zero)==int(self.scalar)
        for q,r in zero:
            if q==1:
                assert self.scalar and r[-1]==1 and all(c==0 for c in r[:-1])
            else:
                assert sum(any(c!=0 for c in r[1+4*b:5+4*b]) for b in range(self.blocks))==1
                assert not self.scalar or r[-1]==0
        walls=[tuple(-O if i==j else Z for i in range(self.n+1)) for j in range(1,self.n+1)]
        controls=sorted((self.alpha*x*x/q,q,x) for q in self.norms
            for x in bounded_golden_integers((height*q/self.alpha,-q/self.alpha))
            if x.sign()>0 and (2*self.alpha*x/q).is_integral)
        levels=[]
        payload={'host':'<-a>+H4^'+str(self.blocks)+('+<1>' if self.scalar else ''),
                 'dimension':self.n,'alpha':self.alpha.to_json(),'blocks':self.blocks,'scalar':self.scalar,
                 'requested_height':height.to_json(),'norm_pool':enc(self.norms),
                 'gram':[enc(row) for row in self.gram],
                 'zero_height_dominant_census':[{'norm':q.to_json(),'vector':enc(r)} for q,r in zero],
                 'levels':levels,'scope':'Exact bounded search; no global nonreflectivity inference.'}
        print(json.dumps({'phase':'start','dimension':self.n,'control_levels':len(controls)}),flush=True)
        for h,packet in groupby(controls,key=lambda row:row[0]):
            accepted=[];records=[]
            for _,q,x in packet:
                candidates,nodes,block_count=self.candidates(x,q)
                kept=[r for r in candidates if all(pair(r,self.gram,w).sign()<=0 for w in walls)]
                for r in kept:
                    assert all((2*z/q).is_integral for z in mv(self.gram,r))
                accepted.extend(kept)
                records.append({'height':h.to_json(),'norm':q.to_json(),'control':x.to_json(),
                                'candidates':len(candidates),'roots':len(kept),'nodes':nodes,'block_vectors':block_count})
            assert all(pair(r,self.gram,s).sign()<=0 for i,r in enumerate(accepted) for s in accepted[:i])
            walls.extend(accepted);levels.extend(records)
            payload.update({'all_walls':[enc(r) for r in walls],'wall_count':len(walls),
                            'completed':False,'completed_height':h.to_json()})
            output.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
            print(json.dumps({'height':str(h),'new_roots':len(accepted),'walls':len(walls),
                              'nodes':sum(row['nodes'] for row in records)}),flush=True)
        spatial=tuple(c for _ in range(self.blocks) for c in mv(H4_INV,(O,)*4))
        epsilon=G(Fraction(1,2))
        while True:
            inside=(O,*(epsilon*c for c in spatial),*((epsilon,) if self.scalar else ()))
            if pair(inside,self.gram,inside).sign()<0 and all(pair(inside,self.gram,w).sign()<0 for w in walls):break
            epsilon/=2
        payload.update({'completed':True,'completed_height':height.to_json(),
                        'all_walls':[enc(r) for r in walls],'wall_count':len(walls),
                        'strict_interior':enc(inside),'root_seconds':round(time.monotonic()-start,3)})
        output.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
        payload['geometry'],payload['hyperideal_witness']=self.geometry(tuple(walls),inside)
        payload['seconds']=round(time.monotonic()-start,3)
        output.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
        print(json.dumps({'phase':'done','walls':len(walls),'geometry':payload['geometry'],'seconds':payload['seconds']}),flush=True)


def enc(vector):return [v.to_json() for v in vector]


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--alpha',default='0,1')
    parser.add_argument('--blocks',type=int,choices=(1,2,3),default=2)
    parser.add_argument('--no-scalar',action='store_true')
    parser.add_argument('--height',type=int,default=30)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    Host(G(*(int(x) for x in args.alpha.split(','))),args.blocks,not args.no_scalar).run(G(args.height),args.output)
