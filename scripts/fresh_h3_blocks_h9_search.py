#!/usr/bin/env python3
"""Exact bounded exploration of <-alpha>+H3^3 over Z[phi].

The primitive norm alphabet follows from the discriminant exponent, and the
complete zero-height dominant census is checked before positive heights.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import product,groupby
import json
from pathlib import Path
import time
from coxeter_polytopes.exact.golden import Golden as G,PHI,ZERO as Z,ONE as O
from coxeter_polytopes.exact.matrix import determinant,inertia
from coxeter_polytopes.search.golden_path_host import bounded_golden_integers
from scripts.fresh_h4_blocks_h9_search_20261013 import inverse,mv,pair,enc
from scripts.fresh_rank3_h9_cones import escape,compatible

ROOT=Path(__file__).resolve().parents[1]
H3=tuple(tuple(G(2) if i==j else -PHI if {i,j}=={0,1} else -O if abs(i-j)==1 else Z for j in range(3)) for i in range(3))
HI=inverse(H3);SCHUR=tuple(inverse(tuple(row[:k] for row in H3[:k])) for k in range(1,4))

@lru_cache(None)
def block_vectors(target,q):
    alphabet=tuple(d for d in bounded_golden_integers((2*target,2*target)) if d.sign()>=0 and (2*d/q).is_integral)
    found=[];nodes=0
    def recur(prefix):
        nonlocal nodes
        nodes+=1;k=len(prefix)
        if k:
            energy=pair(prefix,SCHUR[k-1],prefix)
            if any((target-energy).sign(place)<0 for place in (0,1)):return
            if (target-pair(prefix,tuple(row[:k] for row in HI[:k]),prefix)).sign()<0:return
        if k==3:
            v=mv(HI,prefix)
            if all(x.is_integral for x in v):found.append((energy,v))
            return
        for d in alphabet:recur((*prefix,d))
    recur(())
    return tuple(found),nodes

def host(alpha):
    return tuple(tuple(-alpha if i==j==0 else H3[(i-1)%3][(j-1)%3] if i*j and (i-1)//3==(j-1)//3 else Z for j in range(10)) for i in range(10))

def norm_pool(alpha):
    assert alpha in (PHI,2*PHI)
    return (O,G(2),G(4)) if alpha==PHI else (G(2),G(4))

def candidates(alpha,x,q):
    target=q+alpha*x*x;vectors,nodes=block_vectors(target,q)
    buckets={}
    for energy,v in vectors:buckets.setdefault(energy,[]).append(v)
    out=[]
    for e1,v1 in vectors:
        for e2,v2 in vectors:
            for v3 in buckets.get(target-e1-e2,[]):out.append((x,*v1,*v2,*v3))
    return tuple(out),nodes,len(vectors)

def zero_census(alpha):
    out=[]
    for q in norm_pool(alpha):
        rs,_,_=candidates(alpha,Z,q)
        out.extend({'norm':q.to_json(),'root':enc(r)} for r in rs)
    assert len(out)==3 and all(G.from_json(r['norm'])==2 for r in out)
    assert all(sum(any(x!=Z for x in tuple(map(G.from_json,r['root']))[1+3*b:4+3*b]) for b in range(3))==1 for r in out)
    return out

def run(alpha,height,output):
    start=time.monotonic();gram=host(alpha);norms=norm_pool(alpha)
    assert inertia(gram).positive==9 and inertia(gram).negative==1 and inertia(gram,embedding=1).positive==10
    assert all((2*x).is_integral for row in inverse(gram) for x in row)
    zero=zero_census(alpha)
    walls=[tuple(-O if i==j else Z for i in range(10)) for j in range(1,10)]
    controls=sorted((alpha*x*x/q,q,x) for q in norms for x in bounded_golden_integers((height*q/alpha,-q/alpha))
        if x.sign()>0 and (2*alpha*x/q).is_integral)
    out={'host':'<-alpha>+H3^3','alpha':alpha.to_json(),'gram':[enc(row) for row in gram],
        'norm_pool':enc(norms),'height':height.to_json(),'zero_height_dominant_census':zero,'levels':[],'completed':False,
        'scope':'Bounded full-primitive-norm search. Field-space duplication is known; no worldwide or conjugacy novelty is asserted.'}
    print(json.dumps({'controls':len(controls),'zero_roots':len(zero)}),flush=True)
    for h,packet in groupby(controls,key=lambda c:c[0]):
        accepted=[]
        for _,q,x in packet:
            roots,nodes,blocks=candidates(alpha,x,q);new=[]
            for r in roots:
                assert pair(r,gram,r)==q and all((2*y/q).is_integral for y in mv(gram,r))
                if q in (4,28) and all((c/2).is_integral for c in r):continue
                if all(pair(r,gram,w).sign()<=0 for w in walls):new.append(r)
            accepted.extend(new);out['levels'].append({'height':h.to_json(),'norm':q.to_json(),'control':x.to_json(),
                'candidates':len(roots),'roots':len(new),'block_vectors':blocks,'nodes':nodes})
        assert all(compatible(r,s,gram) for i,r in enumerate(accepted) for s in (*walls,*accepted[:i]))
        walls.extend(accepted);out.update(all_walls=[enc(r) for r in walls],wall_count=len(walls),completed_height=h.to_json())
        output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
        print(json.dumps({'height':str(h),'new':len(accepted),'walls':len(walls)}),flush=True)
    rho=tuple(v for _ in range(3) for v in mv(HI,(O,O,O)));eps=G(F(1,2))
    while True:
        inside=(O,*(eps*x for x in rho))
        if pair(inside,gram,inside).sign()<0 and all(pair(inside,gram,r).sign()<0 for r in walls):break
        eps/=2
    witness=escape(walls,gram,20261013)
    out.update(completed=True,all_walls=[enc(r) for r in walls],wall_count=len(walls),completed_height=height.to_json(),
        strict_interior=enc(inside),hyperideal_witness=witness,status='EXACT_SPACELIKE_ESCAPE' if witness else 'VOLUME_UNDECIDED',
        seconds=round(time.monotonic()-start,3))
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps({'walls':len(walls),'status':out['status'],'seconds':out['seconds']}),flush=True)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--alpha',choices=('phi','2phi'),default='phi');p.add_argument('--height',type=int,default=30);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(PHI if a.alpha=='phi' else 2*PHI,G(a.height),a.output)
