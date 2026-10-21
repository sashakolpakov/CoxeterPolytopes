#!/usr/bin/env python3
"""Independent exact trace-ellipsoid enumeration of cubic-nine root prefixes.

For a totally positive cap c, Tr(x^2/c)<=3 is necessary. Rational LDL
decomposition enumerates every integral coefficient vector in that ellipsoid;
the three exact embedding tests then give the requested alphabet. This uses
neither the producer's trace-dual boxes nor its interval pruning.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import groupby
import json
from math import isqrt
from pathlib import Path
import time

import fresh_cubic9_dihedral_independent_audit_20261013 as f
import fresh_cubic9_enumeration_audit_20261013 as census

box_alphabet=census.alphabet

@lru_cache(None)
def ellipsoid(cap):
    if any(f.sign(cap,e)<=0 for e in range(3)):
        assert cap==f.ZERO
        return (f.ZERO,)
    basis=(f.ONE,f.THETA,f.mul(f.THETA,f.THETA));reciprocal=f.inverse(cap)
    trace=lambda x:3*x[0]+6*x[2]
    Q=[[trace(f.mul(f.mul(a,b),reciprocal)) for b in basis] for a in basis]
    A=[r[:] for r in Q];U=[[F(i==j) for j in range(3)] for i in range(3)];delta=[]
    for i in range(3):
        d=A[i][i];assert d>0;delta.append(d)
        for j in range(i+1,3):U[i][j]=A[i][j]/d
        for j in range(i+1,3):
            for k in range(i+1,3):A[j][k]-=d*U[i][j]*U[i][k]
    assert all(Q[i][j]==sum(delta[k]*U[k][i]*U[k][j] for k in range(3)) for i in range(3) for j in range(3))
    x=[0]*3;out=[]
    def descend(i,remaining):
        if i<0:
            value=tuple(map(F,x));square=f.mul(value,value)
            if all(f.sign(f.sub(cap,square),e)>=0 for e in range(3)):out.append(value)
            return
        offset=sum(U[i][j]*x[j] for j in range(i+1,3));a,b=offset.numerator,offset.denominator
        radius2=remaining*b*b/delta[i];radius=isqrt(radius2.numerator//radius2.denominator)
        low=-((radius+a)//b);high=(radius-a)//b
        for value in range(low,high+1):
            x[i]=value;left=remaining-delta[i]*(value+offset)**2
            assert left>=0;descend(i-1,left)
    descend(2,F(3))
    assert len(set(out))==len(out)
    return tuple(out)

def alphabet(caps):
    if caps[0]==caps[1]==caps[2]:return ellipsoid(caps[0])
    return box_alphabet(caps)

def verify(source,output):
    started=time.monotonic();d=json.loads(source.read_text())
    assert d['completed'];dominant,_=f.binary_census();geometry=f.verify(d,dominant)
    # Shared consumer is the independent Cartesian enumerator, not producer code.
    census.alphabet=alphabet;census.block_list.cache_clear()
    alpha=f.decode(d['alpha']);height=f.decode(d['height_bound']);norms=f.vecdecode(d['norm_pool'])
    G=f.make_gram(alpha,d['blocks'],d['scalar']);walls=tuple(f.vecdecode(w) for w in d['all_walls'])
    controls=set()
    for q in norms:
        caps=(f.divide(f.mul(height,q),alpha),f.neg(f.divide(q,alpha)),f.neg(f.divide(q,alpha)))
        for x in alphabet(caps):
            if f.sign(x)>0 and f.integral(f.divide(f.scale(f.mul(alpha,x),2),q)):
                controls.add((f.divide(f.mul(alpha,f.mul(x,x)),q),q,x))
    saved={(f.decode(r['height']),f.decode(r['norm']),f.decode(r['time'])) for r in d['levels']}
    assert controls==saved and len(saved)==len(d['levels'])
    zero=[]
    for q in norms:
        roots,_=census.candidates(d,f.ZERO,q);zero.extend((q,r) for r in roots)
    assert set(zero)=={(f.decode(r['norm']),f.vecdecode(r['vector'])) for r in d['zero_height_dominant_census']}
    selected=list(walls[:d['dimension']]);keys={f.ray_key(r) for r in selected};counts=[]
    for _,packet in groupby(d['levels'],key=lambda row:row['height']):
        accepted=[];cov=tuple(f.mv(G,r) for r in selected)
        for row in packet:
            q=f.decode(row['norm']);x=f.decode(row['time']);roots,bc=census.candidates(d,x,q)
            assert len(roots)==row['candidate_roots'] and bc==row['block_vectors']
            added=[]
            for r in roots:
                if any(f.sign(f.dot(r,n))>0 for n in cov):continue
                key=f.ray_key(r)
                if key in keys:continue
                keys.add(key);added.append(r)
            assert len(added)==row['accepted'];accepted.extend(added);counts.append(len(roots))
        selected.extend(accepted)
        print(json.dumps({'controls_replayed':len(counts),'walls':len(selected)}),flush=True)
    assert set(selected)==set(walls)
    out={'status':'PASS_FULL_INDEPENDENT_PREFIX_ENUMERATION','source':str(source),
         'arithmetic':'Standard-library Fraction polynomial arithmetic; spatial trace ellipsoids with exact rational LDL; unpruned boxes for time; Cartesian block product.',
         'complete_time_norm_census':True,'complete_spatial_census':True,'complete_zero_census':True,
         'controls':len(controls),'walls':len(walls),'candidate_counts':counts,'geometry':geometry,
         'seconds':round(time.monotonic()-started,3)}
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps({k:out[k] for k in ('status','controls','walls','seconds')}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();verify(a.source,a.output)
