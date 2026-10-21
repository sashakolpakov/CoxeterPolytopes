#!/usr/bin/env python3
"""Exact bounded Vinberg search on <-2T> + B^2 + A4 + A1.

B has diagonal 2 and off-diagonal -T, T^3-3T-1=0. The safe complete
primitive norm-class pool follows from evenness and dual exponent 10 beta.
Only individually imprimitive vectors are removed from that pool.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations_with_replacement,groupby,permutations
import json
from pathlib import Path
import time
from fresh_cubic9_dihedral_h9_20261013 import (
    K,T,Z,O,C,D,BETA,BLOCK,ALLOWED,inverse,mv,dot,pair,enc,encode_matrix,
    ray_key,finite_angle_counts,geometry,cubic9_inertia,le_all)
from fresh_cubic9_alphabet_20261013 import alphabet

ROOT=Path(__file__).resolve().parents[1]
A4=tuple(tuple(K(2 if i==j else -1 if abs(i-j)==1 else 0) for j in range(4)) for i in range(4))
AINV=inverse(A4)
PREFIX_INV=tuple(inverse(tuple(row[:k] for row in A4[:k])) for k in range(1,5))
NORMS=tuple(K(n)*b for b in (O,BETA) for n in (2,4,10,20))
G=tuple(tuple(-2*T if i==j==0 else BLOCK[(i-1)%2][(j-1)%2]
    if 1<=i<=4 and 1<=j<=4 and (i-1)//2==(j-1)//2 else A4[i-5][j-5]
    if 5<=i<=8 and 5<=j<=8 else K(2) if i==j==9 else Z for j in range(10)) for i in range(10))

def initial():
    out=[]
    for k in range(2):
        r=[Z]*10;r[1+2*k]=-O;out.append(tuple(r))
        r=[Z]*10;r[1+2*k]=D;r[2+2*k]=-D;out.append(tuple(r))
    for i in range(5,10):
        r=[Z]*10;r[i]=-O;out.append(tuple(r))
    return out

@lru_cache(None)
def spatial_lists(target,q):
    # Cauchy: each Cartan covector w_i has w_i^2 <= 2*target in all
    # embeddings. Reflection integrality writes w_i=(q/2) k_i, k_i in O.
    weights=tuple(q*k/2 for k in alphabet((8*target/(q*q),)*3) if k.sign()>=0)
    bs=[]
    for i,u in enumerate(weights):
        for v in weights[i:]:
            y=((2*u+C*v)/BETA,(C*u+2*v)/BETA)
            if not all(x.is_integral for x in y):continue
            e=dot(y,(u,v))
            if le_all(e,target):bs.append((e,y))
    aa=[];prefix_counts=[0]*4
    def extend(ws):
        k=len(ws)
        for w in weights:
            ww=(*ws,w);low=dot(ww,mv(PREFIX_INV[k],ww))
            if not le_all(low,target):continue
            prefix_counts[k]+=1
            if k<3:extend(ww)
            else:
                y=mv(AINV,ww)
                if all(t.is_integral for t in y):aa.append((low,y))
    extend(())
    ss=tuple((2*z*z,z) for z in alphabet((target/2,)*3)
             if z.sign()>=0 and (4*z/q).is_integral)
    return tuple(bs),tuple(aa),ss,len(weights),tuple(prefix_counts)

def candidates(x,q):
    target=q+2*T*x*x;bs,aa,ss,nw,pc=spatial_lists(target,q)
    tails={}
    for e,y in aa:
        for s,z in ss:
            total=e+s
            if le_all(total,target):tails.setdefault(total,[]).append((*y,z))
    out=[];imprimitive=0
    for i,j in combinations_with_replacement(range(len(bs)),2):
        leftover=target-bs[i][0]-bs[j][0]
        for tail in tails.get(leftover,()):
            for a,b in sorted(set(permutations((i,j)))):
                r=(x,*bs[a][1],*bs[b][1],*tail)
                assert pair(r,G,r)==q and all(t.is_integral for t in r)
                assert all((2*t/q).is_integral for t in mv(G,r))
                if all((t/2).is_integral for t in r):imprimitive+=1;continue
                out.append(r)
    return tuple(out),{'weights':nw,'block_vectors':len(bs),'a4_vectors':len(aa),
        'scalar_vectors':len(ss),'a4_prefix_counts':pc,'imprimitive':imprimitive}

def compatible(r,q,cov,norms):
    for n,pq in zip(cov,norms):
        p=dot(r,n)
        if p.sign()>0:return False
        c=4*p*p/(q*pq)
        if (c-4).sign()<0 and c not in ALLOWED:return False
    return True

def run(height,output):
    start=time.monotonic();zero=[]
    assert cubic9_inertia(G).negative==1
    assert all(cubic9_inertia(G,embedding=e).positive==10 for e in (1,2))
    for q in NORMS:
        roots,counts=candidates(Z,q)
        zero.extend({'norm':q.to_json(),'vector':enc(r)} for r in roots)
        print(json.dumps({'zero_norm':str(q),'roots':len(roots),**counts}),flush=True)
    walls=initial();keys={ray_key(r) for r in walls}
    # The full zero census must be checked before using these simple walls.
    assert len(zero)==6,[(r['norm'],r['vector']) for r in zero]
    controls=sorted((2*T*x*x/q,q,x) for q in NORMS
        for x in alphabet((height*q/(2*T),-q/(2*T),-q/(2*T)))
        if x.sign()>0 and (4*T*x/q).is_integral)
    out={'status':'RUNNING','field':'Q(2cos(pi/9))','host':'<-2T>+B^2+A4+A1',
        'dimension':9,'gram':encode_matrix(G),'alpha':(2*T).to_json(),
        'safe_primitive_norm_pool':enc(NORMS),'zero_height_dominant_census':zero,
        'initial_walls':encode_matrix(walls),'height_bound':height,'control_count':len(controls),
        'levels':[],'completed':False,
        'scope':'Complete safe primitive norm pool; exact bounded full positive-root search. No worldwide novelty or full-host reflectivity claim.'}
    print(json.dumps({'phase':'start','controls':len(controls),'zero_roots':len(zero)}),flush=True)
    for h,packet in groupby(controls,key=lambda row:row[0]):
        accepted=[];cov=tuple(mv(G,r) for r in walls);norms=tuple(dot(r,c) for r,c in zip(walls,cov))
        for _,q,x in packet:
            roots,counts=candidates(x,q);kept=0;rejected=0
            for r in roots:
                if any(dot(r,c).sign()>0 for c in cov):continue
                if not compatible(r,q,cov,norms):rejected+=1;continue
                key=ray_key(r)
                if key in keys:continue
                keys.add(key);accepted.append(r);kept+=1
            out['levels'].append({'height':h.to_json(),'norm':q.to_json(),'time':x.to_json(),
                'candidates':len(roots),'accepted':kept,'acute_noncoxeter_rejected':rejected,**counts})
        finite_angle_counts(accepted,G);walls.extend(accepted)
        out.update(all_walls=encode_matrix(walls),wall_count=len(walls),completed_height=h.to_json())
        output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
        print(json.dumps({'height':str(h),'added':len(accepted),'walls':len(walls)}),flush=True)
    base=mv(inverse(BLOCK),(O,K(2)));tail=mv(AINV,(O,O,O,O));eps=K(F(1,2))
    while True:
        inside=(O,*(eps*y for _ in range(2) for y in base),*(eps*y for y in tail),eps)
        if pair(inside,G,inside).sign()<0 and all(pair(inside,G,r).sign()<0 for r in walls):break
        eps/=2
    out.update(status='COMPLETE_BOUNDED_MIXED_BLOCK_SEARCH',completed=True,
        completed_height=K(height).to_json(),strict_interior=enc(inside),
        finite_angle_counts=finite_angle_counts(walls,G))
    out['geometry']=geometry(walls,G,9,inside=inside)
    out['seconds']=round(time.monotonic()-start,3)
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'phase':'done','walls':len(walls),'geometry':out['geometry']['status'],'seconds':out['seconds']}),flush=True)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--height',type=int,default=30)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.height,a.output)
