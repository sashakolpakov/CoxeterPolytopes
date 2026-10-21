#!/usr/bin/env python3
"""Changed-endpoint scout for coherent covers of one complete radial segment.

Composition, pair, and baseline gates use integers. Interval endpoint ordering
and the bounded DFS are numerical proposals, not an exact volume certificate.
"""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import time
import numpy as np
from fresh_h20_balanced_integer_obstruction_20261014 import factor

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_integer_segment_20261014'


def near_segment(N,Z,V):
    if Z <= 0:
        return N <= 112
    if Z >= V:
        return N-2*Z+V <= 112
    return N*V-Z*Z <= 112*V


def prepare():
    ff=[factor(n) for n in (6,6,7)]
    pp=[[r for r in rows if near_segment(r['norm28'],r['dot_vertex28'],vn)]
        for rows,_,vn in ff]
    ab=[]
    for i,a in enumerate(pp[0]):
        for j,b in enumerate(pp[1]):
            N=a['norm28']+b['norm28'];Z=a['dot_vertex28']+b['dot_vertex28']
            if near_segment(N,Z,ff[0][2]+ff[1][2]):ab.append((i,j,N,Z))
    rows=[]
    for i,j,N,Z in ab:
        for k,c in enumerate(pp[2]):
            n=N+c['norm28'];z=Z+c['dot_vertex28']
            if not near_segment(n,z,196*28):continue
            parts=[pp[0][i],pp[1][j],c]
            rows.append(dict(gaps=sum([list(r['gaps']) for r in parts],[]),
                             scalar=sum([r['scalar'] for r in parts],[]),
                             coeff28=sum([r['coeff28'] for r in parts],[]),
                             norm28=n,dot_vertex28=z))
    vertices=np.array([sum(parts,[]) for parts in product(*[f[1] for f in ff])],dtype=np.int64)
    N=np.array([r['norm28'] for r in rows],dtype=np.int64)
    Z=np.array([r['dot_vertex28'] for r in rows],dtype=np.int64)
    S=np.array([r['scalar'] for r in rows],dtype=np.int64)
    K=np.array([r['coeff28'] for r in rows],dtype=np.int64)
    central=(N==224)|(N==336)|(N>=448)
    distance=N[:,None]+5488-56*S@vertices.T
    projection=Z/5488
    perpendicular=N/28-Z.astype(float)**2/(28*5488)
    half=np.sqrt(np.maximum(0,4-perpendicular)/196)
    intervals=np.column_stack((projection-half,projection+half))
    assert len(rows)==8066
    return rows,N,S,K,central,distance,intervals


def solve(q,N,S,K,central,distance,intervals,limit):
    v=distance-28*(4+q)
    sq=v*v
    gate=(v==0)|((v>0)&((sq==4*q*28**2)|(sq==8*q*28**2)|
                         (sq==12*q*28**2)|(sq>=16*q*28**2)))
    live=central&gate.all(axis=1)
    ids=np.where(live)[0]
    a,b=intervals.T
    start=1/7;target=1-np.sqrt(q)/14
    cache={};nodes=0;deepest=start;best=[];finished=None
    def compatible(i):
        if i not in cache:
            dd=N+N[i]-2*S@K[i]
            cache[i]=(dd==0)|(dd==224)|(dd==336)|(dd>=448)
        return cache[i]
    def walk(bound,possible,chosen):
        nonlocal nodes,deepest,best,finished
        nodes+=1
        if nodes>limit:return None
        if bound>deepest:deepest=bound;best=chosen
        if bound>=target-1e-10:
            finished=chosen;return True
        pool=np.where(possible&(a<=bound+1e-10)&(b>bound+1e-10))[0]
        if len(pool)==0:return False
        pool=pool[np.argsort(-b[pool],kind='stable')]
        for i in pool:
            nxt=possible&compatible(int(i))&(b>b[i]+1e-10)
            verdict=walk(float(b[i]),nxt,chosen+[int(i)])
            if verdict is not False:return verdict
        return False
    began=time.time();verdict=walk(start,live.copy(),[])
    return dict(endpoint_radius_square=q,baseline_compatible_caps=len(ids),
                initial_extensions=int(np.sum(live&(a<=start+1e-10)&(b>start+1e-10))),
                verdict='NUMERICAL_SEGMENT_COVER'if verdict else('NODE_LIMIT'if verdict is None else'NO_SEGMENT_COVER_IN_NUMERICAL_DFS'),
                evaluated_states=nodes,limit=limit,target=target,
                furthest=deepest,selected=finished if verdict else best,
                compatibility_rows=len(cache),elapsed_seconds=time.time()-began)


def run(args):
    rows,N,S,K,central,distance,intervals=prepare()
    results=[]
    for q in args.radius_squares:
        r=solve(q,N,S,K,central,distance,intervals,args.nodes)
        results.append(r);print({k:v for k,v in r.items() if k!='selected'},flush=True)
    array_path=f'data/enumerations/{TAG}.npz'
    np.savez_compressed(ROOT/array_path,gaps=np.array([r['gaps']for r in rows],dtype=np.int8),
                        scalar=S,coeff28=K,norm28=N,intervals=intervals,
                        vertex_distances28=distance,central_compatible=central)
    paths=[f'scripts/{TAG}.py','scripts/fresh_h20_balanced_integer_obstruction_20261014.py',array_path]
    out=dict(status='CHANGED_ENDPOINT_SEGMENT_SCOUT',
             sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
             family='Radius2, nonnegative integer seed gaps7,7,8; all centers whose ball meets segment[0,z_000]',
             candidate_count=len(rows),results=results,
             selected_recipes={str(i):rows[i] for i in sorted({i for r in results for i in r['selected']})},
             scope='Exact integer enumeration and contact gates; numerical interval endpoints and bounded path-cover DFS. Covers only one segment, not the full19-dimensional ideal cross-section. No finite-volume or exact infeasibility claim.')
    (ROOT/f'results/{TAG}.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--nodes',type=int,default=3000)
    ap.add_argument('--radius-squares',type=int,nargs='+',default=[1,4,9,16,21])
    run(ap.parse_args())
