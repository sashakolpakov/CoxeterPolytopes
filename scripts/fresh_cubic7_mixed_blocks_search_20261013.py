#!/usr/bin/env python3
"""Exact bounded Vinberg search on <-2theta> + I2(7)^3 + A3."""
import argparse
from functools import lru_cache
from fractions import Fraction as F
from itertools import combinations_with_replacement, groupby, permutations
import json
from pathlib import Path
import time

from fresh_cubic7_dihedral_h9_20261013 import (
    K,T,Z,O,C,D,BETA,BLOCK,ALLOWED,inverse,mv,dot,pair,enc,
    encode_matrix,ray_key,finite_angle_counts,geometry,cubic7_inertia,le_all)
from fresh_cubic7_target_alphabet_20261013 import thin_alphabet

ROOT=Path(__file__).resolve().parents[1]
A3=tuple(tuple(K(2 if i==j else -1 if abs(i-j)==1 else 0) for j in range(3)) for i in range(3))
AINV=inverse(A3)
NORMS=(K(2),K(4),2*BETA,4*BETA)

def gram(blocks=3):
    n=2*blocks+3
    return tuple(tuple(-2*T if i==j==0 else
        BLOCK[(i-1)%2][(j-1)%2] if 1<=i<=2*blocks and 1<=j<=2*blocks and (i-1)//2==(j-1)//2 else
        A3[i-2*blocks-1][j-2*blocks-1] if i>2*blocks and j>2*blocks else Z
        for j in range(n+1)) for i in range(n+1))

def initial(blocks=3):
    n=2*blocks+4;out=[]
    for k in range(blocks):
        r=[Z]*n;r[1+2*k]=-O;out.append(tuple(r))
        r=[Z]*n;r[1+2*k]=D;r[2+2*k]=-D;out.append(tuple(r))
    for tail in ((-O,Z,Z),(Z,-O,Z),(O,Z,-O)):
        out.append((Z,*(Z for _ in range(2*blocks)),*tail))
    return out

@lru_cache(None)
def spatial_lists(target,q):
    weights=tuple(w for w in thin_alphabet((2*target,)*3)
                  if w.sign()>=0 and (2*w/q).is_integral)
    blocks=[];a3=[]
    for i,u in enumerate(weights):
        for w in weights[i:]:
            y=((2*u+C*w)/BETA,(C*u+2*w)/BETA)
            if all(t.is_integral for t in y):
                e=u*y[0]+w*y[1]
                if le_all(e,target):blocks.append((e,y,(u,w)))
            for v in weights:
                y=mv(AINV,(u,v,w))
                if not all(t.is_integral for t in y):continue
                e=dot(y,(u,v,w))
                if le_all(e,target):a3.append((e,y,(u,v,w)))
    return tuple(blocks),tuple(a3),len(weights)

def candidates(x,q,g,blocks=3):
    target=q+2*T*x*x;bs,aa,nw=spatial_lists(target,q)
    energies={}
    for e,y,w in aa:energies.setdefault(e,[]).append(y)
    out=[];imprimitive=0
    for ids in combinations_with_replacement(range(len(bs)),blocks):
        e=sum((bs[i][0] for i in ids),Z)
        if not le_all(e,target):continue
        for tail in energies.get(target-e,()):
            for order in sorted(set(permutations(ids))):
                r=(x,*(y for i in order for y in bs[i][1]),*tail)
                assert pair(r,g,r)==q and all(t.is_integral for t in r)
                assert all((2*t/q).is_integral for t in mv(g,r))
                if q in (K(4),4*BETA) and all((t/2).is_integral for t in r):
                    imprimitive+=1;continue
                out.append(r)
    return tuple(out),{'weights':nw,'block_vectors':len(bs),'a3_vectors':len(aa),'imprimitive':imprimitive}

def compatible(r,q,cov,norms):
    for n,pq in zip(cov,norms):
        p=dot(r,n)
        if p.sign()>0:return False
        c=4*p*p/(q*pq)
        if (c-4).sign()<0 and c not in ALLOWED:return False
    return True

def run(height,output,blocks=3):
    started=time.monotonic();g=gram(blocks);n=len(g)-1
    assert cubic7_inertia(g).negative==1
    assert all(cubic7_inertia(g,embedding=e).positive==n+1 for e in (1,2))
    zero=[]
    for q in NORMS:
        roots,counts=candidates(Z,q,g,blocks)
        zero.extend({'norm':q.to_json(),'vector':enc(r)} for r in roots)
        print(json.dumps({'zero_norm':str(q),'roots':len(roots),**counts}),flush=True)
    walls=initial(blocks);keys={ray_key(r) for r in walls}
    # Complete normalized primitive zero-height dominant census: one norm-two
    # and one norm-2beta orbit per B, plus the two B3 root-length orbits.
    assert len(zero)==2*blocks+2
    controls=sorted((2*T*x*x/q,q,x) for q in NORMS
        for x in thin_alphabet((height*q/(2*T),-q/(2*T),-q/(2*T)))
        if x.sign()>0 and (4*T*x/q).is_integral)
    out={'status':'RUNNING','host':'<-2theta>+B^'+str(blocks)+'+A3',
         'blocks':blocks,'dimension':n,'gram':encode_matrix(g),'alpha':(2*T).to_json(),
         'primitive_norm_pool':enc(NORMS),'zero_height_dominant_census':zero,
         'initial_walls':encode_matrix(walls),'height_bound':height,'control_count':len(controls),
         'levels':[],'completed':False,
         'scope':'Complete primitive norm pool and complete zero-height chamber; bounded full positive-root search. No full-host reflectivity or worldwide novelty claim.'}
    print(json.dumps({'phase':'start','controls':len(controls),'zero_roots':len(zero)}),flush=True)
    for h,packet in groupby(controls,key=lambda row:row[0]):
        accepted=[];cov=tuple(mv(g,r) for r in walls);norms=tuple(dot(r,c) for r,c in zip(walls,cov))
        for _,q,x in packet:
            roots,counts=candidates(x,q,g,blocks);kept=0;angle_rejected=0
            for r in roots:
                if any(dot(r,c).sign()>0 for c in cov):continue
                if not compatible(r,q,cov,norms):angle_rejected+=1;continue
                key=ray_key(r)
                if key in keys:continue
                keys.add(key);accepted.append(r);kept+=1
            out['levels'].append({'height':h.to_json(),'norm':q.to_json(),'time':x.to_json(),
                'candidates':len(roots),'accepted':kept,'acute_noncoxeter_rejected':angle_rejected,**counts})
        # A failure here would invalidate chronological chamber selection.
        finite_angle_counts(accepted,g)
        walls.extend(accepted)
        out.update(all_walls=encode_matrix(walls),wall_count=len(walls),completed_height=h.to_json())
        output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
        print(json.dumps({'height':str(h),'added':len(accepted),'walls':len(walls)}),flush=True)
    block_inside=mv(inverse(BLOCK),(O,K(2)));tail_inside=mv(AINV,(O,O,K(2)));eps=K(F(1,2))
    while True:
        inside=(O,*(eps*y for _ in range(blocks) for y in block_inside),*(eps*y for y in tail_inside))
        if pair(inside,g,inside).sign()<0 and all(pair(inside,g,r).sign()<0 for r in walls):break
        eps/=2
    out.update(status='COMPLETE_BOUNDED_MIXED_BLOCK_SEARCH',completed=True,
        completed_height=K(height).to_json(),strict_interior=enc(inside),
        finite_angle_counts=finite_angle_counts(walls,g))
    out['geometry']=geometry(walls,g,n,inside=inside)
    out['seconds']=round(time.monotonic()-started,3)
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'phase':'done','walls':len(walls),'geometry':out['geometry']['status'],'seconds':out['seconds']}),flush=True)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--height',type=int,default=30)
    p.add_argument('--blocks',type=int,choices=(1,3),default=3);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.height,a.output,a.blocks)
