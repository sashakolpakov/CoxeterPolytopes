#!/usr/bin/env python3
"""Delete wall 3 before a two-height lift; search a bounded exact cap bank.

Numerics propose roots/escapes. Every archived algebraic claim is replayed
using Fraction arithmetic. No completeness claim for cap generation.
"""
import argparse
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np
from scripts.fresh_h8_alternate_caps_20261013 import (
    Z,O,P,parse,enc,add,sub,neg,scale,mul,div,sign,total,dot,mv,pair,
    matrix,inverse,rref,real,vec,mat,strict_escape,
)
from scripts.fresh_h4_blocks_independent_verify import determinant

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/enumerations/fresh_h8_twoheight_lift_20261013.json'
IDS=(0,1,2,4,5,6,7,8,9,10)
LEFT=(0,1,2,8)
RIGHT=tuple(i for i in range(10) if i not in LEFT)
Q=add(scale(O,3),scale(P,2))
LABELS=((2,Z),(3,O),(4,scale(O,2)),(5,add(O,P)),(6,scale(O,3)),(10,add(scale(O,2),P)))

def angle(b,n,m):
    if sign(b)>0:return None
    c=div(scale(mul(b,b),4),mul(n,m))
    if sign(sub(c,scale(O,4)))>0:return 'dotted'
    if c==scale(O,4):return 'parallel'
    return next((label for label,val in LABELS if val==c),None)

def inertia(A,conjugate=False):
    A=[list(row) for row in A];positive=negative=0;pivots=[]
    for k in range(len(A)):
        p=A[k][k];assert p!=Z
        val=(add((p[0],0),(p[1],-p[1])) if conjugate else p)
        s=sign(val);positive+=s>0;negative+=s<0;pivots.append(p)
        for i in range(k+1,len(A)):
            for j in range(k+1,len(A)):
                A[i][j]=sub(A[i][j],div(mul(A[i][k],A[k][j]),p))
    return (positive,negative),pivots

def seed():
    source=json.loads((ROOT/'data/known_examples/bugaenko_compact_h8.json').read_text())
    old=matrix(source['wall_gram'])
    G=[[old[i][j] for j in IDS] for i in IDS]
    for i in LEFT:
        for j in LEFT:G[i][j]=mul(Q,G[i][j])
    G[0][1]=G[1][0]=neg(Q)
    G[0][8]=G[8][0]=neg(Q)
    G[8][9]=G[9][8]=neg(Q)
    G=tuple(map(tuple,G));I=inverse(G)
    assert inertia(G)[0]==(9,1)
    assert inertia(G,True)[0]==(9,1)
    assert len(rref(G)[1])==10
    # Exact square identities deriving this Gram from the actual two-height lift.
    x0=scale(P,2);x9=div(add(scale(O,2),scale(P,8)),scale(O,11))
    n0=add(scale(O,2),x0);n9=add(scale(O,2),x9)
    assert div(scale(mul(P,P),4),scale(n0,2))==O
    assert div(scale(mul(x0,x9),4),mul(n0,n9))==O
    assert div(scale(mul(scale(P,2),scale(P,2)),4),scale(n9,2))==Q
    a=np.array([[real(x) for x in row] for row in G]);_,v=np.linalg.eigh(a)
    center=v[:,0];center*=np.sign(sum(center));assert np.all(center>0)
    from fractions import Fraction
    center=tuple((Fraction(round(float(x)*10**8),10**8),Z[1]) for x in center)
    assert all(sign(x)>0 for x in center) and sign(pair(center,I,center))<0
    assert all(angle(G[i][j],G[i][i],G[j][j]) not in (None,'parallel') for i in range(10) for j in range(i))
    return G,I,center,x0,x9

def run(max_support):
    G,I,center,x0,x9=seed();In=np.array([[real(x) for x in r] for r in I])
    caps=[];seen=set();stats=[]
    for typ in (0,1):
        norm=scale(O if typ==0 else Q,2)
        same=RIGHT if typ==0 else LEFT
        fixed=(O,P,add(O,P),scale(P,2))
        samevals=tuple(mul(Q,x) for x in fixed) if typ else fixed
        cross=(Q,mul(Q,P),scale(Q,2))
        alphabets=[samevals if i in same else cross for i in range(10)]
        for k in range(1,max_support+1):
            count=0;accepted=0
            for support in combinations(range(10),k):
                choices=list(product(*(range(len(alphabets[i])) for i in support)))
                pp=np.array([[real(alphabets[i][j]) for i,j in zip(support,row)] for row in choices])
                ns=np.einsum('bi,ij,bj->b',pp,In[np.ix_(support,support)],pp)
                for ci in np.flatnonzero(abs(ns-real(norm))<1e-7):
                    count+=1;p=[Z]*10
                    for i,j in zip(support,choices[ci]):p[i]=neg(alphabets[i][j])
                    p=tuple(p);c=mv(I,p)
                    if dot(p,c)!=norm or sign(dot(c,center))<=0:continue
                    labels=[angle(p[i],norm,G[i][i]) for i in range(10)]
                    if any(x in (None,'parallel') for x in labels):continue
                    if c in seen:continue
                    seen.add(c);accepted+=1
                    caps.append({'norm':enc(norm),'pairings':vec(p),'coefficients':vec(c),'seed_angles':labels,'support':k,'norm_type':typ})
            row={'norm_type':typ,'support':k,'numeric_proposals':count,'accepted':accepted};stats.append(row);print(row,flush=True)
    compatible=[];bad={}
    for i,a in enumerate(caps):
        c=tuple(map(parse,a['coefficients']));n=parse(a['norm'])
        for j,b in enumerate(caps[:i]):
            p=tuple(map(parse,b['pairings']));m=parse(b['norm']);relation=angle(dot(c,p),n,m)
            if relation not in (None,'parallel'):compatible.append([j,i,relation])
            else:bad[f'{j},{i}']=enc(dot(c,p))
    # Deterministic maximal cliques provide diverse finite configurations.
    adjacency=[set() for _ in caps]
    for i,j,_ in compatible:adjacency[i].add(j);adjacency[j].add(i)
    selected=[];used=set()
    orders=[sorted(range(len(caps)),key=lambda i:(-len(adjacency[i]),i))]
    rng=np.random.default_rng(20261013)
    orders.extend(rng.permutation(len(caps)).tolist() for _ in range(30))
    for order in orders:
        chosen=[]
        for i in order:
            if all(j in adjacency[i] for j in chosen):chosen.append(i)
        chosen=tuple(sorted(chosen))
        if chosen in used:continue
        used.add(chosen);A=[tuple(map(parse,caps[i]['coefficients'])) for i in chosen]
        gate=strict_escape(I,A,center,tries=300)
        selected.append({'caps':list(chosen),'facets':10+len(chosen),'volume_gate':gate})
        print({'selected':len(chosen),'volume':gate['status']},flush=True)
    out={'schema':'deleted-wall-two-height-H9-seed-v1','source':'data/known_examples/bugaenko_compact_h8.json','deleted_original_wall':3,'retained_original_walls':IDS,
         'lifted_original_walls':[0,9],'opposite_height_signs':True,'height_squares':vec([x0,x9]),
         'left_rescaled_original_walls':[IDS[i] for i in LEFT],'left_norm_factor':enc(Q),'seed_gram':mat(G),'seed_inverse':mat(I),
         'seed_signature':inertia(G)[0],'conjugate_signature':inertia(G,True)[0],'seed_determinant':enc(determinant(G)),
         'ldl_pivots':vec(inertia(G)[1]),'strict_center_y':vec(center),'strict_center_norm':enc(pair(center,I,center)),
         'dual_norms':vec([I[i][i] for i in range(10)]),'maximum_support':max_support,'caps':caps,'generation_counts':stats,
         'compatible_pairs':compatible,'other_pairings':bad,'selections':selected,
         'all_caps_escape':strict_escape(I,[tuple(map(parse,c['coefficients'])) for c in caps],center,tries=1000),
         'scope':'Bounded cap magnitudes and maximum support, numerical proposal filters, exact accepted roots/angles/escapes. No complete cap-family or full-host verdict.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print({'caps':len(caps),'configurations':len(selected),'output':str(OUT)},flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--support',type=int,default=5);args=p.parse_args();run(args.support)
