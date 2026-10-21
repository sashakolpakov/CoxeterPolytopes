#!/usr/bin/env python3
"""Exact S4 x S4 seed, exhaustive finite-leaf cap census, clique escapes.

No numerical search module is imported. All coefficients lie in Q(phi).
"""
import argparse, copy, json
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement, product
from pathlib import Path
from scripts.fresh_h9_symmetric_balanced_caps_verify_20261013 import (
    Z,O,T,add,sub,mul,sc,sign,conj,dot,mv)
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/enumerations/fresh_h9_split_leaf_caps_20261013_exact.json'

def inv(x):
    y=conj(x);n=mul(x,y);assert n[1]==0 and n[0]!=0
    return sc(y,1/n[0])
def div(x,y):return mul(x,inv(y))
def total(xs):
    r=Z
    for x in xs:r=add(r,x)
    return r
def enc(x):return [str(t) for t in x]
def ev(xs):return list(map(enc,xs))
A=add(T,sc(O,Q(1,2)));B=sc(div(A,add(sc(A,12),O)),-1)
E=(1,)*4+(-1,)*4

def features(x,y):
    d=tuple(sub(u,v) for u,v in zip(x,y))
    t=sc(add(add(x[0],y[0]),add(x[4],y[4])),Q(1,2))
    h=sc(sub(add(x[0],y[0]),add(x[4],y[4])),Q(1,2))
    return t,h,d,total(d),total(sc(v,e) for v,e in zip(d,E))

def pair(f,g):
    t,h,d,D,F=f;u,k,e,J,L=g
    r=sc(div(mul(t,u),A),Q(-1,2))
    r=sub(r,sc(div(mul(h,k),B),Q(1,2)))
    r=add(r,sc(dot(d,e),Q(1,8)))
    r=sub(r,div(mul(A,mul(D,J)),sc(add(O,sc(A,2)),32)))
    r=sub(r,div(mul(B,mul(F,L)),sc(add(O,sc(B,2)),32)))
    return r

def normal(f):
    t,h,d,D,F=f
    mean=sc(div(t,A),Q(1,16));bias=sc(div(h,B),Q(1,16))
    da=div(mul(A,D),sc(add(O,sc(A,2)),32))
    db=div(mul(B,F),sc(add(O,sc(B,2)),32))
    w=tuple(add(add(sc(v,Q(-1,8)),da),sc(db,e)) for v,e in zip(d,E))
    c=tuple(add(mean,sc(bias,e)) for e in E)
    return tuple(add(x,y) for x,y in zip(c,w))+tuple(sub(x,y) for x,y in zip(c,w))

def finite_pair(g):
    # 2cos(pi/m) in Q(phi) permits precisely m=2,3,5.
    return g in (Z,sc(O,-1),sc(T,-1)) or sign(add(g,sc(O,2)))<=0

def assignments():
    vals=(Z,O,T)
    for l,h in combinations_with_replacement(range(3),2):
        for k in ([4] if l==h else range(5)):
            for inds in combinations(range(4),k):
                ix=set(inds)
                yield tuple(vals[l if i in ix else h] for i in range(4)),tuple(vals[h if i in ix else l] for i in range(4))

def maximal_cliques(adj):
    answer=[]
    def go(r,p,x):
        if not p and not x:answer.append(tuple(r));return
        u=max((i for i in range(len(adj)) if ((p|x)>>i)&1),key=lambda i:(p&adj[i]).bit_count())
        choices=p&~adj[u]
        while choices:
            bit=choices&-choices;v=bit.bit_length()-1;choices^=bit
            go(r+[v],p&adj[v],x&adj[v]);p^=bit;x|=bit
    go([], (1<<len(adj))-1,0)
    return sorted(tuple(sorted(c)) for c in answer)

def construct():
    assert sign(A)>0 and sign(add(B,sc(O,Q(1,2))))>0 and sign(B)<0
    assert sign(sub(add(A,B),sc(O,2)))>0
    # Explicit invariant decomposition: six zero sum modes; sum directions
    # -128a,-128b; difference directions8 (six),8+16a,8+16b.
    assert sign(add(sc(O,8),sc(B,16)))>0
    G=[]
    for i in range(16):
        row=[]
        for j in range(16):
            if i==j:v=sc(O,2)
            elif i//8==j//8:v=Z
            else:v=sc(add(add(A,sc(B,E[i%8]*E[j%8])),sc(O,2 if i%8==j%8 else 0)),-1)
            row.append(v)
        G.append(tuple(row))
    G=tuple(G);S=(O,)*16;GS=mv(G,S)
    assert all(x==sc(A,-8) for x in GS) and dot(S,GS)==sc(A,-128)
    half=list(assignments());assert len(half)==51
    caps=[];fs=[];rs=[];cov=[]
    for (x1,y1),(x2,y2) in product(half,repeat=2):
        x=x1+x2;y=y1+y2;f=features(x,y)
        if pair(f,f)!=sc(O,2):continue
        r=normal(f);gr=mv(G,r)
        assert gr==tuple(sc(v,-1) for v in x+y)
        assert dot(r,gr)==sc(O,2) and sign(dot(S,gr))<0
        caps.append({'x':ev(x),'y':ev(y),'normal':ev(r)})
        fs.append(f);rs.append(r);cov.append(gr)
    assert len(caps)==192
    adj=[0]*len(caps);edges=[]
    for i,j in combinations(range(len(caps)),2):
        p=pair(fs[i],fs[j]);assert p==dot(rs[i],cov[j])
        if finite_pair(p):adj[i]|=1<<j;adj[j]|=1<<i;edges.append([i,j,enc(p)])
    cliques=maximal_cliques(adj)
    witnesses=[];witness_cov=[];cap_witness=[]
    # Each eta has two plus signs in each four-block.
    for a,c in product(combinations(range(4),2),combinations(range(4,8),2)):
        inds=set(a+c);eta=tuple(1 if i in inds else -1 for i in range(8))
        scalar=add(inv(sc(A,2)),sc(O,Q(1,1024)))
        w=tuple(add(scalar,sc(O,e)) for e in eta)+tuple(sub(scalar,sc(O,e)) for e in eta)
        gw=mv(G,w)
        assert sign(dot(w,gw))>0 and sign(dot(S,gw))<0
        assert all(sign(x)<0 for x in gw)
        witnesses.append({'plus_indices':sorted(inds),'coefficients':ev(w),'square':enc(dot(w,gw))})
        witness_cov.append(gw)
        cap_witness.append([sign(dot(r,gw))<0 for r in rs])
    proof=[]
    for clique in cliques:
        wi=next((i for i,row in enumerate(cap_witness) if all(row[j] for j in clique)),None)
        assert wi is not None
        proof.append({'caps':list(clique),'strict_escape':wi})
    assert len(edges)==4986 and len(cliques)==8478 and max(map(len,cliques))==6
    chosen=next(c for c in cliques if len(c)==6)
    return {'status':'PASS_EXACT_SPLIT_SEED_AND_FINITE_LEAF_BANK',
        'field':'phi^2=phi+1, phi>0','a':enc(A),'b':enc(B),
        'gram':list(map(ev,G)),'seed_signature':[9,1,6],
        'strict_center':ev(S),'strict_center_square':enc(dot(S,GS)),
        'leaf_pair_alphabet':ev((Z,O,T)),
        'ordered_assignments_exhausted':51**2,'cap_count':len(caps),'caps':caps,
        'compatibility_edges':edges,'maximal_clique_count':len(cliques),
        'maximum_compatible_caps':6,'selected_six_caps':list(chosen),
        'strict_spacelike_witnesses':witnesses,'all_maximal_clique_escapes':proof,
        'scope':'No finite-volume completion of this literal16-wall seed using only normalized norm2 normals with every leaf pairing in {0,-1,-phi}. Arbitrary dotted leaf pairings and other finite labels remain open. The selected22-wall Coxeter configuration has an exact strict spacelike escape.'}

def verify(saved, expected):
    assert saved==expected

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    expected=construct()
    if args.verify:verify(json.loads(OUT.read_text()),expected)
    else:OUT.write_text(json.dumps(expected,indent=2)+'\n')
    bad=copy.deepcopy(expected);bad['caps'][0]['normal'][0][0]='100'
    try:verify(bad,expected)
    except AssertionError:pass
    else:raise AssertionError('Forged cap accepted')
    print({'status':expected['status'],'caps':len(expected['caps']),
           'compatible_edges':len(expected['compatibility_edges']),
           'maximal_cliques':expected['maximal_clique_count'],'max_caps':expected['maximum_compatible_caps'],
           'forged_cap_rejected':True})
if __name__=='__main__':main()
