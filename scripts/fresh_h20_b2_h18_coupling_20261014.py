#!/usr/bin/env python3
"""Exact direct-B2 coupling test for the fixed Vinberg H18 source.

One common strict escape covers every possible inducing-wall profile. No
claim is made about additional walls disjoint from the source face.
"""
import argparse
import hashlib
import json
import sys
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from coxeter_polytopes.exact.silver import Silver as K, silver_inertia

TAG='fresh_h20_b2_h18_coupling_20261014'
SOURCE='data/enumerations/vinberg_h18_mixture_facets_2026_10_13.json'
Z=K(0);ONE=K(1);RT=K(0,1)
def add(a,b):return[x+y for x,y in zip(a,b)]
def scale(a,c):return[x*c for x in a]
def mv(A,x):return[sum((a*b for a,b in zip(r,x)),Z)for r in A]
def dot(a,b,Q):return sum((x*y for x,y in zip(a,mv(Q,b))),Z)
def rref(A):
    A=[list(r)for r in A];p=[];k=0
    for j in range(len(A[0])):
        i=next((i for i in range(k,len(A))if A[i][j]!=Z),None)
        if i is None:continue
        A[k],A[i]=A[i],A[k];v=A[k][j];A[k]=[a/v for a in A[k]]
        for i in range(len(A)):
            if i!=k:
                v=A[i][j]
                if v!=Z:A[i]=[a-v*b for a,b in zip(A[i],A[k])]
        p.append(j);k+=1
        if k==len(A):break
    return A,p
def inverse(A):
    n=len(A);R,p=rref([r+[K(i==j)for j in range(n)]for i,r in enumerate(A)])
    assert p[:n]==list(range(n));return[r[n:]for r in R]
def enc(x):
    if isinstance(x,K):return x.to_json()
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
    if isinstance(x,(list,tuple)):return[enc(v)for v in x]
    return x
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def gate(g):
    assert g<=0
    if g==0:return'right'
    if g==-1:return'pi3'
    if g==-RT:return'pi4'
    assert g<=-2,g
    return'parallel'if g==-2 else'dotted'

def run():
    d=json.loads((ROOT/SOURCE).read_text());names=d['labels']
    G=[[K.from_json(a)for a in row]for row in d['gram']]
    assert len(G)==40
    _,basis=rref(G);assert len(basis)==19
    Q0=[[G[i][j]for j in basis]for i in basis];J=inverse(Q0)
    old=[mv(J,[G[j][i]for j in basis])for i in range(40)]
    assert all(dot(old[i],old[j],Q0)==G[i][j]for i in range(40)for j in range(40))
    Q=[r+[Z,Z]for r in Q0]+[[Z]*19+[K(2),-RT],[Z]*19+[-RT,K(2)]]
    inertia=silver_inertia(Q);assert(inertia.positive,inertia.negative,inertia.zero)==(20,1,0)
    support=[[Z]*19+[ONE,Z],[Z]*19+[Z,ONE]]
    weights=[[Z]*19+[ONE,RT/2],[Z]*19+[RT/2,ONE]]
    assert all(dot(weights[i],support[j],Q)==K(i==j)for i,j in product(range(2),repeat=2))
    assert all(dot(w,w,Q)==1 for w in weights)
    assert dot(weights[0],weights[1],Q)==RT/2
    # Complete inducing profile: positive residual square requires zero
    # contacts or exactly one pi3 contact to one B2 support node.
    assert 2+RT>2 and K(2)==2
    pi3_edges=[(i,j)for i,j in combinations(range(40),2)
               if G[i][j]<0 and 4*G[i][j]**2==G[i][i]*G[j][j]]
    forced=sorted(set(i for p in pi3_edges for i in p))
    eligible=sorted(set(range(40))-set(forced))
    assert eligible==[17]+list(range(23,36))+[38,39]
    assert len(forced)==24 and all(G[i][i]==2 for i in forced)
    forcing=[dict(source_wall=i,source_neighbor=next(j if i==k else k for k,j in pi3_edges if i in(k,j)))for i in forced]
    # A source pi3 pair has ambient products -1/sqrt2 (one special)
    # or -1/2+1, -1/2+1/sqrt2 (two special), all forbidden.
    assert -1 < -RT/2 < 0 and -K(F(1,2))+RT/2>0
    profiles=[]
    for i,r in enumerate(old):
        ordinary=scale(r,RT if G[i][i]==1 else ONE)+[Z,Z]
        prof=[ordinary]
        if i in eligible:
            y=scale(r,ONE if G[i][i]==1 else RT/2)+[Z,Z]
            prof += [add(y,scale(w,-1))for w in weights]
        assert all(dot(n,n,Q)==2 for n in prof)
        profiles.append(prof)

    # Literal exact point recovered from a bounded proposal; no LP dependency.
    x=[K(a,b)/23 for a,b in [(142,-56),(282,-84),(422,-112),(351,-84),
       (280,-56),(233,-42),(186,-28),(139,-14),(92,0),(45,14),(-2,28),
       (-49,42),(-96,56),(-119,56),(-142,56),(-165,56),(-188,56),
       (-211,56),(211,-56)]]
    assert dot(x,x,Q0)==K(F(3857,529),F(-2800,529))<0
    xp=[dot(x,r,Q0)for r in old];assert all(v<=0 for v in xp)
    assert all(xp[i]<=(-ONE if G[i][i]==1 else -RT)for i in eligible)
    h=[sum((r[j]for r in old),Z)for j in range(19)]
    assert dot(h,h,Q0)==K(-2952,-1280)<0
    assert dot(x,h,Q0)==K(-16,-40)<0
    assert all(dot(h,r,Q0)<0 for r in old)
    w=scale(add(weights[0],weights[1]),F(-1,2))
    assert dot(w,w,Q)==K(F(1,2),F(1,4))
    assert all(dot(w,z,Q)==-K(F(1,2),F(1,4))for z in weights)
    center=add(h+[Z,Z],w)
    escape=add(add(x,scale(h,F(1,1000)))+[Z,Z],w)
    bank=support+[n for prof in profiles for n in prof]
    assert len(bank)==74
    assert dot(center,center,Q)<0 and dot(escape,center,Q)<0
    assert dot(escape,escape,Q)==K(F(512876299,66125000),F(-33884339,6612500))>0
    hp=[dot(center,n,Q)for n in bank];ep=[dot(escape,n,Q)for n in bank]
    assert all(v<0 for v in hp+ep)

    # One coherent representative: all14 source norm-one walls special in
    # the same color, all26 norm-two walls ordinary. It preserves the old
    # non-quasi-arithmetic conditional witness.
    chosen=[1 if G[i][i]==1 else 0 for i in range(40)]
    normals=support+[profiles[i][a]for i,a in enumerate(chosen)]
    counts=Counter(gate(dot(a,b,Q))for a,b in combinations(normals,2))
    assert len(normals)==42 and sum(counts.values())==861
    facet_points=[]
    for i,n in enumerate(normals):
        p=add(center,scale(n,-dot(center,n,Q)/2))
        signs=[dot(p,m,Q)for m in normals]
        assert dot(p,p,Q)<0 and signs[i]==0 and all(v<0 for j,v in enumerate(signs)if j!=i)
        facet_points.append(p)
    assert all(dot(h+[Z,Z],n,Q)<0 for n in normals[2:])
    a,b=names.index('L1'),names.index('R10')
    assert G[a][b]**2==K(12,8)
    a,b=names.index('R23'),names.index('R24')
    assert G[a][b]==-6 and G[a][a]*G[b][b]-G[a][b]**2==-32
    paths=[SOURCE,f'scripts/{TAG}.py','src/coxeter_polytopes/exact/silver.py',
           'src/coxeter_polytopes/exact/matrix.py']
    return enc(dict(status='PASS_ALL_DIRECT_B2_H18_COUPLINGS_HAVE_COMMON_STRICT_ESCAPE',
      sources={p:sha(p)for p in paths},source_dimension=18,target_dimension=20,
      source_basis=basis,source_names=names,source_gram=G,form=Q,
      support=support,weights=weights,forced_ordinary=forcing,eligible_special=eligible,
      inducing_profiles=profiles,bank_normals=bank,profile_count=72,bank_count=74,
      source_point=x,source_point_norm=dot(x,x,Q0),source_pairings=xp,
      strict_center=center,center_norm=dot(center,center,Q),all_bank_center_pairings=hp,
      strict_escape=escape,escape_norm=dot(escape,escape,Q),escape_center_pairing=dot(escape,center,Q),
      all_bank_escape_pairings=ep,
      coherent_example=dict(walls=42,choices=chosen,normals=normals,pair_counts=dict(counts),
         actual_facet_points=facet_points,source_face_interior=h+[Z,Z],
         conditional_nonquasi=dict(cycle=K(12,8),rational_indefinite_pair=['R23','R24'],determinant=-32)),
      scope='Every direct reverse-B2 coupling of the fixed complete H18 source, using only its40 inducing facets and the2 support walls. All inducing contacts/all finite labels are covered. Extra walls not inducing source facets, changed source metrics and generalH20 completions remain outside. No finite-volume example or complete vertex census claimed.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status']);print(d['coherent_example']['pair_counts']);print(d['escape_norm'])
