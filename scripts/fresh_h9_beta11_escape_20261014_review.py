#!/usr/bin/env python3
"""Standalone Fraction review of the two finite-prefix volume witnesses.

No producer imports and no full cone or vertex census repetition.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_beta11_escape_20261014_review'
Z,O=(F(0),F(0)),(F(1),F(0))
def dec(x):return F(x['a']),F(x['b'])
def enc(x):return list(map(str,x))
def add(x,y):return x[0]+y[0],x[1]+y[1]
def mul(x,y):return x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
def conj(x):return x[0]+x[1],-x[1]
def sign(x):
    a,b=2*x[0]+x[1],x[1]
    sg=lambda q:(q>0)-(q<0)
    if not a:return sg(b)
    if not b or sg(a)==sg(b):return sg(a)
    return sg(a)*sg(a*a-5*b*b)
def total(xs):
    s=Z
    for x in xs:s=add(s,x)
    return s
def determinant(A):
    M=[list(row)for row in A];answer=F(1)
    for k in range(len(M)):
        p=next((i for i in range(k,len(M))if M[i][k]),None)
        if p is None:return F(0)
        if p!=k:M[k],M[p]=M[p],M[k];answer=-answer
        q=M[k][k];answer*=q
        for i in range(k+1,len(M)):
            t=M[i][k]/q
            for j in range(k+1,len(M)):M[i][j]-=t*M[k][j]
    return answer
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()


def run():
    seen={}
    def read(p):
        if p in seen:return seen[p]
        d=json.loads((ROOT/p).read_text());seen[p]=d
        for s,h in d.get('sources',{}).items():
            assert sha(s)==h,s
            if s.endswith('.json'):read(s)
        return d
    cases=[(100,'results/fresh_h9_beta11_vinberg_20261014.json',
                 'results/fresh_h9_beta11_geometry_20261014.json',21),
           (300,'results/fresh_h9_beta11_vinberg_h300_20261014.json',
                 'results/fresh_h9_beta11_geometry_h300_20261014.json',35)]
    records=[];data={}
    for height,vpath,gpath,wallcount in cases:
        v,g=read(vpath),read(gpath)
        G=[[dec(x)for x in row]for row in v['form']]
        assert len(G)==10 and all(len(row)==10 for row in G)
        assert G==list(map(list,zip(*G)))
        assert G[0][0]==(F(0),F(-2)) and G[9][9]==(F(2),F(0))
        assert G[0][9]==G[9][0]==(F(-1),F(0))
        assert all(G[i][j]==Z for i in (0,9)for j in range(1,9))
        assert all(x[1]==0 for row in G[1:9]for x in row[1:9])
        E8=[[x[0]for x in row[1:9]]for row in G[1:9]]
        minors=[determinant([row[:k]for row in E8[:k]])for k in range(1,9)]
        assert all(x>0 for x in minors)and minors[-1]==1
        binary_det=(F(-1),F(-4))
        assert sign(binary_det)<0 and sign(conj(binary_det))>0
        assert sign(conj(G[0][0]))>0
        # Thus real signature(9,1), conjugate signature(10,0).
        entries=[(i,j,x)for i,row in enumerate(G)for j,x in enumerate(row)if x!=Z]
        def pair(a,b):return total(mul(mul(a[i],x),b[j])for i,j,x in entries)
        N=[list(map(dec,r))for r in v['roots']]
        assert len(N)==wallcount==g['wall_count']
        norms=[pair(r,r)for r in N]
        assert all(q in((F(2),F(0)),(F(6),F(2)))for q in norms)
        assert all(sign(q)>0 for q in norms)
        H=list(map(dec,v['basepoint']))
        assert sign(pair(H,H))<0
        P=list(map(dec,g['strict_center']['point']))
        E=list(map(dec,g['strict_escape']['point']))
        qP,qE=pair(P,P),pair(E,E)
        assert qP==dec(g['strict_center']['norm'])and sign(qP)<0
        assert qE==dec(g['strict_escape']['norm'])and sign(qE)>0
        assert sign(pair(P,H))<0
        p_pairs=[pair(P,r)for r in N];e_pairs=[pair(E,r)for r in N]
        assert all(sign(x)<0 for x in p_pairs+e_pairs)
        assert e_pairs==list(map(dec,g['strict_escape']['pairings']))
        ep=pair(E,P)
        assert ep==dec(g['strict_escape']['center_pairing'])and sign(ep)<0
        Q=list(map(dec,v['strict_center']))
        qQ=pair(Q,Q);q_pairs=[pair(Q,r)for r in N]
        assert sign(qQ)<0 and sign(pair(Q,H))<0
        assert all(sign(x)<0 for x in q_pairs)
        records.append(dict(height=height,walls=wallcount,real_signature=[9,1],
          conjugate_signature=[10,0],E8_leading_minors=list(map(str,minors)),
          geometry_center_norm=enc(qP),geometry_center_pairings=list(map(enc,p_pairs)),
          vinberg_center_norm=enc(qQ),vinberg_center_pairings=list(map(enc,q_pairs)),
          strict_escape_norm=enc(qE),strict_escape_pairings=list(map(enc,e_pairs)),
          escape_center_pairing=enc(ep)))
        data[height]=(G,N,E,entries)
    G,N,E100,entries=data[100]
    G300,N300,E300,_=data[300]
    assert G300==G and N300[:len(N)]==N
    old_pairs=[total(mul(mul(E100[i],x),r[j])for i,j,x in entries)for r in N300]
    positive=[i for i,x in enumerate(old_pairs)if sign(x)>0]
    zero=[i for i,x in enumerate(old_pairs)if sign(x)==0]
    assert positive==[25] and not zero
    sources=[p for _,vp,gp,_ in cases for p in(vp,gp)]+[f'scripts/{TAG}.py']
    return dict(status='PASS_INDEPENDENT_HEIGHT100_300_STRICT_ESCAPE_REVIEW',
        sources={p:sha(p)for p in sources},recursive_manifests=len(seen),
        arithmetic='Standalone Fraction pairs in1,phi; all products/signs reconstructed from literal form and vectors. No producer imports.',
        cases=records,checked_escape_pairings=56,checked_center_pairings=112,
        old100_escape_against35=dict(pairings=list(map(enc,old_pairs)),positive_indices=positive,zero_indices=zero),
        finite_vertex_preservation_proof_review='PASS: for a finite vertex with complete independent PD acute support B, write H=a v+sum b_i n_i, a>0,b=B^-1(H,n)<=0. Any later actual positive-height simple root has(r,H)<0 and(r,n_i)<=0, hence(r,v)<0. Complete zero sector and true accepted-prefix premises are required; no near-height gate is needed.',
        vertex_map_corollary='A future-preserving integral isometry mapping the complete incident outward root rays positively onto those at another certified full-chamber finite vertex preserves the full mirror arrangement and the chamber. Equal abstract diagrams alone do not suffice.',
        scope='Independent infinite-volume witnesses for the finite21-wall/35-wall prefixes only. Full cone/extreme-ray counts and individual finite vertices are inherited, not recounted. No reflectivity or eventual-volume conclusion.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
