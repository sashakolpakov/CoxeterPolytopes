#!/usr/bin/env python3
"""Bounded exact graph census; not an all-rank extension theorem."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import fresh_sqrt5_forest_sensitive_weights_20261015_review as C
import fresh_large_affine_white_contact_counterexamples_20261017 as E

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_large_affine_white_continuation_scout_20261017'
Q=E.field(2)


def affine(kind,n):
    if kind=='A':es=[(i,(i+1)%n,3)for i in range(n)]
    elif kind=='C':es=[(i,i+1,4 if i in(0,n-2)else 3)for i in range(n-1)]
    elif kind=='B':
        es=[(i,i+1,4 if i==n-4 else 3)for i in range(n-3)]
        es +=[(0,n-2,3),(0,n-1,3)]
    else:
        assert kind=='D'
        es=[(i,i+1,3)for i in range(n-5)]
        es +=[(0,n-4,3),(0,n-3,3),(n-5,n-2,3),(n-5,n-1,3)]
    return E.graph(n,es),es


FINITE={C.code(g):name for name,g in C.models()}


def finite(g):
    ans=[]
    for h in C.pieces(g):
        if sum(map(len,h.values()))!=2*len(h)-2:return None
        name=FINITE.get(C.code(h))
        if name is None:return None
        ans.append(name)
    return sorted(ans)


def affine_models():
    out={}
    def add(name,g):out[C.code(g)]=name
    for kind,lo in [('B',4),('C',3),('D',5)]:
        for n in range(lo,10):add(kind+'_affine_order_'+str(n),affine(kind,n)[0])
    for name,arms in [('E6',(2,2,2)),('E7',(1,3,3)),('E8',(1,2,5))]:
        es=[];n=1
        for length in arms:
            p=0
            for _ in range(length):es.append((p,n,3));p=n;n+=1
        add(name+'_affine',E.graph(n,es))
    add('F4_affine',E.graph(5,[(0,1,3),(1,2,4),(2,3,3),(3,4,3)]))
    return out


AFFINE=affine_models()


def positive_kernel(g):
    ids=sorted(g);n=len(ids)
    G=[[Q(2)if i==j else -Q(1)if g[i].get(j)==3 else -Q(0,1)if g[i].get(j)==4 else Q()for j in ids]for i in ids]
    # A positive-definite (n-1)-block and its exact Schur zero recover the kernel.
    E.pd_pivots([r[:-1]for r in G[:-1]],Q)
    A=[G[i][:-1]+[-G[i][-1]]for i in range(n-1)]
    for k in range(n-1):
        t=A[k][k];assert t.sign()>0
        for j in range(k,n):A[k][j]/=t
        for i in range(k+1,n-1):
            t=A[i][k]
            for j in range(k,n):A[i][j]-=t*A[k][j]
    x=[Q()for _ in range(n-1)]
    for i in reversed(range(n-1)):x[i]=A[i][-1]-sum((A[i][j]*x[j]for j in range(i+1,n-1)),Q())
    x.append(Q(1))
    cert=E.minimal_affine_certificate(G,x,Q)
    return dict(roots=ids,edges=[[i,j,m]for i in ids for j,m in g[i].items()if i<j],**cert)


def small_affine_witness(g,z):
    # Every visited mask is connected and contains z.  A returned witness is
    # independently checked by a positive kernel and weighted-Laplacian identity.
    layer={1<<z};visited=0
    for size in range(1,10):
        nxt=set()
        for mask in sorted(layer):
            visited+=1;ids=[i for i in g if mask>>i&1];h=E.induced(g,ids)
            degrees=[len(h[i])for i in ids];name=None
            if size>=3 and all(d==2 for d in degrees)and all(m==3 for i in h for m in h[i].values()):
                name='A_affine_order_'+str(size)
            elif sum(degrees)==2*size-2:name=AFFINE.get(C.code(h))
            if name:
                return dict(type=name,visited_connected_subsets=visited,**positive_kernel(h))
            if size<9:
                for i in ids:
                    for j in g[i]:
                        if not mask>>j&1:nxt.add(mask|1<<j)
        layer=nxt
    return None


def census(kind,n):
    old,es=affine(kind,n);z=n;y=n+1;records=[]
    stats=Counter();witnesses=Counter()
    for size in(1,2):
        for support in combinations(range(n),size):
            for labels in product((3,4,5),repeat=size):
                stats['tested_contact_profiles']+=1
                g=E.graph(n+1,es+[(z,i,m)for i,m in zip(support,labels)])
                deletions={}
                for b in range(n):
                    f=finite(E.induced(g,set(g)-{b}))
                    if f is not None:deletions[b]=f
                if len(deletions)<2:continue
                stats['profiles_two_spherical_deletions']+=1
                markings=[bs for bs in combinations(deletions,2)if not set(support)<=set(bs)]
                if not markings:continue
                stats['white_contact_profiles']+=1;stats['white_contact_markings']+=len(markings)
                h=E.graph(n+2,es+[(z,i,m)for i,m in zip(support,labels)]+[(y,z,3)])
                extended={}
                for b in deletions:
                    f=finite(E.induced(h,set(h)-{b}))
                    if f is not None:extended[b]=f
                continued=[bs for bs in markings if all(b in extended for b in bs)]
                rec=dict(support=support,labels=labels,spherical_deletions=[[b,deletions[b]]for b in deletions],
                    white_markings=markings,extended_spherical_deletions=[[b,extended[b]]for b in extended],continued_markings=continued)
                if continued:
                    stats['continued_profiles']+=1;stats['continued_markings']+=len(continued)
                    witness=small_affine_witness(g,z)
                    rec['small_affine_witness']=witness
                    if witness:
                        stats['continued_markings_excluded_by_small_affine']+=len(continued)
                        witnesses[witness['type']]+=len(continued)
                    else:stats['open_continued_markings']+=len(continued)
                records.append(rec)
    stats['open_continued_markings']+=0
    assert stats['open_continued_markings']==0
    print('CASE',kind,n,dict(stats),flush=True)
    return dict(kind=kind,order=n,old_edges=es,counts=dict(stats),affine_witness_types=dict(sorted(witnesses.items())),records=records)


def run():
    cases=[census(kind,n)for n in(18,19,23)for kind in 'ABCD']
    paths=[Path(__file__),Path(C.__file__),Path(E.__file__)]
    totals=Counter()
    for c in cases:totals.update(c['counts'])
    return dict(status='PASS_BOUNDED_WHITE_CONTINUATION_SCOUT',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        domain=dict(orders=[18,19,23],affine_kinds=list('ABCD'),support_sizes=[1,2],contact_labels=[3,4,5],new_leaf_label=3),
        counts=dict(totals),cases=cases,
        scope='Complete bounded graph enumeration only for the saved orders, supports and labels. Every profile admitting two spherical deletions and the additional leaf continuation has an explicit induced affine witness of order<=9. This is not an all-rank theorem and does not independently establish the contact-support/label reductions for actual quadrilaterals. No ambient chamber realization or dimension bound is asserted.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(json.dumps(run()));f=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(f.read_text())
    else:f.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],flush=True)
