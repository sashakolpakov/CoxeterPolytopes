#!/usr/bin/env python3
"""Exact coupled D4-supported Makarov P2 face, complete seed rays and polars.

This constructs and tests geometry, not a finite-volume assertion. The compact
source and all seed extreme rays are enumerated in Q(sqrt5).
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from fresh_h9_makarov_face_extension_20261013 import K, Z, O, PHI, kk, inv, null1, pair, ident, inertia

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_d4_prism_coupling_20261014'
ENGINE='scripts/fresh_h9_makarov_face_extension_20261013.py'

def enc(x):
    if isinstance(x,K):return x.enc()
    if isinstance(x,(list,tuple)):return [enc(a)for a in x]
    if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
    return x

def dotrow(v,Q):
    return [sum((v[i]*Q[i][j]for i in range(len(Q))),Z)for j in range(len(Q))]

def rays(Q,N):
    dim=len(Q);dual=[sum((v[j]for v in N),Z)for j in range(dim)]
    assert pair(dual,dual,Q).sign()<0
    rows=[dotrow(v,Q)for v in N];found={};tried=0;rankdrops=0
    for support in combinations(range(len(N)),dim-1):
        tried+=1;v=null1([rows[i]for i in support])
        if v is None:rankdrops+=1;continue
        s=pair(v,dual,Q).sign()
        if not s:continue
        if s>0:v=[-x for x in v]
        vals=[sum((a*b for a,b in zip(row,v)),Z)for row in rows]
        if any(x.sign()>0 for x in vals):continue
        first=next(x for x in v if x);scale=first if first.sign()>0 else -first
        v=[x/scale for x in v];key=tuple((x.a,x.b)for x in v)
        q=pair(v,v,Q)
        found[key]=dict(vector=v,norm=q,kind={-1:'finite',0:'ideal',1:'spacelike'}[q.sign()],
                        incident=[i for i,row in enumerate(rows)if not sum((a*b for a,b in zip(row,v)),Z)])
    out=list(found.values())
    return dict(positive_dual_sum=dual,dual_sum_norm=pair(dual,dual,Q),tested_supports=tried,
                rank_deficient_supports=rankdrops,counts=dict(Counter(v['kind']for v in out)),rays=out)

COS2={3:K(1)/4,4:K(1)/2,5:(K(3)+K(0,1))/8,6:K(3)/4,10:(K(5)+K(0,1))/8}

def relation(g,q,r):
    if not g:return dict(type='right',label=2)
    if g.sign()>0:return dict(type='obtuse',product=g)
    c=g*g/(q*r)
    if (c-O).sign()>=0:return dict(type='parallel'if c==O else'dotted',squared_weight=c)
    for m,v in COS2.items():
        if c==v:return dict(type='finite',label=m,squared_weight=c)
    return dict(type='forbidden_finite',squared_weight=c)

def strict_escape(Q,N,h,census):
    for row in census['rays']:
        if row['kind']!='spacelike':continue
        for k in range(1,10):
            v=[x+y/F(10**k) for x,y in zip(row['vector'],h)]
            if pair(v,v,Q).sign()<=0 or pair(v,h,Q).sign()>=0:continue
            v=[K(F(round(x.approx()*10**9),10**9)) for x in v]
            vals=[pair(v,n,Q) for n in N]
            if pair(v,v,Q).sign()>0 and pair(v,h,Q).sign()<0 and all(x.sign()<0 for x in vals):
                return dict(vector=v,norm=pair(v,v,Q),center_pairing=pair(v,h,Q),all_pairings=vals)
    raise AssertionError('No inward rational escape proposal certified')

def host_certificate(Q,N,omega):
    # D4 plus an outer weight is an odd integral unimodular overlattice.
    B=ident(10)
    for i in range(4):B[6+i][9]=omega[i]
    columns=[list(x) for x in zip(*B)];Bi=inv(B)
    G=[[pair(a,b,Q) for b in columns] for a in columns]
    integral=lambda x:(2*x.b).denominator==1 and (x.a-x.b).denominator==1
    assert all(integral(x) for row in G for x in row)
    A=[row[:] for row in G];det=O
    for k in range(10):
        pivot=A[k][k];assert pivot;det*=pivot
        for i in range(k+1,10):
            for j in range(k+1,10):A[i][j]-=A[i][k]*A[k][j]/pivot
    assert integral(det) and integral(O/det)
    conjugate=[[K(x.a,-x.b) for x in row] for row in G]
    assert inertia(G)==[9,1] and inertia(conjugate)==[10,0] and G[5][5]==O
    reflections=[]
    for n in N:
        c=[sum((x*y for x,y in zip(row,n)),Z) for row in Bi]
        dual=dotrow(c,G);q=pair(c,c,G)
        T=[[K(i==j)-2*c[i]*dual[j]/q for j in range(10)] for i in range(10)]
        reflections.append(dict(coordinates=c,matrix=T,integral=all(integral(x) for row in T for x in row)))
    assert all(r['integral'] for r in reflections)
    return dict(basis_columns=columns,gram=G,determinant=det,signature=[9,1],
                conjugate_signature=[10,0],norm_one_basis_index=5,reflections=reflections,
                scope='All listed reflections preserve this literal odd unimodular Z[phi] host. Bugaenko nonreflectivity excludes finite-volume completions within its integral reflection group; it does not exclude arbitrary nonintegral or larger-field caps.')

def run():
    J=ident(6)
    for i in range(5):J[i][i]=K(2)
    for i,x in enumerate([PHI,O,O,O,O]):J[i][i+1]=J[i+1][i]=-x
    Ji=inv(J)
    source=ident(6)+[[-PHI*Ji[i][5]for i in range(6)]]
    assert pair(source[5],source[5],J)==O
    assert pair(source[6],source[6],J)==K(2)
    assert pair(source[5],source[6],J)==-PHI
    source_rays=rays(J,source)
    assert source_rays['counts']=={'finite':10}
    D=ident(4)
    for i in range(4):D[i][i]=K(2)
    for j in [1,2,3]:D[0][j]=D[j][0]=-O
    Di=inv(D);omega=[Di[i][1]for i in range(4)]
    assert pair(omega,omega,D)==O
    Q=[row+[Z]*4 for row in J]+[[Z]*6+row for row in D]
    assert inertia(Q)==[9,1]
    N=[v+[Z]*4 for v in source]
    N[5]=source[5]+[-x for x in omega]
    N +=[[Z]*6+v for v in ident(4)]
    assert len(N)==11 and all(pair(n,n,Q)==K(2)for n in N)
    h=[-sum(row,Z)for row in Ji]+[-sum(row,Z)/10 for row in Di]
    assert pair(h,h,Q).sign()<0 and all(pair(h,n,Q).sign()<0 for n in N)
    G=[[pair(n,m,Q)for m in N]for n in N]
    pairs=[]
    for i,j in combinations(range(11),2):
        r=relation(G[i][j],K(2),K(2));assert r['type']in('finite','right')
        pairs.append(dict(pair=[i,j],**r))
    seed=rays(Q,N)
    facets=[]
    for i,n in enumerate(N):
        p=[a-pair(h,n,Q)*b/2 for a,b in zip(h,n)]
        assert pair(p,p,Q).sign()<0
        assert all(pair(p,m,Q).sign()==(0 if i==j else -1)for j,m in enumerate(N))
        facets.append(p)
    candidates=[]
    for ri,row in enumerate(seed['rays']):
        if row['kind']!='spacelike':continue
        n=row['vector'];q=row['norm'];rels=[relation(pair(n,m,Q),q,K(2))for m in N]
        bad=[i for i,r in enumerate(rels)if r['type']in('obtuse','forbidden_finite')]
        candidates.append(dict(ray=ri,normal=n,norm=q,relations=rels,bad_old_walls=bad,
                               center_pairing=pair(n,h,Q),passes_all_old_coxeter_pairs=not bad))
    passed=[i for i,row in enumerate(candidates)if row['passes_all_old_coxeter_pairs']]
    mutual=[]
    for i,j in combinations(passed,2):
        a,b=candidates[i],candidates[j]
        mutual.append(dict(candidates=[i,j],**relation(pair(a['normal'],b['normal'],Q),a['norm'],b['norm'])))
    assert len(passed)==4 and all(r['type']=='dotted'for r in mutual)
    extension=N+[candidates[i]['normal']for i in passed]
    assert all(pair(h,n,Q).sign()<0 for n in extension)
    extended=rays(Q,extension)
    print('extended complete rays',extended['counts'],flush=True)
    extension_facets=[]
    for i,n in enumerate(extension):
        p=[a-pair(h,n,Q)*b/pair(n,n,Q) for a,b in zip(h,n)]
        assert pair(p,p,Q).sign()<0
        assert all(pair(p,m,Q).sign()==(0 if i==j else -1) for j,m in enumerate(extension))
        extension_facets.append(p)
    escape=strict_escape(Q,extension,h,extended)
    host=host_certificate(Q,extension,omega)
    files=[ENGINE,f'scripts/{TAG}.py']
    return enc(dict(status='PASS_EXACT_COUPLED_D4_P2_SEED_AND_COMPLETE_RAY_DIAGNOSTIC',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in files},
        source_form=J,source_normals=source,source_complete_ray_census=source_rays,
        support_gram=D,support_weight=omega,form=Q,normals=N,gram=G,
        center=h,center_norm=pair(h,h,Q),actual_facet_points=facets,seed_pairs=pairs,
        complete_ray_census=seed,polar_candidates=candidates,passing_polar_indices=passed,
        passing_polar_mutual_relations=mutual,
        extension_normals=extension,extension_complete_ray_census=extended,
        extension_actual_facet_points=extension_facets,extension_strict_escape=escape,
        literal_arithmetic_host=host,
        scope='Complete exact compact P2 source, coupled eleven-wall H9 seed, and coherent fifteen-wall polar extension. Both ambient chambers have infinite volume, certified by complete ray enumeration and a strict rational escape for the extension. All fifteen reflections preserve a classified nonreflective odd unimodular arithmetic host; arbitrary nonintegral caps remain outside that exclusion.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['complete_ray_census']['counts'],'passing_polars',d['passing_polar_indices'],flush=True)
