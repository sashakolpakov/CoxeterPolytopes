#!/usr/bin/env python3
"""Exact verification of121172 actual H22 facets and a strict common escape.

All geometric arithmetic is integer/rational. NumPy integer products are
used only after explicit absolute overflow bounds. No scout is imported.
"""
import argparse,hashlib,json
from collections import Counter,deque
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
SOURCE='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
SYMS='data/enumerations/fresh_h20_global_interior_sections_20261013.json'
ORBIT='data/enumerations/fresh_h22_symmetry_orbit_20261013.npz'
SCOUT='data/enumerations/fresh_h22_symmetry_orbit_20261013.json'
OUT='results/fresh_h22_symmetry_orbit_verify_20261013.json'
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def classify(p,q):
    assert p<=0
    z=p*p/q
    if z>1:return 'dotted'
    if z==1:return 'parallel'
    allowed={F(0):'2',F(1,4):'3',F(1,2):'4',F(3,4):'6'}
    assert z in allowed,(str(p),str(q),str(z));return allowed[z]

def run(saved=None):
    raw=read(SOURCE)['coxeter_input'];Q=-np.array(raw['basis_gram'],dtype=np.int64);W=np.array(raw['wall_coordinates_as_columns_transposed'],dtype=np.int64)
    U=np.load(ROOT/ORBIT)['twice_old_projections'];assert U.shape==(60480,22)and U.dtype==np.int64
    syms=read(SYMS)['ambient_symmetry_generators'];gens=[np.array(x['matrix'],dtype=np.int64)for x in syms]
    # Crude absolute bounds cover every native-integer matrix product below.
    C=max(int(np.max(abs(x)))for x in[Q,W,U]+gens);assert 22**3*C**4<2**63
    h=np.sum(W,axis=0);h//=gcd(*map(int,h));assert h@Q@h==-14
    assert max(map(abs,map(int,h)))<=C
    assert Counter(map(int,W@Q@h))=={-1:42,-8:168}
    hc=h@Q;j=next(i for i in range(22)if hc[i])
    Z=np.array([[int(hc[j])*(k==i)-int(hc[i])*(k==j)for i in range(22)if i!=j]for k in range(22)],dtype=object)
    assert np.all(hc.astype(object)@Z==0)
    positive=[[F(x)for x in row]for row in (Z.T@Q.astype(object)@Z).tolist()];pivots=[]
    for i in range(21):
        p=positive[i][i];assert p>0;pivots.append(str(p))
        for k in range(i+1,21):
            for ell in range(i+1,21):positive[k][ell]-=positive[k][i]*positive[i][ell]/p
    G=W@Q@W.T;assert Counter(map(int,G.diagonal()))=={2:42,4:168}
    assert np.all(np.einsum('ij,jk,ik->i',U,Q,U)==-6)and np.all(U@Q@h==-33)
    lookup={tuple(map(int,u)):i for i,u in enumerate(U)};assert len(lookup)==len(U)
    edges=[]
    for gen,meta in zip(gens,syms,strict=True):
        assert np.array_equal(gen.T@Q@gen,Q)and np.array_equal(gen@h,h)
        assert np.array_equal(W@gen.T,W[meta['wall_permutation']])
        image=U@gen.T;edges.append([lookup[tuple(map(int,u))]for u in image])
    seen={0};queue=deque([0])
    while queue:
        i=queue.popleft()
        for permutation in edges:
            j=permutation[i]
            if j not in seen:seen.add(j);queue.append(j)
    assert len(seen)==60480
    old=Counter(classify(F(int(G[i,j])),F(int(G[i,i]*G[j,j])))for i in range(210)for j in range(i))
    ps=U[0]@Q@W.T;cap=Counter(classify(F(int(p),2),F(int(n)))for p,n in zip(ps,G.diagonal(),strict=True))
    same=U@(Q@U[0])+10;assert same[0]==4
    orbitpairs=Counter(classify(F(int(p),4),F(1))for p in same[1:])
    opposite=U@(Q@U[0])-10;assert np.max(opposite)<-4
    # Center cap(h,T), T² in(29/2,16], has norm in(1/2,2]. Every
    # orbit cap pairs at most-33/2+2sqrt10<-2, hence strictly dotted.
    assert F(29,2)**2>40
    den=10**8
    if saved is None:
        proposal=read(SCOUT)['runs'][0]['escape']['point'];assert len(proposal)==23
        v=[round(float(x)*den) for x in proposal]
        v[:22]=[int(x)+10000*int(y)for x,y in zip(v[:22],h,strict=True)]
    else:
        assert saved['strict_escape']['denominator']==str(den);v=list(map(int,saved['strict_escape']['numerator']))
    # Integer-object arithmetic avoids any dependence on native overflow in
    # the large rational escape. Both height signs reduce to one squared test.
    spatial=np.array(v[:22],dtype=object);cov=Q.astype(object)@spatial;z=int(v[-1]);q=int(spatial@cov)+z*z;hp=int(h.astype(object)@cov)
    oldp=W.astype(object)@cov;assert q>0 and hp<0 and max(oldp)<0
    up=U.astype(object)@cov;assert max(up)<0
    margins=[int(x)**2-10*z*z for x in up];assert min(margins)>0
    assert hp+4*abs(z)<0
    pair_total=121172*121171//2
    counts=Counter(old);N=len(U)
    for k,n in cap.items():counts[k]+=2*N*n
    for k,n in orbitpairs.items():counts[k]+=N*n
    counts['dotted']+=N*N+4*N+1+336
    counts['label_m']+=84
    assert sum(counts.values())==pair_total
    return {'status':'PASS_EXACT_121172_FACET_H22_FAMILY_WITH_STRICT_ESCAPE','source_sha256':{p:sha(p)for p in(SOURCE,SYMS,ORBIT)},
        'dimension':22,'facets':121172,'ambient_inertia':[22,1],'old_center_orthogonal_gram_pivots':pivots,'orbit_size':N,'signed_orbit_caps':2*N,'orbit_representative_twice_projection':U[0].tolist(),
        'orbit_cap_norm':'1','orbit_height_squared':'5/2','orbit_cap_center_pairing':'-33/2',
        'old_pair_counts':dict(old),'one_orbit_cap_old_pair_counts':dict(cap),'same_height_pair_counts_with_representative':dict(orbitpairs),
        'total_pair_counts_parameterized':dict(counts),'total_unordered_pairs':pair_total,
        'parameter':'Any single integer m>=3 for center caps(h,+/-sqrt(14+1/(2cos²(pi/m)))); label_m counts merge with ordinary labels when applicable.',
        'strict_escape':{'denominator':str(den),'numerator':list(map(str,v)),'norm':str(F(q,den*den)),'center_pairing':str(F(hp,den)),
            'maximum_old_pairing':str(F(int(max(oldp)),den)),'maximum_twice_orbit_spatial_pairing':str(F(int(max(up)),den)),
            'minimum_orbit_squared_margin_numerator':str(min(margins)),'center_cap_pairing_upper_bound':str(F(hp+4*abs(z),den))},
        'proof':'Generator closure and connectivity certify a single exact orbit. Every old/orbit and orbit/orbit pair transfers to one representative; opposite height signs are all dotted. The strict timelike center proves actual facets by projection. One rational strict future spacelike escape works for every m>=3 and proves infinite volume.',
        'scope':'An exact family of infinite-volume chambers inH22, not a new finite-volume dimension or a global novelty claim. No finite-volume completion is excluded.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();saved=read(OUT)if a.verify else None;out=run(saved)
    if a.verify:assert saved==out
    else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['total_unordered_pairs'])
