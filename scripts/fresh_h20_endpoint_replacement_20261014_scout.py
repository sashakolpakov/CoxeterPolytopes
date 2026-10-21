#!/usr/bin/env python3
"""Bounded changed-endpoint/face-center experiment; numerical proposal only."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import numpy as np
from fresh_h20_balanced_20261014 import data

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_endpoint_replacement_20261014_scout'

def run():
    d=data();Z=np.array([v['point'] for v in d['vertices']]);V=d['vectors']
    def idx(a):return a[0]*56+a[1]*8+a[2]
    centers=[np.zeros(19)]+list(Z);r2=[4.]+[16.]*392
    rows=[dict(kind='central')]+[dict(kind='endpoint',vertex=i) for i in range(392)]
    for kind in ('nearest_edge','consecutive_triangle'):
        for factor in (0,1):
            for a in product(range(7),range(7),range(8)):
                b=list(a);b[factor]=(b[factor]+1)%7
                if kind=='nearest_edge':
                    ids=[idx(a),idx(b)];c=(Z[ids[0]]+Z[ids[1]])/2;q=5.
                else:
                    cc=list(a);cc[factor]=(cc[factor]+2)%7
                    ids=[idx(a),idx(b),idx(cc)]
                    c=(16*Z[ids[0]]+17*Z[ids[1]]+16*Z[ids[2]])/49;q=240/49
                centers.append(c);r2.append(q)
                rows.append(dict(kind=kind,factor=factor,vertices=ids))
    C=np.array(centers);Q=np.array(r2);N=np.sum(C*C,axis=1)
    assert len(C)==1961
    gaps=1-C@V.T;cos=gaps/np.sqrt(Q[:,None])
    def label(a):
        if a < -1e-8:return 'positive'
        if a > 1+1e-8:return 'dotted'
        if abs(a-1)<1e-8:return 'parallel'
        if abs(a)<1e-8:return '2'
        m=np.pi/np.arccos(np.clip(a,-1,1));k=int(round(m))
        return str(k) if k>=3 and abs(m-k)<1e-7 else 'bad'
    seed_bad=[(i,j,float(x)) for i,row in enumerate(cos) for j,x in enumerate(row) if label(x) in ('bad','positive')]
    counts=Counter();bad=[];minimum_pair_by_types={}
    for i in range(len(C)):
        dist=N[i]+N[:i]-2*C[:i]@C[i]
        cs=(dist-Q[i]-Q[:i])/(2*np.sqrt(Q[i]*Q[:i]))
        for j,x in enumerate(cs):
            typ=label(x);counts[typ]+=1
            kt=tuple(sorted((rows[i]['kind'],rows[j]['kind'])))
            minimum_pair_by_types[kt]=min(minimum_pair_by_types.get(kt,1e99),float(x))
            if typ in ('bad','positive'):bad.append(dict(i=i,j=j,cosine=float(x)))
    # Coverage probes include full source edges and explicit two/three-factor
    # face centers. These are not a complete 19D power-cell enumeration.
    probes=[]
    for lam in (0.,.15,.25,.4,.6,.8,1.):
        for i,z in enumerate(Z):probes.append((f'ray{lam}:{i}',lam*z))
    for factor,n in enumerate((7,7,8)):
        for i in range(n):
            for j in range(i+1,n):
                a=[0,0,0];b=[0,0,0];a[factor]=i;b[factor]=j
                probes.append((f'edge{factor}:{i},{j}',.99*(Z[idx(a)]+Z[idx(b)])/2))
    for a in product((1,2,3),(1,2,3),(1,2,3)):
        ids=[idx(v) for v in product(range(a[0]),range(a[1]),range(a[2]))]
        probes.append((f'facebox{a}',.99*np.mean(Z[ids],axis=0)))
    best=[]
    for name,x in probes:
        p=np.sum((C-x)**2,axis=1)-Q
        best.append(dict(name=name,minimum_power=float(p.min()),point=x.tolist(),nearest=int(p.argmin())))
    best.sort(key=lambda a:a['minimum_power'],reverse=True)
    # Frozen old .4-middle layer is optional; measure whether it cuts these gaps.
    old=json.loads((ROOT/'results/fresh_h20_balanced_integer_exact_20261014.json').read_text())
    L=d['cholesky'];middle=np.array([[float(__import__('fractions').Fraction(x)) for x in r['center']] for r in old['new_caps']])@L
    for p in best[:10]:
        x=np.array(p['point']);p['old_middle_min_power']=float((np.sum((middle-x)**2,axis=1)-4).min())
    paths=[f'scripts/{TAG}.py','scripts/fresh_h20_balanced_20261014.py','results/fresh_h20_balanced_integer_exact_20261014.json']
    return dict(status='NUMERICAL_FACE_HIERARCHY_SCOUT',sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        endpoint_radius=4,cap_count=len(C),rows=rows,seed_bad=seed_bad,pair_types=dict(counts),
        bad_count=len(bad),bad_examples=bad[:30],minimum_pair_by_types={str(k):v for k,v in minimum_pair_by_types.items()},
        probe_count=len(probes),uncovered_probes=sum(p['minimum_power']>1e-8 for p in best),worst_probes=best[:10],
        scope='Numerical bounded proposal. No finite-volume claim; probe list is not a complete coverage census.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run();f=ROOT/f'results/{TAG}.json'
    if a.verify:assert json.loads(f.read_text())==out
    else:f.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['pair_types'],out['bad_count'],out['uncovered_probes'])
    print(out['worst_probes'][:2])
