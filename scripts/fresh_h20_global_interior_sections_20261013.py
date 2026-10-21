#!/usr/bin/env python3
"""Balanced interior H20 sections of the Borcherds H21 polytope.

Certifies ambient symmetries and actual forbidden-ridge witnesses. No floating
point geometric decision enters the certificates. This is a bounded section
family, not a nonexistence result for other H20 sections or polytopes.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import json
from math import gcd
from pathlib import Path
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
OUTPUT=ROOT/'data/enumerations/fresh_h20_global_interior_sections_20261013.json'
ALLOWED={F(0),F(1,4),F(1,2),F(3,4)}

def primitive(v):
    v=tuple(map(int,v));d=gcd(*v)
    assert d>0
    return np.array([x//d for x in v],dtype=object)

def load():
    raw=json.loads(SOURCE.read_text())['coxeter_input']
    Q=-np.array([[int(x) for x in row] for row in raw['basis_gram']],dtype=object)
    W=np.array([[int(x) for x in row] for row in raw['wall_coordinates_as_columns_transposed']],dtype=object)
    G=W@Q@W.T;h=primitive(W.sum(axis=0));cp=W@Q@h
    assert int(h@Q@h)==-14 and set(map(int,cp))=={-1,-8}
    return Q,W,G,h

def generators(Q,W,G):
    import pynauty as pn
    short=[i for i in range(210) if G[i,i]==2]
    graph=pn.Graph(42,False,adjacency_dict={i:[j for j in range(42) if G[short[i],short[j]]==-1] for i in range(42)})
    aut=pn.autgrp(graph);assert int(aut[1])*10**int(aut[2])==241920
    basis=sp.Matrix(W[short].T).rref()[1]
    B=sp.Matrix(W[[short[k] for k in basis]].T);Bi=B.inv()
    lookup={tuple(v):i for i,v in enumerate(W)};result=[]
    for g in aut[0]:
        A=sp.Matrix(W[[short[g[k]] for k in basis]].T)*Bi
        assert all(x.q==1 for x in A)
        A=np.array(A,dtype=object)
        assert np.array_equal(A.T@Q@A,Q)
        perm=[lookup[tuple(row)] for row in W@A.T]
        assert len(set(perm))==210
        result.append({'matrix':[[int(x) for x in row] for row in A], 'wall_permutation':perm})
    return result

def pair_orbits(perms):
    remaining=set(combinations(range(210),2));result=[]
    while remaining:
        root=min(remaining);seen={root};queue=[root]
        while queue:
            i,j=queue.pop()
            for g in perms:
                nxt=tuple(sorted((g[i],g[j])))
                if nxt not in seen:seen.add(nxt);queue.append(nxt)
        assert seen<=remaining
        remaining-=seen;result.append(sorted(seen))
    return result

def normals(Q,W,h,e):
    qe=int(e@Q@e);assert qe>0 and int(h@Q@e)==0
    V=qe*W-np.outer(W@Q@e,e)
    lookup={};removed=Counter()
    for i,row in enumerate(V):
        assert int(row@Q@h)<0
        if row@Q@row<=0:
            removed['nonpositive_norm']+=1;continue
        key=tuple(primitive(row))
        if key in lookup:lookup[key].append(i)
        else:lookup[key]=[i]
    R=np.array(list(lookup),dtype=object)
    return R,list(lookup.values()),removed

def actual_bad_ridge(Q,R,h):
    H=R@Q@R.T;cp=R@Q@h;counts=Counter();witness=None
    for i,j in combinations(range(len(R)),2):
        aa,bb,ab=map(int,(H[i,i],H[j,j],H[i,j]));D=aa*bb-ab*ab
        if D<=0:
            if ab>0:counts['positive_nonelliptic']+=1
            else:counts['nonelliptic']+=1
            continue
        ratio=F(ab*ab,aa*bb)
        if ab<=0 and ratio in ALLOWED:counts['allowed_elliptic']+=1;continue
        counts['forbidden_elliptic']+=1
        if witness is not None:continue
        p=primitive(D*h-(bb*int(cp[i])-ab*int(cp[j]))*R[i]-(aa*int(cp[j])-ab*int(cp[i]))*R[j])
        vals=R@Q@p
        if any(x>0 for x in vals):continue
        zeros=[k for k,x in enumerate(vals) if x==0]
        if zeros!=[i,j]:continue
        assert p@Q@p<0 and p@Q@h<0
        witness={'normal_indices':[i,j], 'normal_ambient_wall_labels_one_based':None,
                 'gram':[[aa,ab],[ab,bb]],'cosine_squared':str(ratio),
                 'outward_pairing_positive':ab>0,'point':list(map(int,p)),
                 'point_square':int(p@Q@p),'center_pairing':int(p@Q@h),
                 'maximum_other_wall_pairing':max(int(x) for x in vals if x!=0),
                 'zero_normal_indices':zeros}
    return dict(counts),witness

def build():
    Q,W,G,h=load();gens=generators(Q,W,G);orbits=pair_orbits([x['wall_permutation'] for x in gens]);rows=[]
    cp=W@Q@h
    for orbit in orbits:
        i,j=orbit[0];e=primitive(int(cp[j])*W[i]-int(cp[i])*W[j])
        R,labels,removed=normals(Q,W,h,e);counts,witness=actual_bad_ridge(Q,R,h)
        if witness is not None:witness['normal_ambient_wall_labels_one_based']=[[x+1 for x in labels[k]] for k in witness['normal_indices']]
        row={'pair_one_based':[i+1,j+1],'ambient_pair_gram':[[int(G[i,i]),int(G[i,j])],[int(G[i,j]),int(G[j,j])]],
             'orbit_size':len(orbit),'section_normal':list(map(int,e)),'section_normal_square':int(e@Q@e),
             'positive_projected_normal_rays':len(R),'discarded':dict(removed),
             'projected_duplicate_count':sum(len(x)-1 for x in labels),'pair_gate_counts':counts,
             'actual_forbidden_ridge':witness,'status':'NON_COXETER_ACTUAL_RIDGE' if witness else 'UNRESOLVED'}
        rows.append(row);print(row['pair_one_based'],len(orbit),len(R),row['status'],witness['cosine_squared'] if witness else '',flush=True)
    out={'schema':'balanced-interior-H20-sections-v1','source':str(SOURCE.relative_to(ROOT)),
         'center':list(map(int,h)),'center_square':-14,'ambient_symmetry_generators':gens,
         'pair_orbit_count':len(orbits),'pair_count':sum(map(len,orbits)),
         'finite_volume_reason':'Each section contains a strict interior timelike point and is an intersection of a finite-volume finite-vertex Klein polytope with a linear hyperplane. Its ideal vertices are inherited from the ambient polytope; finite/ideal simplex triangulation gives finite volume.',
         'sections':rows,'all_representatives_non_coxeter':all(x['actual_forbidden_ridge'] for x in rows)}
    OUTPUT.write_text(json.dumps(out,indent=2)+'\n');verify(out)

def verify(out):
    Q,W,G,h=load();assert out['source']==str(SOURCE.relative_to(ROOT))
    assert out['center']==list(map(int,h)) and int(h@Q@h)==out['center_square']<0
    perms=[]
    for row in out['ambient_symmetry_generators']:
        A=np.array(row['matrix'],dtype=object);g=row['wall_permutation']
        assert sorted(g)==list(range(210)) and np.array_equal(A.T@Q@A,Q)
        assert np.array_equal(W@A.T,W[g]) and np.array_equal(A@h,h)
        perms.append(g)
    orbits=pair_orbits(perms);assert len(orbits)==out['pair_orbit_count']
    assert sum(map(len,orbits))==out['pair_count']==21945
    assert len(out['sections'])==len(orbits)
    cp=W@Q@h
    for orbit,row in zip(orbits,out['sections']):
        i,j=orbit[0];assert row['pair_one_based']==[i+1,j+1] and row['orbit_size']==len(orbit)
        e=primitive(int(cp[j])*W[i]-int(cp[i])*W[j]);assert list(map(int,e))==row['section_normal']
        assert int(e@Q@e)==row['section_normal_square']>0 and int(e@Q@h)==0
        R,labels,removed=normals(Q,W,h,e)
        assert len(R)==row['positive_projected_normal_rays'] and dict(removed)==row['discarded']
        counts,w=actual_bad_ridge(Q,R,h)
        assert counts==row['pair_gate_counts']
        if w is not None:w['normal_ambient_wall_labels_one_based']=[[x+1 for x in labels[k]] for k in w['normal_indices']]
        assert w==row['actual_forbidden_ridge']
        assert row['status']==('NON_COXETER_ACTUAL_RIDGE' if w else 'UNRESOLVED')
    assert out['all_representatives_non_coxeter']==all(x['actual_forbidden_ridge'] is not None for x in out['sections'])
    print('PASS: exact source, isometries, all 21945 pair choices, and actual ridge witnesses',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    if args.verify:verify(json.loads(OUTPUT.read_text()))
    else:build()
