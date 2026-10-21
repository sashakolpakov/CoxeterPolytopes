#!/usr/bin/env python3
"""Exact Q(sqrt2) mixed short/long gluing around a B3 face.

All stored vectors are pairs of integer coefficient vectors (a,b), meaning
a+sqrt(2)b. Positive projective rescaling is permitted for wall normals.
The eight-sector angle cycle is gamma,beta,pi/4,pi/4,beta,gamma,pi/4,pi/4.
Global convexity is a separate gate, never inferred from local closure.
"""
from __future__ import annotations
from collections import Counter
from itertools import combinations
import json
from math import gcd
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
OUTPUT=ROOT/'data/enumerations/fresh_h20_mixed_b3_20261013.json'

def sign(a,b=0):
    a,b=int(a),int(b)
    if a==b==0:return 0
    if a>=0 and b>=0:return 1
    if a<=0 and b<=0:return -1
    return (1 if a>0 else -1)*(1 if a*a>2*b*b else -1)

def canonical(a,b):
    a,b=tuple(map(int,a)),tuple(map(int,b));p=next(i for i in range(len(a)) if a[i] or b[i])
    d=a[p]**2-2*b[p]**2;s=sign(a[p],b[p])*(1 if d>0 else -1)
    A=[s*(x*a[p]-2*y*b[p]) for x,y in zip(a,b)]
    B=[s*(y*a[p]-x*b[p]) for x,y in zip(a,b)]
    g=gcd(*(A+B));return tuple(x//g for x in A),tuple(x//g for x in B)

def pairing(x,y,Q):
    a,b=np.array(x,dtype=object);c,d=np.array(y,dtype=object)
    return int(a@Q@c+2*b@Q@d),int(a@Q@d+b@Q@c)

def multiply(x,y):
    a,b=x;c,d=y;return a*c+2*b*d,a*d+b*c

def scaled_sum(s,x,t,y):
    a,b=np.array(x,dtype=object);c,d=np.array(y,dtype=object)
    A,B=multiply(s,(a,b));C,D=multiply(t,(c,d));return A+C,B+D

def reflection(v,Q):
    q=int(v@Q@v);A=q*np.eye(22,dtype=object)-2*np.outer(v,v@Q)
    return reduce_matrix(A,np.zeros((22,22),dtype=object),q)

def reduce_matrix(A,B,den):
    g=gcd(int(den),*(map(int,A.flat)),*(map(int,B.flat)))
    return A//g,B//g,int(den)//g

def compose(T,U):
    A,B,d=T;C,D,e=U
    return reduce_matrix(A@C+2*B@D,A@D+B@C,d*e)

def transport(T,x):
    A,B,den=T;a,b=np.array(x,dtype=object)
    return A@a+2*B@b,A@b+B@a

def mixed(short,long,Q):
    s=int(short@Q@long);I=np.eye(22,dtype=object)
    gg=np.outer(long,long@Q);ff=np.outer(short,short@Q)
    cross=np.outer(short,long@Q)+np.outer(long,short@Q)
    assert int(short@Q@short)==2 and int(long@Q@long)==4
    if s==0:return reduce_matrix(4*I-gg-2*ff,-cross,4)
    assert s==-2
    return reduce_matrix(4*I-2*gg-4*ff-2*cross,-gg-2*ff-2*cross,4)

def data():
    raw=json.loads(INPUT.read_text())['coxeter_input']
    Q=-np.array(raw['basis_gram'],dtype=int).astype(object)
    W=np.array(raw['wall_coordinates_as_columns_transposed'],dtype=int).astype(object)
    G=W@Q@W.T;return Q,W,G

def section(Q,W,G,parent):
    f=W[parent];q=int(G[parent,parent]);V=q*W-np.outer(G[:,parent],f)
    ids=[i for i in range(210) if V[i]@Q@V[i]>0]
    h=W.sum(axis=0);center=q*h-(h@Q@f)*f
    assert center@Q@center<0 and max(V[ids]@Q@center)<0
    return ids,V[ids],center

def build_tiles():
    Q,W,G=data();a,b=0,21
    c=next(i for i in range(210) if G[i,i]==4 and G[a,i]==0 and G[b,i]==-2)
    I=np.eye(22,dtype=object);Z=np.zeros((22,22),dtype=object)
    Mab=reflection(W[a]+W[b],Q);Mbc=mixed(W[b],W[c],Q);Mac=mixed(W[a],W[c],Q)
    Rac=reflection(W[a],Q);Rbc=reflection(2*W[b]+W[c],Q)
    parents=[a,b,c,c,b,a,c,c]
    steps=[Mab,Mbc,Rac,Mbc,Mab,Mac,Rbc,Mac]
    T=(I,Z,1);tiles=[]
    for index,(parent,U) in enumerate(zip(parents,steps)):
        A,B,den=T
        assert np.array_equal(A.T@Q@A+2*B.T@Q@B,den*den*Q)
        assert not np.any(A.T@Q@B+B.T@Q@A)
        ids,V,center=section(Q,W,G,parent)
        roots=[canonical(*transport(T,(v,np.zeros(22,dtype=object)))) for v in V]
        image_center=transport(T,(center,np.zeros(22,dtype=object)))
        assert sign(*pairing(image_center,image_center,Q))<0
        assert all(sign(*pairing(r,image_center,Q))<0 for r in roots)
        assert all(pairing(r,(W[a],np.zeros(22,dtype=object)),Q)==(0,0) for r in roots)
        tiles.append({'name':chr(65+index),'parent':parent,'ids':ids,'roots':roots,'center':image_center,'T':T})
        T=compose(T,U)
    # Final attachment returns the first complete tile and the seam exactly.
    ids,V,center=section(Q,W,G,a)
    closure={canonical(*transport(T,(v,np.zeros(22,dtype=object)))) for v in V}==set(tiles[0]['roots'])
    return Q,W,G,tiles,{'ambient_B3_labels_one_based':[a+1,b+1,c+1],
                         'full_tile_cycle_closure':closure,
                         'transport_return_numerator_a':T[0].tolist(),
                         'transport_return_numerator_b':T[1].tolist(),'transport_return_denominator':T[2]}

def face_roots(roots,s,Q):
    q=pairing(s,s,Q);result=set()
    for r in roots:
        p=pairing(r,s,Q);z=scaled_sum(q,r,(-p[0],-p[1]),s)
        if sign(*pairing(z,z,Q))>0:result.add(canonical(*z))
    return result

def analyze(tiles,Q):
    internal=set();matches=[];coincidences=[]
    for i,j in combinations(range(len(tiles)),2):
        A,B=tiles[i],tiles[j];lookup={r:k for k,r in enumerate(B['roots'])}
        for k,r in enumerate(A['roots']):
            neg=tuple(tuple(-v for v in row) for row in r)
            l=lookup.get(neg)
            if l is None:continue
            same=face_roots(A['roots'],r,Q)==face_roots(B['roots'],r,Q)
            row={'tiles':[i,j],'local_indices':[k,l],'ambient_labels':[A['ids'][k]+1,B['ids'][l]+1],'full_facet_match':same}
            coincidences.append(row)
            if same:matches.append(row);internal.update(((i,k),(j,l)))
    boundary={}
    for i,tile in enumerate(tiles):
        for j,r in enumerate(tile['roots']):
            if (i,j) not in internal:boundary.setdefault(r,[]).append([i,j,tile['ids'][j]+1])
    violations=[]
    for r,owners in boundary.items():
        for i,tile in enumerate(tiles):
            value=pairing(r,tile['center'],Q)
            if sign(*value)>0:
                violations.append({'wall_owners':owners,'violated_by_tile_center':i,'pairing':list(value)});break
        if len(violations)>=8:break
    return {'tiles':len(tiles),'full_facet_matches':matches,'opposite_normal_coincidences':coincidences,
            'distinct_uncancelled_boundary_planes':len(boundary),
            'boundary_normal_pairs':[[list(x) for x in r] for r in sorted(boundary)],
            'boundary_sources':[boundary[r] for r in sorted(boundary)],
            'boundary_planes_with_nonrational_normals':sum(any(r[1]) for r in boundary),
            'center_support_violations':violations,
            'global_convexity_status':'DISPROVED_BY_EXACT_CENTER' if violations else 'PENDING_EXACT_FULL_SUPPORT'},boundary

def main():
    Q,W,G,tiles,cycle=build_tiles()
    result={'source':str(INPUT.relative_to(ROOT)),'number_field':'Q(sqrt2)',
            'scope':'Mixed short/long B3 cycle and its two-tile seams; no Coxeter or global convexity assertion without separate certificates.',
            'cycle':cycle,'tile_records':[]}
    for t in tiles:
        A,B,d=t['T'];result['tile_records'].append({'name':t['name'],'parent_one_based':t['parent']+1,
          'wall_ambient_labels_one_based':[i+1 for i in t['ids']],
          'transport_a':A.tolist(),'transport_b':B.tolist(),'transport_denominator':d,
          'wall_normal_pairs':[[list(x) for x in r] for r in t['roots']],
          'center_pair':[list(map(int,x)) for x in t['center']]})
    result['eight_tile'],unused=analyze(tiles,Q)
    result['two_tile_B2'],unused=analyze(tiles[1:3],Q)
    result['two_tile_orthogonal'],unused=analyze(tiles[5:7],Q)
    OUTPUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    for label in ['cycle','eight_tile','two_tile_B2','two_tile_orthogonal']:
        print(label,{k:v for k,v in result[label].items() if k not in ('boundary_normal_pairs','boundary_sources','transport_return_numerator_a','transport_return_numerator_b')},flush=True)

if __name__=='__main__':main()
