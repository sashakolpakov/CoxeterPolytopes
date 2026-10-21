#!/usr/bin/env python3
"""Exact algebraic cone certificates for the mixed B3 eight-tile patch.

LP selects candidate cone generators only. Small rational linear systems
recover both coefficients of Q(sqrt2), and every identity and coefficient
sign is checked exactly. No approximate feasibility proves a result.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import json
from math import gcd,lcm
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
import sympy as sp
import fresh_h20_mixed_b3_20261013 as mixed

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/enumerations/fresh_h20_mixed_b3_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h20_mixed_support_20261013.json'

def sign(a,b):
    d=lcm(a.denominator,b.denominator);return mixed.sign(int(a*d),int(b*d))

def setup():
    source=json.loads(SOURCE.read_text());Q,W,G=mixed.data()
    Qi=sp.Matrix(Q).inv();d=lcm(*(int(x.q) for x in Qi));Qi=np.array(d*Qi,dtype=object)
    tiles=[]
    for row in source['tile_records']:
        A=np.array(row['transport_a'],dtype=object);B=np.array(row['transport_b'],dtype=object)
        inverse=(Qi@A.T@Q,Qi@B.T@Q,d*row['transport_denominator'])
        ids,V,center=mixed.section(Q,W,G,row['parent_one_based']-1)
        V=np.array([tuple(int(x)//gcd(*map(int,v)) for x in v) for v in V],dtype=object)
        roots=[tuple(tuple(v) for v in r) for r in row['wall_normal_pairs']]
        tiles.append({'name':row['name'],'V':V,'inverse':inverse,'roots':roots,
                      'floatV':np.array(V.T,dtype=float),'lookup':{tuple(v):i for i,v in enumerate(V)},'cache':{}})
    return source,Q,tiles

def pullback(normal,tile):
    return mixed.canonical(*mixed.transport(tile['inverse'],normal))

def check_terms(normal,tile,terms):
    a,b=pullback(normal,tile);V=tile['V']
    rows=[(int(i),F(x),F(y)) for i,x,y in terms]
    assert rows and len({r[0] for r in rows})==len(rows)
    assert all(0<=i<len(V) and sign(x,y)>0 for i,x,y in rows)
    for j in range(22):
        assert sum(x*int(V[i,j]) for i,x,y in rows)==a[j]
        assert sum(y*int(V[i,j]) for i,x,y in rows)==b[j]

def find_terms(normal,tile):
    a,b=pullback(normal,tile);key=(a,b)
    if key in tile['cache']:return tile['cache'][key]
    if not any(b) and a in tile['lookup']:
        result=[[tile['lookup'][a],'1','0']];tile['cache'][key]=result;return result
    y=np.array(a,dtype=float)+np.sqrt(2)*np.array(b,dtype=float);y/=max(abs(y))
    proposal=linprog(np.ones(len(tile['V'])),A_eq=tile['floatV'],b_eq=y,bounds=(0,None),method='highs')
    if not proposal.success:return None
    indices=[int(i) for i in np.where(proposal.x>1e-8)[0]]
    M=sp.Matrix(tile['V'][indices].T);target=sp.Matrix.hstack(sp.Matrix(a),sp.Matrix(b))
    try:coeff,parameters=M.gauss_jordan_solve(target)
    except ValueError:return None
    assert not parameters.rows
    terms=[]
    for i,x,y in zip(indices,coeff[:,0],coeff[:,1]):
        if x or y:terms.append([i,str(x),str(y)])
    check_terms(normal,tile,terms);tile['cache'][key]=terms;return terms

def noncoxeter_ridge(source,Q):
    tile=source['tile_records'][0];roots=np.array(tile['wall_normal_pairs'],dtype=object)
    ids=tile['wall_ambient_labels_one_based']
    r,s=roots[ids.index(24),0],roots[ids.index(49),0]
    h=np.array(tile['center_pair'][0],dtype=object)
    assert r@Q@r==s@Q@s==6 and r@Q@s==-2
    a,b=h@Q@r,h@Q@s;p=32*h-(6*a+2*b)*r-(2*a+6*b)*s;p=p//gcd(*map(int,p))
    boundary=np.array(source['eight_tile']['boundary_normal_pairs'],dtype=object)
    pa,pb=boundary[:,0]@Q@p,boundary[:,1]@Q@p
    assert p@Q@p<0 and all(mixed.sign(x,y)<=0 for x,y in zip(pa,pb))
    zero=[i for i,(a,b) in enumerate(zip(pa,pb)) if a==b==0];assert len(zero)==2
    return {'ambient_tile_A_labels':[24,49],'boundary_indices_zero_based':zero,
            'normal_gram':[[6,-2],[-2,6]],'cosine_squared':'1/9',
            'relative_interior_timelike_point':list(map(int,p)),'square':int(p@Q@p),
            'all_other_boundary_pairings_strictly_negative':True}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    source,Q,tiles=setup();boundary=[tuple(tuple(x) for x in r) for r in source['eight_tile']['boundary_normal_pairs']]
    internal={}
    for row in source['eight_tile']['full_facet_matches']:
        for tile,wall in zip(row['tiles'],row['local_indices']):internal.setdefault(tile,set()).add(wall)
    if args.verify:
        result=json.loads(OUTPUT.read_text())
    else:
        result={'source':str(SOURCE.relative_to(ROOT)),'support_proofs':[],'disjointness':[],
                'number_field':'Q(sqrt2)','scope':'Only the literal mixed eight-tile B3 patch.'}
        for index,r in enumerate(boundary):
            proofs=[find_terms(r,tile) for tile in tiles]
            if any(p is None for p in proofs):
                result['status']='FULL_SUPPORT_NOT_CERTIFIED';result['failed_boundary_index']=index
                result['failed_tiles']=[i for i,p in enumerate(proofs) if p is None]
                OUTPUT.write_text(json.dumps(result,indent=2)+'\n');print(result['status'],index,result['failed_tiles'],flush=True);return
            result['support_proofs'].append(proofs)
            if index%100==0:print('Certified boundary',index,'of',len(boundary),flush=True)
        for a,b in combinations(range(8),2):
            for i in sorted(internal[a]):
                r=tiles[a]['roots'][i];negative=tuple(tuple(-v for v in row) for row in r)
                terms=find_terms(negative,tiles[b])
                if terms is not None:
                    result['disjointness'].append({'tiles':[a,b],'wall_in_first':i,'terms_in_second':terms});break
            else:raise AssertionError(('No separator',a,b))
        result['status']='EXACT_CONVEX_FINITE_VOLUME_MIXED_PATCH'
        result['facets']=len(boundary);result['support_certificate_count']=len(boundary)*8
        result['non_coxeter_ridge']=noncoxeter_ridge(source,Q)
        OUTPUT.write_text(json.dumps(result,indent=2)+'\n')
    assert result['status']=='EXACT_CONVEX_FINITE_VOLUME_MIXED_PATCH'
    assert len(result['support_proofs'])==len(boundary)==767
    for r,proofs in zip(boundary,result['support_proofs']):
        assert len(proofs)==8
        for tile,terms in zip(tiles,proofs):check_terms(r,tile,terms)
    assert {tuple(r['tiles']) for r in result['disjointness']}==set(combinations(range(8),2))
    for row in result['disjointness']:
        a,b=row['tiles'];i=row['wall_in_first'];assert i in internal[a]
        r=tiles[a]['roots'][i];negative=tuple(tuple(-v for v in x) for x in r)
        check_terms(negative,tiles[b],row['terms_in_second'])
    assert result['non_coxeter_ridge']==noncoxeter_ridge(source,Q)
    print(result['status'],'facets',len(boundary),'support',len(boundary)*8,'disjoint pairs',len(result['disjointness']),flush=True)
    print('Cone-generator count histogram',dict(Counter(len(p) for row in result['support_proofs'] for p in row)),flush=True)

if __name__=='__main__':main()
