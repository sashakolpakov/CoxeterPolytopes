#!/usr/bin/env python3
"""Exact actual quadrilaterals for the fifteen D7/E8 weight-barrier witnesses.

Uses the separately reviewed Fraction Qsqrt5 arithmetic; no root producer
imports and no numerical point construction.
"""
import argparse,json
from hashlib import sha256
from itertools import combinations,permutations,product
from pathlib import Path
import fresh_sqrt5_quad_filters_20261015_review as A
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_sealed_weights_20261015_geometry'
Z=A.Z;ONE=A.ONE

def dot(a,b):
    s=Z
    for x,y in zip(a,b):s=A.add(s,A.mul(x,y))
    return s

def determinant(M):
    n=len(M);s=Z
    for p in permutations(range(n)):
        v=ONE
        for i,j in enumerate(p):v=A.mul(v,M[i][j])
        if sum(p[i]>p[j]for i,j in combinations(range(n),2))%2:v=A.neg(v)
        s=A.add(s,v)
    return s

def encode(M):return [[A.enc(x)for x in r]for r in M]

def geometry(row):
    M=A.matrix(row,A.coloring(row));S=A.residual(M,row['white'])
    assert len(S)==4 and A.inertia(S)==[2,1,1]
    B=next(list(ids)for ids in combinations(range(4),3)if determinant([[S[i][j]for j in ids]for i in ids])!=Z)
    G=[[S[i][j]for j in B]for i in B];corners=[];rejected=[]
    for a,b in combinations(range(4),2):
        v=[S[a][i]for i in B];w=[S[b][i]for i in B]
        x=[A.sub(A.mul(v[1],w[2]),A.mul(v[2],w[1])),
           A.sub(A.mul(v[2],w[0]),A.mul(v[0],w[2])),
           A.sub(A.mul(v[0],w[1]),A.mul(v[1],w[0]))]
        assert any(z!=Z for z in x)
        products=[dot([S[i][j]for j in B],x)for i in range(4)]
        sign=next(A.sign(z)for z in products if z!=Z)
        if sign>0:x=list(map(A.neg,x));products=list(map(A.neg,products))
        norm=dot(x,[products[i]for i in B])
        record=dict(pair=[a,b],coordinates=list(map(A.enc,x)),products=list(map(A.enc,products)),norm=A.enc(norm))
        if all(A.sign(z)<=0 for z in products):
            assert A.sign(norm)<0
            assert [i for i,z in enumerate(products)if z==Z]==[a,b]
            corners.append(record)
        else:rejected.append(record)
    assert {tuple(r['pair'])for r in corners}=={(0,1),(0,3),(1,2),(2,3)}
    assert len(rejected)==2
    vectors=[[tuple(map(A.F,z))for z in r['coordinates']]for r in corners]
    center=[Z]*3
    for v in vectors:center=[A.add(a,b)for a,b in zip(center,v)]
    products=[dot([S[i][j]for j in B],center)for i in range(4)]
    norm=dot(center,[products[i]for i in B]);assert A.sign(norm)<0 and all(A.sign(z)<0 for z in products)
    assert all(A.sign(dot(v,[products[i]for i in B]))<0 for v in vectors)
    return dict(white_rank=len(row['white']),full_physical_inertia=A.inertia(M),
                full_conjugate_inertia=A.inertia([[A.conj(z)for z in r]for r in M]),
                projected_gram=encode(S),basis=B,basis_gram=encode(G),corners=corners,
                rejected_opposites=rejected,strict_center=dict(coordinates=list(map(A.enc,center)),products=list(map(A.enc,products)),norm=A.enc(norm)))

def run():
    lp=json.loads((ROOT/'results/fresh_sqrt5_sealed_weights_20261015.json').read_text())
    union=json.loads((ROOT/'results/fresh_sqrt5_quad_exhaustion_20261015.json').read_text())
    ids=sorted({r['class_index']for c in lp['component_minima']for r in c['dual']})
    assert ids==[535,563,843,940,958,967,982,984,993,995,1032,1052,1107,1143,1145]
    records=[]
    for j in ids:
        ui=union['marked_isomorphism_classes'][j]['representative'];row=union['labelled_union_supports'][ui]
        rec=geometry(row);records.append(dict(class_index=j,union_index=ui,**rec))
    files=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_quad_filters_20261015_review.py','results/fresh_sqrt5_quad_filters_20261015_review.json',
           'scripts/fresh_sqrt5_sealed_weights_20261015.py','results/fresh_sqrt5_sealed_weights_20261015.json',
           'results/fresh_sqrt5_quad_exhaustion_20261015.json']
    return dict(status='PASS_FIFTEEN_EXACT_COMPACT_QUADRILATERAL_WEIGHT_WITNESSES',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in files},
                records=records,pair_intersections=6*len(records),finite_vertices=4*len(records),
                scope='Actual compact quadrilateral white faces of local arithmetic-compatible Coxeter chambers. No compactness/finite volume of the full local chamber or global extension is claimed.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['pair_intersections'],out['finite_vertices'])
