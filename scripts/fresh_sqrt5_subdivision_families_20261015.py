#!/usr/bin/env python3
"""Exact all-k principal-star subdivision families; no rank extrapolation.

Ten literal sealed starting graphs yield nine unbounded necessary families
and one family with the exact conjugate bound k<=10. Sealing is not asserted
for arbitrary k. The finite rank<=29 table is illustrative weighting input.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_quad_filters_20261015_review as R
import fresh_sqrt5_quad_contacts_20261015_review as Q
import fresh_sqrt5_singular_attachment_closure_20261015 as P

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_subdivision_families_20261015'
ARCHIVE='data/enumerations/fresh_sqrt5_nonsingular_attachment_search_20261015_batch610.json'
BASES=[333,343,344,353,367,368,377,395,396,430]
K=Q.K;U,V=7,8

def subdivide(row,k):
    assert k>=0
    n=row['nodes'];g=Q.edges(row);assert g[U,V]==3
    if k:
        g[U,V]=2
        route=[U]+list(range(n,n+k))+[V]
        for a,b in zip(route,route[1:]):g[Q.pair(a,b)]=3
    return dict(nodes=n+k,white=list(range(4,n+k)),labels=[g.get(p,2)for p in combinations(range(n+k),2)])

def qform(J,ids,d):
    z=[K(int(i==U)-int(i==V))for i in ids]
    return K(*d)/2*sum((x*K(*J[i][j])*y for i,x in enumerate(z)for j,y in enumerate(z)),K())

def run():
    D=json.loads((ROOT/ARCHIVE).read_text());assert D['state']['next_core']==610
    weight=[json.loads((ROOT/'results'/f'{p}.json').read_text())for p in('fresh_sqrt5_forest_sensitive_weights_20261015','fresh_sqrt5_sealed_weights_20261015','fresh_sqrt5_extended_forest_weights_20261015')]
    out=[];patterns=[];counts=Counter()
    for index in BASES:
        row=D['records'][index]['row'];n=row['nodes'];assert n==12
        assert D['records'][index]['result']['status']=='UNPAID_SEALED_FRONTIER'
        M=R.matrix(row,R.coloring(row));B=[[K(*x)for x in r]for r in M]
        assert R.inertia(M)==[10,1,1]and R.inertia([[R.conj(x)for x in r]for r in M])==[11,0,1]
        kernel=Q.kernel(B);assert len(kernel)==1;v=kernel[0];assert v[U]==v[V]
        ground=next(i for i,x in enumerate(v)if x);ids=[i for i in range(n)if i!=ground]
        H=[[M[i][j]for j in ids]for i in ids];J=R.inverse(H)
        assert R.inertia(H)==[10,1,0]and R.inertia([[R.conj(x)for x in r]for r in H])==[11,0,0]
        d=M[U][U];assert d==M[V][V]and M[U][V]==R.div(R.neg(d),R.TWO)
        q=qform(J,ids,d);corners=[]
        for a,b in product((0,2),(1,3)):
            ii=row['white']+[a,b];C=[[M[i][j]for j in ii]for i in ii]
            assert R.inertia(C)==[10,0,0]and R.inertia([[R.conj(x)for x in r]for r in C])==[10,0,0]
            c=qform(R.inverse(C),ii,d)
            assert (1-c).sign()>=0 and (1-c.conj()).sign()>=0
            corners.append(dict(marked=[a,b],cost=c.enc()))
        if index==333:
            assert q==K(F(37,38),F(-1,19))
            cutoff=K(1)/(q.conj()-1)
            assert cutoff==K(2,4)and(cutoff-10).sign()>0 and(11-cutoff).sign()>0
            permitted='0<=k<=10';maxk=10
        else:
            assert q in(K(1),K(F(67,76),F(-1,76)))
            assert(1-q).sign()>=0 and(1-q.conj()).sign()>=0
            permitted='all integers k>=0';maxk=19
        # Direct literal Schur complements verify the universal formula on
        # small cases; its all-k proof is the exact Cartan(A_k) inverse.
        checks=[]
        for k in(1,2,3):
            ext=subdivide(row,k);G=R.matrix(ext,R.coloring(ext));inside=list(range(n,n+k))
            I=R.inverse([[G[i][j]for j in inside]for i in inside])
            assert K(*I[0][0])==K(2*F(k,k+1))/K(*d)
            assert K(*I[0][-1])==K(2*F(1,k+1))/K(*d)
            for i,j in product(range(n),repeat=2):
                correction=sum((K(*G[i][u])*K(*I[x][y])*K(*G[w][j])for x,u in enumerate(inside)for y,w in enumerate(inside)),K())
                val=K(*G[i][j])-correction
                target=B[i][j]-K(*d)/2*F(k,k+1)*(int(i==U)-int(i==V))*(int(j==U)-int(j==V))
                assert val==target
            assert R.inertia(G)==[n+k-2,1,1]
            assert R.inertia([[R.conj(x)for x in r]for r in G])==[n+k-1,0,1]
            checks.append(k)
        out.append(dict(source_index=index,base_row=row,edge=[U,V],kernel=[x.enc()for x in v],ground=ground,norm=R.enc(d),core_cost=q.enc(),corner_costs=corners,permitted_subdivisions=permitted,direct_schur_checks=checks))
        for k in range(maxk+1):
            ext=subdivide(row,k);g=Q.edges(ext)
            assert Q.no_short(ext['nodes'],g)and not Q.paid(ext)
            pay=P.payment(ext,*weight);assert pay['sums_without_with_external_a1']==['0','0']
            # This exact finite table checks the diagrams used for target
            # dimensions<=29. The all-rank PD statement follows separately.
            assert all(Q.finite(ext['white']+[a,b],g)for a,b in product((0,2),(1,3)))
            patterns.append(dict(source_index=index,subdivisions=k,corner_rank=10+k,row=ext,payment=pay));counts[10+k]+=1
    assert len(out)==10 and len(patterns)==191
    paths=[f'scripts/{TAG}.py',ARCHIVE,'results/fresh_sqrt5_nonsingular_attachment_search_20261015_batch610.json',
      'scripts/fresh_sqrt5_quad_filters_20261015_review.py','scripts/fresh_sqrt5_quad_contacts_20261015_review.py',
      'scripts/fresh_sqrt5_singular_attachment_closure_20261015.py','results/fresh_sqrt5_extended_forest_weights_20261015.json']
    manifests=Q.manifest_check(paths)
    return dict(status='PASS_EXACT_SUBDIVISION_PRINCIPAL_STAR_FAMILIES',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},recursive_manifests=manifests,families=out,rank_at_most29_patterns=patterns,rank_counts={str(k):v for k,v in sorted(counts.items())},
      all_rank_identity='For k inserted norm-d roots on a pi3 edge U,V, eliminate (d/2)Cartan(A_k). The original core becomes H-[d*k/(2(k+1))](eU-eV)(eU-eV)^T. The old kernel remains because its U,V coordinates agree. The grounded determinant multiplier is1-[k/(k+1)]q, q=(d/2)(eU-eV)^T H_ground^-1(eU-eV). Same in each corner and conjugate embedding.',
      scope='Nine all-integer-k necessary principal-star families and one exact finite k<=10 family, with four PD corners, physical one-negative corank1 and conjugate PSD corank1. No rank cutoff or extrapolation supplies these identities. The191 explicit instances describe corner ranks10..29 for weight exploration. Arbitrary-k sealing, exhaustive affine absence, future quadrilateral realization and global compact realizability are not asserted. The active BFS is not modified or declared complete.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();out=run();path=ROOT/'results'/(TAG+'.json')
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],len(out['families']),len(out['rank_at_most29_patterns']))
