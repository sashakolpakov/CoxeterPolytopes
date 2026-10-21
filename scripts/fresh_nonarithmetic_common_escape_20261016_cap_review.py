#!/usr/bin/env python3
"""Independent exact linear reconstruction of the49th nonpolar cap."""
import argparse
from collections import Counter
import json
import fresh_nonarithmetic_joint_regions_20261016_batch_review as B

ROOT,T,K,Z,O=B.ROOT,B.T,B.K,B.Z,B.O
TAG='fresh_nonarithmetic_common_escape_20261016_cap_review'
PRIMARY='results/fresh_nonarithmetic_common_escape_20261016_cap.json'
OLD_REVIEW='results/fresh_nonarithmetic_joint_regions_20261016_batch_review.json'


def geometry():
    d=B.read(PRIMARY);review=B.read(OLD_REVIEW)
    assert review['status']=='PASS_INDEPENDENT_SIXTEEN_POLARS_AND_SIMULTANEOUS_48_FACETS'
    assert review['sources'][B.BATCH]==B.digest(B.BATCH)
    G=B.geometry();N,cov,Q,H=[G[k] for k in ('N','cov','Q','H')]
    ids=[27,17,10,8,6,5,4,2,1,0]
    phi=(O+B.A.R.BASIS[1])/2
    targets=[-O if i in (8,10) else -phi if i==17 else Z for i in ids]
    rows=[[B.A.base(x) for x in cov[i]]+[v] for i,v in zip(ids,targets)]
    assert all(B.pair(N[i],N[i],Q)==T(2) for i in ids)
    pivots=[]
    for col in reversed(range(10)):
        k=len(pivots);j=next(j for j in range(k,10) if rows[j][col])
        rows[k],rows[j]=rows[j],rows[k]
        den=rows[k][col];rows[k]=[x/den for x in rows[k]]
        for j in range(10):
            if j!=k and rows[j][col]:
                factor=rows[j][col];rows[j]=[x-factor*y for x,y in zip(rows[j],rows[k])]
        pivots.append(col)
    v=[Z]*10
    for row,col in zip(rows,pivots):v[col]=row[-1]
    r=list(map(T,v));rcov=B.mv(Q,r)
    assert B.dot(r,rcov)==T(2)
    assert r==list(map(B.dec,d['cap']['normal']))
    assert N==[list(map(B.dec,row)) for row in d['normals'][:48]]
    assert all(B.dot(r,cov[i])==T(t) for i,t in zip(ids,targets))
    G['N']=N+[r];G['cov']=cov+[rcov]
    G['NI'] += [[B.iv(x) for x in r]]
    G['covI'] += [[B.iv(x) for x in rcov]]
    G['cap']=d
    return G


def run():
    G=geometry();d,Q,H,N,cov=[G[k] for k in ('cap','Q','H','N','cov')]
    norm=B.dot(N[-1],cov[-1]);assert norm==T(2)
    relations=[]
    for i in range(48):
        pairing=B.dot(N[-1],cov[i])
        relation=B.relation(pairing,norm*B.dot(N[i],cov[i]))
        assert relation['kind'] in ('right','finite','dotted')
        original=d['cap']['relations'][i]
        assert relation=={k:v for k,v in original.items() if k not in ('wall','pairing')}
        assert pairing==B.dec(original['pairing'])
        relations.append(dict(wall=i,**relation,pair_bounds=B.iv(pairing).saved()))
    assert Counter(r['kind'] for r in relations)==dict(right=8,finite=3,dotted=37)
    hh=B.pair(H,H,Q);hp=[B.dot(H,c) for c in cov]
    assert B.sg(hh)<0 and all(B.sg(x)<0 for x in hp)
    facets=[]
    for i in range(49):
        q=B.dot(N[i],cov[i]);assert B.sg(q)>0
        point=[q*h-hp[i]*n for h,n in zip(H,N[i])]
        square=q*q*hh-q*hp[i]*hp[i]
        assert B.sg(square)<0 and B.sg(q*hh)<0 and B.dot(point,cov[i])==0
        assert point==list(map(B.dec,d['facet_witnesses'][i]['point']))
        for j in range(49):
            if i!=j:assert B.sg(B.dot(N[i],cov[j]))<=0
        facets.append(dict(wall=i,norm_bounds=B.iv(square).saved()))
    escape=list(map(B.dec,G['batch']['common_strict_escape']['point']))
    assert escape==list(map(B.dec,d['cut_escape']['point']))
    cut=B.dot(escape,cov[-1]);assert B.sg(cut)>0
    assert cut==B.dec(d['cut_escape']['new_pairing'])
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',PRIMARY,OLD_REVIEW,B.BATCH,B.OLD_REVIEW,B.CONE_REVIEW]
    return dict(status='PASS_INDEPENDENT_49TH_NONPOLAR_CAP_AND_ALL_FACETS',
                sources={p:B.digest(p) for p in paths},recursive_integrity=B.R.integrity(paths),
                linear_rank=10,relations=relations,actual_facets=49,facets=facets,
                cut_pairing_bounds=B.iv(cut).saved(),
                proof='Reverse-column exact elimination over the independently instantiated octic field reconstructs the unique normal. All48 contacts and49 actual facet witnesses are checked with independent nested radical arithmetic. No new producer imports or numerical coordinates.',
                scope='Exact coherent49wall geometry and oldescape cut. Complete cone and finitevolume status are separate; no ambient incidence claim here.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==B.read(str(out.relative_to(ROOT)))
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
