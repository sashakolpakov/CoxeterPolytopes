#!/usr/bin/env python3
"""Independent construction of the three filtered caps from reviewed49 geometry.

Old-field projection identities are exact. New square roots are defined
algebraically and only strict signs use independent rational intervals.
"""
import argparse
import json
from fractions import Fraction as F
from itertools import combinations
import fresh_nonarithmetic_common_escape_20261016_cap_review as C

B=C.B
ROOT,T,I,K,Z,O=B.ROOT,B.T,B.I,B.K,B.Z,B.O
TAG='fresh_nonarithmetic_affine_filtered_20261016_cluster_review'
PRIMARY='results/fresh_nonarithmetic_affine_filtered_20261016_cluster.json'
OLD_CONE='results/fresh_nonarithmetic_common_escape_20261016_cone.json'
OLD_CONE_REVIEW='results/fresh_nonarithmetic_common_escape_20261016_cone_review.json'
SUPPORT=[0,1,2,3,5,6,7,30]
PROFILES=[{30:5,7:3,33:10},{3:3,7:5,33:10},{5:6,6:4,33:10}]


def invert_nested(x):
    if isinstance(x,T):
        d=B.A.INNER(B.A.D)
    elif isinstance(x,B.A.P.Y):
        d=B.A.P.X(B.A.P.D)
    elif isinstance(x,B.A.P.P.Radical):
        d=B.A.P.P.D
    else:
        assert x
        return O/x
    q=invert_nested(x.a*x.a-d*x.b*x.b)
    result=type(x)(x.a*q,-x.b*q)
    assert x*result==1
    return result


def inverse_base(matrix):
    n=len(matrix)
    a=[row[:]+[O if i==j else Z for j in range(n)] for i,row in enumerate(matrix)]
    for col in reversed(range(n)):
        p=next(i for i in range(col+1) if a[i][col])
        a[col],a[p]=a[p],a[col]
        q=a[col][col];a[col]=[x/q for x in a[col]]
        for i in range(n):
            if i!=col:
                q=a[i][col];a[i]=[x-q*y for x,y in zip(a[i],a[col])]
    assert [row[:n] for row in a]==[[O if i==j else Z for j in range(n)] for i in range(n)]
    return [row[n:] for row in a]


def geometry():
    g=C.geometry();target=B.read(PRIMARY);oldcone=B.read(OLD_CONE)
    review=B.read(OLD_CONE_REVIEW)
    assert review['sources'][OLD_CONE]==B.digest(OLD_CONE)
    assert target['spherical_support']==SUPPORT and target['anchor']==33
    N,cov,Q,H=[g[k] for k in ('N','cov','Q','H')]
    gram=[[B.A.base(B.dot(N[i],cov[j])) for j in SUPPORT] for i in SUPPORT]
    positive=[row[:] for row in gram]
    for i in range(8):
        assert B.sg(T(positive[i][i]))>0
        for j in range(i+1,8):
            for k in range(j,8):
                positive[j][k]=positive[k][j]=positive[j][k]-positive[j][i]*positive[i][k]/positive[i][i]
    J=inverse_base(gram)
    def project(v):
        c=[B.dot(v,cov[j]) for j in SUPPORT]
        weights=[B.dot(row,c) for row in J]
        return [B.dot(weights,[N[j][k] for j in SUPPORT]) for k in range(10)]
    a=[x-y for x,y in zip(N[33],project(N[33]))]
    qa=B.pair(a,a,Q)
    assert B.sg(qa)<0 and all(B.dot(a,cov[i])==0 for i in SUPPORT)
    t=[x-y for x,y in zip(H,project(H))]
    alpha=B.pair(t,a,Q)*invert_nested(qa)
    kernel=[x-alpha*y for x,y in zip(t,a)]
    assert all(B.dot(kernel,cov[i])==0 for i in SUPPORT+[33])
    qk=B.pair(kernel,kernel,Q);assert B.sg(qk)>0
    E=list(map(B.dec,oldcone['strict_escape']['point']))
    if B.sg(B.pair(kernel,E,Q))<0:kernel=[-x for x in kernel]
    assert B.sg(B.pair(kernel,E,Q))>0
    assert a==list(map(B.dec,target['projected_anchor']))
    assert kernel==list(map(B.dec,target['positive_kernel']))
    assert qa==B.dec(target['projected_anchor_norm']) and qk==B.dec(target['kernel_norm'])
    assert all(B.dot(N[i],cov[i])==2 for i in SUPPORT+[33])
    JI=[[B.ib(x) for x in row] for row in J]
    ai,ki=[B.iv(x) for x in a],[B.iv(x) for x in kernel]
    qai,qki=B.iv(qa),B.iv(qk)
    theta={2:I(),3:I(1),4:I(2).sqrt(),5:(I(1)+I(5).sqrt())/2,
           6:I(3).sqrt(),10:((I(5)+I(5).sqrt())/2).sqrt()}
    caps=[];records=[]
    for offset,profile in enumerate(PROFILES):
        old=target['caps'][offset]
        assert dict(old['profile'])==profile and old['norm_exact']==2
        prescribed=[-theta[profile.get(i,2)] for i in SUPPORT]
        weights=B.imv(JI,prescribed)
        pS=[B.idot(weights,[g['NI'][j][k] for j in SUPPORT]) for k in range(10)]
        factor=(-theta[profile[33]]-B.idot(pS,g['covI'][33]))/qai
        p=[x+factor*y for x,y in zip(pS,ai)]
        pp=B.idot(prescribed,weights)+factor*factor*qai
        delta=(2-pp)/qki
        assert delta.lo>0
        cap=[x+delta.sqrt()*y for x,y in zip(p,ki)]
        for x,raw in zip(cap,old['normal_interval']):
            lo,hi=map(F,raw);assert x.lo<=hi and lo<=x.hi
        # Exact norm and contacts follow from the projection decomposition:
        # p.k=0, k.S=k.n33=0, r=p+sqrt((2-p^2)/k^2)*k.
        contacts=[]
        for i,c in enumerate(g['covI']):
            value=B.idot(cap,c)
            if i in SUPPORT+[33]:
                label=profile.get(i,2);kind='right' if label==2 else 'finite'
                expected=-theta[label]
                assert value.lo<=expected.hi and expected.lo<=value.hi
            else:
                label=None;kind='dotted';norm=B.iv(B.dot(N[i],cov[i]))
                assert value.hi<0 and value.hi*value.hi>2*norm.hi,(offset,i)
            assert old['relations'][i]['kind']==kind
            if label is not None:assert old['relations'][i]['label']==label
            contacts.append(dict(wall=i,kind=kind,label=label,bounds=value.saved()))
        hp=B.idot(cap,g['hI']);cut=B.idot(cap,B.imv(g['QI'],[B.iv(x) for x in E]))
        assert hp.hi<0<cut.lo
        records.append(dict(wall=49+offset,profile=sorted(profile.items()),norm_exact=2,
                            radical_square=delta.saved(),contacts=contacts,
                            center_pair=hp.saved(),previous_escape_cut=cut.saved()))
        caps.append(cap)
    mutual=[]
    for i,j in combinations(range(3),2):
        value=B.idot(caps[i],B.imv(g['QI'],caps[j]));assert value.hi<-2
        mutual.append(dict(walls=[49+i,49+j],kind='dotted',bounds=value.saved()))
    full=g['NI']+caps
    norms=[B.iv(B.dot(n,c)) for n,c in zip(N,cov)]+[I(2)]*3
    hh=B.idot(g['HI'],g['hI']);assert hh.hi<0
    facets=[]
    for i,(normal,norm) in enumerate(zip(full,norms)):
        hp=B.idot(normal,g['hI']);assert hp.hi<0
        value=hh-hp*hp/norm;assert value.hi<0
        facets.append(dict(wall=i,projected_center_norm=value.saved()))
    g['NI']=full;g['covI']+=[B.imv(g['QI'],cap) for cap in caps]
    g.update(new_caps=caps,new_records=records,new_mutual=mutual,facets=facets)
    return g


def run():
    g=geometry()
    paths=[f'scripts/{TAG}.py',PRIMARY,f'scripts/{C.TAG}.py',f'scripts/{B.TAG}.py',
           C.PRIMARY,f'results/{C.TAG}.json',OLD_CONE,OLD_CONE_REVIEW]
    return dict(status='PASS_INDEPENDENT_FILTERED_THREE_CAPS_AND_52_FACETS',
                sources={p:B.digest(p) for p in paths},caps=g['new_records'],mutual=g['new_mutual'],
                actual_facets=52,facet_witnesses=g['facets'],
                normals=[[x.saved() for x in v] for v in g['NI']],
                proof='Independent49geometry, reverse base-field inversion and exact nested-field anchor/kernel projection. New positive-root formulas prove all prescribed norm/contact identities. Independent rational bounds prove120 unspecified oldcontacts and3mutualpairs dotted, all3oldescape cuts and52actualfacets.',
                scope='Exact local52wall chamber and previousescape cuts. Completecone and right33gate are separate; no finitevolume conclusion.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)
