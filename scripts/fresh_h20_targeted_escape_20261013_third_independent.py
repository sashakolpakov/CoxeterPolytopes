#!/usr/bin/env python3
"""Independent exact audit of the 18-short H20 extension and finite vertex."""
import argparse
from collections import Counter
from copy import deepcopy
from fractions import Fraction as Q
import json
from fresh_h20_targeted_escape_20261013_second_independent import (
    Field,ROOT,BASE,POLAR,CAPS,FIRST,read,sha,add,sub,enc)

SECOND='data/enumerations/fresh_h20_targeted_escape_20261013_second.json'
SOURCE='data/enumerations/fresh_h20_targeted_escape_20261013_third.json'
VOLUME='data/enumerations/fresh_h20_targeted_escape_20261013_third_volume.json'
OUT='results/fresh_h20_targeted_escape_20261013_third_independent.json'


def check(data):
    assert all(sha(p)==s for p,s in data['source_sha256'].items())
    assert all(p in data['source_sha256']for p in(BASE,POLAR,CAPS,FIRST,SECOND))
    b,p,c,first,second=map(read,(BASE,POLAR,CAPS,FIRST,SECOND))
    f=Field(data['field']['radicand']);zero=f.rational(0)
    M=tuple(tuple(f.base(x)for x in row)for row in b['form'])
    roots=[tuple(map(f.base,r))for r in b['normals']]+[tuple(map(f.base,r['vector']))for r in p['cuts']]
    cov=[f.mv(M,r)for r in roots];norms=[f.dot(r,v)for r,v in zip(roots,cov,strict=True)]
    h=tuple(map(f.base,b['strict_center']));hc=f.mv(M,h);assert f.sign(f.dot(h,hc))<0
    r=tuple(map(f.parse,data['normal']));rc=f.mv(M,r);q=f.dot(r,rc)
    assert len(r)==21 and q==f.rational(2)==f.parse(data['norm'])
    center=f.dot(r,hc);assert f.sign(center)<0 and center==f.parse(data['strict_center_pairing'])
    records=[]
    for i,(v,n)in enumerate(zip(cov,norms,strict=True)):
        pp=f.dot(r,v);assert f.sign(pp)<=0
        square=f.mul(pp,pp);qn=f.mul(q,n)
        if pp==zero:kind='2'
        elif f.sign(sub(square,qn))>0:kind='dotted'
        elif square==qn:kind='parallel'
        else:
            kinds=[str(m)for m,scale in[(3,Q(1,4)),(4,Q(1,2)),(6,Q(3,4))]if square==f.mul(f.rational(scale),qn)]
            assert len(kinds)==1;kind=kinds[0]
        saved=data['old_pair_checks'][i]
        assert saved['old']==i and saved['kind']==kind and pp==f.parse(saved['pairing'])
        records.append({'old':i,'kind':kind,'pairing':enc(pp)})
    raw=[{'field':x['field'],'normal':x['vector'],'norm':[x['norm'],['0','0']]}for x in c['caps']]+[first,second]
    extras=[]
    for i,cap in enumerate(raw,82):
        d=f.base(cap['field']['radicand']);other=tuple((f.base(x[0]),f.base(x[1]))for x in cap['normal'])
        pp=(f.dot(rc,tuple(x[0]for x in other)),f.dot(rc,tuple(x[1]for x in other)))
        oldq=tuple(f.base(x)for x in cap['norm']);square=f.extmul(pp,pp,d)
        delta=tuple(sub(x,f.mul(q,y))for x,y in zip(square,oldq,strict=True))
        assert f.extsign(pp,d)<0 and f.extsign(delta,d)>0
        saved=data['old_pair_checks'][i]
        assert saved['old']==i and saved['kind']=='dotted'and pp==tuple(map(f.parse,saved['pairing']))
        records.append({'old':i,'kind':'dotted','pairing':[enc(x)for x in pp]});extras.append((d,other))
    assert len(records)==len(data['old_pair_checks'])==87
    counts=dict(sorted(Counter(x['kind']for x in records).items()));assert counts==data['pair_counts']
    S=data['short_orthogonal_support'];assert len(S)==len(set(S))==18 and all(0<=i<24 for i in S)
    assert all(records[i]['kind']=='2'for i in S)
    assert all(records[i]['kind']==str(m)for i,m in zip(data['anchor_indices'],data['anchor_labels'],strict=True))
    cusp_pairs=[f.dot(r,f.mv(M,tuple(map(f.base,x['ray']))))for x in b['cusp_rows']]
    assert len(cusp_pairs)==107 and all(f.sign(x)<0 for x in cusp_pairs)

    def old_signs(x,allow_zero=()):
        xc=f.mv(M,x);pairs=[]
        for i,v in enumerate(cov):
            pp=f.dot(x,v);assert f.sign(pp)==(0 if i in allow_zero else -1);pairs.append(enc(pp))
        for d,other in extras:
            pp=(f.dot(xc,tuple(y[0]for y in other)),f.dot(xc,tuple(y[1]for y in other)))
            assert f.extsign(pp,d)<0;pairs.append([enc(y)for y in pp])
        return xc,pairs

    target=data['strict_target'];t=tuple(map(f.base,target['vector']));tc,_=old_signs(t)
    tq=f.dot(t,tc);th=f.dot(t,hc);cut=f.dot(r,tc)
    assert tq==f.base(target['norm'])and f.sign(tq)>0
    assert th==f.base(target['center_pairing'])and f.sign(th)<0
    assert cut==f.parse(target['new_cap_pairing'])and f.sign(cut)>0

    vertex=data['actual_finite_vertex'];J=vertex['incident_walls'];assert len(J)==len(set(J))==20 and J[-1]==87
    assert set(S).issubset(J)and all(0<=i<82 for i in J[:-1])
    incident=[roots[i]for i in J[:-1]]+[r];gr=[[f.dot(a,f.mv(M,b))for b in incident]for a in incident]
    expected=[[f.rational(Q(x))for x in row]for row in vertex['gram']];assert gr==expected
    # Exact positive definite Gaussian congruence, over the rational Gram.
    reduced=[[Q(x)for x in row]for row in vertex['gram']];pivots=[]
    for j in range(20):
        pivot=reduced[j][j];assert pivot>0;pivots.append(str(pivot))
        for k in range(j+1,20):
            for l in range(j+1,20):reduced[k][l]-=reduced[k][j]*reduced[j][l]/pivot
    x=tuple(map(f.parse,vertex['vector']));xc,pairs=old_signs(x,J[:-1]);xq=f.dot(x,xc)
    assert xq==f.parse(vertex['norm'])and f.sign(xq)<0 and f.sign(f.dot(x,hc))<0 and f.dot(x,rc)==zero
    saved_pairs=vertex['old_pairings'];assert len(saved_pairs)==87
    for i,saved in enumerate(saved_pairs):
        expected=enc(f.parse(saved))if i<82 else[enc(f.parse(y))for y in saved]
        assert pairs[i]==expected

    volume=read(VOLUME);assert all(sha(path)==digest for path,digest in volume['source_hashes'].items())
    assert volume['source_hashes'][SOURCE]==sha(SOURCE)
    escape=tuple(map(f.base,volume['certificate']['vector']));ec,_=old_signs(escape)
    eq=f.dot(escape,ec);eh=f.dot(escape,hc)
    assert eq==f.base(volume['certificate']['norm'])and f.sign(eq)>0
    assert eh==f.base(volume['certificate']['center_pairing'])and f.sign(eh)<0
    assert f.sign(f.dot(escape,rc))<0
    return {'status':'PASS_INDEPENDENT_88_WALL_H20_EXTENSION_VERTEX_AND_ESCAPE',
            'source':SOURCE,'source_sha256':sha(SOURCE),'dependencies':data['source_sha256'],
            'new_norm':enc(q),'new_pair_counts':counts,'new_pairs':records,
            'strict_center_pairing':enc(center),'inherited_cusps_strictly_preserved':107,
            'strict_target_norm':enc(tq),'target_cut_pairing':enc(cut),
            'actual_finite_vertex':{'incident_walls':J,'norm':enc(xq),'positive_gram_pivots':pivots},
            'remaining_escape':{'source':VOLUME,'source_sha256':sha(VOLUME),'norm':enc(eq),'center_pairing':enc(eh),'strict_inequalities':88},
            'scope':'Independent Fraction arithmetic checks87 new pairs,norm2,strict center,107 inherited cusp signs,an exact strict87-wall target cut,and an actual finite vertex with PD20x20 Gram. A new strict future spacelike escape proves infinite volume. The numerical1593-bank is not certified by this replay.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();data=read(SOURCE);out=check(data)
    for field in('normal','norm','strict_target'):
        bad=deepcopy(data)
        if field=='normal':bad[field][0][0][0]=str(Q(bad[field][0][0][0])+1)
        elif field=='norm':bad[field][0][0]='0'
        else:bad[field]['new_cap_pairing'][0][0]=str(Q(bad[field]['new_cap_pairing'][0][0])+1)
        try:check(bad)
        except AssertionError:pass
        else:raise AssertionError('forged '+field+' accepted')
    out['forged_normal_norm_cut_rejected']=True
    if a.verify:assert read(OUT)==out
    else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_pair_counts'])
