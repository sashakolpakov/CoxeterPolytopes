#!/usr/bin/env python3
"""Standalone Fraction replay of the 86th H20 wall and its finite vertex.

No project arithmetic or construction modules are imported. The base algebra
uses coefficients in 1,sqrt6,sqrt10,sqrt60, with exact sign comparisons.
"""
import argparse
from collections import Counter
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE='data/enumerations/fresh_h20_targeted_escape_20261013_exact.json'
OUT='results/fresh_h20_targeted_escape_20261013_independent.json'
BASE='data/enumerations/fresh_dimension_bridge_20261013.json'
POLAR='data/enumerations/fresh_dimension_bridge_polar_20261013.json'
CAPS='data/enumerations/fresh_dimension_bridge_two_angle_20261013_exact.json'
VOLUME='data/enumerations/fresh_h20_targeted_escape_20261013_volume.json'

def read(path):return json.loads((ROOT/path).read_text())
def sha(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
def rational(x):return (F(x),F(0),F(0),F(0))
ZERO=rational(0)
def base(x):return(F(x[0]),F(x[1]),F(0),F(0))
def radical(x):return(F(x[0][0]),F(x[0][1]),F(x[1][0]),F(x[1][1]))
def add(a,b):return tuple(x+y for x,y in zip(a,b,strict=True))
def neg(a):return tuple(-x for x in a)
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    out=[F(0)]*4
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            common=i&j
            out[i^j]+=x*y*(6 if common&1 else 1)*(10 if common&2 else 1)
    return tuple(out)
def signq(a):
    x,y=a
    sx=(x>0)-(x<0);sy=(y>0)-(y<0)
    if not sy:return sx
    if not sx or sx==sy:return sy
    delta=x*x-6*y*y
    return sx*((delta>0)-(delta<0))
def sign(a):
    sa,sb=signq(a[:2]),signq(a[2:])
    if not sb:return sa
    if not sa or sa==sb:return sb
    aa=(a[0],a[1],F(0),F(0));bb=(a[2],a[3],F(0),F(0))
    return sa*signq(sub(mul(aa,aa),mul(rational(10),mul(bb,bb)))[:2])
def extmul(x,y,d):return(add(mul(x[0],y[0]),mul(d,mul(x[1],y[1]))),add(mul(x[0],y[1]),mul(x[1],y[0])))
def extsign(x,d):
    sa,sb=sign(x[0]),sign(x[1])
    if not sb:return sa
    if not sa or sa==sb:return sb
    return sa*sign(sub(mul(x[0],x[0]),mul(d,mul(x[1],x[1]))))
def dot(a,b):
    out=ZERO
    for x,y in zip(a,b,strict=True):out=add(out,mul(x,y))
    return out
def mv(M,x):return tuple(dot(row,x)for row in M)
def encode(x):return list(map(str,x))
def rankq(M):
    a=[list(row)for row in M];rank=0
    for j in range(len(a[0])):
        k=next((i for i in range(rank,len(a))if a[i][j]),None)
        if k is None:continue
        a[rank],a[k]=a[k],a[rank];p=a[rank][j];a[rank]=[x/p for x in a[rank]]
        for i in range(rank+1,len(a)):
            p=a[i][j];a[i]=[x-p*y for x,y in zip(a[i],a[rank],strict=True)]
        rank+=1
        if rank==len(a):break
    return rank
def pair_label(p,q,n):
    assert sign(p)<=0
    if p==ZERO:return '2'
    square=mul(p,p);qn=mul(q,n);delta=sign(sub(square,qn))
    if delta>0:return 'dotted'
    if delta==0:return 'parallel'
    for m,c in [(3,F(1,4)),(4,F(1,2)),(6,F(3,4))]:
        if square==mul(rational(c),qn):return str(m)
    raise AssertionError('forbidden angle')

def check(data):
    b,p,c=read(BASE),read(POLAR),read(CAPS)
    assert data['source_sha256']=={s:sha(s)for s in(BASE,POLAR,CAPS)}
    assert p['source_sha256']==sha(BASE)and c['base_sha256']==sha(POLAR)
    Q=tuple(tuple(base(v)for v in row)for row in b['form'])
    assert len(Q)==21 and all(len(row)==21 for row in Q)
    assert all(x[1:]==(0,0,0)for row in Q for x in row)
    old=[tuple(map(base,v))for v in b['normals']]+[tuple(map(base,v['vector']))for v in p['cuts']]
    cov=[mv(Q,v)for v in old];norms=[dot(r,x)for r,x in zip(old,cov,strict=True)]
    h=tuple(map(base,b['strict_center']));hc=mv(Q,h);assert sign(dot(h,hc))<0
    w=tuple(map(base,p['strict_escape']['vector']));wc=mv(Q,w)
    assert sign(dot(w,wc))>0 and sign(dot(w,hc))<0
    r=tuple(map(radical,data['normal']));assert len(r)==21
    v=tuple(map(base,p['rejected_polar_proposals'][27]['vector']))
    assert r[:-1]==v[:-1]and r[-1]==(F(-4),F(0),F(1),F(0))
    rc=mv(Q,r);q=dot(r,rc)
    assert q==(F(76),F(0),F(-24),F(0))==radical(data['norm'])and sign(q)>0
    center=dot(r,hc);cut=dot(r,wc)
    assert center==rational(-110)==radical(data['strict_center_pairing'])
    assert cut==(F(1225,128),F(0),F(-3),F(0))==radical(data['old_escape_pairing'])and sign(cut)>0
    records=[]
    for i,(cc,n)in enumerate(zip(cov,norms,strict=True)):
        pp=dot(r,cc);kind=pair_label(pp,q,n);saved=data['old_pair_checks'][i]
        assert saved['old']==i and saved['kind']==kind and radical(saved['pairing'])==pp
        assert sign(dot(old[i],hc))<0
        records.append({'old':i,'kind':kind,'pairing':encode(pp)})
    extra=[]
    for i,cap in enumerate(c['caps'],82):
        d=base(cap['field']['radicand']);assert sign(d)>0
        root=tuple((base(x[0]),base(x[1]))for x in cap['vector'])
        pp=(dot(rc,tuple(x[0]for x in root)),dot(rc,tuple(x[1]for x in root)))
        assert extsign(pp,d)<0
        pp2=extmul(pp,pp,d);delta=(sub(pp2[0],mul(q,base(cap['norm']))),pp2[1])
        assert extsign(delta,d)>0
        saved=data['old_pair_checks'][i]
        assert saved['old']==i and saved['kind']=='dotted'
        assert tuple(map(radical,saved['pairing']))==pp
        records.append({'old':i,'kind':'dotted','pairing':[encode(x)for x in pp]})
        extra.append((d,root))
    assert len(records)==len(data['old_pair_checks'])==85
    counts=dict(sorted(Counter(x['kind']for x in records).items()))
    assert counts==data['pair_counts']=={'2':19,'4':3,'dotted':63}
    cusp_pairs=[dot(r,mv(Q,tuple(map(base,x['ray']))))for x in b['cusp_rows']]
    assert len(cusp_pairs)==107 and all(sign(x)<0 for x in cusp_pairs)
    assert cusp_pairs==list(map(radical,data['inherited_cusp_strict_pairings']))
    # A finite vertex from the two-dimensional orthogonal complement.
    t=(F(4),F(0),F(-1),F(0))
    x=tuple(mul(mul(rational(3),t),a)for a in v[:-1])+(rational(-2),)
    xc=mv(Q,x);xq=dot(x,xc)
    assert xq==mul(rational(-6),q)and sign(xq)<0 and sign(dot(x,hc))<0
    incidence=[]
    for i,cc in enumerate(cov):
        pp=dot(x,cc);assert sign(pp)<=0
        if pp==ZERO:incidence.append(i)
    expected=[0,1,2,3,5,6,7,9,10,12,14,15,16,17,18,19,20,22,23]
    assert incidence==expected and all(y[1:]==(0,0,0)for i in incidence for y in old[i])
    assert rankq([[y[0]for y in old[i]]for i in incidence])==19
    for d,root in extra:
        pp=(dot(xc,tuple(y[0]for y in root)),dot(xc,tuple(y[1]for y in root)))
        assert extsign(pp,d)<0
    assert dot(x,rc)==ZERO
    # The nineteen normals lie in x-perp, so their independent span is
    # positive definite. The new positive normal is orthogonal to that span.
    assert all(dot(r,cov[i])==ZERO for i in incidence)
    volume=read(VOLUME)
    assert volume['source_hashes']=={path:sha(path)for path in(BASE,POLAR,CAPS,SOURCE)}
    escape=tuple(map(base,volume['certificate']['vector']));ec=mv(Q,escape)
    eq=dot(escape,ec);eh=dot(escape,hc)
    assert sign(eq)>0 and sign(eh)<0
    assert eq==base(volume['certificate']['norm'])and eh==base(volume['certificate']['center_pairing'])
    assert all(sign(dot(escape,cc))<0 for cc in cov)and sign(dot(escape,rc))<0
    for d,root in extra:
        pp=(dot(ec,tuple(y[0]for y in root)),dot(ec,tuple(y[1]for y in root)))
        assert extsign(pp,d)<0
    return {'status':'PASS_STANDALONE_86_WALL_H20_CONSTRUCTION_AND_ESCAPE',
            'source':SOURCE,'source_sha256':sha(SOURCE),'dependencies':data['source_sha256'],
            'wall_count':86,'new_pair_counts':counts,'new_norm':encode(q),
            'strict_center_pairing':encode(center),'old_escape_cut_pairing':encode(cut),
            'inherited_cusps_strictly_preserved':107,'new_pairs':records,
            'finite_vertex':{'vector':list(map(encode,x)),'norm':encode(xq),
                             'incident_walls':incidence+[85],'incident_normal_rank':20},
            'strict_escape':{'source':VOLUME,'source_sha256':sha(VOLUME),'norm':encode(eq),
                             'center_pairing':encode(eh),'strict_inequalities':86},
            'scope':'Independent Fraction multiquadratic arithmetic and one-extra-radical sign tests verify the new facet,85 new pairs,107 cusp signs,old escape cut,and a genuine finite vertex. A new strict future spacelike escape proves infinite volume. Existing85-wall certificate remains a dependency; no exhaustive search conclusion.'}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    data=read(SOURCE);result=check(data)
    for field in('normal','norm','old_escape_pairing'):
        bad=deepcopy(data)
        if field=='normal':bad[field][-1][0][0]='-5'
        else:bad[field][0][0]=str(F(bad[field][0][0])+1)
        try:check(bad)
        except AssertionError:pass
        else:raise AssertionError('forged '+field+' accepted')
    result['forged_normal_norm_cut_rejected']=True
    if args.verify:assert read(OUT)==result
    else:(ROOT/OUT).write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['new_pair_counts'])
