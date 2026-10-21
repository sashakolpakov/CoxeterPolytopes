#!/usr/bin/env python3
"""Exact radical reconstruction of proposed non-right H20 polar repairs."""
import json
import hashlib
from pathlib import Path
from fractions import Fraction as Q
import math
from collections import Counter
import fresh_dimension_bridge_20261013 as K
from fresh_dimension_bridge_polar_20261013 import div
from scripts.fresh_h9_quartic_completion_20261013_arithmetic import Quadratic,dot,matvec,finite_label

ROOT=Path(__file__).resolve().parents[1]
SOURCE='data/enumerations/fresh_dimension_bridge_angle_repairs_20261013_scout.json'
OUT=ROOT/'data/enumerations/fresh_dimension_bridge_angle_repairs_20261013_exact.json'


class Base:
    zero,one=K.Z,K.O
    add=staticmethod(K.add);sub=staticmethod(K.sub);neg=staticmethod(K.neg)
    mul=staticmethod(K.mul);div=staticmethod(div);sign=staticmethod(K.sign)
    rational=staticmethod(K.k);encode=staticmethod(K.enc);parse=staticmethod(K.dec)

B=Base()
ANGLES={5:Quadratic(B,K.k(5)),8:Quadratic(B,K.k(2)),10:Quadratic(B,K.k(5)),12:Quadratic(B,K.k(3))}
POOL=[(2,[1,0],Q(0),Q(0)),(3,[1,-1],Q(1),Q(1)),(4,[1,-2],Q(2),Q(2)),(6,[1,-3],Q(3),Q(3)),
      (5,[1,-3,1],Q(26,10),Q(27,10)),(8,[1,-4,2],Q(34,10),Q(35,10)),
      (10,[1,-5,5],Q(36,10),Q(37,10)),(12,[1,-4,1],Q(37,10),Q(38,10))]


def lift(f,x):return x if f is B else f.lift(lift(f.base,x))
def value(f,x):return K.value(x)if f is B else value(f.base,x[0])+math.sqrt(value(f.base,f.d))*value(f.base,x[1])
def field_data(f):return {'kind':'sqrt6'}if f is B else{'kind':'quadratic','base':field_data(f.base),'radicand':f.base.encode(f.d)}


def line(rows,rhs=None):
    if rhs is None:rhs=[K.Z]*19+[K.k(-1)]
    a=[list(r)+[rhs[i]]for i,r in enumerate(rows)];piv=[];level=0
    for j in range(21):
        k=next((i for i in range(level,20)if a[i][j]!=K.Z),None)
        if k is None:continue
        a[level],a[k]=a[k],a[level];p=a[level][j];a[level]=[div(x,p)for x in a[level]]
        for i in range(20):
            if i!=level:
                p=a[i][j]
                if p!=K.Z:a[i]=[K.sub(x,K.mul(p,y))for x,y in zip(a[i],a[level])]
        piv.append(j);level+=1
        if level==20:break
    assert level==20
    free,=set(range(21))-set(piv);u=[K.Z]*21;v=[K.Z]*21;v[free]=K.O
    for i,j in enumerate(piv):u[j]=a[i][-1];v[j]=K.neg(a[i][free])
    assert all(K.dot(r,u)==rhs[i]and K.dot(r,v)==K.Z for i,r in enumerate(rows))
    return tuple(u),tuple(v)


def relation(f,p,q,n):
    assert f.sign(p)<=0
    ratio=f.div(f.mul(f.rational(4),f.mul(p,p)),f.mul(q,n))
    cmp=f.sign(f.sub(ratio,f.rational(4)))
    if cmp>0:return 'dotted'
    if cmp==0:return 'parallel'
    label=finite_label(f,ratio,POOL);assert label is not None,('forbidden exact angle',f.encode(ratio))
    return str(label)


def run():
    data=json.loads((ROOT/SOURCE).read_text());polar=json.loads((ROOT/data['source']).read_text())
    assert hashlib.sha256((ROOT/data['source']).read_bytes()).hexdigest()==data['source_sha256']
    base=json.loads((ROOT/polar['source']).read_text());form=tuple(tuple(map(K.dec,r))for r in base['form'])
    roots=list(map(K.parsev,base['normals']))+[K.parsev(r['vector'])for r in polar['cuts']]
    cov=[K.mv(form,r)for r in roots];norms=[K.dot(r,c)for r,c in zip(roots,cov)]
    h=K.parsev(base['strict_center']);hc=K.mv(form,h);cache={};caps=[];counts=Counter()
    for index,row in enumerate(data['candidates']):
        target=polar['rejected_polar_proposals'][row['target_index']]
        ids=[i for i in target['active_indices']if i!=row['released']]+[row['anchor']]
        key=tuple(ids)
        if key not in cache:cache[key]=line([cov[i]for i in ids])
        u,v=cache[key];m=row['label'];targetcov=K.mv(form,K.parsev(target['vector']))
        if m==2:
            q=K.pair(v,v,form);assert K.sign(q)>0
            f=B;candidates=[tuple(v),tuple(map(K.neg,v))]
        else:
            e=ANGLES.get(m,B)
            if e is B:cosq=B.rational({3:Q(1,4),4:Q(1,2),6:Q(3,4)}[m])
            else:
                a,b,den={5:(3,1,8),8:(2,1,4),10:(5,1,8),12:(2,1,4)}[m]
                cosq=(K.k(Q(a,den)),K.k(Q(b,den)))
            qe=e.div(e.one,e.mul(lift(e,norms[row['anchor']]),cosq))
            aa=lift(e,K.pair(v,v,form));bb=lift(e,K.sc(K.pair(u,v,form),2));cc=e.sub(lift(e,K.pair(u,u,form)),qe)
            dd=e.sub(e.mul(bb,bb),e.mul(e.rational(4),e.mul(aa,cc)));assert e.sign(dd)>=0 and e.sign(aa)
            f=Quadratic(e,dd);q=f.lift(qe);candidates=[]
            for s in(-1,1):
                t=(e.div(e.neg(bb),e.mul(e.rational(2),aa)),e.div(e.rational(s),e.mul(e.rational(2),aa)))
                candidates.append(tuple(f.add(lift(f,x),f.mul(t,lift(f,y)))for x,y in zip(u,v)))
        success=[]
        for r in candidates:
            if f.sign(dot(f,r,tuple(lift(f,x)for x in hc)))>=0:continue
            if f.sign(dot(f,r,tuple(lift(f,x)for x in targetcov)))<=0:continue
            qq=q
            cv=tuple(dot(f,tuple(lift(f,x)for x in rr),r)for rr in form)
            assert f.sign(f.sub(dot(f,r,cv),qq))==0
            pairs=[dot(f,r,tuple(lift(f,x)for x in c))for c in cov]
            try:kinds=[relation(f,p,qq,lift(f,n))for p,n in zip(pairs,norms)]
            except AssertionError:continue
            normalized=[value(f,x)*math.sqrt(2/value(f,qq))for x in r]
            error=max(abs(x-y)for x,y in zip(normalized,row['normal']))
            if error>1e-5:continue
            success.append({'scout_index':index,'field':field_data(f),'vector':[f.encode(x)for x in r],
                'norm':f.encode(qq),'active_indices':ids,'released':row['released'],'anchor':row['anchor'],'anchor_label':m,
                'pair_kinds':kinds,'center_pairing':f.encode(dot(f,r,tuple(lift(f,x)for x in hc))),
                'target_pairing':f.encode(dot(f,r,tuple(lift(f,x)for x in targetcov))),
                'normalized_approximation_error':error})
        assert len(success)==1,(index,len(success))
        caps.extend(success);counts[str(m)]+=1
        if len(caps)%20==0:print({'certified':len(caps)},flush=True)
    OUT.write_text(json.dumps({'status':'EXACT_H20_ANGLE_REPAIR_CAP_BANK','source':SOURCE,
        'source_sha256':hashlib.sha256((ROOT/SOURCE).read_bytes()).hexdigest(),
        'base_source':data['source'],'base_sha256':data['source_sha256'],
        'caps':caps,'counts_by_anchor_label':dict(counts),'linear_systems':len(cache),
        'scope':'Every saved cap is exactly positive norm, oriented, cuts its indicated rejected polar ray, and is Coxeter-compatible with all82 retained walls. Branch equations and signs use exact radical towers. The cap entries may repeat mirrors and have not yet been required mutually compatible. No finite-volume or exhaustive-search claim.'},indent=2)+'\n')
    print({'done':len(caps),'counts':dict(counts)},flush=True)


if __name__=='__main__':run()
