#!/usr/bin/env python3
"""Three mutually compatible exact two-angle caps for the H20 cone."""
import json
import hashlib
import math
from pathlib import Path
from fractions import Fraction as Q
import fresh_dimension_bridge_20261013 as K
from fresh_dimension_bridge_angle_repairs_20261013_exact import B,Quadratic,line,div,lift,dot,relation,value,field_data

ROOT=Path(__file__).resolve().parents[1]
SOURCE='data/enumerations/fresh_dimension_bridge_two_angle_20261013_scout.json'
OUT=ROOT/'data/enumerations/fresh_dimension_bridge_two_angle_20261013_exact.json'


def run():
    data=json.loads((ROOT/SOURCE).read_text());polar=json.loads((ROOT/data['base_source']).read_text())
    base=json.loads((ROOT/polar['source']).read_text());form=tuple(tuple(map(K.dec,r))for r in base['form'])
    roots=list(map(K.parsev,base['normals']))+[K.parsev(r['vector'])for r in polar['cuts']]
    cov=[K.mv(form,r)for r in roots];norms=[K.dot(r,c)for r,c in zip(roots,cov)]
    h=K.parsev(base['strict_center']);hc=K.mv(form,h);escape=K.parsev(polar['strict_escape']['vector']);ec=K.mv(form,escape)
    chosen=[0,580,914];caps=[];live=[]
    for index in chosen:
        row=data['candidates'][index];assert row['labels']==[6,6]
        target=polar['rejected_polar_proposals'][row['target_index']]
        ids=[i for i in target['active_indices']if i not in row['released']]+row['anchors']
        a,b=row['anchors'];kappa2=div(norms[b],norms[a]);assert kappa2[1]==0
        n,d=math.isqrt(kappa2[0].numerator),math.isqrt(kappa2[0].denominator)
        kappa=Q(n,d);assert kappa*kappa==kappa2[0]
        rhs=[K.Z]*18+[K.k(-1),K.k(-kappa)];u,v=line([cov[i]for i in ids],rhs)
        q=div(K.k(Q(4,3)),norms[a]);aa=K.pair(v,v,form);bb=K.sc(K.pair(u,v,form),2);cc=K.sub(K.pair(u,u,form),q)
        dd=K.sub(K.mul(bb,bb),K.sc(K.mul(aa,cc),4));assert K.sign(dd)>0 and K.sign(aa)
        f=Quadratic(B,dd);candidates=[]
        for s in(-1,1):
            t=(div(K.neg(bb),K.sc(aa,2)),div(K.k(s),K.sc(aa,2)))
            r=tuple(f.add(f.lift(x),f.mul(t,f.lift(y)))for x,y in zip(u,v))
            if f.sign(dot(f,r,tuple(map(f.lift,hc))))>=0:continue
            rc=tuple(dot(f,tuple(map(f.lift,rr)),r)for rr in form)
            assert f.sign(f.sub(dot(f,r,rc),f.lift(q)))==0
            pp=[dot(f,r,tuple(map(f.lift,c)))for c in cov]
            try:kinds=[relation(f,p,f.lift(q),f.lift(n))for p,n in zip(pp,norms)]
            except AssertionError:continue
            targetpair=dot(f,r,tuple(map(f.lift,K.mv(form,K.parsev(target['vector'])))))
            if f.sign(targetpair)<=0:continue
            assert f.sign(dot(f,r,tuple(map(f.lift,ec))))<0
            normalized=[value(f,x)*math.sqrt(2/K.value(q))for x in r]
            assert max(abs(x-y)for x,y in zip(normalized,row['normal']))<1e-5
            candidates.append((r,rc,kinds,targetpair))
        assert len(candidates)==1
        r,rc,kinds,targetpair=candidates[0];live.append((f,r,rc,q))
        caps.append({'scout_index':index,'target_index':row['target_index'],'field':field_data(f),
             'vector':list(map(f.encode,r)),'norm':K.enc(q),'active_indices':ids,'prescribed_pairings':list(map(K.enc,rhs)),
             'old_pair_kinds':kinds,'center_pairing':f.encode(dot(f,r,tuple(map(f.lift,hc)))),
             'target_pairing':f.encode(targetpair),'escape_pairing':f.encode(dot(f,r,tuple(map(f.lift,ec))))})
    cross=[]
    for i,(f,r,rc,q)in enumerate(live):
        for j,(g,s,sc,qq)in enumerate(live[:i]):
            common=Quadratic(f,f.lift(g.d))
            se=tuple((f.lift(x[0]),f.lift(x[1]))for x in s)
            ce=tuple(common.lift(x)for x in rc)
            pp=dot(common,se,ce);kind=relation(common,pp,common.lift(f.lift(q)),common.lift(f.lift(qq)))
            assert kind=='dotted';cross.append({'indices':[j,i],'kind':kind,'pairing':common.encode(pp)})
    OUT.write_text(json.dumps({'status':'EXACT_COHERENT_85_WALL_H20_CONE_WITH_STRICT_ESCAPE',
        'source':SOURCE,'source_sha256':hashlib.sha256((ROOT/SOURCE).read_bytes()).hexdigest(),
        'base_source':data['base_source'],'base_sha256':hashlib.sha256((ROOT/data['base_source']).read_bytes()).hexdigest(),
        'caps':caps,'new_cross_pairs':cross,'wall_count':85,'strict_escape_from_base':polar['strict_escape'],
        'scope':'Three exact caps retain all82 old walls, cut three specified hyperideal rays, and are mutually dotted. Strict interior and surviving strict escape are exact. Each has18 prescribed right angles and two pi/6 anchor angles. This is an infinite-volume85-facet H20 construction; neither full cap enumeration nor finite volume is claimed.'},indent=2)+'\n')
    print({'wall_count':85,'new_caps':len(caps),'new_cross_pairs':len(cross),'norms':[K.enc(t[3])for t in live]},flush=True)


if __name__=='__main__':run()
