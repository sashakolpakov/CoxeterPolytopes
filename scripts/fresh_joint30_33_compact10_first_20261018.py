#!/usr/bin/env python3
"""Reverse the joint hinge construction, prescribing compact10's missing corner."""
import argparse,json
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_joint30_33_hinge_native_pilot_20261018 as H

T,C,ROOT=H.T,H.C,H.ROOT
TAG='fresh_joint30_33_compact10_first_20261018'


def added_frame_checks(s):
    """Original native interval identities for the newly used old44 contacts."""
    p=s.primary;dec=lambda z:tuple(map(Q,z))
    add=lambda a,b:(a[0]+b[0],a[1]+b[1])
    def mul(a,b):
        v=[x*y for x in a for y in b];return min(v),max(v)
    ns={};h=list(map(dec,p['common_center']));g=[[dec(x)for x in r]for r in p['computational_gram']]
    for r in p['facet_witnesses']:
        hp=dec(p['center_products'][str(r['wall'])]);assert hp[1]<0
        ns[r['wall']]=[mul((2*(a[0]-b[1]),2*(a[1]-b[0])),(1/hp[1],1/hp[0]))for a,b in zip(h,map(dec,r['point']))]
    records=[]
    for i,j in[(7,44),(10,44),(44,55),(44,66)]:
        if j==66:q=dec(next(r for r in p['native_gates']if r['wall']==i)['product'])
        else:
            q=(Q(0),Q(0))
            for a,x in enumerate(ns[i]):
                for b,y in enumerate(ns[j]):q=add(q,mul(mul(x,g[a][b]),y))
        val=s.dot(s.old[i],s.ps(i,{}),s.old[j],s.ps(j,{}))
        assert float(q[0])-2e-5<=val<=float(q[1])+2e-5
        records.append(dict(walls=[i,j],profiles=[s.ps(i,{}),s.ps(j,{})],native_product=val,
            native_relation=T.G.relation(val),exact_native_enclosure=list(map(str,q))))
    assert records[1]['native_relation']==2 and records[0]['native_relation']=='dotted'
    return records


def bank(s,base,ids,frame,duties,extra,newid,stage,oldwall):
    bonds=[2]*9;frame_pairs=[]
    for a,b in combinations(range(3),2):
        val=s.dot(base[frame[a]],s.ps(frame[a],{}),base[frame[b]],s.ps(frame[b],{}))
        rel=T.G.relation(val);frame_pairs.append(dict(walls=[frame[a],frame[b]],native_product=val,relation=rel))
        if(a,b)in duties:assert isinstance(rel,int);bonds[3*a+b]=rel
    dom=T.E.recipes(tuple(s.ps(i,{})for i in frame),tuple(bonds),tuple(duties))
    recipes=[(pi,list(ms))for pi,ms,high in dom if not high and all(2<=m<=5 for m in ms)]
    direct=[(pi,list(ms))for pi in range(7)for ms in product(range(2,6),repeat=3)if all(T.P.corner3(s.ps(frame[a],{}),s.ps(frame[b],{}),pi,bonds[3*a+b],ms[a],ms[b])is not None for a,b in duties)]
    assert sorted(recipes)==sorted(direct)
    tests=[];caps=[];seen=set()
    for pi,ms in recipes:
        rhs=[-T.B.theta(m)-s.profiles[pi]['offsets'][s.index[i]]for i,m in zip(frame,ms)]
        roots,info=T.N.norm_roots([base[i]for i in frame],rhs,s.profiles[pi]['square'])
        row=dict(index=len(tests),profile=pi,labels=ms,norm=info,branches=[]);tests.append(row);s.counts['norm_equations']+=1
        if info['status'].startswith('unresolved'):
            s.unresolved.append(dict(stage=stage+'_norm',recipe=row['index'],info=info));continue
        for sign,v in roots:
            s.counts['norm_branches']+=1
            err=max(abs(float(v@T.M@v)-s.profiles[pi]['square']),max(abs(float(v@T.M@base[i])-q)for i,q in zip(frame,rhs)))
            vals,bad,cp=s.native(v,pi,base,{})
            br=dict(branch=sign,norm_residual=err,center_product=cp,blockers=[[i,vals[i]]for i in bad]);row['branches'].append(br)
            if frame==[10,44,7]and pi==0 and ms==[2,2,5]:
                oldid=min((8,30),key=lambda i:float(np.max(abs(v-s.old[i]))))
                difference=float(np.max(abs(v-s.old[oldid])))
                oldcp=float(s.old[oldid]@T.M@s.center)+s.profiles[s.ps(oldid,{})]['center_offset']
                assert s.ps(oldid,{})==pi and difference<2e-6 and abs(cp-oldcp)<2e-6 and cp<0
                br['retained_control_comparison']=dict(old_wall=oldid,profile=pi,normal=v.tolist(),old_normal=s.old[oldid].tolist(),
                    maximum_coordinate_difference=difference,old_center_product=oldcp,
                    scope='Numerical full intrinsic-vector equality plus identical zero white profile; no new physical interval certificate.')
            if err>2e-5:
                s.unresolved.append(dict(stage=stage+'_residual',recipe=row['index'],branch=sign,error=err));continue
            if bad or cp>=-1e-8:continue
            nn={**base,newid:v};tt={newid:pi}
            assert all(s.corner([frame[a],frame[b],newid],nn,tt)is not None for a,b in duties)
            failed=[tri for tri in extra if s.corner(tri,nn,tt)is None];br['extra_corner_failures']=failed
            if failed:continue
            key=(pi,)+tuple(np.round(v,8))
            if key in seen:br['duplicate']=True;continue
            seen.add(key);before=len(s.unresolved);sid=s.score(nn,tt,ids+[newid],stage)
            br['cap']=len(caps)
            caps.append(dict(index=len(caps),recipe=row['index'],profile=pi,labels=ms,branch=sign,normal=v.tolist(),
                native_products=sorted(vals.items()),center_product=cp,norm_residual=err,score=sid,
                cone_validated=sid is not None,cone_diagnostic=None if sid is not None else s.unresolved[before:],
                restores_old_numerically=pi==s.ps(oldwall,{})and float(np.max(abs(v-s.old[oldwall])))<2e-6))
    return dict(new_wall=newid,replaced_old_wall=oldwall,frame=frame,frame_pairs=frame_pairs,duties=duties,extra_corners=extra,
        complete_low_recipe_count=len(recipes),symbolic_high_families_not_tested=sum(bool(h)for _,_,h in dom),recipes=tests,caps=caps)


def run():
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    interface=H.old_frame_interface(s);added=added_frame_checks(s)
    debit=T.H.deletion(s.primary,s.cur,[30,33]);assert debit['lost_finite_count']==15 and debit['compact_faces_by_dimension']=={2:10,3:1}
    assert all(8 not in v for v in debit['lost_finite_supports'])
    base={i:v for i,v in s.old.items()if i not in(30,33)};ids=[i for i in s.cur['actual_sides']if i not in(30,33)]
    control=s.score(s.old,{},s.cur['actual_sides'],'control63');deleted=s.score(base,{},ids,'delete30_33');assert control is not None
    b30=bank(s,base,ids,[10,44,7],[(0,1)],[[7,55,67]],67,'new30',30)
    print('COMPACT10_NEW30',len(b30['recipes']),len(b30['caps']),flush=True)
    b33=bank(s,base,ids,[7,10,66],[(0,1),(0,2)],[],68,'new33',33)
    print('COMPACT10_NEW33',len(b33['recipes']),len(b33['caps']),flush=True)
    assert sum(c['restores_old_numerically']for c in b30['caps'])==1
    assert sum(c['restores_old_numerically']for c in b33['caps'])==1
    pairs=[];digest=sha256();hist=Counter();coherent=[]
    required=[[10,44,67],[7,55,67],[7,10,68],[7,66,68],[10,67,68]]
    for x,y in product(b30['caps'],b33['caps']):
        val=s.dot(np.array(x['normal']),x['profile'],np.array(y['normal']),y['profile']);rel=T.G.relation(val)
        row=dict(index=len(pairs),new30=x['index'],new33=y['index'],native_product=val,relation=rel,score=None);pairs.append(row)
        hist[str(rel)]+=1;digest.update(json.dumps([x['index'],y['index'],val,rel],separators=(',',':')).encode())
        if rel in(None,'parallel_or_unresolved'):continue
        nn={**base,67:np.array(x['normal']),68:np.array(y['normal'])};tt={67:x['profile'],68:y['profile']}
        for i in(67,68):
            values,bad,cp=s.native(nn[i],tt[i],{j:v for j,v in nn.items()if j!=i},tt);assert not bad and cp<-1e-8
        sid=s.score(nn,tt,ids+[67,68],'joint_pair');row['score']=sid
        if sid is None:continue
        l=s.ledgers[sid];failed=[f for f in required if s.corner(f,nn,tt)is None]
        faces=[H.face(l,i)for i in(10,55,67,68)]
        coherent.append(dict(pair=row['index'],score=sid,old30_restored=x['restores_old_numerically'],old33_restored=y['restores_old_numerically'],
            finite_rays=l['counts'].get('finite',0),open_ends=l['open_ends'],failed_corners=failed,
            compact10_preserved=faces[0]['all_rays_finite'],compact55_preserved=faces[1]['all_rays_finite'],whole_face_records=faces,
            exterior_supports=[r['active']for r in l['rays']if r['kind']=='exterior'],open_edges=[e for e in l['edges']if e['kind']in('one_open','two_open')]))
    paths=[Path(__file__),Path(H.__file__),Path(C.__file__),Path(T.__file__),Path(T.H.__file__),Path(T.E.__file__),Path(T.P.__file__),Path(T.N.__file__),Path(T.G.__file__),
        Path(T.H.A.__file__),Path(T.H.A.C.__file__),T.H.A.INPUT,T.H.A.C.PARENT,T.H.PARENT,T.H.BOUNDARY,Path(T.S.__file__),Path(T.S.C.__file__)]
    return dict(status='BOUNDED_JOINT30_33_COMPACT10_FIRST',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        deleted_walls=[30,33],retained_old_walls=61,exact_deletion_debit=debit,old_frame_interface=interface,added_old44_frame_checks=added,
        control_score=control,deletion_score=deleted,low_labels=[2,3,4,5],banks=[b30,b33],required_corners=required,pairs=pairs,
        complete_pair_count=len(pairs),complete_pair_digest=digest.hexdigest(),pair_histogram=dict(hist),coherent_pairs=coherent,
        boundary_ledgers=s.ledgers,counts=dict(s.counts),histograms=[dict(stage=a,finite=b,exterior=c,open_ends=d,count=n)for(a,b,c,d),n in sorted(s.hist.items())],unresolved=s.unresolved,
        scope='Changed reverse-order native reconstruction: new30 first in10/44/7, preserving finite10/44/new30 and7/55/new30; new33 complete low7/10/66 bank joined on all pairs. No artificial7/44 incidence. All seven inherited profiles and2..5 labels, both norm branches, all61 retained old walls/allmutual native products/strictcenter. Intermediate nonsimple exterior diagnostics stay provisional; every coherent final pair has clean complete numericalS6 cone/PD support agreement and explicit whole10/55/newwall records. Exact old15/10/1 support debit; no exact new normal, high/continuous completeness, ambient finite census/volume or newdimension claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:assert not p.exists();p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],'banks',[len(b['caps'])for b in d['banks']],'pairs',len(d['pairs']),'coherent',len(d['coherent_pairs']),'unresolved',len(d['unresolved']),flush=True)
