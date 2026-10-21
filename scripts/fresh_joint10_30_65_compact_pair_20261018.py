#!/usr/bin/env python3
"""Six selected positive65 prefixes; a new N/J quadrilateral/pentagon duty."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import numpy as np
import fresh_joint10_30_native_replacement_20261018 as Z

T,H,ROOT=Z.T,Z.H,Z.ROOT
TAG='fresh_joint10_30_65_compact_pair_20261018'
SEED=ROOT/f'results/{Z.TAG}.json'


def topology():
    boundary=[8,44,63,64,52]
    ts=[[8,44,68],[44,63,68],[8,68,69],[63,68,69],[8,52,69],[52,64,69],[63,64,69]]
    edges=Counter(tuple(sorted(e)) for t in ts for e in combinations(t,2));vertices=set(sum(ts,[]))
    assert (len(vertices),len(edges),len(ts))==(7,13,7) and len(vertices)-len(edges)+len(ts)==1
    assert {e for e,n in edges.items() if n==1}=={tuple(sorted(e)) for e in zip(boundary,boundary[1:]+boundary[:1])}
    assert set(edges.values())=={1,2};links=[]
    for v in sorted(vertices):
        es=[tuple(i for i in t if i!=v) for t in ts if v in t];degrees=Counter(i for e in es for i in e)
        seen={next(iter(degrees))}
        while True:
            new=seen|{b for a,b in es if a in seen}|{a for a,b in es if b in seen}
            if new==seen:break
            seen=new
        assert seen==set(degrees)
        assert sorted(degrees.values())==([2]*len(degrees) if v in (68,69) else [1,1]+[2]*(len(degrees)-2))
        links.append(dict(vertex=v,edges=es))
    return dict(boundary=boundary,triangles=ts,vertices=7,edges=13,faces=7,links=links)


def bank(s,nd,types,ids):
    frame=[8,63,68];duties=[(0,2),(1,2)];bonds=[2]*9;checks=[]
    for a,b in combinations(range(3),2):
        v=s.dot(nd[frame[a]],s.ps(frame[a],types),nd[frame[b]],s.ps(frame[b],types));rel=T.G.relation(v)
        checks.append(dict(walls=[frame[a],frame[b]],native_product=v,relation=rel))
        if (a,b) in duties:assert isinstance(rel,int);bonds[3*a+b]=rel
    dom=T.E.recipes(tuple(s.ps(i,types) for i in frame),tuple(bonds),tuple(duties))
    recipes=[(pi,list(ms)) for pi,ms,high in dom if not high and all(2<=m<=5 for m in ms)]
    direct=[(pi,list(ms)) for pi in range(7) for ms in product(range(2,6),repeat=3)
        if all(T.P.corner3(s.ps(frame[a],types),s.ps(frame[b],types),pi,bonds[3*a+b],ms[a],ms[b]) is not None for a,b in duties)]
    assert sorted(recipes)==sorted(direct)
    tests=[];caps=[];seen=set()
    for pi,ms in recipes:
        rhs=[-T.B.theta(m)-(s.profiles[pi]['offsets'][s.index[i]] if i in s.index else s.off(pi,types[i])) for i,m in zip(frame,ms)]
        roots,info=T.N.norm_roots([nd[i] for i in frame],rhs,s.profiles[pi]['square'])
        row=dict(index=len(tests),profile=pi,labels=ms,norm=info,branches=[]);tests.append(row);s.counts['J_equations']+=1
        if info['status'].startswith('unresolved'):
            s.unresolved.append(dict(stage='J_norm',recipe=row['index'],info=info));continue
        for sign,v in roots:
            s.counts['J_branches']+=1
            err=max(abs(float(v@T.M@v)-s.profiles[pi]['square']),max(abs(float(v@T.M@nd[i])-q) for i,q in zip(frame,rhs)))
            vals,bad,cp=s.native(v,pi,nd,types)
            br=dict(branch=sign,normal=v.tolist(),norm_residual=err,center_product=cp,native_products=sorted(vals.items()),blockers=[[i,vals[i]] for i in bad])
            row['branches'].append(br)
            if err>2e-5:
                s.unresolved.append(dict(stage='J_residual',recipe=row['index'],branch=sign,error=err));continue
            if bad or cp>=-1e-8:continue
            key=(pi,)+tuple(np.round(v,8))
            if key in seen:br['duplicate']=True;continue
            seen.add(key);nn={**nd,69:v};tt={**types,69:pi}
            assert all(s.corner([frame[a],frame[b],69],nn,tt) is not None for a,b in duties)
            for new in (67,68,69):
                vv,bb,cc=s.native(nn[new],tt[new],{i:x for i,x in nn.items() if i!=new},tt);assert not bb and cc<-1e-8
            failed=[tri for tri in topology()['triangles'] if s.corner(tri,nn,tt) is None]
            before=len(s.unresolved);sid=s.score(nn,tt,ids+[69],'complete_J')
            raw=s.ledgers[sid] if sid is not None else s.unresolved[before]['boundary']
            rawfaces=[H.face(raw,i) for i in (8,11,33,44,55,63,64,67,68,69)]
            def rawwhole(v):
                f=next(f for f in rawfaces if f['wall']==v)
                wanted={tuple(sorted(t)) for t in topology()['triangles'] if v in t}
                return f['all_rays_finite'] and {tuple(sorted(r['active'])) for r in f['rays']}==wanted
            br['cap']=len(caps);caps.append(dict(index=len(caps),recipe=row['index'],profile=pi,labels=ms,branch=sign,
                normal=v.tolist(),native_products=sorted(vals.items()),center_product=cp,norm_residual=err,
                failed_disk_corners=failed,score=sid,cone_diagnostic=None if sid is not None else s.unresolved[before:],
                raw_whole_face_records=rawfaces,raw_whole_N=rawwhole(68),raw_whole_J=rawwhole(69),
                raw_face_scope='Whole-face records on the saved numerical cone; a nonsimple warning elsewhere prevents strict global acceptance. Local N closure is not an accepted chamber or finite-volume certificate.'))
    return dict(frame=frame,duties=duties,frame_native_checks=checks,symbolic_high_families_not_tested=sum(bool(h) for _,_,h in dom),
        complete_low_recipe_count=len(recipes),recipes=tests,caps=caps)


def run():
    old=json.loads(SEED.read_text())
    for p,h in old['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    interface=H.old_frame_interface(s);checks=Z.changed_frame_checks(s)
    debit=T.H.deletion(s.primary,s.cur,[10,30,65]);prior=old['exact_deletion_debit']
    assert debit['lost_finite_count']==80 and debit['compact_faces_by_dimension']=={2:314,3:325,4:244,5:129,6:46,7:10,8:1}
    assert sum(8 in r for r in debit['lost_finite_supports'])==46
    added=sorted(set(map(tuple,debit['lost_finite_supports']))-set(map(tuple,prior['lost_finite_supports'])))
    assert len(added)==5 and all(65 in v for v in added)
    base={i:v for i,v in s.old.items() if i not in (10,30,65)};ids=[i for i in s.cur['actual_sides'] if i not in (10,30,65)]
    control=s.score(s.old,{},s.cur['actual_sides'],'control63');assert control is not None
    deleted=s.score(base,{},ids,'delete10_30_65')
    # Positive65 candidates are selected by least conflict separately in each
    # of the seven physical white profiles, then only0/3/6 are used for banks.
    inventory=[];selected=[]
    for bi in (0,1):
        rows=[]
        for r in old['new10_banks'][bi]['recipes']:
            for b in r['branches']:
                if len(b['blockers'])==1 and b['blockers'][0][0]==65 and b['center_product']<-1e-8:
                    rows.append((b['blockers'][0][1],r['profile'],r['labels'],b['branch'],r['index'],b))
        assert len(rows)==94
        minimum=[]
        for pi in range(7):
            v,_,ms,sign,ri,b=min(z for z in rows if z[1]==pi)
            assert ms==[2,2,2] and sign==-1
            rec=dict(R_prefix=bi,N_profile=pi,source_recipe=ri,branch=sign,labels=ms,old65_product=v,selected=pi in (0,3,6))
            minimum.append(rec)
            if rec['selected']:selected.append(rec)
        inventory.append(dict(R_prefix=bi,sole65_branches=94,profile_minima=minimum))
    assert len(selected)==6
    records=[];coherent=[];digest=sha256()
    for choose in selected:
        bi,pi=choose['R_prefix'],choose['N_profile'];seed=old['prefixes'][bi]
        rhs=[-T.B.theta(m)-s.profiles[0]['offsets'][s.index[i]] for i,m in zip([7,55,33],seed['labels'])]
        roots,info=T.N.norm_roots([s.old[i] for i in [7,55,33]],rhs,s.profiles[0]['square'])
        R=next(v for sign,v in roots if sign==seed['branch']);assert np.max(abs(R-np.array(seed['normal'])))<2e-6
        rhs=[-s.profiles[pi]['offsets'][s.index[i]] for i in [8,44,63]]
        roots,infoN=T.N.norm_roots([s.old[i] for i in [8,44,63]],rhs,s.profiles[pi]['square'])
        N=next(v for sign,v in roots if sign==-1)
        saved=old['new10_banks'][bi]['recipes'][choose['source_recipe']]
        saved=next(b for b in saved['branches'] if b['branch']==-1)
        assert np.max(abs(N-np.array(saved['normal'])))<2e-6
        before={i:v for i,v in s.old.items() if i not in (10,30)};before[67]=R
        vals,bad,cp=s.native(N,pi,before,{67:0});assert bad==[65] and abs(vals[65]-choose['old65_product'])<2e-7
        nd={**base,67:R,68:N};types={67:0,68:pi};actual=ids+[67,68]
        native=[]
        for new in (67,68):
            vv,bb,cc=s.native(nd[new],types[new],{i:v for i,v in nd.items() if i!=new},types)
            assert not bb and cc<-1e-8
            native.append(dict(wall=new,profile=types[new],normal=nd[new].tolist(),native_products=sorted(vv.items()),center_product=cc))
        assert s.corner([8,44,68],nd,types) is not None and s.corner([44,63,68],nd,types) is not None
        before=len(s.unresolved);sid=s.score(nd,types,actual,'selected_prefix')
        rec=dict(index=len(records),selection=choose,R_reconstruction=info,N_reconstruction=infoN,native=native,
            score=sid,cone_diagnostic=None if sid is not None else s.unresolved[before:])
        b=bank(s,nd,types,actual);rec['J_bank']=b;records.append(rec)
        print('NEW65',rec['index'],bi,pi,'recipes',len(b['recipes']),'caps',len(b['caps']),flush=True)
        for cap in b['caps']:
            digest.update(json.dumps([rec['index'],cap['index'],cap['profile'],cap['labels'],cap['branch'],cap['native_products']],separators=(',',':')).encode())
            if cap['score'] is None:continue
            bd=s.ledgers[cap['score']];fs=[H.face(bd,i) for i in (8,11,33,44,55,63,64,67,68,69)]
            def whole(v):
                f=next(f for f in fs if f['wall']==v);wanted={tuple(sorted(t)) for t in topology()['triangles'] if v in t}
                return f['all_rays_finite'] and {tuple(sorted(r['active'])) for r in f['rays']}==wanted
            coherent.append(dict(record=rec['index'],cap=cap['index'],score=cap['score'],finite_rays=bd['counts'].get('finite',0),
                open_ends=bd['open_ends'],whole_faces=fs,whole_N=whole(68),whole_J=whole(69),failed_disk_corners=cap['failed_disk_corners'],
                affected_old_compact_faces_finite=all(f['all_rays_finite'] for f in fs if f['wall'] in (8,11,44,55,63,64)),
                exterior_supports=[r['active'] for r in bd['rays'] if r['kind']=='exterior'],open_edges=[e for e in bd['edges'] if e['kind'] in ('one_open','two_open')]))
    # Separate prescribed original65 control. Its44 label15 is recovered as
    # a native gate, not quietly omitted because the new banks use labels<=5.
    oldpi=s.ps(65,{});frame=[8,63,64]
    labels=[T.G.relation(s.dot(s.old[65],oldpi,s.old[i],s.ps(i,{}))) for i in frame];assert labels==[2,2,2]
    rhs=[-s.profiles[oldpi]['offsets'][s.index[i]] for i in frame]
    roots,info=T.N.norm_roots([s.old[i] for i in frame],rhs,s.profiles[oldpi]['square']);ctrl=[]
    for sign,v in roots:
        nn={i:w for i,w in s.old.items() if i!=65};vv,bb,cc=s.native(v,oldpi,nn,{})
        row=dict(branch=sign,normal=v.tolist(),profile=oldpi,native_products=sorted(vv.items()),blockers=[[i,vv[i]] for i in bb],center_product=cc,score=None)
        if not bb and cc<-1e-8:
            assert np.max(abs(v-s.old[65]))<2e-6 and T.G.relation(vv[44])==15
            nn[69]=v;row['score']=s.score(nn,{69:oldpi},[i for i in s.cur['actual_sides'] if i!=65]+[69],'old65_restoration')
            assert row['score'] is not None
        ctrl.append(row)
    assert sum(r['score'] is not None for r in ctrl)==1
    paths=[Path(__file__),Path(Z.__file__),SEED,ROOT/f'results/{Z.TAG}.md',Path(H.__file__),Path(T.__file__),Path(T.H.__file__),
        Path(T.E.__file__),Path(T.P.__file__),Path(T.N.__file__),Path(T.G.__file__),Path(T.S.__file__),Path(T.S.C.__file__),
        Path(T.H.A.__file__),Path(T.H.A.C.__file__),T.H.A.INPUT,T.H.A.C.PARENT,T.H.PARENT,T.H.BOUNDARY]
    return dict(status='BOUNDED_JOINT10_30_65_COMPACT_PAIR',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        deleted_walls=[10,30,65],retained_old_walls=60,exact_deletion_debit=debit,additional_lost_supports=added,
        old_frame_interface=interface,old8_44_63_frame_checks=checks,numerical_replay_snapshot=old['numerical_replay_snapshot'],
        topology=topology(),selection_inventory=inventory,selected_prefixes=selected,control_score=control,deletion_score=deleted,
        records=records,complete_native_triple_count=sum(len(r['J_bank']['caps']) for r in records),complete_native_triple_digest=digest.hexdigest(),
        coherent_triples=coherent,old65_control=dict(frame=frame,labels=labels,norm=info,branches=ctrl),
        boundary_ledgers=s.ledgers,counts=dict(s.counts),unresolved=s.unresolved,
        scope='Six fixed selected prefixes: two frozen R branches times least65-positive-conflict N in profiles0/3/6. Every J bank exhausts seven profiles and frame labels2..5 for the two required finite N-edge duties. All60 retained old walls and all mutual native/strictcenter gates precede complete numerical S6 cone/PD checks. WholeN quadrilateral/J pentagon, all disk corners and affected compact-face ledgers are separate necessary completion checks. Exact80 oldfinite-support debit; provisional nonsimple diagnostics never accepted. High/continuous families, exact new normals and ambient H9 volume are not certified.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=json.loads(json.dumps(run()))
    path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'triples',out['complete_native_triple_count'],'coherent',len(out['coherent_triples']),'unresolved',len(out['unresolved']),flush=True)
