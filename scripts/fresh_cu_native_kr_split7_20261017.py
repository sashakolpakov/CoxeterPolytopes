#!/usr/bin/env python3
"""Native-first KR search after subdividing the old 7--R duty.

Bounded numerical representative domain; no exact or all-label conclusion.
Frozen pentagonal records supply only C/A/L/P/K prefixes, never their
six rejected near-parallel R normals as accepted geometry.
"""
import argparse
import json
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path

import numpy as np
import fresh_cu_release47_pentagonal_cap_20261017 as F

U,E,H,B,M,N,S,P,ROOT=F.U,F.E,F.H,F.B,F.M,F.N,F.S,F.P,F.ROOT
TAG='fresh_cu_native_kr_split7_20261017'
HIGH=[6,8,12,24]


def relation(x):
    """Floating native gate with spacing-relative finite-label tolerance."""
    tol=2e-7
    if x < -2-tol:return 'dotted'
    if abs(x+2)<=tol:return 'parallel_or_unresolved'
    if abs(x)<=tol:return 2
    if x>0 or x < -2:return None
    m=int(round(np.pi/np.arccos(-x/2)))
    if m<3:return None
    spacing=min(abs(B.theta(m)-B.theta(m-1)),abs(B.theta(m+1)-B.theta(m)))
    return m if abs(x+B.theta(m))<=min(tol,spacing/100) else None


def disk():
    old=F.disk();removed=[{7,68,71},{7,70,71}]
    fs=[t for t in old['triangles']if set(t)not in removed]
    assert len(fs)==9
    fs += [[7,68,73],[68,71,73],[7,70,73],[70,71,73]]
    es=Counter(tuple(sorted(e))for t in fs for e in combinations(t,2))
    vs=set(sum(fs,[]));bc=old['boundary']
    assert (len(vs),len(es),len(fs))==(10,22,13)
    assert len(vs)-len(es)+len(fs)==1 and (7,71)not in es
    assert set(es.values())=={1,2}
    assert {e for e,n in es.items()if n==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
    links=[]
    for v in sorted(vs):
        ee=[tuple(i for i in t if i!=v)for t in fs if v in t]
        dd=Counter(i for e in ee for i in e);seen={next(iter(dd))}
        while True:
            new=seen|{b for a,b in ee if a in seen}|{a for a,b in ee if b in seen}
            if new==seen:break
            seen=new
        assert seen==set(dd)
        assert sorted(dd.values())==([1,1]+[2]*(len(dd)-2)if v in bc else[2]*len(dd))
        links.append(dict(vertex=v,edges=ee))
    return dict(boundary=bc,triangles=fs,vertices=10,edges=22,faces=13,
                links=links,removed_edge=[7,71],J_link=[7,68,71,70,7])


def select(parent):
    events={e['index']:e for e in parent['events']};groups=defaultdict(list)
    for r in parent['R_caps']:
        e=events[r['event']]
        key=(e['family'],parent['P_seeds'][e['P']]['labels'][0],
             parent['L_seeds'][e['L']]['profile'],parent['K_seeds'][e['K']]['profile'])
        groups[key].append(r)
    rows=[]
    for key,rs in sorted(groups.items()):
        # The product only ranks prefix geometry. No old R is reused here.
        rr=max(rs,key=lambda r:(dict(r['native_products'])[72],-r['event'],-r['index']))
        rows.append(dict(group=list(key),event=rr['event'],ranking_R=rr['index'],
                         ranking_KR=dict(rr['native_products'])[72],
                         available_events=len({r['event']for r in rs})))
    assert len(rows)==134 and len({r['group'][0]for r in rows})==24
    return rows,events


def run():
    path=ROOT/f'results/{F.TAG}.json';parent=json.loads(path.read_text())
    selections,events=select(parent);topology=disk()
    moving=json.loads((ROOT/f'results/{U.TAG}.json').read_text())
    cur,labels,base,h,profiles,*_=U.setup();base.pop(47);base.pop(67)
    oldindex={i:j for j,i in enumerate(labels)}
    actual=[i for i in cur['actual_sides']if i!=47]
    off=lambda a,b:float(H.A.native(np.zeros(4),a,np.zeros(4),b,profiles))
    def native(v,pi,nd,types):
        vals={i:float(v@M@w)+(profiles[pi]['offsets'][oldindex[i]]if i in oldindex else off(pi,types.get(i,0)))for i,w in nd.items()}
        bad=[i for i,x in vals.items()if relation(x)in[None,'parallel_or_unresolved']]
        return vals,bad,float(v@M@h)+profiles[pi]['center_offset']
    ledgers=[];ledgerkeys={}
    def score(nd,types):
        ids=actual+[i for i in sorted(nd)if i not in oldindex]
        bd=S.C.detailed_cone(nd,ids)
        sq=[profiles[types.get(i,0)]['square']if i not in oldindex else float(nd[i]@M@nd[i])for i in ids]
        ic=H.incidence_score([nd[i]for i in ids],ids,sq)
        good=not bd['warnings']and ic['status']=='numerically_separated_acute_incidence'
        good=good and bd['open_ends']==ic['open_ends']and bd['counts'].get('finite',0)==ic['finite']
        good=good and {tuple(sorted(v['active']))for v in bd['rays']if v['kind']=='finite'}=={tuple(sorted(v))for v in ic['finite_supports']}
        if not good:return None,dict(boundary=bd,PD_incidence=ic)
        rec=F.D.compact_score(dict(boundary=bd,PD_incidence=ic));key=json.dumps(rec,sort_keys=True)
        if key not in ledgerkeys:ledgerkeys[key]=len(ledgers);ledgers.append(rec)
        return ledgerkeys[key],None
    counts=Counter();blockers=Counter();hist=Counter();unresolved=[];records=[];improvements=[];near=[]
    for selection in selections:
        eid=selection['event'];ev=events[eid]
        nd=dict(base);nd.update({int(i):np.array(v)for i,v in ev['normals'].items()})
        ty={int(i):p for i,p in ev['profiles'].items()}
        ap,pp,kp=ty[68],ty[70],ty[72]
        kr=parent['K_seeds'][ev['K']];pr=parent['P_seeds'][ev['P']]
        fam=moving['moving_families'][ev['family']]
        assert fam['index']==ev['family'] and pp==0
        prefix_bad=[]
        for i in [67,68,69,70,72]:
            vals,bad,cp=native(nd[i],ty.get(i,0),{j:v for j,v in nd.items()if j!=i},ty)
            if bad or cp>=-1e-8:prefix_bad.append(dict(wall=i,blockers=[[j,vals[j]]for j in bad],center_product=cp))
        if prefix_bad:
            counts['strict_prefix_rejected']+=1;records.append(dict(selection=selection,prefix_rejected=prefix_bad));continue
        # A/P is a nonedge, so neither its label nor rightness enters this classifier.
        rb=[2]*9;rb[2]=ev['AK'];rb[5]=kr['labels'][2]
        rd=E.recipes((ap,pp,kp),tuple(rb),((0,2),(1,2)))
        rr=[];seen=set()
        for rp,low,hi in rd:
            for hh in product(HIGH,repeat=len(hi)):
                ms=list(low)
                for j,m in zip(hi,hh):ms[j]=m
                ids=[68,70,72];rhs=[-B.theta(m)-off(rp,ty[i])for i,m in zip(ids,ms)]
                roots,info=N.norm_roots([nd[i]for i in ids],rhs,profiles[rp]['square']);counts['R_equations']+=1
                if info['status'].startswith('unresolved'):
                    unresolved.append(dict(stage='R_norm',event=eid,profile=rp,labels=ms,norm=info));continue
                for sign,v in roots:
                    counts['R_branches']+=1
                    err=max(abs(float(v@M@v)-profiles[rp]['square']),max(abs(float(v@M@nd[i])-r)for i,r in zip(ids,rhs)))
                    if err>2e-5:unresolved.append(dict(stage='R_residual',event=eid,profile=rp,labels=ms,error=err));continue
                    vals,bad,cp=native(v,rp,nd,ty)
                    if bad or cp>=-1e-8:
                        blockers['R',tuple(bad)]+=1
                        if len(bad)<=2 and len(near)<150:near.append(dict(event=eid,profile=rp,labels=ms,branch=sign,blockers=[[i,vals[i]]for i in bad],center_product=cp))
                        continue
                    key=(rp,)+tuple(np.round(v,8))
                    if key in seen:counts['R_duplicates']+=1;continue
                    seen.add(key);nn=dict(nd);nn[71]=v;tt=dict(ty);tt[71]=rp
                    sid,issue=score(nn,tt)
                    if issue:unresolved.append(dict(stage='R_cone',event=eid,**issue));continue
                    rec=dict(index=len(rr),profile=rp,labels=ms,branch=sign,normal=v.tolist(),native_products=sorted(vals.items()),center_product=cp,norm_residual=err,score=sid)
                    rr.append(rec);counts['native_R']+=1
                    bd=ledgers[sid];hist['R',bd['counts']['finite'],bd['counts']['exterior'],bd['open_ends']]+=1
                    if bd['open_ends']<10:improvements.append(dict(stage='R',event=eid,R=rec))
        # J owns 7/A and7/P. Its bank is independent of R and constructed once.
        jb=[2]*9;jb[1]=fam['labels'][0];jb[2]=pr['labels'][0]
        jd=E.recipes((4,ap,pp),tuple(jb),((0,1),(0,2)))
        jj=[];joins=[];seen=set()
        if rr:
            for jp,low,hi in jd:
                for hh in product(HIGH,repeat=len(hi)):
                    ms=list(low)
                    for j,m in zip(hi,hh):ms[j]=m
                    ids=[7,68,70];rhs=[-B.theta(m)-off(jp,4 if i==7 else ty[i])for i,m in zip(ids,ms)]
                    roots,info=N.norm_roots([nd[i]for i in ids],rhs,profiles[jp]['square']);counts['J_equations']+=1
                    if info['status'].startswith('unresolved'):
                        unresolved.append(dict(stage='J_norm',event=eid,profile=jp,labels=ms,norm=info));continue
                    for sign,v in roots:
                        counts['J_branches']+=1
                        err=max(abs(float(v@M@v)-profiles[jp]['square']),max(abs(float(v@M@nd[i])-r)for i,r in zip(ids,rhs)))
                        if err>2e-5:unresolved.append(dict(stage='J_residual',event=eid,profile=jp,labels=ms,error=err));continue
                        vals,bad,cp=native(v,jp,nd,ty)
                        if bad or cp>=-1e-8:blockers['J',tuple(bad)]+=1;continue
                        key=(jp,)+tuple(np.round(v,8))
                        if key in seen:counts['J_duplicates']+=1;continue
                        seen.add(key)
                        jr=dict(index=len(jj),profile=jp,labels=ms,branch=sign,normal=v.tolist(),center_product=cp,norm_residual=err)
                        jj.append(jr);counts['native_J_prefix']+=1
                        for r in rr:
                            rv=np.array(r['normal']);rp=r['profile'];x=float(H.A.native(v,jp,rv,rp,profiles));rel=relation(x);counts['RJ_tests']+=1
                            if rel in[None,'parallel_or_unresolved']:blockers['RJ',(rel,)]+=1;continue
                            nn=dict(nd);nn.update({71:rv,73:v});tt=dict(ty);tt.update({71:rp,73:jp})
                            for i in [71,73]:
                                _,bad,cp=native(nn[i],tt[i],{j:w for j,w in nn.items()if j!=i},tt)
                                assert not bad and cp<-1e-8
                            sid,issue=score(nn,tt)
                            if issue:unresolved.append(dict(stage='RJ_cone',event=eid,R=r['index'],J=jr['index'],**issue));continue
                            finite=isinstance(rel,int)and P.corner3(ap,rp,jp,r['labels'][0],ms[1],rel)is not None and P.corner3(pp,rp,jp,r['labels'][1],ms[2],rel)is not None
                            bd=ledgers[sid];rayset={tuple(sorted(z['active']))for z in bd['rays']if z['kind']=='finite'}
                            alltriangles=all(tuple(sorted(t))in rayset for t in topology['triangles'])
                            whole={}
                            for i in [69,70,71,72,73]:
                                expected={tuple(sorted(t))for t in topology['triangles']if i in t};rs=[z for z in bd['rays']if i in z['active']]
                                whole[i]=all(z['kind']=='finite'for z in rs)and{tuple(sorted(z['active']))for z in rs}==expected
                            rec=dict(R=r['index'],J=jr['index'],RJ_product=x,RJ_relation=rel,last_corners_finite=finite,whole_faces=whole,completed_disk=finite and alltriangles and all(whole.values()),score=sid)
                            joins.append(rec);counts['coherent_RJ']+=1
                            if rec['completed_disk']:counts['completed_disks']+=1
                            hist['RJ',bd['counts']['finite'],bd['counts']['exterior'],bd['open_ends']]+=1
                            if bd['open_ends']<10:improvements.append(dict(stage='RJ',event=eid,R=r,J=jr,join=rec))
        records.append(dict(selection=selection,R_domain=[list(z)for z in rd],J_domain=[list(z)for z in jd],R_caps=rr,J_caps=jj,joins=joins))
        print('NATIVE_KR_SPLIT7',len(records),eid,'R',len(rr),'J',len(jj),'joins',len(joins),'improvements',len(improvements),'unresolved',len(unresolved),flush=True)
        assert counts['native_R']<10000 and counts['coherent_RJ']<20000,'bounded output guard'
    paths=[Path(__file__),Path(F.__file__),path,Path(U.__file__),ROOT/f'results/{U.TAG}.json',Path(E.__file__),Path(N.__file__),Path(P.__file__),Path(H.__file__),Path(B.__file__),Path(S.__file__)]
    return dict(status='BOUNDED_NATIVE_KR_SPLIT7_SCOUT',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        topology=topology,deleted_walls=[47],retained_control_walls=62,exact_deletion_ledger=parent['exact_deletion_ledger'],
        selection_rule='One prefix per nonempty(Afamily,P7label,Lprofile,Kprofile) group, maximizing prior native KR only as a geometric ranking; no old R normal reused.',
        selected_groups=len(selections),selected_events=len({r['event']for r in selections}),selected_A_families=sorted({r['group'][0]for r in selections}),selected_P_labels=sorted({r['group'][1]for r in selections}),
        high_label_samples=HIGH,finite_match_tolerance='min(2e-7, one hundredth of adjacent cosine spacing)',
        counts=dict(counts),records=records,boundary_ledgers=ledgers,improvements=improvements,unresolved=unresolved,near_native_failures=near,
        histograms=[dict(stage=s,finite=f,exterior=e,open_ends=i,count=n)for(s,f,e,i),n in sorted(hist.items())],
        blockers=[dict(stage=k[0],walls=list(k[1]),count=n)for k,n in sorted(blockers.items(),key=lambda x:str(x[0]))],
        scope='New incidence removes7/R, retainsKR with prescribed native finite label, addsJ.134 representative groups across24 A families andP3/P4; P5 has no saved pentagonal prefix in the chosen input. Seven A2profiles and both numerical norm branches, complete low finite-corner classifier with only6/8/12/24 high representatives. All62 retained old and all new native gates plus strict commoncenter; fullS6cone and independentPD support comparison for every coherentR/RJ chamber. No exact native proof, continuous/all-label exclusion, new ambient finite census or finite-volume assertion. The exact63wall199finite10end infinite-volume control is unchanged.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'improvements',len(out['improvements']),'unresolved',len(out['unresolved']),flush=True)
