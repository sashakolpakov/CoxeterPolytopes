#!/usr/bin/env python3
"""Changed A-K events for a K-pentagon replacing the old AL/LR adjacency.

Numerical discovery on all116 saved moving-A branches and98 fixed-L seeds.
High labels are only6/8/12/24; no all-label or exact chamber assertion.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_cu_release47_quadratic_events_20261017 as D
U,E,H,B,M,N,S,P,ROOT=D.U,D.E,D.H,D.B,D.M,D.N,D.S,D.P,D.ROOT
TAG='fresh_cu_release47_pentagonal_cap_20261017';HIGH=[6,8,12,24]

def disk():
    old=D.Q.disk_check();removed=[{30,68,69},{68,69,71},{70,69,71}]
    fs=[t for t in old['triangles']if set(t)not in removed]
    assert len(fs)==6
    fs += [[30,68,72],[30,69,72],[69,70,72],[70,71,72],[71,68,72]]
    es=Counter(tuple(sorted(e))for t in fs for e in combinations(t,2));vs=set(sum(fs,[]));bc=old['boundary']
    assert(len(vs),len(es),len(fs))==(9,19,11)and len(vs)-len(es)+len(fs)==1
    assert(68,69)not in es and(69,71)not in es and set(es.values())=={1,2}
    assert{e for e,c in es.items()if c==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
    links=[]
    for v in sorted(vs):
        ee=[tuple(i for i in t if i!=v)for t in fs if v in t];dd=Counter(i for e in ee for i in e);seen={next(iter(dd))}
        while True:
            new=seen|{b for a,b in ee if a in seen}|{a for a,b in ee if b in seen}
            if new==seen:break
            seen=new
        assert seen==set(dd)and sorted(dd.values())==([1,1]+[2]*(len(dd)-2)if v in bc else[2]*len(dd))
        links.append(dict(vertex=v,edges=ee))
    return dict(boundary=bc,triangles=fs,vertices=9,edges=19,faces=11,links=links,removed_edges=[[68,69],[69,71]],K_link=[30,69,70,71,68,30])

def run():
    topology=disk();dp=ROOT/'results'/f'{D.TAG}.json';up=ROOT/'results'/f'{U.TAG}.json'
    old=json.loads(dp.read_text());moving=json.loads(up.read_text())
    cur,labels,base,h,profiles,primary,ef,ei,e0,e1,gamma,t0,cv=U.setup();base.pop(47);base.pop(67)
    oldindex={i:j for j,i in enumerate(labels)};actual=[i for i in cur['actual_sides']if i!=47]
    off=lambda a,b:float(H.A.native(np.zeros(4),a,np.zeros(4),b,profiles))
    def native(v,pi,nd,types):
        vals={i:float(v@M@w)+(profiles[pi]['offsets'][oldindex[i]]if i in oldindex else off(pi,types.get(i,0)))for i,w in nd.items()}
        bad=[i for i,x in vals.items()if B.relation(x)in[None,'parallel_or_unresolved']]
        return vals,bad,float(v@M@h)+profiles[pi]['center_offset']
    ledgers=[];ledgerkeys={}
    def score(nd,types):
        ids=actual+[i for i in sorted(nd)if i not in oldindex];bd=S.C.detailed_cone(nd,ids)
        sq=[profiles[types.get(i,0)]['square']if i not in oldindex else float(nd[i]@M@nd[i])for i in ids]
        ic=H.incidence_score([nd[i]for i in ids],ids,sq)
        good=not bd['warnings']and ic['status']=='numerically_separated_acute_incidence'and bd['open_ends']==ic['open_ends']
        good=good and bd['counts'].get('finite',0)==ic['finite']and{tuple(sorted(v['active']))for v in bd['rays']if v['kind']=='finite'}=={tuple(sorted(v))for v in ic['finite_supports']}
        if not good:return None,dict(boundary=bd,PD_incidence=ic)
        compact=D.compact_score(dict(boundary=bd,PD_incidence=ic));key=json.dumps(compact,sort_keys=True)
        if key not in ledgerkeys:ledgerkeys[key]=len(ledgers);ledgers.append(compact)
        return ledgerkeys[key],None
    counts=Counter();blockers=Counter();unresolved=[];banks=[];domains=[];events=[];fourth=[];improvements=[];hist=Counter()
    pbank=old['P_seeds'];lbank=old['L_seeds'];assert len(lbank)==98
    for lr in lbank:
        li=lr['index'];lp=lr['profile'];l=np.array(lr['normal']);pj=lr['P'];pr=pbank[pj];pv=np.array(pr['normal']);pp=pr['profile'];assert pp==0
        nd=dict(base);nd[69]=l;nd[70]=pv;ty={69:lp,70:pp}
        for i in[69,70]:
            _,bad,cp=native(nd[i],ty[i],{j:v for j,v in nd.items()if j!=i},ty);assert not bad and cp<-1e-8
        bonds=[2]*9;bonds[1]=lr['labels'][0];bonds[5]=lr['labels'][2]
        # 30/P is not an edge required by either finite K corner.
        dom=E.recipes((0,lp,pp),tuple(bonds),((0,1),(1,2)));domains.append(dict(L=li,symbolic_recipes=[list(z)for z in dom]))
        seen=set()
        for kp,low,hi in dom:
            for hh in product(HIGH,repeat=len(hi)):
                ms=list(low)
                for j,m in zip(hi,hh):ms[j]=m
                ids=[30,69,70];rhs=[-B.theta(m)-off(kp,ty.get(i,0))for i,m in zip(ids,ms)]
                vv,info=N.norm_roots([nd[i]for i in ids],rhs,profiles[kp]['square']);counts['K_equations']+=1
                if info['status'].startswith('unresolved'):unresolved.append(dict(stage='K_norm',L=li,profile=kp,labels=ms,norm=info));continue
                for sg,k in vv:
                    counts['K_branches']+=1;vals,bad,cp=native(k,kp,nd,ty)
                    if bad or cp>=-1e-8:blockers['K',tuple(bad)]+=1;continue
                    err=max(abs(float(k@M@k)-profiles[kp]['square']),max(abs(float(k@M@nd[i])-r)for i,r in zip(ids,rhs)))
                    if err>2e-5:unresolved.append(dict(stage='K_residual',L=li,profile=kp,labels=ms,error=err));continue
                    key=(kp,)+tuple(np.round(k,9))
                    if key in seen:counts['K_duplicates']+=1;continue
                    seen.add(key);banks.append(dict(index=len(banks),L=li,P=pj,profile=kp,labels=ms,branch=sg,normal=k.tolist(),native_products=sorted(vals.items()),center_product=cp,norm_residual=err))
        print('PENTAGON_K_BANK',li,'caps',len(banks),'unresolved',len(unresolved),flush=True)
    assert len(banks)<5000,('bounded bank guard',len(banks))
    seen=set()
    for fam in moving['moving_families']:
        fid=fam['index'];ap=fam['profile'];a0,at,ai=map(np.array,[fam['constant'],fam['positive_t'],fam['inverse_t']]);am7,am30,amc=fam['labels']
        for kr in banks:
            ki=kr['index'];kp=kr['profile'];k=np.array(kr['normal']);lr=lbank[kr['L']];lp=lr['profile'];l=np.array(lr['normal']);pj=kr['P'];pr=pbank[pj];pp=pr['profile'];pv=np.array(pr['normal']);pm7=pr['labels'][0]
            allowed=[m for m in range(2,7)if P.corner3(0,ap,kp,am30,kr['labels'][0],m)is not None]
            allowed=[m for m in allowed if m!=6]+(HIGH if 6 in allowed else[])
            ca,cb,cc=float(at@M@k),float(ai@M@k),float(a0@M@k)+off(ap,kp)
            for akm in allowed:
                counts['A_K_quadratics']+=1
                if max(abs(ca),abs(cb),abs(cc+B.theta(akm)))<1e-10:
                    unresolved.append(dict(stage='identical_AK',family=fid,K=ki,label=akm));continue
                for t in U.equation_roots(ca,cb,cc+B.theta(akm)):
                    if t<t0-1e-9:continue
                    counts['positive_parameter_roots']+=1;u=float((t-1/t)/(2*gamma));a=a0+t*at+ai/t;c=cv(t)
                    nd=dict(base);nd.update({67:c,68:a,69:l,70:pv,72:k});types={68:ap,69:lp,70:pp,72:kp}
                    err=max(abs(float(a@M@a)-profiles[ap]['square']),abs(float(c@M@c)-2),abs(float(H.A.native(a,ap,k,kp,profiles))+B.theta(akm)),abs(float(a@M@c)+B.theta(amc)))
                    if err>2e-5:unresolved.append(dict(stage='event_residual',family=fid,K=ki,u=u,error=err));continue
                    gates={};failed=False
                    for name,i in[('C',67),('A',68),('L',69),('P',70),('K',72)]:
                        vals,bad,cp=native(nd[i],types.get(i,0),{j:w for j,w in nd.items()if j!=i},types);gates[name]=dict(native_products=sorted(vals.items()),center_product=cp)
                        if bad or cp>=-1e-8:blockers['event_'+name,tuple(bad)]+=1;failed=True
                    if failed:continue
                    key=(ap,lp,kp)+tuple(np.round(np.r_[c,a,l,pv,k],8))
                    if key in seen:counts['event_duplicates']+=1;continue
                    seen.add(key);scoreid,issue=score(nd,types)
                    if issue:unresolved.append(dict(stage='event_cone',family=fid,K=ki,u=u,**issue));continue
                    ev=dict(index=len(events),family=fid,K=ki,L=kr['L'],P=pj,u=u,t=float(t),AK=akm,normals={str(i):nd[i].tolist()for i in[67,68,69,70,72]},profiles=types,norm_residual=err,gates=gates,score=scoreid)
                    events.append(ev);counts['native_AK_events']+=1;bd=ledgers[scoreid];hist['event',bd['counts']['finite'],bd['counts']['exterior'],bd['open_ends']]+=1
                    if bd['open_ends']<10:improvements.append(dict(stage='event',index=ev['index'],score=scoreid))
                    rb=[2]*9;rb[1]=am7;rb[2]=pm7;rd=E.recipes((4,ap,pp),tuple(rb),((0,1),(0,2)))
                    for rp,low,hi in rd:
                        for hh in product(HIGH,repeat=len(hi)):
                            ms=list(low)
                            for j,m in zip(hi,hh):ms[j]=m
                            ids=[7,68,70];rhs=[-B.theta(m)-off(rp,4 if i==7 else types[i])for i,m in zip(ids,ms)]
                            vv,info=N.norm_roots([nd[i]for i in ids],rhs,profiles[rp]['square']);counts['R_equations']+=1
                            if info['status'].startswith('unresolved'):unresolved.append(dict(stage='R_norm',event=ev['index'],profile=rp,labels=ms,norm=info));continue
                            for sg,r in vv:
                                counts['R_branches']+=1;vals,bad,cp=native(r,rp,nd,types)
                                if bad or cp>=-1e-8:blockers['R',tuple(bad)]+=1;continue
                                er=max(abs(float(r@M@r)-profiles[rp]['square']),max(abs(float(r@M@nd[i])-v)for i,v in zip(ids,rhs)))
                                if er>2e-5:unresolved.append(dict(stage='R_residual',event=ev['index'],error=er));continue
                                nn=dict(nd);nn[71]=r;ty=dict(types);ty[71]=rp;sid,issue=score(nn,ty)
                                if issue:unresolved.append(dict(stage='R_cone',event=ev['index'],**issue));continue
                                rr=B.relation(vals[72]);closed=isinstance(rr,int)and P.corner3(ap,rp,kp,ms[1],akm,rr)is not None and P.corner3(pp,rp,kp,ms[2],kr['labels'][2],rr)is not None
                                rayset={tuple(sorted(v['active']))for v in ledgers[sid]['rays']if v['kind']=='finite'}
                                alltriangles=all(tuple(sorted(t))in rayset for t in topology['triangles'])
                                whole={}
                                for i in[69,70,71,72]:
                                    expected={tuple(sorted(t))for t in topology['triangles']if i in t};rs=[v for v in ledgers[sid]['rays']if i in v['active']]
                                    whole[i]=all(v['kind']=='finite'for v in rs)and{tuple(sorted(v['active']))for v in rs}==expected
                                complete=closed and alltriangles and all(whole.values())
                                rec=dict(index=len(fourth),event=ev['index'],profile=rp,labels=ms,branch=sg,normal=r.tolist(),native_products=sorted(vals.items()),center_product=cp,norm_residual=er,KR_relation=rr,required_last_corners_finite=closed,whole_faces=whole,completed_disk=complete,score=sid)
                                fourth.append(rec);counts['native_R_caps']+=1
                                if complete:counts['completed_disks']+=1
                                bd=ledgers[sid];hist['R',bd['counts']['finite'],bd['counts']['exterior'],bd['open_ends']]+=1
                                if bd['open_ends']<10:improvements.append(dict(stage='R',index=rec['index'],score=sid))
            assert len(events)<20000 and len(fourth)<100000,('bounded output guard',len(events),len(fourth))
        print('PENTAGON_A_FAMILY',fid,'events',len(events),'R',len(fourth),'improved',len(improvements),'unresolved',len(unresolved),flush=True)
    paths=[Path(__file__),Path(D.__file__),dp,Path(U.__file__),up,Path(E.__file__),Path(N.__file__),Path(P.__file__),Path(B.__file__),Path(H.__file__),Path(S.__file__)]
    return dict(status='NUMERICAL_RELEASE47_PENTAGONAL_K_EVENTS',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        topology=topology,deleted_walls=[47],retained_control_walls=62,exact_deletion_ledger=old['exact_deletion_ledger'],
        P_seeds=pbank,L_seeds=lbank,K_domains=domains,K_seeds=banks,moving_A_families=116,high_label_samples=HIGH,
        counts=dict(counts),events=events,R_caps=fourth,boundary_ledgers=ledgers,improvements=improvements,unresolved=unresolved,
        histograms=[dict(stage=s,finite=f,exterior=e,open_ends=i,count=n)for(s,f,e,i),n in sorted(hist.items())],
        blockers=[dict(stage=k[0],walls=list(k[1]),count=n)for k,n in sorted(blockers.items())],
        scope='Changed incidence removes both AL and LR, gives K wholepentagon30/L/P/R/A.98 fixedL seeds, all116 movingA branches; K solved once per(P,L), then nativeAK quadratics foru>=3 without a parametergrid/uppercutoff. Highlabels only6/8/12/24, floating cancellation/tangencies/near-singularities not intervalcertified. Every accepted event/R addition passes all62 retainedold/new native gates and strictcommoncenter, complete numericalS6cone/independentPDsupports. Wholecap faces checked for completion; no exact chamber/ambientfinitecensus/volume claim. Allpredecessors preserved.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'improvements',len(out['improvements']),'unresolved',len(out['unresolved']),flush=True)
