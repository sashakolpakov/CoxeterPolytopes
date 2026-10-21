#!/usr/bin/env python3
"""Flip the B-region cap diagonal: impose native LR and allow PQ dotted.

Bounded numerical moving-contact discovery, full retained native gates and
complete S6 cone/PD-support comparison. No exact new normal is asserted.
"""
import argparse
import json
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import fresh_b_region_dotted7_four_cap_20261017 as F

N,D,ROOT,M=F.N,F.D,F.ROOT,F.M
TAG='fresh_b_region_lr_diagonal_20261017'
GRID=F.GRID;HIGH=F.HIGH;JOIN=F.JOIN


def topology():
    old=F.topology();fs=old['triangles'][:8]+[[67,68,70],[67,69,70]]
    bc=old['boundary'];es=Counter(tuple(sorted(e))for f in fs for e in combinations(f,2));links=[]
    assert len(set(sum(fs,[])))==8 and len(es)==17 and len(fs)==10
    assert (68,69)not in es and es[67,70]==2
    assert {e for e,n in es.items()if n==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
    assert set(es.values())=={1,2}
    for v in sorted(set(sum(fs,[]))):
        ee=[tuple(i for i in f if i!=v)for f in fs if v in f];ds=Counter(i for e in ee for i in e);seen={next(iter(ds))}
        while True:
            new=seen|{j for e in ee if set(e)&seen for j in e}
            if new==seen:break
            seen=new
        assert seen==set(ds)
        assert sorted(ds.values())==([1,1]+[2]*(len(ds)-2)if v in bc else[2]*len(ds))
        links.append(dict(vertex=v,edges=ee))
    return dict(boundary=bc,triangles=fs,vertices=8,edges=17,faces=10,links=links,
                removed_diagonal=[68,69],central_diagonal=[67,70])


def run():
    predecessor=json.loads((ROOT/f'results/{F.TAG}.json').read_text())
    for p,h in predecessor['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    s=N.Scout();disk=topology()
    vertices=[tuple(v)for v in s.primary['retained_old_finite_supports']]+[tuple(v['support'])for v in s.primary['new_finite_vertices']]
    debit=json.loads(json.dumps(D.debit(vertices,D.compact_incidence(vertices),[47])))
    assert debit==predecessor['chosen_debit']
    base={i:v for i,v in s.old.items()if i!=47};ids=[i for i in s.cur['actual_sides']if i!=47]
    deleted=s.score(base,{},ids,'delete47');assert s.ledgers[deleted]['open_ends']==10
    counts=Counter();blockers=Counter();unresolved=[];samples=[];brackets=[];events=[];completions=[];improvements=[];rdomains={};seeds=[]

    def solve(frame,pi,ms,nd,types):
        rhs=[-N.B.theta(m)-(s.profiles[pi]['offsets'][s.index[i]]if i in s.index else s.off(pi,types[i]))for i,m in zip(frame,ms)]
        return F.roots([nd[i]for i in frame],rhs,s.profiles[pi]['square'])

    def L(pi,sg,u):
        rhs=[-x-s.profiles[pi]['offsets'][s.index[i]]for i,x in zip([12,52,7],[0,0,u])]
        rr,status=F.roots([base[i]for i in[12,52,7]],rhs,s.profiles[pi]['square'])
        return next((v for j,v in rr if j==sg),None),status

    def native(v,pi,nd,types,stage):
        vals,bad,cp=s.native(v,pi,nd,types)
        if bad or cp>=-1e-8:blockers[stage,tuple(bad)]+=1;return False
        return True

    def recipes(frame,types,duties,nd):
        bonds=[2]*9
        for a,b in duties:
            rel=N.G.relation(s.dot(nd[frame[a]],s.ps(frame[a],types),nd[frame[b]],s.ps(frame[b],types)))
            assert isinstance(rel,int);bonds[a*3+b]=rel
        out=[]
        for pi,lo,hi in N.E.recipes(tuple(s.ps(i,types)for i in frame),tuple(bonds),tuple(duties)):
            for hh in product(HIGH,repeat=len(hi)):
                ms=list(lo)
                for i,m in zip(hi,hh):ms[i]=m
                for sg in[-1,1]:out.append((pi,tuple(ms),sg))
        return out

    def endpoint(frame,key,nd,ty):
        pi,ms,sg=key;rr,status=solve(frame,pi,ms,nd,ty)
        return next((v for j,v in rr if j==sg),None),status

    for lp in range(7):
        for ls in[-1,1]:
            previous={};prevu=None;pdomain=None
            for u in GRID:
                lv,status=L(lp,ls,u);counts['L_samples']+=1
                if status.startswith('unresolved'):unresolved.append(dict(stage='L',profile=lp,branch=ls,u=u,status=status))
                if lv is None or not native(lv,lp,base,{},'L'):
                    previous={};prevu=u;continue
                nd=dict(base);nd[67]=lv;ty={67:lp};counts['native_L_samples']+=1
                sid=s.score(nd,ty,ids+[67],'L')
                if sid is None:previous={};prevu=u;continue
                if pdomain is None:
                    pdomain=recipes([7,12,67],ty,[(0,1),(1,2)],nd)
                    seeds.append(dict(profile=lp,branch=ls,P_recipes=pdomain))
                current={};pn=rn=0;closest=None
                for pj,pk in enumerate(pdomain):
                    pv,status=endpoint([7,12,67],pk,nd,ty);counts['P_evaluations']+=1
                    if status.startswith('unresolved'):unresolved.append(dict(stage='P',profile=lp,branch=ls,u=u,recipe=pj,status=status))
                    if pv is None or not native(pv,pk[0],nd,ty,'P'):continue
                    pn+=1;counts['native_P_evaluations']+=1
                    # L is deliberately absent from this native screening set:
                    # its product with R is the moving equation being solved.
                    rbase=dict(base);rbase[68]=pv;rty={68:pk[0]};cachekey=(pk[0],pk[1][0])
                    if cachekey not in rdomains:rdomains[cachekey]=recipes([7,30,68],rty,[(0,1),(0,2)],rbase)
                    for rj,rk in enumerate(rdomains[cachekey]):
                        rv,status=endpoint([7,30,68],rk,rbase,rty);counts['R_screen_evaluations']+=1
                        if status.startswith('unresolved'):unresolved.append(dict(stage='R_screen',profile=lp,branch=ls,u=u,P=pj,R=rj,status=status))
                        if rv is None or not native(rv,rk[0],rbase,rty,'R_without_L'):continue
                        rn+=1;counts['native_R_without_L_evaluations']+=1
                        dot=s.dot(lv,lp,rv,rk[0])
                        for m in JOIN:
                            if N.P.corner3(lp,pk[0],rk[0],pk[1][2],m,rk[1][2])is None:continue
                            key=(pj,rj,m);gap=dot+N.B.theta(m);current[key]=gap
                            if closest is None or abs(gap)<abs(closest['gap']):closest=dict(key=key,gap=gap)
                            if key in previous and gap*previous[key]<0:
                                brackets.append(dict(L_profile=lp,L_branch=ls,P=pk,R=rk,LR_label=m,
                                    lower=prevu,upper=u,lower_gap=previous[key],upper_gap=gap))
                samples.append(dict(profile=lp,branch=ls,u=u,score=sid,P=pn,R_without_L=rn,LR_equations=len(current),closest=closest))
                previous=current;prevu=u
            print('MOVING_LR',lp,ls,'brackets',len(brackets),'counts',dict(counts),'unresolved',len(unresolved),flush=True)
    assert len(brackets)<4000,('bounded_event_guard',len(brackets))
    seen=set()
    for bi,br in enumerate(brackets):
        lp,ls=br['L_profile'],br['L_branch'];pk,rk,m=br['P'],br['R'],br['LR_label']
        def state(u):
            lv,status=L(lp,ls,u)
            if lv is None:raise ValueError(status)
            nd=dict(base);nd[67]=lv;ty={67:lp}
            pv,status=endpoint([7,12,67],pk,nd,ty)
            if pv is None:raise ValueError(status)
            nd[68]=pv;ty[68]=pk[0]
            rv,status=endpoint([7,30,68],rk,nd,ty)
            if rv is None:raise ValueError(status)
            return nd,ty,rv
        def gap(u):
            nd,ty,rv=state(u);return s.dot(nd[67],lp,rv,rk[0])+N.B.theta(m)
        try:u=float(brentq(gap,br['lower'],br['upper'],xtol=1e-12,maxiter=100));nd,ty,rv=state(u)
        except(ValueError,RuntimeError)as ex:unresolved.append(dict(stage='LR_event',bracket=bi,error=str(ex)));continue
        counts['LR_event_roots']+=1
        key=(lp,ls,pk[0],tuple(pk[1]),pk[2],rk[0],tuple(rk[1]),rk[2],m,round(u,9))
        if key in seen:continue
        seen.add(key)
        if not native(nd[67],lp,base,{},'event_L'):continue
        if not native(nd[68],pk[0],{**base,67:nd[67]},{67:lp},'event_P'):continue
        if not native(rv,rk[0],nd,ty,'event_R'):continue
        if N.G.relation(s.dot(nd[67],lp,rv,rk[0]))!=m:
            unresolved.append(dict(stage='event_spacing',bracket=bi,u=u,gap=gap(u)));continue
        nd[70]=rv;ty[70]=rk[0]
        sid=s.score(nd,ty,ids+[67,68,70],'LR_event')
        if sid is None:continue
        required=[[7,12,68],[12,52,67],[30,7,70],[12,68,67],[7,70,68],[67,68,70]]
        assert all(s.corner(t,nd,ty)is not None for t in required)
        rec=dict(index=len(events),u=u,L_profile=lp,L_branch=ls,P=pk,R=rk,LR_label=m,bracket=bi,score=sid,
            normals={str(i):nd[i].tolist()for i in[67,68,70]})
        events.append(rec)
        if s.ledgers[sid]['open_ends']<10:improvements.append(dict(stage='LR_event',event=rec['index'],score=sid))
        near_start=len(s.near)
        bank=s.bank([30,52,67],[(0,1),(1,2)],nd,ty,ids+[67,68,70],69,'Q_completion')
        for rejection in s.near[near_start:]:rejection['event']=rec['index']
        for q in bank['caps']:
            full,tt=s.add(nd,ty,[q],69);sid4=q['score'];corners=[s.corner(t,full,tt)for t in disk['triangles']];faces=[]
            for node in[67,68,69,70]:
                got=[x for x in s.ledgers[sid4]['rays']if node in x['active']]
                wanted={tuple(sorted(t))for t in disk['triangles']if node in t}
                faces.append(dict(node=node,finite_exactly_expected=all(x['kind']=='finite'for x in got)and{tuple(sorted(x['active']))for x in got}==wanted))
            done=all(c is not None for c in corners)and all(x['finite_exactly_expected']for x in faces)
            completions.append(dict(event=rec['index'],cap=q,finite_corners=corners,whole_cap_faces=faces,completed_disk=done))
            if s.ledgers[sid4]['open_ends']<10:improvements.append(dict(stage='Q_completion',record=len(completions)-1,score=sid4,completed_disk=done))
        print('LR_EVENT',rec['index'],'u',u,'profile',lp,'label',m,'ends',s.ledgers[sid]['open_ends'],'Q',len(bank['caps']),'improvements',len(improvements),flush=True)
    sources=dict(predecessor['sources'])
    for folder,ext in [('scripts','py'),('results','json'),('results','md')]:
        p=ROOT/f'{folder}/{F.TAG}.{ext}';sources[str(p.relative_to(ROOT))]=sha256(p.read_bytes()).hexdigest()
    sources[f'scripts/{TAG}.py']=sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(status='BOUNDED_B_REGION_LR_DIAGONAL_SCOUT',sources=sources,deleted_walls=[47],retained_old_wall_count=62,
        exact_deletion_debit=debit,disk=disk,control_score=s.control,deletion_score=deleted,grid=GRID,
        high_label_samples=HIGH,LR_labels=JOIN,seed_recipes=seeds,R_recipe_domains=[dict(key=k,recipes=v)for k,v in sorted(rdomains.items())],
        samples=samples,brackets=brackets,events=events,completions=completions,counts=dict(counts),bank_counts=dict(s.counts),
        boundary_ledgers=s.ledgers,histograms=[dict(stage=a,finite=b,exterior=c,open_ends=d,count=n)for(a,b,c,d),n in sorted(s.hist.items())],
        blockers=[dict(stage=a,walls=b,count=n)for(a,b),n in sorted(blockers.items())],
        bank_blockers=[dict(stage=a,walls=b,count=n)for(a,b),n in sorted(s.blockers.items())],
        near_native_rejections=s.near,
        improvements=improvements,unresolved=unresolved,bank_unresolved=s.unresolved,
        scope='New centralLR diagonal and nativeLR moving events; PQ is not a prescribed adjacency. Exact old47 debit and abstract disk only. Numerical grid/sign-bracket search uses all7 profiles, both norm signs, low2..5 and high6/8/12/24 samples. R screening omits only its unresolvedL contact; every accepted event/completion restores all62 old walls and all new native products/common center, completeS6cone and independentPD support agreement. Endpoint-native preselection can miss tangencies/thin windows; no all-parameter/all-label exclusion, exact new wall, ambient finite census or finite-volume proof. Original63/199/10 control unchanged.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],'events',len(out['events']),'Q',len(out['completions']),'improvements',len(out['improvements']),
          'unresolved',len(out['unresolved'])+len(out['bank_unresolved']),flush=True)
