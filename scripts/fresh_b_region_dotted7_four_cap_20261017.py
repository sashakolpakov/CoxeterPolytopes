#!/usr/bin/env python3
"""Changed B-region: a separate12/52 cap is dotted to7.

Bounded continuous native-contact events and full numerical cone scoring.
No exact new normal or all-parameter exclusion is asserted.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import fresh_release55_59_61_prism_replacement_20261017 as D
N=D.N;ROOT=N.ROOT;M=N.M;TAG='fresh_b_region_dotted7_four_cap_20261017'
HIGH=[6,8,12,24];JOIN=[2,3,4,5,6,8,12,24]
GRID=np.unique(np.r_[np.linspace(2.01,3.,21),np.linspace(3.,4.,9),6.]).tolist()


def topology():
    fs=[[7,12,68],[12,52,67],[52,30,69],[30,7,70],[12,68,67],
        [52,67,69],[30,69,70],[7,70,68],[67,68,69],[68,69,70]]
    bc=[7,12,52,30];es=Counter(tuple(sorted(e))for f in fs for e in combinations(f,2));vs=set(sum(fs,[]));links=[]
    assert(len(vs),len(es),len(fs))==(8,17,10)and set(es.values())=={1,2}
    assert {e for e,n in es.items()if n==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
    for v in sorted(vs):
        ee=[tuple(i for i in f if i!=v)for f in fs if v in f];ds=Counter(sum((list(e)for e in ee),[]));seen={next(iter(ds))}
        while True:
            new=seen|{j for i,j in ee if i in seen}|{i for i,j in ee if j in seen}
            if new==seen:break
            seen=new
        assert seen==set(ds)and sorted(ds.values())==([1,1]+[2]*(len(ds)-2)if v in bc else[2]*len(ds))
        links.append(dict(vertex=v,edges=ee))
    return dict(boundary=bc,triangles=fs,vertices=8,edges=17,faces=10,links=links,central_diagonal=[68,69])


def roots(frame,rhs,square):
    f=np.array(frame);g=f@M@f.T
    if abs(np.linalg.det(g))<1e-10:return [],'unresolved_singular'
    part=np.array(rhs)@np.linalg.inv(g)@f
    # Ordered cofactors give a continuous kernel orientation on every
    # nonsingular frame interval; no SVD sign changes enter root brackets.
    a=f@M;k=np.array([(-1.)**i*np.linalg.det(np.delete(a,i,axis=1))for i in range(4)])
    k/=np.linalg.norm(k);kq=float(k@M@k)
    if abs(kq)<1e-10:return [],'unresolved_null'
    ratio=(square-float(part@M@part))/kq
    if abs(ratio)<1e-10:return [],'unresolved_tangent'
    if ratio<0:return [],'no_real'
    out=[]
    for sg in[-1,1]:
        v=part+sg*np.sqrt(ratio)*k
        if max(abs(float(v@M@v)-square),max(abs(float(v@M@w)-r)for w,r in zip(f,rhs)))>2e-5:return [],'unresolved_residual'
        out.append((sg,v))
    return out,'two_real'


def run():
    s=N.Scout();disk=topology();vs=[tuple(v)for v in s.primary['retained_old_finite_supports']]+[tuple(v['support'])for v in s.primary['new_finite_vertices']]
    cs=D.compact_incidence(vs);costs=[D.debit(vs,cs,ds)for ds in[[47],[30],[52],[30,47],[47,53]]]
    cost=costs[0];assert cost['lost_finite_count']==2 and cost['compact_faces_by_dimension']=={2:1}and cost['lost_supports_containing8']==0
    base={i:v for i,v in s.old.items()if i!=47};ids=[i for i in s.cur['actual_sides']if i!=47]
    deleted=s.score(base,{},ids,'delete47');assert s.ledgers[deleted]['counts']=={'finite':24,'exterior':6}and s.ledgers[deleted]['open_ends']==10
    assert N.G.relation(s.dot(base[12],0,base[52],0))==6
    assert N.G.relation(s.dot(base[30],0,base[52],0))==4
    counts=Counter();blockers=Counter();unresolved=[];samples=[];brackets=[];events=[];completions=[];improvements=[];seed_equations=[]
    def solve(frame,pi,ms,nd,types):
        rhs=[-N.B.theta(m)-(s.profiles[pi]['offsets'][s.index[i]]if i in s.index else s.off(pi,types[i]))for i,m in zip(frame,ms)]
        return roots([nd[i]for i in frame],rhs,s.profiles[pi]['square'])
    def L(pi,sg,u):
        ff=[12,52,7];rhs=[-x-s.profiles[pi]['offsets'][s.index[i]]for i,x in zip(ff,[0,0,u])]
        rr,status=roots([base[i]for i in ff],rhs,s.profiles[pi]['square'])
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
            previous={};prevu=None;seed_done=False;domains=None
            for u in GRID:
                lv,status=L(lp,ls,u);counts['L_samples']+=1
                if status.startswith('unresolved'):unresolved.append(dict(stage='L',profile=lp,branch=ls,u=u,status=status))
                if lv is None or not native(lv,lp,base,{},'L'):
                    previous={};prevu=u;continue
                nd=dict(base);nd[67]=lv;ty={67:lp}
                score=s.score(nd,ty,ids+[67],'L');counts['native_L_samples']+=1
                if score is None:previous={};prevu=u;continue
                if domains is None:
                    domains=[recipes(f,ty,[(0,1),(1,2)],nd)for f in[[7,12,67],[30,52,67]]]
                    seed_equations.append(dict(profile=lp,branch=ls,frame=[12,52,7],native_contacts=['0','0','-u'],projected_square=s.profiles[lp]['square'],P_recipes=domains[0],Q_recipes=domains[1]))
                banks=[]
                for side,frame,dom in zip(['P','Q'],[[7,12,67],[30,52,67]],domains):
                    bank={}
                    for j,key in enumerate(dom):
                        v,status=endpoint(frame,key,nd,ty);counts['endpoint_evaluations']+=1
                        if status.startswith('unresolved'):unresolved.append(dict(stage=side,profile=lp,branch=ls,u=u,recipe=j,status=status))
                        if v is not None and native(v,key[0],nd,ty,side):bank[j]=v;counts['native_endpoint_evaluations']+=1
                    banks.append(bank)
                current={};min_gap=None
                for pj,pv in banks[0].items():
                    pp,pms,ps=domains[0][pj]
                    for qj,qv in banks[1].items():
                        qp,qms,qs=domains[1][qj];dot=s.dot(pv,pp,qv,qp)
                        for m in JOIN:
                            if N.P.corner3(lp,pp,qp,pms[2],qms[2],m)is None:continue
                            key=(pj,qj,m);gap=dot+N.B.theta(m);current[key]=gap
                            if min_gap is None or abs(gap)<abs(min_gap['gap']):min_gap=dict(key=key,gap=gap)
                            if key in previous and gap*previous[key]<0:brackets.append(dict(L_profile=lp,L_branch=ls,key=key,lower=prevu,upper=u,lower_gap=previous[key],upper_gap=gap,P=domains[0][pj],Q=domains[1][qj]))
                samples.append(dict(profile=lp,branch=ls,u=u,score=score,P=len(banks[0]),Q=len(banks[1]),join_equations=len(current),closest=min_gap))
                previous=current;prevu=u
            print('MOVING_L',lp,ls,'brackets',len(brackets),'counts',dict(counts),'unresolved',len(unresolved),flush=True)
    assert len(brackets)<4000,('event domain unexpectedly large',len(brackets))
    seen=set()
    for bi,br in enumerate(brackets):
        lp,ls=br['L_profile'],br['L_branch'];pk=br['P'];qk=br['Q'];m=br['key'][2]
        def state(u):
            lv,status=L(lp,ls,u)
            if lv is None:raise ValueError(status)
            nd=dict(base);nd[67]=lv;ty={67:lp}
            pv,ps=endpoint([7,12,67],pk,nd,ty);qv,qs=endpoint([30,52,67],qk,nd,ty)
            if pv is None or qv is None:raise ValueError((ps,qs))
            return nd,ty,pv,qv
        def gap(u):
            nd,ty,pv,qv=state(u);return s.dot(pv,pk[0],qv,qk[0])+N.B.theta(m)
        try:u=float(brentq(gap,br['lower'],br['upper'],xtol=1e-12,maxiter=100));nd,ty,pv,qv=state(u)
        except(ValueError,RuntimeError)as ex:unresolved.append(dict(stage='PQ_event',bracket=bi,error=str(ex)));continue
        counts['PQ_event_roots']+=1
        key=(lp,ls,tuple(pk[1]),pk[0],pk[2],tuple(qk[1]),qk[0],qk[2],m,round(u,9))
        if key in seen:continue
        seen.add(key)
        if not native(nd[67],lp,base,{},'event_L')or not native(pv,pk[0],nd,ty,'event_P')or not native(qv,qk[0],dict(nd,**{}),ty,'event_Q'):continue
        if N.G.relation(s.dot(pv,pk[0],qv,qk[0]))!=m:unresolved.append(dict(stage='event_spacing',bracket=bi,u=u,gap=gap(u)));continue
        nd[68]=pv;nd[69]=qv;ty[68]=pk[0];ty[69]=qk[0]
        sid=s.score(nd,ty,ids+[67,68,69],'PQ_event')
        if sid is None:continue
        corners=[[7,12,68],[12,52,67],[52,30,69],[12,68,67],[52,67,69],[67,68,69]]
        assert all(s.corner(t,nd,ty)is not None for t in corners)
        rec=dict(index=len(events),u=u,L_profile=lp,L_branch=ls,P=pk,Q=qk,PQ_label=m,bracket=bi,score=sid,normals={str(i):nd[i].tolist()for i in[67,68,69]})
        events.append(rec)
        if s.ledgers[sid]['open_ends']<10:improvements.append(dict(stage='PQ_event',event=rec['index'],score=sid))
        bank=s.bank([7,30,68],[(0,1),(0,2)],nd,ty,ids+[67,68,69],70,'R_completion')
        for r in bank['caps']:
            full,tt=s.add(nd,ty,[r],70);sid4=r['score'];required=[s.corner(t,full,tt)for t in disk['triangles']]
            faces=[]
            for node in[67,68,69,70]:
                got=[x for x in s.ledgers[sid4]['rays']if node in x['active']]
                wanted={tuple(sorted(t))for t in disk['triangles']if node in t}
                faces.append(dict(node=node,finite_exactly_expected=all(x['kind']=='finite'for x in got)and{tuple(sorted(x['active']))for x in got}==wanted))
            done=all(x is not None for x in required)and all(x['finite_exactly_expected']for x in faces)
            rr=dict(event=rec['index'],cap=r,finite_corners=required,whole_cap_faces=faces,completed_disk=done);completions.append(rr)
            if s.ledgers[sid4]['open_ends']<10:improvements.append(dict(stage='R_completion',record=len(completions)-1,score=sid4,completed_disk=done))
        print('PQ_EVENT',rec['index'],'u',u,'Lprofile',lp,'score',s.ledgers[sid]['counts'],s.ledgers[sid]['open_ends'],'R',len(bank['caps']),'improvements',len(improvements),flush=True)
    paths=[Path(__file__),Path(D.__file__),Path(N.__file__),Path(N.H.__file__),Path(N.E.__file__),Path(N.N.__file__),Path(N.S.__file__),Path(N.P.__file__),Path(N.G.__file__),N.H.PARENT,N.H.BOUNDARY]
    return dict(status='BOUNDED_DOTTED7_B_REGION_FOUR_CAP_SCOUT',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        deleted_walls=[47],retained_old_wall_count=62,exact_deletion_comparisons=costs,chosen_debit=cost,disk=disk,
        control_score=s.control,deletion_score=deleted,grid=GRID,high_label_samples=HIGH,PQ_labels=JOIN,
        seed_equations=seed_equations,samples=samples,brackets=brackets,events=events,completions=completions,
        counts=dict(counts),bank_counts=dict(s.counts),boundary_ledgers=s.ledgers,
        histograms=[dict(stage=a,finite=b,exterior=c,open_ends=d,count=n)for(a,b,c,d),n in sorted(s.hist.items())],
        blockers=[dict(stage=a,walls=b,count=n)for(a,b),n in sorted(blockers.items())],bank_blockers=[dict(stage=a,walls=b,count=n)for(a,b),n in sorted(s.blockers.items())],
        improvements=improvements,unresolved=unresolved,bank_unresolved=s.unresolved,
        scope='Exact old incidence debit and abstract disk; numerical conic norm reconstruction and sign-bracket search over the stated parameter grid only. All7 L profiles and both oriented conic branches are sampled; P/Q use all7 profiles, low2..5 and isolated high6/8/12/24. Both endpoint recipes must pass native gates at both bracket ends. Tangencies, thin admissible windows, other finite labels and parameters outside2.01..6 remain open. Every accepted joint event/R chamber passes all62 retained native walls, all new pairs, shared center and complete S6 cone/independent PD support agreement. No numerical native normal is exactified, no full ambient finite census or volume proof, and the exact63/199/10 control stays unchanged.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],'events',len(out['events']),'R',len(out['completions']),'improvements',len(out['improvements']),'unresolved',len(out['unresolved'])+len(out['bank_unresolved']),flush=True)
