#!/usr/bin/env python3
"""Changed-A release47 bridges via fixed-L quadratic contact events.

Fix L from30/52/P; A(t).L=-theta is a quadratic in t. This removes a
parameter grid, but floating quadratic roots and sampled high labels do
not certify a continuous exclusion or exact geometry.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_cu_release47_rotated_bridge_20261017 as Q
U,E,H,B,M,N,S,P,ROOT=Q.U,Q.E,Q.H,Q.B,Q.M,Q.N,Q.S,Q.P,Q.ROOT
TAG='fresh_cu_release47_quadratic_events_20261017'
HIGH=[6,8,12,24]

def compact_score(sc):
    b=sc['boundary'];i=sc['PD_incidence']
    return dict(counts=b['counts'],open_ends=b['open_ends'],warnings=b['warnings'],
      rays=[dict(active=r['active'],kind=r['kind'])for r in b['rays']],
      edges=[dict(sides=e['sides'],kind=e['kind'])for e in b['edges']],
      PD_status=i['status'],PD_open_ends=i['open_ends'],PD_finite_supports=i['finite_supports'])

def run():
    src=ROOT/'results'/f'{Q.TAG}.json';seed=json.loads(src.read_text())
    upath=ROOT/'results'/f'{U.TAG}.json';moving=json.loads(upath.read_text())
    cur,oldlabels,base,h,profiles,primary,ef,ei,e0,e1,gamma,t0,cv=U.setup()
    base.pop(47);base.pop(67);actual=[i for i in cur['actual_sides']if i!=47]
    oldindex={i:j for j,i in enumerate(oldlabels)}
    off=lambda a,b:float(H.A.native(np.zeros(4),a,np.zeros(4),b,profiles))
    def native(v,pi,nd,types):
        vals={i:float(v@M@w)+(profiles[pi]['offsets'][oldindex[i]]if i in oldindex else off(pi,types.get(i,0)))for i,w in nd.items()}
        bad=[i for i,x in vals.items()if B.relation(x)in[None,'parallel_or_unresolved']]
        cp=float(v@M@h)+profiles[pi]['center_offset']
        return vals,bad,cp
    def score(nd,types):
        sel=actual+[i for i in sorted(nd)if i not in oldindex]
        bd=S.C.detailed_cone(nd,sel)
        sq=[profiles[types.get(i,0)]['square']if i not in oldindex else float(nd[i]@M@nd[i])for i in sel]
        ic=H.incidence_score([nd[i]for i in sel],sel,sq)
        good=not bd['warnings']and ic['status']=='numerically_separated_acute_incidence'
        good=good and bd['open_ends']==ic['open_ends']and bd['counts'].get('finite',0)==ic['finite']
        good=good and{tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'}=={tuple(sorted(r))for r in ic['finite_supports']}
        return dict(boundary=bd,PD_incidence=ic),good
    counts=Counter();blockers=Counter();unresolved=[];banks=[];events=[];third=[];improvements=[];hist=Counter()
    pbank=seed['P_seeds'];domains=[]
    for pj,pr in enumerate(pbank):
        pv=np.array(pr['normal']);pp=pr['profile'];assert pp==0
        vals,bad,cp=native(pv,pp,base,{})
        assert not bad and cp<0 and abs(pv@M@pv-2)<2e-5
        nd=dict(base);nd[70]=pv
        bonds=[2]*9;bonds[1]=4;bonds[5]=2
        dom=E.recipes((0,0,pp),tuple(bonds),((0,1),(1,2)))
        domains.append(dict(P=pj,symbolic_recipes=[list(z)for z in dom]))
        seen=set()
        for lp,low,hi in dom:
            for hh in product(HIGH,repeat=len(hi)):
                ms=list(low)
                for k,m in zip(hi,hh):ms[k]=m
                ids=[30,52,70];rhs=[-B.theta(m)-off(lp,pp if i==70 else 0)for i,m in zip(ids,ms)]
                vv,info=N.norm_roots([nd[i]for i in ids],rhs,profiles[lp]['square']);counts['L_equations']+=1
                if info['status'].startswith('unresolved'):unresolved.append(dict(stage='L_norm',P=pj,profile=lp,labels=ms,norm=info));continue
                for sg,l in vv:
                    counts['L_branches']+=1
                    vals,bad,cp=native(l,lp,nd,{70:pp})
                    if bad or cp>=-1e-8:blockers['L',tuple(bad)]+=1;continue
                    err=max(abs(float(l@M@l)-profiles[lp]['square']),max(abs(float(l@M@nd[i])-r)for i,r in zip(ids,rhs)))
                    if err>2e-5:unresolved.append(dict(stage='L_residual',P=pj,error=err));continue
                    key=(lp,)+tuple(np.round(l,9))
                    if key in seen:counts['L_duplicates']+=1;continue
                    seen.add(key);banks.append(dict(index=len(banks),P=pj,profile=lp,labels=ms,branch=sg,normal=l.tolist(),norm_residual=err,native_products=sorted(vals.items()),center_product=cp))
        print('RELEASE47_FIXED_L',pj,'banks',len(banks),flush=True)
    assert len(banks)<1000
    # All116 moving A branches are inherited analytic formulas; high A labels
    # in this fixed recipe domain are only the four saved representatives.
    seen=set()
    for fam in moving['moving_families']:
        fid=fam['index'];aprof=fam['profile'];a0,at,ai=map(np.array,[fam['constant'],fam['positive_t'],fam['inverse_t']])
        am7,am30,amc=fam['labels']
        for lr in banks:
            li=lr['index'];lp=lr['profile'];lm30,lm52,lmp=lr['labels'];l=np.array(lr['normal'])
            pj=lr['P'];pr=pbank[pj];pv=np.array(pr['normal']);pp=pr['profile'];pm7=pr['labels'][0]
            allowed=[m for m in range(2,7)if P.corner3(0,aprof,lp,am30,lm30,m)is not None]
            allowed=[m for m in allowed if m!=6]+(HIGH if 6 in allowed else[])
            ca,cb,cc=float(at@M@l),float(ai@M@l),float(a0@M@l)+off(aprof,lp)
            for alm in allowed:
                counts['A_L_quadratics']+=1
                # The all-real identity case needs separate treatment; never
                # silently regard an identically zero equation as exhausted.
                if max(abs(ca),abs(cb),abs(cc+B.theta(alm)))<1e-10:
                    unresolved.append(dict(stage='identical_AL_equation',family=fid,L=li,label=alm));continue
                for t in U.equation_roots(ca,cb,cc+B.theta(alm)):
                    if t<t0-1e-9:continue
                    counts['positive_parameter_roots']+=1
                    u=float((t-1/t)/(2*gamma));a=a0+t*at+ai/t;c=cv(t)
                    nd=dict(base);nd.update({67:c,68:a,69:l,70:pv});types={68:aprof,69:lp,70:pp}
                    err=max(abs(float(a@M@a)-profiles[aprof]['square']),abs(float(c@M@c)-2),abs(float(H.A.native(a,aprof,l,lp,profiles))+B.theta(alm)),abs(float(a@M@c)+B.theta(amc)))
                    if err>2e-5:unresolved.append(dict(stage='event_residual',family=fid,L=li,u=u,error=err));continue
                    gate_data={};failed=False
                    for name,v,pi,i in [('C',c,0,67),('A',a,aprof,68),('L',l,lp,69),('P',pv,pp,70)]:
                        vals,bad,cp=native(v,pi,{j:w for j,w in nd.items()if j!=i},types)
                        gate_data[name]=dict(native_products=sorted(vals.items()),center_product=cp)
                        if bad or cp>=-1e-8:blockers['event_'+name,tuple(bad)]+=1;failed=True
                    if failed:continue
                    key=(aprof,lp,pp)+tuple(np.round(np.r_[c,a,l,pv],8))
                    if key in seen:counts['event_duplicates']+=1;continue
                    seen.add(key)
                    sc,ok=score(nd,types)
                    if not ok:unresolved.append(dict(stage='event_cone',family=fid,L=li,u=u,**sc));continue
                    ev=dict(index=len(events),family=fid,L=li,P=pj,u=u,t=float(t),AL=alm,normals={str(j):nd[j].tolist()for j in[67,68,69,70]},profiles=types,norm_residual=err,gates=gate_data,score=compact_score(sc))
                    events.append(ev);counts['native_AL_events']+=1;bd=sc['boundary'];hist['event',bd['counts']['finite'],bd['counts']['exterior'],bd['open_ends']]+=1
                    if bd['open_ends']<10:improvements.append(dict(kind='event',index=ev['index'],complete=sc))
                    # Required finite7/A/P pair corners are tested by the full
                    # inherited A2 classification before any R frame solve.
                    rb=[2]*9;rb[1]=am7;rb[2]=pm7
                    rd=E.recipes((4,aprof,pp),tuple(rb),((0,1),(0,2)))
                    for rp,low,hi in rd:
                        for hh in product(HIGH,repeat=len(hi)):
                            ms=list(low)
                            for k,m in zip(hi,hh):ms[k]=m
                            ids=[7,68,70];rhs=[-B.theta(m)-off(rp,4 if i==7 else types[i])for i,m in zip(ids,ms)]
                            vv,info=N.norm_roots([nd[i]for i in ids],rhs,profiles[rp]['square']);counts['R_equations']+=1
                            if info['status'].startswith('unresolved'):unresolved.append(dict(stage='R_norm',event=ev['index'],profile=rp,labels=ms,norm=info));continue
                            for sg,r in vv:
                                counts['R_branches']+=1;vals,bad,cp=native(r,rp,nd,types)
                                if bad or cp>=-1e-8:blockers['R',tuple(bad)]+=1;continue
                                er=max(abs(float(r@M@r)-profiles[rp]['square']),max(abs(float(r@M@nd[i])-x)for i,x in zip(ids,rhs)))
                                if er>2e-5:unresolved.append(dict(stage='R_residual',event=ev['index'],error=er));continue
                                rel=B.relation(vals[69]);closed=isinstance(rel,int)and P.corner3(aprof,lp,rp,alm,ms[1],rel)is not None and P.corner3(pp,lp,rp,lmp,ms[2],rel)is not None
                                nn=dict(nd);nn[71]=r;ty=dict(types);ty[71]=rp;sc2,good=score(nn,ty)
                                if not good:unresolved.append(dict(stage='R_cone',event=ev['index'],**sc2));continue
                                rec=dict(index=len(third),event=ev['index'],profile=rp,labels=ms,branch=sg,normal=r.tolist(),native_products=sorted(vals.items()),center_product=cp,norm_residual=er,remaining_required_corners_finite=closed,score=compact_score(sc2))
                                third.append(rec);counts['native_R_caps']+=1
                                if closed:counts['completed_disks']+=1
                                bd=sc2['boundary'];hist['R',bd['counts']['finite'],bd['counts']['exterior'],bd['open_ends']]+=1
                                if bd['open_ends']<10:improvements.append(dict(kind='R',index=rec['index'],complete=sc2))
            assert len(events)<3000 and len(third)<30000,('bounded output guard',len(events),len(third))
        print('RELEASE47_QUADRATIC_FAMILY',fid,'events',len(events),'R',len(third),'improved',len(improvements),'unresolved',len(unresolved),flush=True)
    paths=[Path(__file__),Path(Q.__file__),src,Path(U.__file__),upath,Path(E.__file__),Path(N.__file__),Path(P.__file__),H.PARENT,H.BOUNDARY]
    return dict(status='NUMERICAL_RELEASE47_ALL_U_QUADRATIC_CONTACT_EVENTS',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
      deleted_walls=[47],retained_control_walls=62,exact_deletion_ledger=seed['exact_deletion_ledger'],topology=Q.disk_check(),P_seeds=pbank,L_domains=domains,L_seeds=banks,
      moving_A_family_count=len(moving['moving_families']),u_domain='u>=3; no numerical parameter grid or upper cutoff',high_label_samples=HIGH,
      counts=dict(counts),events=events,third_caps=third,improvements=improvements,unresolved=unresolved,
      histograms=[dict(stage=s,finite=f,exterior=e,open_ends=i,count=n)for(s,f,e,i),n in sorted(hist.items())],
      blockers=[dict(stage=k[0],walls=list(k[1]),count=n)for k,n in sorted(blockers.items())],
      scope='New fixed-L event domain for the rotated release47 disk. P restores12; L is solved in30/52/P with both finite corners; A(t).L is an explicit quadratic, solved for every positive root withu>=3 in116 inherited moving branches. HighA/LP/LA/R labels remain bounded6/8/12/24 representatives. Nearzero polynomial coefficients, tangencies, floating cancellation and inherited numerical norm classifiers are not certified continuous coverage. Every accepted record passes full62old/new native gates, strict commoncenter and complete numericalS6cone/PD support equality. Exact deletion debit2finite certificates andold12quad; H8compactfacet preserved conditionally by acute projection. No exact chamber, ambient finite census, or finite-volume claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'improvements',len(out['improvements']),'unresolved',len(out['unresolved']),flush=True)
