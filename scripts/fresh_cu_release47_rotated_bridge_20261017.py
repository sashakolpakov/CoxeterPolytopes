#!/usr/bin/env python3
"""Numerical changed-duty bridge after releasing old47.

P closes the old12-face, L owns the two30 boundary arcs, and R the7/A arc.
The parameter u of exact-admissible C(u) tunes finite P/L contacts. Neither
sampled root detection nor high-label representatives prove nonexistence.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import product,combinations
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import fresh_cu_three_cap_bridge_20261017 as T
U,E,H,B,M,N,S,P,ROOT=T.U,T.E,T.H,T.B,T.M,T.N,T.S,T.P,T.ROOT
TAG='fresh_cu_release47_rotated_bridge_20261017'
HIGH=[6,8,12,24]
PL_LABELS=[2,3,4,5,6,8,12,24]

def disk_check():
    triangles=[[7,12,70],[12,52,70],[52,70,69],[52,30,69],[30,68,69],
      [68,69,71],[68,7,71],[7,71,70],[70,69,71]]
    boundary=[7,12,52,30,68];edges=Counter(tuple(sorted(e))for f in triangles for e in combinations(f,2))
    assert set(edges.values())=={1,2}
    assert {e for e,m in edges.items()if m==1}=={tuple(sorted(e))for e in zip(boundary,boundary[1:]+boundary[:1])}
    vertices=set(sum(triangles,[]));assert len(vertices)-len(edges)+len(triangles)==1
    links=[]
    for v in sorted(vertices):
        ee=[tuple(x for x in f if x!=v)for f in triangles if v in f]
        dd=Counter(x for e in ee for x in e);seen={next(iter(dd))}
        while True:
            new=seen|{b for a,b in ee if a in seen}|{a for a,b in ee if b in seen}
            if new==seen:break
            seen=new
        assert seen==set(dd)
        assert sorted(dd.values())==([1,1]+[2]*(len(dd)-2)if v in boundary else[2]*len(dd))
        links.append(dict(vertex=v,edges=ee))
    return dict(boundary=boundary,triangles=triangles,vertices=len(vertices),edges=len(edges),faces=len(triangles),links=links)

def run():
    path=ROOT/'results'/f'{U.TAG}.json';parent=json.loads(path.read_text())
    cur,oldlabels,base,h,profiles,primary,ef,ei,e0,e1,gamma,t0,cv=U.setup()
    deletion=H.deletion(primary,cur,[47]);assert deletion['lost_finite_count']==2
    fam=parent['moving_families'][1];assert fam['profile']==0 and fam['labels']==[2,3,2]
    a0,ap,am=map(np.array,[fam['constant'],fam['positive_t'],fam['inverse_t']])
    base.pop(47);retained=[i for i in oldlabels if i!=47]
    actual=[i for i in cur['actual_sides']if i!=47]
    offset=lambda p,q:float(H.A.native(np.zeros(4),p,np.zeros(4),q,profiles))
    counts=Counter();blockers=Counter();unresolved=[];events=[];third=[];improved=[];closest=[]
    def setup_u(u):
        t=gamma*u+np.sqrt(1+gamma*gamma*u*u)
        nd=dict(base);nd[67]=cv(t);nd[68]=a0+ap*t+am/t
        return nd
    def native(v,pi,nd,types):
        vals={i:float(v@M@w)+(profiles[pi]['offsets'][oldlabels.index(i)]if i in oldlabels else offset(pi,types.get(i,0)))for i,w in nd.items()}
        bad=[i for i,x in vals.items()if B.relation(x)in[None,'parallel_or_unresolved']]
        cp=float(v@M@h)+profiles[pi]['center_offset']
        return vals,bad,cp
    def score(nd,types):
        sel=actual+[i for i in sorted(nd)if i not in oldlabels]
        bd=S.C.detailed_cone(nd,sel)
        sq=[profiles[types.get(i,0)]['square']if i not in oldlabels else float(nd[i]@M@nd[i])for i in sel]
        ic=H.incidence_score([nd[i]for i in sel],sel,sq)
        good=not bd['warnings']and ic['status']=='numerically_separated_acute_incidence'
        good=good and bd['open_ends']==ic['open_ends']and bd['counts'].get('finite',0)==ic['finite']
        good=good and{tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'}=={tuple(sorted(r))for r in ic['finite_supports']}
        return dict(boundary=bd,PD_incidence=ic),good
    deletion_only={i:v for i,v in base.items()if i!=67}
    delscore,ok=score(deletion_only,{});assert ok and delscore['boundary']['open_ends']==10
    # Complete mandatory finite-corner classification for the12-restoring P.
    bonds=[2]*9;bonds[1]=2;bonds[5]=6
    pdomain=E.recipes((4,0,0),tuple(bonds),((0,1),(1,2)))
    assert len(pdomain)==8 and all(not z[2]for z in pdomain)
    pbank=[]
    for pi,ms,hi in pdomain:
        ids=[7,12,52];rhs=[-B.theta(m)-offset(pi,4 if i==7 else 0)for i,m in zip(ids,ms)]
        vv,info=N.norm_roots([base[i]for i in ids],rhs,profiles[pi]['square'])
        if info['status'].startswith('unresolved'):unresolved.append(dict(stage='P_frame',profile=pi,labels=ms,norm=info))
        for sg,v in vv:
            counts['P_branches']+=1
            vals,bad,cp=native(v,pi,deletion_only,{})
            if bad or cp>=-1e-8:blockers['P',tuple(bad)]+=1;continue
            err=max(abs(float(v@M@v)-profiles[pi]['square']),max(abs(float(v@M@base[i])-r)for i,r in zip(ids,rhs)))
            if err>2e-5:unresolved.append(dict(stage='P_residual',error=err));continue
            pbank.append(dict(profile=pi,labels=list(ms),branch=sg,normal=v.tolist(),native_products=sorted(vals.items()),center_product=cp,norm_residual=err))
    assert len(pbank)==3 and sorted(z['labels']for z in pbank)==[[3,2,2],[4,2,2],[5,2,2]]
    # The old47 is included only as a control; the other two P normals differ.
    old47=U.setup()[2][47]
    assert np.max(np.abs(np.array(next(z['normal']for z in pbank if z['labels'][0]==3))-old47))<1e-7
    bonds=[2]*9;bonds[1]=4;bonds[2]=3
    ldomain=E.recipes((0,0,0),tuple(bonds),((0,1),(0,2)))
    lr=[]
    for pi,low,hi in ldomain:
        for hh in product(HIGH,repeat=len(hi)):
            ms=list(low)
            for k,m in zip(hi,hh):ms[k]=m
            for sg in[-1,1]:lr.append((pi,tuple(ms),sg))
    def left(j,u):
        pi,ms,sg=lr[j];nd=setup_u(u)
        vv=T.roots([nd[i]for i in[30,52,68]],[-B.theta(m)for m in ms],profiles[pi]['square'])
        return next((v for sign,v in vv if sign==sg),None)
    targets=[]
    for pj,pp in enumerate(pbank):
        for lj,(lp,lms,sg)in enumerate(lr):
            for m in PL_LABELS:
                if P.corner3(0,pp['profile'],lp,2,lms[1],m)is not None:targets.append((pj,lj,m))
    grid=np.unique(np.r_[np.linspace(3.,4.,65),np.geomspace(4.,32.,20)]).tolist()
    previous={};prevu=None;brackets=[];samples=[];seed_scores=[]
    for u in grid:
        nd=setup_u(u);_,ab,ac=native(nd[68],0,{i:v for i,v in nd.items()if i!=68},{})
        counts['A_grid_tests']+=1
        lbank={j:v for j in range(len(lr))if(v:=left(j,u))is not None}
        lnative={}
        for lj,v in lbank.items():
            lp,lms,_=lr[lj];vals,bad,cp=native(v,lp,nd,{})
            counts['L_grid_branches']+=1
            if bad or cp>=-1e-8:blockers['L_grid',tuple(bad)]+=1;continue
            if not ab and ac<-1e-8:lnative[lj]=v
        if lnative:
            # Save one full cone per distinct encountered L recipe, keeping the artifact compact.
            for lj,v in lnative.items():
                if any(r['L_recipe']==lj for r in seed_scores):continue
                nn=dict(nd);nn[69]=v;sc,good=score(nn,{69:lr[lj][0]})
                if not good:unresolved.append(dict(stage='L_cone',u=u,recipe=lj,**sc));continue
                seed_scores.append(dict(u=u,L_recipe=lj,**sc))
        current={};best=None
        for pj,lj,m in targets:
            if lj not in lbank:continue
            pp=pbank[pj];lp,_,_=lr[lj]
            g=float(H.A.native(np.array(pp['normal']),pp['profile'],lbank[lj],lp,profiles))+B.theta(m)
            key=pj,lj,m;current[key]=g
            if key in previous and g*previous[key]<0:brackets.append(dict(key=list(key),lower=prevu,upper=u,lower_gap=previous[key],upper_gap=g))
            if lj in lnative and (best is None or abs(g)<abs(best['gap'])):best=dict(u=u,key=list(key),gap=g)
        if best:closest.append(best)
        samples.append(dict(u=u,A_native=not ab and ac<-1e-8,L_real=len(lbank),L_native=len(lnative),targets=len(current)))
        previous=current;prevu=u
    assert len(brackets)<2000,('bounded root budget',len(brackets))
    print('RELEASE47_GRID',len(grid),'brackets',len(brackets),'L_recipes_native',len(seed_scores),flush=True)
    seen=set()
    for br in brackets:
        pj,lj,m=br['key'];pp=pbank[pj];pv=np.array(pp['normal']);lp,lms,lsg=lr[lj]
        def gap(u):
            l=left(lj,u)
            if l is None:raise ValueError('norm branch absent')
            return float(H.A.native(pv,pp['profile'],l,lp,profiles))+B.theta(m)
        try:u=float(brentq(gap,br['lower'],br['upper'],xtol=1e-12,maxiter=80))
        except(ValueError,RuntimeError)as ex:unresolved.append(dict(stage='PL_root',bracket=br,error=str(ex)));continue
        key=pj,lj,m,round(u,9)
        if key in seen:continue
        seen.add(key);counts['PL_event_roots']+=1
        nd=setup_u(u);l=left(lj,u);assert l is not None
        tests=[('A',nd[68],0,{i:v for i,v in nd.items()if i!=68},{}),('P',pv,pp['profile'],nd,{}),('L',l,lp,dict(nd,**{}),{})]
        rejected=[];products={}
        for name,v,pi,nn,ty in tests:
            vals,bad,cp=native(v,pi,nn,ty);products[name]=dict(native_products=sorted(vals.items()),center_product=cp)
            if bad or cp>=-1e-8:blockers['event_'+name,tuple(bad)]+=1;rejected.append((name,bad))
        if rejected:continue
        err=max(abs(float(l@M@l)-profiles[lp]['square']),abs(gap(u)))
        if err>2e-5:unresolved.append(dict(stage='event_residual',u=u,error=err));continue
        nd[69]=l;nd[70]=pv;types={69:lp,70:pp['profile']}
        sc,good=score(nd,types)
        if not good:unresolved.append(dict(stage='PL_cone',u=u,**sc));continue
        event=dict(index=len(events),u=u,P_index=pj,L_recipe=lj,PL=m,normals={str(i):nd[i].tolist()for i in[67,68,69,70]},products=products,norm_residual=err,**sc)
        events.append(event);counts['native_PL_events']+=1
        if sc['boundary']['open_ends']<10:improved.append(dict(kind='PL_event',index=event['index']))
        # R owns7/A and7/P; both other required corners A/L/R and P/L/R checked literally.
        rb=[2]*9;rb[1]=2;rb[2]=pp['labels'][0]
        rdom=E.recipes((4,0,pp['profile']),tuple(rb),((0,1),(0,2)))
        for rp,low,hi in rdom:
            for hh in product(HIGH,repeat=len(hi)):
                ms=list(low)
                for k,mm in zip(hi,hh):ms[k]=mm
                ids=[7,68,70];rhs=[-B.theta(mm)-offset(rp,4 if i==7 else types.get(i,0))for i,mm in zip(ids,ms)]
                vv,info=N.norm_roots([nd[i]for i in ids],rhs,profiles[rp]['square']);counts['R_equations']+=1
                if info['status'].startswith('unresolved'):unresolved.append(dict(stage='R_norm',event=event['index'],profile=rp,labels=ms,norm=info));continue
                for sg,r in vv:
                    counts['R_branches']+=1
                    vals,bad,cp=native(r,rp,nd,types)
                    if bad or cp>=-1e-8:blockers['R',tuple(bad)]+=1;continue
                    er=max(abs(float(r@M@r)-profiles[rp]['square']),max(abs(float(r@M@nd[i])-z)for i,z in zip(ids,rhs)))
                    if er>2e-5:unresolved.append(dict(stage='R_residual',event=event['index'],error=er));continue
                    rl=B.relation(vals[69]);closed=isinstance(rl,int)and P.corner3(0,lp,rp,lms[2],ms[1],rl)is not None and P.corner3(pp['profile'],lp,rp,m,ms[2],rl)is not None
                    nn=dict(nd);nn[71]=r;ty=dict(types);ty[71]=rp;sc2,ok=score(nn,ty)
                    if not ok:unresolved.append(dict(stage='R_cone',event=event['index'],**sc2));continue
                    rec=dict(event=event['index'],profile=rp,labels=ms,branch=sg,normal=r.tolist(),native_products=sorted(vals.items()),center_product=cp,norm_residual=er,remaining_finite_corners=closed,**sc2)
                    third.append(rec);counts['native_R_caps']+=1
                    if closed:counts['completed_three_cap_disks']+=1
                    if sc2['boundary']['open_ends']<10:improved.append(dict(kind='R_cap',index=len(third)-1))
        print('RELEASE47_EVENT',event['index'],u,'P',pp['labels'],'L',lms,'PL',m,'ends',sc['boundary']['open_ends'],'R',len(third),flush=True)
    paths=[Path(__file__),Path(T.__file__),Path(U.__file__),path,Path(E.__file__),Path(H.__file__),Path(N.__file__),Path(P.__file__),H.PARENT,H.BOUNDARY]
    return dict(status='BOUNDED_RELEASE47_ROTATED_THREE_CAP_BRIDGE',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
      deleted_walls=[47],retained_control_walls=62,exact_deletion_ledger=deletion,deletion_only= delscore,
      topology=disk_check(),A_family=fam,P_symbolic_recipes=[list(x)for x in pdomain],P_seeds=pbank,L_symbolic_recipes=[list(x)for x in ldomain],L_recipes=lr,
      high_label_samples=HIGH,PL_label_samples=PL_LABELS,u_grid=grid,grid_summary=samples,L_sample_cones=seed_scores,PL_brackets=brackets,
      closest_sampled_contacts=sorted(closest,key=lambda r:abs(r['gap']))[:20],counts=dict(counts),events=events,third_caps=third,improvements=improved,unresolved=unresolved,
      blockers=[dict(stage=k[0],walls=list(k[1]),count=n)for k,n in sorted(blockers.items())],
      scope='Numerical changed-duty construction after exact debit of old47. The12-restoring P mandatory finite-corner classification is inherited complete; actual normal/native/cone reconstruction is floating point. One moving A family, u3..32 at84 samples, sign-changing P/L event brackets, bounded high labels6/8/12/24. Tangencies, near-singular norm branches, disconnected admissibility, other A families, other high labels and u>32 remain open. Every accepted event/cap passes all62 retained native gates, allnew contacts, a common strict center, and complete S6 cone/independentPD support agreement. No exact new geometry, ambient finite census or finite-volume claim.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'improvements',len(out['improvements']),'unresolved',len(out['unresolved']),flush=True)
