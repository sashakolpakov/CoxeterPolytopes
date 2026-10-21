#!/usr/bin/env python3
"""Bounded continuous three-cap disk on a moving, deletion-free A companion.

Boundary7-47-52-30-A-7; R owns the two7 arcs, L the two30 arcs,
D owns47-52. A finite R-L event fixes the one free boost parameter.
"""
import argparse,json
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import fresh_cu_companion_contact_events_20261017 as U
import fresh_c3_companion_continuation_20261017 as E
C,H,B,M,N,S,P,ROOT=U.C,U.H,U.B,U.M,U.N,U.S,U.P,U.ROOT
TAG='fresh_cu_three_cap_bridge_20261017'
HIGH=[6,8,12]
RL_LABELS=[2,3,4,5,6,8,12]
MAX_BRACKETS=160

def roots(frame,rhs,q):
    f=np.array(frame);g=f@M@f.T
    if abs(np.linalg.det(g))<1e-10:return []
    part=np.array(rhs)@np.linalg.inv(g)@f;a=f@M
    k=np.array([(-1.)**j*np.linalg.det(np.delete(a,j,axis=1))for j in range(4)])
    k/=np.linalg.norm(k);kq=float(k@M@k)
    if abs(kq)<1e-10:return []
    ratio=(q-float(part@M@part))/kq
    if ratio<=1e-10:return []
    return [(sg,part+sg*np.sqrt(ratio)*k)for sg in[-1,1]]

def run():
    src=ROOT/'results'/f'{U.TAG}.json';parent=json.loads(src.read_text())
    cur,labels,base,h,profiles,primary,ef,ei,e0,e1,gamma,t0,cv=U.setup()
    f=parent['moving_families'][1];assert f['profile']==0 and f['labels']==[2,3,2]
    a0,ap,am=map(np.array,[f['constant'],f['positive_t'],f['inverse_t']])
    a47,b47,c47=float(ap@M@base[47]),float(am@M@base[47]),float(a0@M@base[47])
    thresholds=[t for t in U.equation_roots(a47,b47,c47+2)if t>t0]
    assert len(thresholds)==1
    tmin=thresholds[0];umin=(tmin-1/tmin)/(2*gamma)
    grid=umin+np.geomspace(.0001,32.-umin,48)
    frames={'R':[7,47,68],'L':[30,52,68]};types={7:4}
    off=lambda a,b:float(H.A.native(np.zeros(4),a,np.zeros(4),b,profiles))
    ps=lambda i:types.get(i,0)
    recipes={};domains=[]
    for name,ids in frames.items():
        bonds=[2]*9;bonds[1]=3 if name=='R' else 4;bonds[2]=2 if name=='R' else 3
        rows=E.recipes(tuple(ps(i)for i in ids),tuple(bonds),((0,1),(0,2)))
        rr=[]
        for pi,low,high in rows:
            for hh in product(HIGH,repeat=len(high)):
                ms=list(low)
                for i,m in zip(high,hh):ms[i]=m
                for sg in[-1,1]:rr.append((pi,tuple(ms),sg))
        recipes[name]=rr;domains.append(dict(name=name,frame=ids,symbolic_recipes=len(rows),sampled_branches=len(rr)))
    counts=Counter();blockers=Counter();unresolved=[];events=[];accepted=[];improved=[];closest=[];brackets=[]
    def setup_u(u):
        t=gamma*u+np.sqrt(1+gamma*gamma*u*u)
        nd=dict(base);nd[67]=cv(t);nd[68]=a0+ap*t+am/t
        return nd
    def native(v,pi,nd,extra=()):
        vals={i:float(v@M@nd[i])+(profiles[pi]['offsets'][labels.index(i)]if i in labels else 0.)for i in nd}
        for i,w,pj in extra:vals[i]=float(H.A.native(v,pi,w,pj,profiles))
        good=not any(B.relation(x)in[None,'parallel_or_unresolved']for x in vals.values())
        return vals,good and float(v@M@h)+profiles[pi]['center_offset']<-1e-8
    def single(name,j,u,check=True):
        if check:counts[name+'_norm_branches_tested']+=1
        nd=setup_u(u);pi,ms,sg=recipes[name][j];ids=frames[name]
        rhs=[-B.theta(m)-off(pi,ps(i))for i,m in zip(ids,ms)]
        found=next((v for s,v in roots([nd[i]for i in ids],rhs,profiles[pi]['square'])if s==sg),None)
        if found is None:
            if check:counts[name+'_no_real_branch']+=1
            return None
        err=max(abs(float(found@M@found)-profiles[pi]['square']),max(abs(float(found@M@nd[i])-r)for i,r in zip(ids,rhs)))
        if err>2e-5:
            if check:counts[name+'_residual_rejected']+=1
            return None
        if check:
            vals,good=native(found,pi,nd)
            if not good:
                bad=tuple(i for i,x in vals.items()if B.relation(x)in[None,'parallel_or_unresolved'])
                blockers[name,bad]+=1
                if bad==(67,):counts[name+'_blocked_only_by_C']+=1
                return None
        return found
    def score(nd,more):
        nn=dict(nd);nn.update({69+i:v for i,v in enumerate(more)});sel=cur['actual_sides']+[67,68]+list(range(69,69+len(more)))
        bd=S.C.detailed_cone(nn,sel);ic=H.incidence_score([nn[i]for i in sel],sel,[float(nn[i]@M@nn[i])for i in sel])
        ok=not bd['warnings']and ic['status']=='numerically_separated_acute_incidence'and bd['open_ends']==ic['open_ends']and bd['counts'].get('finite',0)==ic['finite']
        ok=ok and {tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'}=={tuple(sorted(r))for r in ic['finite_supports']}
        return dict(boundary=bd,PD_incidence=ic),ok
    previous={};prevu=None;grid_summary=[]
    for u in grid:
        nd=setup_u(float(u));avals,agood=native(nd[68],0,{i:v for i,v in nd.items()if i!=68})
        assert agood and float(nd[68]@M@nd[47])<-2
        banks={name:{j:v for j in range(len(rr))if(v:=single(name,j,float(u)))is not None}for name,rr in recipes.items()}
        current={};best=-1e99
        for ri,r in banks['R'].items():
            rp,rms,_=recipes['R'][ri]
            for li,l in banks['L'].items():
                lp,lms,_=recipes['L'][li];g=float(H.A.native(r,rp,l,lp,profiles));best=max(best,g)
                for m in RL_LABELS:
                    if P.corner3(0,rp,lp,rms[2],lms[2],m)is None:continue
                    key=ri,li,m;current[key]=g+B.theta(m)
                    if key in previous and current[key]*previous[key]<0:
                        brackets.append(dict(key=list(key),lower=prevu,upper=float(u),lower_gap=previous[key],upper_gap=current[key]))
        grid_summary.append(dict(u=float(u),R=len(banks['R']),L=len(banks['L']),finite_corner_pair_targets=len(current),maximum_mutual_product=best))
        if current:closest.append(dict(u=float(u),key=list(min(current,key=lambda k:abs(current[k]))),gap=min(current.values(),key=abs)))
        counts['grid_pair_targets']+=len(current);previous=current;prevu=float(u)
        print('CU_THREE_CAP_GRID',float(u),'R/L',len(banks['R']),len(banks['L']),'maxRL',best,'brackets',len(brackets),flush=True)
    assert len(brackets)<=MAX_BRACKETS,('bounded root budget',len(brackets))
    seen=set()
    for br in brackets:
        ri,li,m=br['key'];rp,rms,_=recipes['R'][ri];lp,lms,_=recipes['L'][li]
        def gap(u):
            r=single('R',ri,u,False);l=single('L',li,u,False)
            if r is None or l is None:raise ValueError('branch does not continue')
            return float(H.A.native(r,rp,l,lp,profiles))+B.theta(m)
        try:u=float(brentq(gap,br['lower'],br['upper'],xtol=1e-12,maxiter=60))
        except(ValueError,RuntimeError)as ex:unresolved.append(dict(stage='event_root',bracket=br,error=str(ex)));continue
        key=ri,li,m,round(u,9)
        if key in seen:continue
        seen.add(key);nd=setup_u(u);r=single('R',ri,u);l=single('L',li,u)
        if r is None or l is None:counts['event_native_rejected']+=1;continue
        sc,good=score(nd,[r,l])
        if not good:unresolved.append(dict(stage='RL_cone',u=u,key=br['key'],**sc));continue
        event=dict(index=len(events),u=u,R_recipe=recipes['R'][ri],L_recipe=recipes['L'][li],RL=m,normals=[r.tolist(),l.tolist()],**sc)
        events.append(event);counts['native_RL_events']+=1
        if sc['boundary']['open_ends']<10:improved.append(dict(kind='RL_event',event=event))
        # Complete symbolic finite-corner recipes for D at47/52 and47/R;
        # the other two required corners52/L and the centralR/L are checked.
        nd[69]=r;nd[70]=l
        bds=[2]*9;bds[1]=2;bds[2]=rms[1]
        dom=E.recipes((0,0,rp),tuple(bds),((0,1),(0,2)))
        for dp,low,high in dom:
            for hh in product(HIGH,repeat=len(high)):
                ms=list(low)
                for j,mj in zip(high,hh):ms[j]=mj
                rhs=[-B.theta(ms[0]),-B.theta(ms[1]),-B.theta(ms[2])-off(dp,rp)]
                counts['D_equations']+=1
                for sg,d in roots([nd[47],nd[52],r],rhs,profiles[dp]['square']):
                    counts['D_branches']+=1
                    err=max(abs(float(d@M@d)-profiles[dp]['square']),
                      max(abs(float(d@M@w)-z)for w,z in zip([nd[47],nd[52],r],rhs)))
                    if err>2e-5:
                        unresolved.append(dict(stage='D_residual',u=u,profile=dp,labels=ms,error=err));continue
                    # nd's added R/L need their A2 offsets explicitly.
                    native_base={i:v for i,v in nd.items()if i not in[69,70]}
                    vals,ng=native(d,dp,native_base,[(69,r,rp),(70,l,lp)])
                    if not ng:continue
                    dl=B.relation(vals[70])
                    finite=isinstance(dl,int)and P.corner3(0,lp,dp,lms[1],ms[1],dl)is not None and P.corner3(rp,lp,dp,m,ms[2],dl)is not None
                    sc2,ok=score(native_base,[r,l,d])
                    if not ok:unresolved.append(dict(stage='D_cone',u=u,**sc2));continue
                    rec=dict(event=event['index'],profile=dp,labels=ms,branch=sg,normal=d.tolist(),norm_residual=err,native_products=sorted(vals.items()),remaining_required_finite_corners=finite,**sc2)
                    accepted.append(rec);counts['native_D_caps']+=1
                    if finite:counts['completed_three_cap_disks']+=1
                    if sc2['boundary']['open_ends']<10:improved.append(dict(kind='D_cap',cap=rec))
        print('CU_THREE_CAP_EVENT',len(events),u,'D',len(accepted),'improved',len(improved),flush=True)
    paths=[Path(__file__),Path(U.__file__),src,Path(E.__file__),Path(P.__file__),H.PARENT,H.BOUNDARY]
    return dict(status='BOUNDED_CONTINUOUS_CU_THREE_CAP_BRIDGE',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
      retained_control_walls=63,deleted_walls=[],A_family=f,parallel_threshold_u=float(umin),u_grid=grid.tolist(),u_search_upper=32,
      boundary_cycle=[7,47,52,30,68],triangles=[[7,47,69],[7,69,68],[68,69,70],[68,70,30],[30,70,52],[52,70,71],[52,71,47],[47,71,69],[69,71,70]],
      domains=domains,recipes={k:[list(z)for z in v]for k,v in recipes.items()},high_label_samples=HIGH,RL_labels=RL_LABELS,
      counts=dict(counts),blockers=[dict(stage=k[0],walls=list(k[1]),count=n)for k,n in sorted(blockers.items())],grid_summary=grid_summary,closest_sampled_contacts=sorted(closest,key=lambda r:abs(r['gap']))[:20],
      brackets=brackets,events=events,accepted_third_caps=accepted,improvements=improved,unresolved=unresolved,
      scope='One genuinely moving A family, beyond A47 parallel threshold, with the C triangle retained. Three caps have disjoint boundary duties in a five-edge disk. Only48 u samples through32 detect sign-changing R-L finite-contact brackets, with at most160 bounded scalar solves. Tangencies, disconnected admissible intervals, u>32, other A families and other high labels remain open. All accepted events/caps use full native/common-center and complete numerical target cone/PD support agreement. No exact new geometry or finite-volume claim; no old deletions.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'events',len(out['events']),'improved',len(out['improvements']),'unresolved',len(out['unresolved']),flush=True)
