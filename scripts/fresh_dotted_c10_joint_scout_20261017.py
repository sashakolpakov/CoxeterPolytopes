#!/usr/bin/env python3
"""Bounded discovery on a genuinely continuous dotted-C10 prefix domain.

Samples locate promising one-variable families; bounded scalar minimization
then varies C10 continuously. Neither step is an interval exclusion proof.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar,brentq
import fresh_joint_rz_circle_scout_20261017 as J

D,O,P,F,T,N,S,H,B,M,ROOT=J.D,J.O,J.P,J.F,J.T,J.N,J.S,J.H,J.B,J.M,J.ROOT
TAG='fresh_dotted_c10_joint_scout_20261017'
UGRID=[2.0001,2.01,2.04,2.1,2.25,2.5,3.,4.,6.,8.,12.,16.]
ACS=[2,3,4,5]

def run():
    cur,labels,ns,h,pro,primary=H.setup();actual=[i for i in cur['actual_sides']if i not in S.DELETE]
    counts=Counter();blockers=Counter();unresolved=[];prefixes=[];events=[];best={}
    f=np.array([ns[i]for i in[7,30,10]]);inv=np.linalg.inv(f@M@f.T)
    v=(inv@np.array([0.,0.,-1.]))@f;alpha=float(-v@M@v)
    roots,info=N.norm_roots(f,[0.,0.,-2.],2.)
    k=next((c-2*v)/np.sqrt(1+2*alpha)for _,c in roots if (c-2*v)@M@ns[55]<0)
    assert abs(k@M@k-2)<1e-7 and abs(v@M@k)<1e-7
    circle=json.loads((ROOT/'results/fresh_r_circle_admissibility_20261017.json').read_text())
    rs={r['label']:np.array(r['normal'])for r in circle['numerical_R_candidates']if r['label']<=15 and r['native_and_center_pass']}
    zcases=[]
    for zp,m8,mr,mx in product(range(7),range(2,6),range(2,6),range(2,6)):
        if P.corner3(0,0,zp,2,m8,mx)is None:continue
        if P.corner3(0,0,zp,2,m8,2)is None or P.corner3(0,0,zp,2,mr,2)is None:continue
        rem=pro[zp]['square']-(B.theta(m8)**2+B.theta(mr)**2)/2
        if rem>=-1e-9:continue
        zcases.append((zp,m8,mr,mx))
    assert len(zcases)==28 and len(rs)==14
    def C(u):return u*v+np.sqrt(1+alpha*u*u/2)*k
    def native(v,pi,assigned=(),whole=False):
        values={i:float(v@M@ns[i]+pro[pi]['offsets'][labels.index(i)])for i in labels if whole or i not in S.DELETE}
        return T.gates(values,assigned)
    def npair(a,pi,b,pj):return float(H.A.native(a,pi,b,pj,pro))
    def dotted(a,pi,b,pj):return float(a@M@b)<-np.sqrt(pro[pi]['square']*pro[pj]['square'])-1e-7
    def nonincident(v,pi,allow):return all(float(v@M@ns[i])<-np.sqrt(pro[pi]['square']*float(ns[i]@M@ns[i]))-1e-7 for i in actual if i not in allow)
    def solve(stage,frame,rhs,pi,assigned=()):
        counts[stage+'_equations']+=1;roots,info=N.norm_roots(frame,rhs,pro[pi]['square'])
        if info['status'].startswith('unresolved'):
            unresolved.append(dict(stage=stage,info=info));return []
        out=[]
        for sg,r in roots:
            counts[stage+'_branches']+=1
            err=max(abs(float(r@M@r)-pro[pi]['square']),max(abs(float(r@M@w)-q)for w,q in zip(frame,rhs)))
            if err>1e-5:unresolved.append(dict(stage=stage,error=err));continue
            bad=native(r,pi,assigned)
            if bad:blockers[stage,str(bad[0][0])]+=1;continue
            if float(r@M@h)+pro[pi]['center_offset']>=-1e-8:blockers[stage,'center']+=1;continue
            out.append((sg,r));counts[stage+'_native']+=1
        return out
    def acx(u,mac):
        c=C(u);out=[]
        if native(c,0,[7,30])or not nonincident(c,0,[7,30]):return []
        for sa,a in solve('A',[ns[7],ns[8],c],[0.,0.,-B.theta(mac)],0,[7,8]):
            if not nonincident(a,0,[7,8]):continue
            for sx,x in solve('X',[ns[8],a,c],[0.,-1.,0.],0,[8]):
                if nonincident(x,0,[8]):out.append((sa,sx,a,c,x))
        return out
    def score(vs,ps,whole=False):
        nd=dict(ns);nd.update({67+i:w for i,w in enumerate(vs)})
        sel=(cur['actual_sides']if whole else actual)+list(range(67,67+len(vs)))
        bd=S.C.detailed_cone(nd,sel);ic=H.incidence_score([nd[i]for i in sel],sel,[float(nd[i]@M@nd[i])for i in sel])
        if bd['warnings']or ic['status']!='numerically_separated_acute_incidence'or bd['open_ends']!=ic['open_ends']or bd['counts'].get('finite',0)!=ic['finite']:
            unresolved.append(dict(stage='boundary',bd=bd,ic=ic));return None
        if {tuple(sorted(q['active']))for q in bd['rays']if q['kind']=='finite'}!={tuple(sorted(q))for q in ic['finite_supports']}:
            unresolved.append(dict(stage='support_equality'));return None
        return dict(boundary=bd,PD_incidence=ic)
    # Whole-control one-wall events are a diagnostic, including deletion-free
    # finite C47 contacts. They are not mistaken for patch completion.
    for m in[2,3,4,5,6,8,12,24]:
        u=brentq(lambda u:float(C(u)@M@ns[47])+B.theta(m),2.,3.)
        c=C(u);bad=native(c,0,[7,30,47],whole=True)
        assert not bad
        sc=score([c],[0],whole=True);assert sc is not None
        events.append(dict(C47=m,u=u,normal=c.tolist(),native_products=[[i,float(c@M@ns[i])]for i in labels],**sc))
    tail=score([C(3.)],[0],whole=True);assert tail is not None
    # Isolated finite old53 contacts are not detected by a u grid. Solve
    # those events directly. AC>=6 changes the finite ACX duty: X must
    # be right to both A and C. This is a new incidence specialization.
    high_ac_events=[]
    for m53 in[3,4,5,6,8,12]:
        for sa,a in solve('event_A',[ns[7],ns[8],ns[53]],[0.,0.,-B.theta(m53)],0,[7,8,53]):
            for mac in[6,8,12,24]:
                for sc,c in solve('event_C',[ns[7],ns[30],a],[0.,0.,-B.theta(mac)],0,[7,30]):
                    u=float(-c@M@ns[10])
                    if u<=2+1e-7:continue
                    xs=[];trials=0
                    for xp,m8 in product(range(7),range(2,6)):
                        if P.corner3(0,0,xp,mac,2,2)is None or P.corner3(0,0,xp,2,m8,2)is None:continue
                        trials+=1
                        for sx,x in solve('event_X',[ns[8],a,c],[-B.theta(m8),0.,0.],xp,[8]):
                            xs.append(dict(profile=xp,X8=m8,branch=sx,normal=x.tolist()))
                    high_ac_events.append(dict(A53=m53,AC=mac,u=u,branches=[sa,sc],normals=[a.tolist(),c.tolist()],X_trials=trials,native_X=xs))
    def zroots(u,mac,mr,case,selected=None):
        zp,m8,mzr,mx=case;out=[]
        for sa,sx,a,c,x in acx(u,mac):
            if selected is not None and (sa,sx)!=selected:continue
            r=rs[mr]
            if not all(dotted(r,0,w,0)for w in[a,c,x]):continue
            roots,info=N.norm_roots([ns[8],r,x],[-B.theta(m8),-B.theta(mzr),-B.theta(mx)],pro[zp]['square'])
            if info['status'].startswith('unresolved'):
                unresolved.append(dict(stage='Z',u=u,AC=mac,R44=mr,case=case,info=info));continue
            for sz,z in roots:
                counts['Z_branches']+=1
                err=max(abs(float(z@M@z)-pro[zp]['square']),abs(float(z@M@ns[8])+B.theta(m8)),abs(float(z@M@r)+B.theta(mzr)),abs(float(z@M@x)+B.theta(mx)))
                if err>1e-5:unresolved.append(dict(stage='Z',error=err));continue
                values=[(i,float(z@M@ns[i]+pro[zp]['offsets'][labels.index(i)]))for i in labels if i not in S.DELETE+[8]]
                worst=max(values,key=lambda p:p[1]);margin=worst[1]+2
                key=(mac,mr,*case,sa,sx,sz)
                if key not in best or margin<best[key]['margin']:
                    best[key]=dict(key=list(key),u=u,margin=margin,worst_wall=worst[0],worst_product=worst[1])
                if selected is not None:out.append((sz,margin))
                if native(z,zp,[8])or float(z@M@h)+pro[zp]['center_offset']>=-1e-8:continue
                if not nonincident(z,zp,[8])or not all(dotted(z,zp,w,0)for w in[a,c]):continue
                for sy,y in solve('Y',[ns[8],r,z],[0.,0.,0.],0,[8]):
                    if not nonincident(y,0,[8])or not all(dotted(y,0,w,0)for w in[a,c,x]):continue
                    vs=[a,r,c,x,y,z];ps=[0,0,0,0,0,zp]
                    pairs=[dict(pair=[i,j],product=npair(vs[i],ps[i],vs[j],ps[j]))for i,j in combinations(range(6),2)]
                    if any(B.relation(p['product'])in[None,'parallel_or_unresolved']for p in pairs):continue
                    sc=score(vs,ps)
                    if sc is None:continue
                    yr=[q for q in sc['boundary']['rays']if 71 in q['active']]
                    if len(yr)!=3 or any(q['kind']!='finite'for q in yr):unresolved.append(dict(stage='Y_triangle'));continue
                    prefixes.append(dict(u=u,AC=mac,R44=mr,Z_profile=zp,Z_labels=[m8,mzr,mx],branches=[sa,sx,sz,sy],normals=[w.tolist()for w in vs],profiles=ps,native_pairs=pairs,**sc))
        return out
    for u,mac in product(UGRID,ACS):
        for mr,case in product(rs,zcases):zroots(u,mac,mr,case)
        print('DOTTED_C10_JOINT',u,mac,'prefixes',len(prefixes),'best',min([q['margin']for q in best.values()],default=None),flush=True)
    # Continuous local searches in the 24 best distinct literal branches.
    # The finite evaluation budget is part of the discovery scope.
    starts=sorted(best.values(),key=lambda q:q['margin'])[:24];optimizations=[]
    for st in starts:
        mac,mr,zp,m8,mzr,mx,sa,sx,sz=st['key'];case=(zp,m8,mzr,mx)
        def objective(t):
            rows=zroots(2+np.exp(t),mac,mr,case,(sa,sx))
            matches=[margin for sg,margin in rows if sg==sz]
            return min(matches)if matches else 1e6
        opt=minimize_scalar(objective,bounds=(np.log(.0001),np.log(14.)),method='bounded',options=dict(maxiter=60,xatol=1e-8))
        optimizations.append(dict(key=st['key'],initial_u=st['u'],initial_margin=st['margin'],u=float(2+np.exp(opt.x)),margin=float(opt.fun),success=bool(opt.success),evaluations=int(opt.nfev)))
    # Repeated optimization evaluations can rediscover exactly the same
    # prefix; do not count those as distinct chambers.
    unique={}
    for p in prefixes:
        key=(p['AC'],p['R44'],p['Z_profile'],tuple(p['Z_labels']),tuple(p['branches']),round(p['u'],9))
        unique[key]=p
    paths=[Path(__file__),Path(J.__file__),Path(D.__file__),H.PARENT,H.BOUNDARY,
           ROOT/'results/fresh_r_circle_admissibility_20261017.json',
           ROOT/'results/fresh_dotted_c10_domain_20261017.json',
           ROOT/'scripts/fresh_dotted_c10_domain_20261017.py']
    return dict(status='BOUNDED_CONTINUOUS_DOTTED_C10_JOINT_DISCOVERY',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        domain=dict(u_samples=UGRID,u_continuous_optimization=[2.0001,16.],AC=ACS,X_profile=0,X_labels=[2,3,2],R44=list(rs),Z_cases=[list(c)for c in zcases],optimized_literal_branches=len(optimizations),max_scalar_iterations=60),
        counts=dict(counts),blockers=[[list(k),v]for k,v in sorted(blockers.items())],best_branches=sorted(best.values(),key=lambda q:q['margin'])[:50],optimizations=optimizations,
        coherent_six_prefixes=list(unique.values()),whole_control_C47_events=events,whole_control_C3=tail,high_AC_finite_A53_events=high_ac_events,
        high_AC_event_scope='OnlyAC6,8,12,24 andA53=3,4,5,6,8,12 are tested. For a finiteACX corner withAC>=6, X must be rightA/C. Within the later polarY topology, X8>=6 is impossible: finite8XZ forcesZ right8/X; hence the8/R/Z frame is the orthogonal sum of8 and its positive-definiteR/Z pair, leaving a timelike polar line. This justifiesX8<=5 here, not an all-AC or all-incidence exclusion.',
        exact_old_deletion_cost=H.deletion(primary,cur,S.DELETE),unresolved=unresolved,
        scope='Numerical discovery on a continuously varying dotted-C10 family certified separately. Finite C47 events retain all63 control walls; the six-wall prefix domain uses the explicit five-wall deletion. Complete numerical intrinsic cones are cross-checked against independent PD support enumeration. No numerical result is promoted to exact geometry or finite volume. Bounded scalar searches and low-Z labels do not exclude otheru, high-Z tails, otherA/X recipes, or changed incidences.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],'prefixes',len(out['coherent_six_prefixes']),'unresolved',len(out['unresolved']),flush=True)
