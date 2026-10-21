#!/usr/bin/env python3
"""Numerical circle/isolated-bond necessary domains at frozen83 K events.

Low K-A/K-P: K-L varies continuously. High K-A or K-P forces K right
to L and R, hence a separate three-frame solve. No interval certificate.
"""
import argparse,json,math
from collections import Counter
from hashlib import sha256
from pathlib import Path
import numpy as np
import fresh_cu_release47_interior_cap_20261017 as F
D,U,E,H,B,M,N,S,P,ROOT=F.D,F.U,F.E,F.H,F.B,F.M,F.N,F.S,F.P,F.ROOT
TAG='fresh_cu_release47_ap_circle_20261017'
TAU=2*math.pi

def intersect(a,b):
    return [(max(x,u),min(y,v))for x,y in a for u,v in b if max(x,u)<=min(y,v)]

def cut(arcs,a,b,c):
    """Closed necessary half-plane a cos(t)+b sin(t)+c<=0, enlarged."""
    rho=math.hypot(a,b);eps=1e-9*max(1.,rho,abs(c));c-=eps
    if rho<1e-13:return arcs if c<=0 else[]
    v=-c/rho
    if v>=1:return arcs
    if v<-1:return[]
    phi=math.atan2(b,a);alpha=math.acos(max(-1.,min(1.,v)))
    lo=(phi+alpha)%TAU;hi=lo+TAU-2*alpha
    good=[(lo,min(hi,TAU))]
    if hi>TAU:good.append((0.,hi-TAU))
    return intersect(arcs,good)

def value(coeff,t):return coeff[0]+coeff[1]*math.cos(t)+coeff[2]*math.sin(t)

def ranges(coeff,arcs):
    phi=math.atan2(coeff[2],coeff[1])%TAU;out=[]
    for lo,hi in arcs:
        ts=[lo,hi]+[t for t in[phi,(phi+math.pi)%TAU]if lo<=t<=hi]
        vv=[value(coeff,t)for t in ts];out.append((min(vv),max(vv)))
    return out

def level_roots(coeff,target,arcs):
    c,a,b=coeff;rho=math.hypot(a,b)
    if rho<1e-12:return None if abs(c-target)<1e-8 else[]
    z=(target-c)/rho
    if z<-1-1e-9 or z>1+1e-9:return[]
    alpha=math.acos(max(-1.,min(1.,z)));phi=math.atan2(b,a)
    ts=[(phi+alpha)%TAU,(phi-alpha)%TAU]
    return sorted(set(t for t in ts if any(lo-1e-8<=t<=hi+1e-8 for lo,hi in arcs)))

def pd_interval(q0,q1,qk,g01,g0k):
    """Closed Schur interval for g1k in the three-root PD completion."""
    d=q0*q1-g01*g01;e=q0*qk-g0k*g0k
    if d<=0 or e< -1e-8:return None
    center=g01*g0k/q0;radius=math.sqrt(max(0.,d*e))/q0
    return center-radius,center+radius

def run():
    fp=ROOT/'results'/f'{F.TAG}.json';pj=ROOT/'results'/f'{D.TAG}.json'
    frozen=json.loads(fp.read_text());parent=json.loads(pj.read_text())
    events={e['index']:e for e in parent['events']};rbank={r['index']:r for r in parent['third_caps']}
    cur,labels,base,h,profiles,*_=U.setup();base.pop(47);base.pop(67)
    oldindex={i:j for j,i in enumerate(labels)};actual=[i for i in cur['actual_sides']if i!=47]
    off=lambda a,b:float(H.A.native(np.zeros(4),a,np.zeros(4),b,profiles))
    counts=Counter();blockers=Counter();records=[];candidates=[];unresolved=[];recover=0
    for row in frozen['records']:
        sel=row['selection'];eid=sel['selected_event'];ev=events[eid]
        types={int(i):v for i,v in ev['profiles'].items()};nd=dict(base)
        nd.update({int(i):np.array(v)for i,v in ev['normals'].items()})
        av,pv,lv=[nd[i]for i in[68,70,69]];ap,pp,lp=[types[i]for i in[68,70,69]]
        qa,qp,ql=[profiles[i]['square']for i in[ap,pp,lp]]
        af=np.array([av,pv]);gg=af@M@af.T;gi=np.linalg.inv(gg)
        assert np.linalg.det(gg)<-1e-7
        # Euclidean nullspace gives a basis of the Lorentz-orthogonal
        # complement, whose induced form must be positive definite.
        _,_,vt=np.linalg.svd(af@M);ker=vt[2:];eig,vec=np.linalg.eigh(ker@M@ker.T)
        assert eig.min()>1e-8
        eb=np.diag(1/np.sqrt(eig))@vec.T@ker
        assert np.max(np.abs(eb@M@eb.T-np.eye(2)))<1e-7 and np.max(np.abs(af@M@eb.T))<1e-7
        lows=set();high=[]
        for kp,ms,hi in row['symbolic_domain']:
            if 0 in hi or 2 in hi:high.append([kp,ms,hi])
            else:lows.add((kp,ms[0],ms[2]))
        eventrows=[]
        for kp,ma,mp in sorted(lows):
            qk=profiles[kp]['square'];rhs=np.array([-B.theta(ma)-off(kp,ap),-B.theta(mp)-off(kp,pp)])
            k0=rhs@gi@af;rr=float(qk-k0@M@k0);assert rr>0
            ee=math.sqrt(rr)*eb
            make=lambda t:k0+math.cos(t)*ee[0]+math.sin(t)*ee[1]
            def coeff(v,offset=0.):return[float(k0@M@v)+offset,float(ee[0]@M@v),float(ee[1]@M@v)]
            constraints=[]
            def upper(name,v,bound,offset=0.):
                z=coeff(v,offset);constraints.append((name,z[1],z[2],z[0]-bound))
            def lower(name,v,bound,offset=0.):
                z=coeff(v,offset);constraints.append((name,-z[1],-z[2],bound-z[0]))
            for i,w in nd.items():
                if i in[68,70]:continue
                o=profiles[kp]['offsets'][oldindex[i]]if i in oldindex else off(kp,types.get(i,0))
                upper('native:'+str(i),w,0.,o)
            upper('center',h,0.,profiles[kp]['center_offset'])
            # Intended K has no old-wall or C side. A PD pair of actual
            # acute facets would be an actual edge, so this is necessary.
            for i in actual+[67]:
                q=float(nd[i]@M@nd[i]);assert q>0
                upper('nonincident:'+str(i),nd[i],-math.sqrt(q*qk))
            il=[pd_interval(qa,ql,qk,float(av@M@lv),rhs[0]),pd_interval(qp,ql,qk,float(pv@M@lv),rhs[1])]
            counts['low_AP_recipes']+=1
            if any(x is None for x in il):
                eventrows.append(dict(profile=kp,labels=[ma,mp],status='nonpositive_required_pair'));counts['empty_required_pair']+=1;continue
            lo=max(-2-off(kp,lp),il[0][0],il[1][0]);hi=min(-off(kp,lp),il[0][1],il[1][1])
            lower('KL_lower',lv,lo);upper('KL_upper',lv,hi)
            arcs=[(0.,TAU)];stopped=None
            for name,a,b,c in constraints:
                arcs=cut(arcs,a,b,c)
                if not arcs:stopped=name;blockers[name]+=1;break
            counts['nonempty_AP_circles'if arcs else'empty_AP_circles']+=1
            # Reproduce all earlier sampled low-A/P K roots as a guard
            # against a mistaken circle sheet, offset or relaxed bound.
            for old in row['caps']:
                if (old['profile'],old['labels'][0],old['labels'][2])!=(kp,ma,mp):continue
                k=np.array(old['normal']);xy=(k-k0)@M@eb.T/math.sqrt(rr);t=math.atan2(xy[1],xy[0])%TAU
                err=max(float(np.max(np.abs(make(t)-k))),abs(float(xy@xy)-1.))
                assert err<2e-5 and all(a*math.cos(t)+b*math.sin(t)+c<=2e-7*max(1.,abs(a),abs(b),abs(c))for _,a,b,c in constraints)
                assert any(a-1e-7<=t<=b+1e-7 for a,b in arcs);recover+=1
            tests=[]
            for rid in sel['R_indices']:
                r=rbank[rid];rp=r['profile'];rv=np.array(r['normal']);ar,pr=r['labels'][1:]
                allowed=[m for m in range(2,7)if P.corner3(ap,rp,kp,ar,ma,m)is not None and P.corner3(pp,rp,kp,pr,mp,m)is not None]
                if not allowed:counts['no_KR_finite_recipe']+=1;continue
                counts['R_circle_tests']+=1
                cr=coeff(rv,off(kp,rp));rg=ranges(cr,arcs)
                jr=[pd_interval(qa,profiles[rp]['square'],qk,float(av@M@rv),rhs[0]),pd_interval(qp,profiles[rp]['square'],qk,float(pv@M@rv),rhs[1])]
                assert all(j is not None for j in jr)
                rlo=max(-2.,*(j[0]+off(kp,rp)for j in jr));rhi=min(0.,*(j[1]+off(kp,rp)for j in jr))
                rg=[(max(x,rlo),min(y,rhi))for x,y in rg if max(x,rlo)<=min(y,rhi)]
                hits=[m for m in allowed if m<6 and any(x-1e-8<=-B.theta(m)<=y+1e-8 for x,y in rg)]
                tail=6 in allowed and any(max(x,-2.)<=min(y,-math.sqrt(3))for x,y in rg)
                counts['KR_discrete_hits']+=len(hits);counts['KR_high_tail_hits']+=int(tail)
                angle_rows=[]
                for m in hits:
                    ts=level_roots(cr,-B.theta(m),arcs)
                    if ts is None:unresolved.append(dict(stage='constant_KR',event=eid,profile=kp,AP=[ma,mp],R=rid,KR=m));continue
                    for t in ts:
                        k=make(t);vals={i:float(k@M@w)+(profiles[kp]['offsets'][oldindex[i]]if i in oldindex else off(kp,types.get(i,0)))for i,w in nd.items()}
                        bad=[i for i,x in vals.items()if B.relation(x)in[None,'parallel_or_unresolved']]
                        cp=float(k@M@h)+profiles[kp]['center_offset']
                        err=max(abs(float(k@M@k)-qk),abs(vals[68]+B.theta(ma)),abs(vals[70]+B.theta(mp)),abs(float(H.A.native(k,kp,rv,rp,profiles))+B.theta(m)))
                        ni=[i for i in actual+[67]if float(k@M@nd[i])>-math.sqrt(qk*float(nd[i]@M@nd[i]))+1e-7]
                        if err>2e-5:unresolved.append(dict(stage='low_AP_residual',event=eid,profile=kp,AP=[ma,mp],R=rid,KR=m,residual=err))
                        ml=B.relation(vals[69]);complete=not bad and cp<-1e-8 and err<=2e-5 and not ni and isinstance(ml,int)and P.corner3(ap,lp,kp,row['AL'],ma,ml)is not None and P.corner3(pp,lp,kp,row['LP'],mp,ml)is not None
                        angle_rows.append(dict(KR=m,angle=t,KL_product=vals[69],KL_relation=ml,bad_native=[[i,vals[i]]for i in bad],center_product=cp,residual=err,nonincidence_failures=ni,complete_native=complete))
                        if complete:
                            both=dict(nd);both[71]=rv;both[72]=k;ids=actual+[67,68,69,70,71,72]
                            ty=dict(types);ty[71]=rp;ty[72]=kp
                            bd=S.C.detailed_cone(both,ids);ic=H.incidence_score([both[i]for i in ids],ids,[profiles[ty.get(i,0)]['square']if i not in oldindex else float(both[i]@M@both[i])for i in ids])
                            good=not bd['warnings']and ic['status']=='numerically_separated_acute_incidence'and bd['open_ends']==ic['open_ends']and{tuple(sorted(v['active']))for v in bd['rays']if v['kind']=='finite'}=={tuple(sorted(v))for v in ic['finite_supports']}
                            rec=dict(event=eid,profile=kp,AP=[ma,mp],R=rid,KR=m,normal=k.tolist(),native_products=sorted(vals.items()),cone=bd,PD=ic)
                            if good:candidates.append(rec)
                            else:unresolved.append(dict(stage='candidate_cone',**rec))
                tests.append(dict(R=rid,allowed_KR=allowed,native_KR_ranges=rg,finite_hits=hits,high_tail_hit=tail,angles=angle_rows))
            eventrows.append(dict(profile=kp,labels=[ma,mp],projection_square=float(k0@M@k0),circle_radius_square=rr,
                KL_projected_interval=[lo,hi],arcs=arcs,first_empty_constraint=stopped,R_tests=tests))
        highrows=[]
        for kp,ms,hi in high:
            # A high finite bond is an isolated I2(m) component in every
            # required spherical corner. In particular KL=KR=2 and K has
            # zero A2 profile. Each saved family has only one high AP bond.
            assert kp==0 and ms[1]==2 and len(hi)==1 and hi[0]in[0,2]
            highpos=hi[0];fixed=70 if highpos==0 else 68;target=68 if highpos==0 else 70
            fixedm=ms[2]if highpos==0 else ms[0]
            assert fixedm<=5 and types[target]==0
            for rid in sel['R_indices']:
                r=rbank[rid];rp=r['profile'];rv=np.array(r['normal']);ar,pr=r['labels'][1:]
                if P.corner3(ap,rp,0,ar,ms[0],2)is None or P.corner3(pp,rp,0,pr,ms[2],2)is None:
                    counts['high_AP_R_classifier_reject']+=1;continue
                counts['high_AP_frame_equations']+=1
                vv,info=N.norm_roots([lv,rv,nd[fixed]],[0.,0.,-B.theta(fixedm)],2.)
                hr=dict(profile=kp,labels=ms,high_position=highpos,R=rid,norm=info,branches=[])
                if info['status'].startswith('unresolved'):
                    unresolved.append(dict(stage='high_AP_norm',event=eid,**hr));highrows.append(hr);continue
                for sg,k in vv:
                    counts['high_AP_norm_branches']+=1
                    vals={i:float(k@M@w)for i,w in nd.items()};vals[71]=float(k@M@rv)
                    product=vals[target];in_tail=-2-1e-8<=product<=-math.sqrt(3)+1e-8
                    counts['high_AP_in_tail']+=int(in_tail)
                    bad=[i for i,x in vals.items()if B.relation(x)in[None,'parallel_or_unresolved']]
                    cp=float(k@M@h);rrr=max(abs(float(k@M@k)-2.),abs(vals[69]),abs(vals[71]),abs(vals[fixed]+B.theta(fixedm)))
                    assert rrr<2e-5
                    nonincident=[i for i in actual+[67]if vals[i]>-math.sqrt(2*float(nd[i]@M@nd[i]))+1e-7]
                    hr['branches'].append(dict(branch=sg,high_product=product,in_tail=in_tail,high_relation=B.relation(product),
                        bad_native=[[i,vals[i]]for i in bad],center_product=cp,nonincidence_failures=nonincident,residual=rrr))
                    if in_tail and not bad and cp<-1e-8 and not nonincident:
                        both=dict(nd);both[71]=rv;both[72]=k;ids=actual+[67,68,69,70,71,72];ty=dict(types);ty[71]=rp;ty[72]=0
                        bd=S.C.detailed_cone(both,ids);ic=H.incidence_score([both[i]for i in ids],ids,[profiles[ty.get(i,0)]['square']if i not in oldindex else float(both[i]@M@both[i])for i in ids])
                        good=not bd['warnings']and ic['status']=='numerically_separated_acute_incidence'and bd['open_ends']==ic['open_ends']and{tuple(sorted(v['active']))for v in bd['rays']if v['kind']=='finite'}=={tuple(sorted(v))for v in ic['finite_supports']}
                        rec=dict(stage='high_AP',event=eid,profile=0,labels=ms,high_position=highpos,R=rid,normal=k.tolist(),native_products=sorted(vals.items()),cone=bd,PD=ic)
                        if good:candidates.append(rec)
                        else:unresolved.append(dict(stage='high_AP_cone',candidate=rec))
                highrows.append(hr)
        records.append(dict(event=eid,AP_gram=gg.tolist(),complement_eigenvalues=eig.tolist(),symbolic_high_AP_families=high,high_AP_tests=highrows,recipes=eventrows))
        print('AP_CIRCLE',eid,'low',len(lows),'high',len(high),'candidates',len(candidates),'unresolved',len(unresolved),flush=True)
    paths=[Path(__file__),Path(F.__file__),fp,Path(D.__file__),pj,Path(U.__file__),Path(P.__file__),Path(N.__file__),Path(B.__file__),Path(H.__file__),Path(S.__file__)]
    return dict(status='NUMERICAL_AP_CIRCLE_NECESSARY_DOMAINS',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        counts=dict(counts),recovered_sampled_low_AP_caps=recover,records=records,candidates=candidates,unresolved=unresolved,empty_constraints=dict(blockers),
        scope='Fixed83 frozen event frames and516 savedR. Low K-A/K-P labels2..5 use all7 inherited A2profiles where finite corner recipes allow them; K-L varies continuously in its closed native/PD necessary interval. Native acuteness, commoncenter and intended intrinsic nonincidence with all retainedold/C are relaxed closed halfplanes on the complete AP circle. Symbolic high K-A/K-P families force KL=KR=2 and reduce to both norm branches of(L,R,other A/P); no high-angle cutoff there. NativeKR labels2..5 tested and the entire symbolic>=6 tail retained as a necessary interval. Numerical arc/range/norm arithmetic only: analytic finite-label reduction, but no interval-certified all-label obstruction and no allparameter claim. Old63/199/10 exact infinite-volume control unchanged.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'recovered',out['recovered_sampled_low_AP_caps'],'candidates',len(out['candidates']),'unresolved',len(out['unresolved']),flush=True)
