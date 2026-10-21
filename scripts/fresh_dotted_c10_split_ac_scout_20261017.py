#!/usr/bin/env python3
"""Continuous mutual-angle events after splitting A--C by old7 connector J.

Bounded numerical discovery; no exact chamber promotion.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import product,combinations
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import fresh_c10_endpoint_disk_20261017 as D
P,T,N,S,H,B,M,ROOT=D.P,D.T,D.N,D.S,D.H,D.B,D.M,D.ROOT
TAG='fresh_dotted_c10_split_ac_scout_20261017'
UGRID=[2.0001,2.005,2.01,2.02,2.03,2.04,2.06,2.08,2.1,2.12,2.16,2.25,2.5]

def run():
    cur,labels,ns,h,pro,primary=H.setup();actual=[i for i in cur['actual_sides']if i not in S.DELETE]
    pp=ROOT/'results/fresh_old8_triangle_double_flip_20261017.json';parent=json.loads(pp.read_text())
    a=np.array(parent['prefixes'][0]['normals'][0]);f=np.array([ns[i]for i in[7,30,10]])
    v=(np.linalg.inv(f@M@f.T)@np.array([0.,0.,-1.]))@f;alpha=float(-v@M@v)
    roots,_=N.norm_roots(f,[0.,0.,-2.],2.)
    k=next((c-2*v)/np.sqrt(1+2*alpha)for _,c in roots if(c-2*v)@M@ns[55]<0)
    def C(u):return u*v+np.sqrt(1+alpha*u*u/2)*k
    counts=Counter();blocks=Counter();unresolved=[];samples=[];events=[];extensions=[]
    def offset(pi,pj):return float(H.A.native(np.zeros(4),pi,np.zeros(4),pj,pro))
    def native(v,pi,assigned=()):
        bad=T.gates(T.native_values(v,pi,labels,ns,pro),assigned)
        if bad:return ('old'+str(bad[0][0]),bad[0][1])
        if float(v@M@h)+pro[pi]['center_offset']>=-1e-8:return ('center',0.)
        return None
    def dotted(v,pi,w,pj):return float(v@M@w)<-np.sqrt(pro[pi]['square']*pro[pj]['square'])-1e-7
    def nonincident(v,pi,allow):return all(float(v@M@ns[i])<-np.sqrt(pro[pi]['square']*float(ns[i]@M@ns[i]))-1e-7 for i in actual if i not in allow)
    assert abs(float(ns[7]@M@ns[7])-pro[4]['square'])<1e-7
    assert native(a,0,[7,8])is None and nonincident(a,0,[7,8])
    def checkC(c,u):
        assert max(abs(float(c@M@c)-2),abs(float(c@M@ns[7])),abs(float(c@M@ns[30])),abs(float(c@M@ns[10])+u))<1e-7
        assert native(c,0,[7,30])is None and nonincident(c,0,[7,30])
    def norm(stage,frame,rhs,pi):
        counts[stage+'_equations']+=1;roots,info=N.norm_roots(frame,rhs,pro[pi]['square'])
        if info['status'].startswith('unresolved'):unresolved.append(dict(stage=stage,info=info));return []
        out=[]
        for sg,w in roots:
            err=max(abs(float(w@M@w)-pro[pi]['square']),max(abs(float(w@M@z)-q)for z,q in zip(frame,rhs)))
            if err>1e-5:unresolved.append(dict(stage=stage,error=err));continue
            out.append((sg,w));counts[stage+'_branches']+=1
        return out
    jdomain=[];xdomain=[]
    for pi,m7,ma,mc in product(range(7),range(2,6),range(2,7),range(2,7)):
        if P.corner3(4,0,pi,2,m7,ma)is not None and P.corner3(4,0,pi,2,m7,mc)is not None:jdomain.append((pi,m7,ma,mc))
    for pi,m8,ma,mc in product(range(7),range(2,6),range(2,6),range(2,6)):
        if P.corner3(0,0,pi,2,m8,ma)is not None:xdomain.append((pi,m8,ma,mc))
    def bank(u,stage):
        c=C(u);checkC(c,u);out={};domain=jdomain if stage=='J'else xdomain
        for pi,mo,ma,mc in domain:
            old=7 if stage=='J'else 8;op=4 if old==7 else 0
            rhs=[-B.theta(mo)-offset(pi,op),-B.theta(ma),-B.theta(mc)]
            for sg,w in norm(stage,[ns[old],a,c],rhs,pi):
                bad=native(w,pi,[old])
                if bad:blocks[stage,bad[0]]+=1;continue
                if not nonincident(w,pi,[old]):blocks[stage,'intrinsic_nonincidence']+=1;continue
                out[pi,mo,ma,mc,sg]=w
        return out
    def root(u,stage,key):
        pi,mo,ma,mc,sg=key;old=7 if stage=='J'else 8;op=4 if old==7 else 0
        rows=norm('event_'+stage,[ns[old],a,C(u)],[-B.theta(mo)-offset(pi,op),-B.theta(ma),-B.theta(mc)],pi)
        return next((w for s,w in rows if s==sg),None)
    def score(vs,ps):
        nd=dict(ns);nd.update({67+i:w for i,w in enumerate(vs)});sel=actual+list(range(67,67+len(vs)))
        bd=S.C.detailed_cone(nd,sel);ic=H.incidence_score([nd[i]for i in sel],sel,[float(nd[i]@M@nd[i])for i in sel])
        if bd['warnings']or ic['status']!='numerically_separated_acute_incidence'or bd['open_ends']!=ic['open_ends']or bd['counts'].get('finite',0)!=ic['finite']:
            unresolved.append(dict(stage='boundary',bd=bd,ic=ic));return None
        if {tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'}!={tuple(sorted(r))for r in ic['finite_supports']}:
            unresolved.append(dict(stage='support_equality'));return None
        return dict(boundary=bd,PD_incidence=ic)
    previous=None;seen=set()
    for u in UGRID:
        jb,xb=bank(u,'J'),bank(u,'X');samples.append(dict(u=u,J_count=len(jb),X_count=len(xb),AC=float(a@M@C(u))))
        if previous is not None:
            lo,pj,px=previous
            for jk,xk in product(set(jb)&set(pj),set(xb)&set(px)):
                jp,j7,ja,jc,js=jk;xp,x8,xa,xc,xs=xk
                for m in range(2,6):
                    if P.corner3(0,jp,xp,ja,xa,m)is None or P.corner3(0,jp,xp,jc,xc,m)is None:continue
                    target=-B.theta(m)-offset(jp,xp)
                    f0=float(pj[jk]@M@px[xk])-target;f1=float(jb[jk]@M@xb[xk])-target
                    if f0*f1>0:continue
                    def fun(t):
                        j,x=root(t,'J',jk),root(t,'X',xk)
                        if j is None or x is None:raise ValueError('missing norm branch')
                        return float(j@M@x)-target
                    try:ue=brentq(fun,lo,u,xtol=1e-13)
                    except ValueError as e:unresolved.append(dict(stage='event_root',lo=lo,hi=u,J=jk,X=xk,message=str(e)));continue
                    key=(jk,xk,m,round(ue,10))
                    if key in seen:continue
                    seen.add(key);j,x=root(ue,'J',jk),root(ue,'X',xk);c=C(ue);checkC(c,ue)
                    if native(j,jp,[7])or native(x,xp,[8])or not nonincident(j,jp,[7])or not nonincident(x,xp,[8]):continue
                    if abs(float(j@M@x)-target)>1e-7:unresolved.append(dict(stage='event_residual'));continue
                    vs=[a,c,j,x];ps=[0,0,jp,xp]
                    pairdata=[dict(pair=[i,l],product=float(H.A.native(vs[i],ps[i],vs[l],ps[l],pro)))for i,l in combinations(range(4),2)]
                    if any(B.relation(p['product'])in[None,'parallel_or_unresolved']for p in pairdata):continue
                    sc=score(vs,ps)
                    if sc is None:continue
                    jr=[r for r in sc['boundary']['rays']if 69 in r['active']]
                    if len(jr)!=4 or any(r['kind']!='finite'for r in jr):unresolved.append(dict(stage='J_quad',rays=jr));continue
                    events.append(dict(u=ue,J_recipe=list(jk),X_recipe=list(xk),JX=m,normals=[w.tolist()for w in vs],profiles=ps,native_pairs=pairdata,**sc))
        previous=(u,jb,xb)
        print('SPLIT_AC',u,'J',len(jb),'X',len(xb),'events',len(events),flush=True)
    # Test the changed four-prefix against the complete necessary R circle
    # and the bounded low-Z polar-Y duty, including all extra J contacts.
    circle=json.loads((ROOT/'results/fresh_r_circle_admissibility_20261017.json').read_text())
    rs={r['label']:np.array(r['normal'])for r in circle['numerical_R_candidates']if r['label']<=15 and r['native_and_center_pass']}
    for ei,event in enumerate(events):
        aa,c,j,x=map(np.array,event['normals']);jp,xp=event['profiles'][2:];m8x=event['X_recipe'][1]
        for mr,r in rs.items():
            if not all(dotted(r,0,w,pi)for w,pi in zip([aa,c,j,x],[0,0,jp,xp])):continue
            for zp,m8z,mrz,mxz in product(range(7),range(2,6),range(2,6),range(2,6)):
                if P.corner3(0,xp,zp,m8x,m8z,mxz)is None:continue
                if P.corner3(0,0,zp,2,m8z,2)is None or P.corner3(0,0,zp,2,mrz,2)is None:continue
                if pro[zp]['square']-(B.theta(m8z)**2+B.theta(mrz)**2)/2>=-1e-9:continue
                for sz,z in norm('Z',[ns[8],r,x],[-B.theta(m8z),-B.theta(mrz),-B.theta(mxz)-offset(zp,xp)],zp):
                    bad=native(z,zp,[8])
                    if bad:blocks['Z',bad[0]]+=1;continue
                    if not nonincident(z,zp,[8])or not all(dotted(z,zp,w,pi)for w,pi in zip([aa,c,j],[0,0,jp])):continue
                    for sy,y in norm('Y',[ns[8],r,z],[0.,0.,0.],0):
                        if native(y,0,[8])or not nonincident(y,0,[8]):continue
                        if not all(dotted(y,0,w,pi)for w,pi in zip([aa,c,j,x],[0,0,jp,xp])):continue
                        vs=[aa,c,j,x,r,z,y];ps=[0,0,jp,xp,0,zp,0]
                        pairs=[dict(pair=[i,l],product=float(H.A.native(vs[i],ps[i],vs[l],ps[l],pro)))for i,l in combinations(range(7),2)]
                        if any(B.relation(q['product'])in[None,'parallel_or_unresolved']for q in pairs):continue
                        sc=score(vs,ps)
                        if sc is not None:extensions.append(dict(event=ei,R44=mr,Z_profile=zp,Z_labels=[m8z,mrz,mxz],branches=[sz,sy],normals=[w.tolist()for w in vs],profiles=ps,native_pairs=pairs,**sc))
        print('SPLIT_AC_EXTENSION',ei+1,'of',len(events),'extensions',len(extensions),flush=True)
    paths=[Path(__file__),Path(D.__file__),pp,H.PARENT,H.BOUNDARY,ROOT/'results/fresh_dotted_c10_domain_20261017.json',ROOT/'results/fresh_r_circle_admissibility_20261017.json']
    return dict(status='BOUNDED_DOTTED_C10_SPLIT_AC_CONTINUOUS_EVENT_DISCOVERY',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
      topology=dict(move='Subdivide interiorAC shared by7AC andACX usingJ',removed_triangles=[['7','A','C'],['A','C','X']],added_triangles=[['7','A','J'],['7','J','C'],['A','J','X'],['J','C','X']],new_nonincidence='A/C dotted',local_goal='whole compactJ quadrilateral with sides7,A,X,C'),
      domain=dict(u_samples=UGRID,event_intervals='adjacent samples where both literal branches pass native gates',J_recipes=len(jdomain),X_recipes=len(xdomain),JX_labels=list(range(2,6)),low_Z_labels=list(range(2,6))),
      samples=samples,counts=dict(counts),blockers=[[list(k),v]for k,v in sorted(blocks.items())],four_prefix_events=events,seven_prefix_extensions=extensions,unresolved=unresolved,
      exact_old_deletion_cost=H.deletion(primary,cur,S.DELETE),
      scope='Bounded numerical continuous-angle discovery after a new edge subdivision. Every reported four-prefix passes all retained native contacts, common center, intrinsic nonincidences, all six new pair gates, complete numerical S6 cone, independent PD support census, and exactly four finiteJ-incident vertices. No exact chamber or volume claim. Unbracketed events, high labels, otherA choices and other patch incidences remain open.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],'events',len(out['four_prefix_events']),'extensions',len(out['seven_prefix_extensions']),'unresolved',len(out['unresolved']),flush=True)
