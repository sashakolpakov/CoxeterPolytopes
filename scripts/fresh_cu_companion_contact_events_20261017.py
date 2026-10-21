#!/usr/bin/env python3
"""Numerical moving-frame companion events on the exact all63 C(u) tail.

The moving-frame identities are analytic. Roots, interval cells, old-wall
tests and complete target cones below use floating point, not interval proof.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_c3_companion_caps_20261017 as C
H,B,M,N,S,P,ROOT=C.H,C.B,C.M,C.N,C.S,C.P,C.ROOT
TAG='fresh_cu_companion_contact_events_20261017'
HIGH=[6,8,12,24]

def setup():
    cur,labels,ns,h,profiles,primary=C.setup()
    f=np.array([ns[i]for i in [7,30,10]])
    v=(np.linalg.inv(f@M@f.T)@np.array([0.,0.,-1.]))@f
    alpha=float(-v@M@v);gamma=np.sqrt(alpha/2)
    k=(ns[67]-3*v)/np.sqrt(1+4.5*alpha)
    e0=v/np.sqrt(alpha);e1=k/np.sqrt(2)
    assert max(abs(e0@M@e0+1),abs(e1@M@e1-1),abs(e0@M@e1))<1e-7
    ef=np.array([ns[i]for i in [7,30]]);eg=ef@M@ef.T
    assert np.linalg.eigvalsh(eg).min()>0
    t0=3*gamma+np.sqrt(1+9*gamma*gamma)
    def cv(t):return (t-1/t)/(2*gamma)*v+(t+1/t)/2*k
    return cur,labels,ns,h,profiles,primary,ef,np.linalg.inv(eg),e0,e1,gamma,t0,cv

def equation_roots(a,b,c):
    # a*t + b/t + c = 0. Retain every positive real quadratic root.
    scale=max(abs(a),abs(b),abs(c),1.)
    if abs(a)<1e-13*scale:
        return [] if abs(c)<1e-13*scale else [-b/c]
    disc=c*c-4*a*b
    if disc<0:return []
    sd=np.sqrt(disc);q=-.5*(c+np.copysign(sd,c))
    rr=[q/a]
    if q:rr.append(b/q)
    return sorted(set(float(t)for t in rr if np.isfinite(t)and t>0))

def run():
    cur,labels,ns,h,profiles,primary,ef,ei,e0,e1,gamma,t0,cv=setup()
    dom=C.domains(ns)[1];assert dom['frame']==[7,30,67]
    counts=Counter();blockers=Counter();families=[];events=[];cells=[];unresolved=[];scores=Counter();improved=[]
    offsets=np.array([p['offsets']for p in profiles]);old=np.array([ns[i]for i in labels])
    targets=[i for i in cur['actual_sides']if i not in[7,30]]
    def score(v,t,pi):
        nd=dict(ns);nd[67]=cv(t);nd[68]=v;sel=cur['actual_sides']+[67,68]
        bd=S.C.detailed_cone(nd,sel)
        ic=H.incidence_score([nd[i]for i in sel],sel,[float(nd[i]@M@nd[i])for i in sel])
        good=not bd['warnings']and ic['status']=='numerically_separated_acute_incidence'
        good=good and bd['open_ends']==ic['open_ends']and bd['counts'].get('finite',0)==ic['finite']
        good=good and {tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'}=={tuple(sorted(r))for r in ic['finite_supports']}
        if not good:return None,dict(boundary=bd,PD_incidence=ic)
        cf=[r for r in bd['rays']if 67 in r['active']]
        return dict(boundary=bd,PD_incidence=ic,complete_C_face=cf),None
    for rec in dom['recipes']:
        for hh in product(HIGH,repeat=len(rec['high_positions'])):
            ms=list(rec['labels'])
            for i,m in zip(rec['high_positions'],hh):ms[i]=m
            pi=rec['profile'];q=profiles[pi]['square']
            rhs=np.array([-B.theta(ms[j])-offsets[pi,labels.index(i)]for j,i in enumerate([7,30])])
            a0=rhs@ei@ef;rc=-B.theta(ms[2])
            lam=float(a0@M@a0+rc*rc/2-q)
            counts['literal_recipes']+=1
            if lam<-1e-8:counts['no_real_branches']+=1;continue
            if abs(lam)<=1e-8:
                unresolved.append(dict(stage='tangent_recipe',profile=pi,labels=ms,lambda_value=lam));continue
            for sign in[-1,1]:
                # A=a0+rc*C/2+sign*sqrt(lam)*T, T=s*e0+gamma*u*e1.
                # t=gamma*u+s turns every native product into at+b/t+c.
                b=sign*np.sqrt(lam)
                ap=(rc*np.sqrt(2)/4+b/2)*(e0+e1)
                am=(-rc*np.sqrt(2)/4+b/2)*(e0-e1)
                # C=t/sqrt2*(e0+e1)+1/(sqrt2*t)*(-e0+e1).
                assert max(abs(ap@M@ap),abs(am@M@am),abs(a0@M@ap),abs(a0@M@am))<1e-6
                fid=len(families);counts['moving_branches']+=1
                aa=old@M@ap;bb=old@M@am;cc=old@M@a0+offsets[pi]
                ca=float(ap@M@h);cb=float(am@M@h);c0=float(a0@M@h)+profiles[pi]['center_offset']
                f=dict(index=fid,profile=pi,labels=ms,sign=sign,lambda_value=lam,
                    constant=a0.tolist(),positive_t=ap.tolist(),inverse_t=am.tolist())
                families.append(f)
                def accept(t,origin):
                    counts['parameter_tests']+=1
                    if not np.isfinite(t)or t<t0-1e-9:return
                    v=a0+t*ap+am/t;c=cv(t)
                    err=max(abs(float(v@M@v)-q),abs(float(v@M@c)-rc),max(abs(float(v@M@w)-r)for w,r in zip(ef,rhs)))
                    if err>2e-5:
                        unresolved.append(dict(family=fid,t=t,stage='residual',error=err));return
                    vals=aa*t+bb/t+cc;bad=[labels[j]for j,x in enumerate(vals)if B.relation(float(x))in[None,'parallel_or_unresolved']]
                    if bad:blockers[tuple(bad)]+=1;return
                    center=ca*t+cb/t+c0
                    if center>=-1e-8:counts['center_rejected']+=1;return
                    # The exact C tail is valid, but verify its numerical realization here too.
                    cvals=old@M@c
                    if any(B.relation(float(x))in[None,'parallel_or_unresolved']for x in cvals):
                        unresolved.append(dict(family=fid,t=t,stage='C_tail_numerics'));return
                    sc,error=score(v,t,pi)
                    if error:
                        unresolved.append(dict(family=fid,t=t,stage='cone',**error));return
                    u=(t-1/t)/(2*gamma)
                    out=dict(family=fid,u=float(u),t=float(t),origin=origin,normal=v.tolist(),C_normal=c.tolist(),
                      norm_residual=err,center_product=float(center),native_products=[[i,float(x)]for i,x in zip(labels,vals)],
                      finite_contacts=[[i,B.relation(float(x))]for i,x in zip(labels,vals)if isinstance(B.relation(float(x)),int)],**sc)
                    key=(sc['boundary']['counts'].get('finite',0),sc['boundary']['counts'].get('exterior',0),sc['boundary']['open_ends'])
                    scores[key]+=1;counts['coherent_parameters']+=1
                    if key[2]<10:improved.append(out)
                    return out
                # All finite2..5 events to a different actual old side; no u cutoff.
                seen=set()
                for target,m in product(targets,range(2,6)):
                    j=labels.index(target)
                    for t in equation_roots(aa[j],bb[j],cc[j]+B.theta(m)):
                        if t<t0-1e-9:continue
                        key=round(t,9)
                        if key in seen:continue
                        seen.add(key);out=accept(t,dict(kind='finite_contact_event',wall=target,label=m))
                        if out:events.append(out)
                # Necessary dotted/acute cell endpoints are algebraic roots.
                # One sample per cell is discovery only, not certified coverage.
                cuts=[t0]
                for a,bb0,c in zip(aa,bb,cc):
                    cuts.extend(t for target in[-2.,0.]for t in equation_roots(a,bb0,c-target)if t>t0+1e-9)
                cuts.extend(t for t in equation_roots(ca,cb,c0)if t>t0+1e-9)
                cuts=sorted(set(round(t,12)for t in cuts));ends=cuts+[None]
                for lo,hi in zip(ends,ends[1:]):
                    t=lo*2 if hi is None else np.sqrt(lo*hi)
                    out=accept(t,dict(kind='open_cell_sample',t_interval=[lo,hi]))
                    if out:cells.append(out)
                if t0 not in seen:
                    out=accept(t0,dict(kind='tail_start',u=3))
                    if out:cells.append(out)
        print('CU_MOVING_RECIPE',rec['profile'],rec['labels'],'families',len(families),'events',len(events),'cells',len(cells),'improved',len(improved),flush=True)
    paths=[Path(__file__),Path(C.__file__),Path(H.__file__),Path(B.__file__),Path(S.C.__file__),
      H.PARENT,H.BOUNDARY,ROOT/'results/fresh_dotted_c10_domain_20261017.json']
    return dict(status='NUMERICAL_CONTINUOUS_CU_COMPANION_CONTACT_EVENTS',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
      retained_control_walls=63,deleted_walls=[],parameter_domain='u>=3',t0=t0,gamma=gamma,
      finite_frame_domain=dom,high_label_representatives=HIGH,event_labels=[2,3,4,5],event_target_walls=targets,
      counts=dict(counts),moving_families=families,coherent_events=events,coherent_cell_samples=cells,improvements=improved,
      full_cone_histogram=[dict(finite=f,exterior=e,open_ends=i,count=n)for(f,e,i),n in sorted(scores.items())],
      blockers=[dict(walls=list(k),count=n)for k,n in sorted(blockers.items())],unresolved=unresolved,
      scope='Analytic moving-frame parameterization covers both branches for each fixed finite recipe, over u>=3 with no parameter cutoff. Recipe labels>=6 are only6/8/12/24 representatives; isolated old contact events are tested only at labels2..5. All roots, interval sampling and geometry are floating point. Every accepted sample retains all63 native gates, the strict common center, and has complete target cone/independent PD support agreement. No numerical negative is an exact continuous exclusion, no ambient cone or finite-volume result, and no old wall deletion. Exact old199finite vertices/compact faces would persist in any exact acute realization.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'events',len(out['coherent_events']),'improved',len(out['improvements']),'unresolved',len(out['unresolved']),flush=True)
