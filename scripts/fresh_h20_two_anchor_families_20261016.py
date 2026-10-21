#!/usr/bin/env python3
"""Numerical one-parameter families with three prescribed old contacts.

All other variable contacts are required dotted; isolated extra finite
contacts are outside this diagnostic. No all-cap or exact gate claim.
"""
import argparse, hashlib, itertools, json, math
from collections import Counter
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
from fresh_h20_quad_adjacent_h3_closure_20261016 import geometry158, face_data, poly, gate_mask
from fresh_h20_release144_cap_20261016_residual_scout import close

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_two_anchor_families_20261016'
SIDES=[31,95,132,142,145,150,152,156,157]
LABELS=[2,3,4,5,6,8,10]

def real_roots(c):
    a=np.array(c,dtype=float);scale=max(abs(a));
    if scale<1e-14:return []
    a=a/scale
    if abs(a[2])<1e-12:
        return [] if abs(a[1])<1e-12 else [-a[0]/a[1]]
    disc=a[1]**2-4*a[2]*a[0]
    if disc < -1e-12:return []
    disc=math.sqrt(max(0.,disc))
    q=-.5*(a[1]+math.copysign(disc,a[1]))
    if abs(q)<1e-15:return [-a[1]/(2*a[2])]
    return sorted(set([q/a[2],a[0]/q]))

def point(f,t):
    p=np.array([1.,t,t*t]);return p@np.array(f['numerator'])/(p@np.array(f['denominator']))

def middle(a,b):
    if a is None and b is None:return 0.
    if a is None:return b-max(1.,abs(b))
    if b is None:return a+max(1.,abs(a))
    return (a+b)/2

def allowed_intervals(P,positive):
    # P[0]+P[1]*t+P[2]*t^2 >= 0; denominators are positive.
    roots=sorted(r for c in P for r in real_roots(c) if np.isfinite(r) and (not positive or r>0))
    distinct=[]
    for r in roots:
        if not distinct or abs(r-distinct[-1])>2e-9*max(1,abs(r)):distinct.append(r)
    edges=([0.] if positive else [None])+distinct+[None]
    intervals=[]
    scale=np.maximum(np.max(abs(P),axis=1),1e-30)
    for a,b in zip(edges,edges[1:]):
        t=middle(a,b)
        if positive and t<=0:continue
        value=P@np.array([1.,t,t*t])
        if min(value/scale)>-2e-9:
            intervals.append([a,b])
    return intervals

def constant_gate(c):
    if abs(c)<2e-6:return 'right'
    if c>=1-2e-6:return 'dotted_or_parallel'
    if c<=0:return None
    m=np.pi/np.arccos(c)
    return int(round(m)) if m>=2.9999 and abs(m-round(m))<2e-4 else None

def charts(p,K,M):
    D,Q=np.linalg.eigh(K.T@M@K);E=K@Q
    for j in range(2):
        E[:,j]/=np.sqrt(abs(D[j]))
        k=np.argmax(abs(E[:,j]));E[:,j]*=1 if E[k,j]>0 else -1
    w=2-p@M@p;out=[]
    if min(D)>0:
        if w < -1e-7:return out
        if abs(w)<1e-7:return [dict(kind='point',numerator=[p.tolist(),[0.]*5,[0.]*5],denominator=[1.,0.,0.],positive=False)]
        e,f=np.sqrt(w)*E.T
        out.append(dict(kind='circle',numerator=np.array([p+e,2*f,p-e]).tolist(),denominator=[1.,0.,1.],positive=False))
    else:
        assert D[0]<0<D[1]
        et,es=E.T
        if abs(w)<1e-7:
            for sign in [-1,1]:out.append(dict(kind='null_line',branch=sign,numerator=np.array([p,es+sign*et,np.zeros(5)]).tolist(),denominator=[1.,0.,0.],positive=False))
        elif w>0:
            for sign in [-1,1]:
                a=np.sqrt(w)/2*(sign*es+et);b=np.sqrt(w)/2*(sign*es-et)
                out.append(dict(kind='spacelike_hyperbola',branch=sign,numerator=np.array([b,p,a]).tolist(),denominator=[0.,1.,0.],positive=True))
        else:
            for sign in [-1,1]:
                a=np.sqrt(-w)/2*(sign*et+es);b=np.sqrt(-w)/2*(sign*et-es)
                out.append(dict(kind='timelike_hyperbola',branch=sign,numerator=np.array([b,p,a]).tolist(),denominator=[0.,1.,0.],positive=True))
    return out

def samples(interval):
    a,b=interval;out=[middle(a,b)]
    if a is not None and b is not None:
        out += [a+(b-a)*i/16 for i in range(1,16)]
    elif a is not None:out += [a+2.**j*max(1.,abs(a)) for j in range(-8,9)]
    elif b is not None:out += [b-2.**j*max(1.,abs(b)) for j in range(-8,9)]
    else:out += [s*2.**j for s in [-1,1] for j in range(-8,9)]
    result=[]
    for t in out:
        t=float(Fraction(float(t)).limit_denominator(10000))
        if (a is None or t>a+1e-10)and(b is None or t<b-1e-10)and t not in result:result.append(t)
    return result

def run():
    G=geometry158();X,M=G['X'],G['metric'];q=G['norms'];H=(G['H']@G['g']@G['F'].T)@M
    base,V,A,b,ids=face_data(G);old=poly(base,V,A,b,ids)
    outside=np.array([p['z']for p in old['vertices']if p['radius2']>1+1e-6])
    joint=json.loads((ROOT/'results/fresh_h20_nonright151_h3_joint_cover_20261016.json').read_text())
    ideal=np.array([p['z']for p in joint['samples']])
    stats=Counter();families=[];pool=[];seen=set();frame_log=[]
    for pair in itertools.combinations(SIDES,2):
        frame=[151,*pair];U=X[frame];B=U@M@U.T
        if min(abs(np.linalg.eigvalsh(B)))<1e-8:stats['singular_frames']+=1;continue
        inv=np.linalg.inv(B);K=null_space(U@M)
        assert K.shape==(5,2)
        start=len(families)
        for labs in itertools.product([2,3,4,5,6],LABELS,LABELS):
            cs=np.cos(np.pi/np.array(labs));cs[np.array(labs)==2]=0
            f=-np.sqrt(2*q[frame])*cs;p=(inv@f)@U
            stats['patterns']+=1
            for model in charts(p,K,M):
                num=np.array(model['numerator']);den=np.array(model['denominator'])
                co=-(num@M@X.T).T/np.sqrt(2*q)[:,None]
                fixed=[];constraints=[];who=[];bad=False
                # Constant contacts are detected from polynomial proportionality.
                di=np.argmax(abs(den))
                for j,c in enumerate(co):
                    value=c[di]/den[di]
                    if j in frame or max(abs(c-value*den))<1e-7:
                        kind=constant_gate(value)
                        if kind is None:bad=True;break
                        fixed.append(dict(wall=j,cosine=float(value),kind=kind))
                    else:
                        constraints.append(c-den);who.append(j)
                if bad:stats['fixed_bad_contact']+=1;continue
                constraints.extend([-num@M@H,num[:,0]]);who+=['center','positive_time']
                P=np.array(constraints)
                ranges=allowed_intervals(P,model['positive'])
                if not ranges:stats['empty_all_dotted_intervals']+=1;continue
                model.update(id=len(families),frame=frame,labels=list(labs),intervals=ranges,
                    fixed_contacts=fixed,inequality_rows=who,inequality_coefficients=P.tolist(),
                    particular_norm=float(p@M@p),frame_inertia=list(map(int,np.sign(np.linalg.eigvalsh(B)))))
                families.append(model)
                for k,interval in enumerate(ranges):
                    for t in samples(interval):
                        r=point(model,t)
                        if abs(r@M@r-2)>2e-4:stats['sample_numeric_norm_failure']+=1;continue
                        co0=-(r@M@X.T)/np.sqrt(2*q)
                        if not gate_mask(co0[None,:])[0]or r@M@H>=0 or r[0]<=0:
                            stats['sample_gate_roundoff_rejection']+=1;continue
                        ar=-r[1:]/r[0];cut=np.where(outside@ar>1+1e-7)[0]
                        if not len(cut):continue
                        key=tuple(np.round(r,7))
                        if key in seen:continue
                        seen.add(key)
                        pool.append(dict(id=len(pool),family=model['id'],interval=k,parameter=t,
                            parameter_rational=str(Fraction(t).limit_denominator(10000)),
                            orthonormal_chart=r.tolist(),residual_constraint=ar.tolist(),
                            old_exterior_cuts=cut.tolist(),joint_sample_cuts=np.where(ideal@ar>1+1e-7)[0].tolist()))
        frame_log.append(dict(frame=frame,real_families=len(families)-start))
        print('FRAME',frame,'families',len(families)-start,'total',len(families),'samples',len(pool),flush=True)
    sources=[f'scripts/{TAG}.py','scripts/fresh_h20_quad_adjacent_h3_closure_20261016.py',
       'scripts/fresh_h20_pentagon_closure_scout_20261016.py','results/fresh_h20_compact_quad_cap_20261016.json',
       'results/fresh_h20_nonright151_h3_joint_cover_20261016.json']
    return dict(status='PASS_NUMERICAL_TWO_ANCHOR_CONTINUOUS_FAMILY_INTERVALS',
       sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sources},
       side_domain=SIDES,side_labels=LABELS,label151=[2,3,4,5,6],old_face=old,
       stats=dict(stats),frames=frame_log,families=families,candidates=pool,joint_samples=joint['samples'],
       scope='Numerical interval isolation for the stated finite prescribed-label domains. All unprescribed variable old contacts are required dotted; isolated extra finite contacts and other anchor/label patterns are not covered. No exact wall, all-cap exclusion or finite-volume claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();path=ROOT/f'results/{TAG}.json'
    if a.verify:assert close(d,json.loads(path.read_text()))
    else:path.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['families']),len(d['candidates']),flush=True)
