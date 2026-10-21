#!/usr/bin/env python3
"""Bounded numerical recovery preserving intrinsic finite neighbor graphs."""
import argparse, ast, hashlib, heapq, itertools, json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from fresh_h20_two_anchor_families_20261016 import point
from fresh_h20_quad_adjacent_h3_closure_20261016 import geometry158, face_data, poly, gate_mask
from fresh_h20_release144_cap_20261016_residual_scout import close

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_near_closed_cluster_recovery_20261016'
BASE='results/fresh_h20_two_anchor_families_20261016.json'
ACUTE='results/fresh_h20_two_anchor_families_20261016_acute_cover.json'

def adjacency(R,G):
    M=G['metric']; p=G['X'][151]
    V=R-(R@M@p)[:,None]*p/2
    q=np.einsum('ij,jk,ik->i',V,M,V)
    # X already projects off S16. Native ambient norms would incorrectly
    # count walls which do not meet this H3 face.
    oldq=np.einsum('ij,jk,ik->i',G['X'],M,G['X'])-(G['X']@M@p)**2/2
    old=np.where(oldq>1e-6)[0]
    assert old.tolist()==[95,142,145,150,152,156,157]
    U=G['X'][old]-(G['X'][old]@M@p)[:,None]*p/2
    gg=V@M@V.T
    adj=q[:,None]*q[None,:]-gg*gg>1e-7
    np.fill_diagonal(adj,False)
    oa=q[:,None]*oldq[old][None,:]-(V@M@U.T)**2>1e-7
    return adj,oa,old,q

def best_products(options,limit):
    """Complete first limit products ordered by sum of independent costs."""
    if not options:
        yield ()
        return
    zero=tuple(0 for _ in options);seen={zero}
    heap=[(sum(x[0][1] for x in options),zero)]
    for _ in range(limit):
        if not heap:return
        _,ix=heapq.heappop(heap)
        yield tuple(options[j][i][0] for j,i in enumerate(ix))
        for j in range(len(ix)):
            if ix[j]+1>=len(options[j]):continue
            yy=list(ix);yy[j]+=1;yy=tuple(yy)
            if yy in seen:continue
            seen.add(yy)
            heapq.heappush(heap,(sum(options[k][i][1] for k,i in enumerate(yy)),yy))

def run(limit,nfev):
    D=json.loads((ROOT/BASE).read_text());A=json.loads((ROOT/ACUTE).read_text())
    G=geometry158();M=G['metric'];base,V,AA,bb,ids=face_data(G)
    Z=np.array(A['samples']);scan=[];eligible=[]
    for row in A['steps']:
        R=np.array([x['orthonormal_chart'] for x in row['roots']])
        adj,oa,old,q=adjacency(R,G);degrees=adj.sum(1)+oa.sum(1)
        scan.append(dict(iteration=row['iteration'],potential_degrees=degrees.tolist(),
            intrinsic_new_pairs=[[i,j] for i in range(len(R)) for j in range(i+1,len(R)) if adj[i,j]],
            old_neighbors=[[int(old[j]) for j in np.where(oa[i])[0]] for i in range(len(R))],
            maximum_squared_radius=row['face']['maximum_squared_radius']))
        if min(degrees)>=3:eligible.append(row['iteration'])
    selected=[22,26,28,63,78]
    assert all(i in eligible for i in selected)
    labels=[2,3,4,5,6,7,8,9,10,12,16,20,32,0]
    cs=np.array([0 if m==2 else 1 if m==0 else np.cos(np.pi/m) for m in labels])
    records=[];success=[]
    for ii in selected:
        row=A['steps'][ii];ff=[D['families'][x['family']] for x in row['roots']]
        t0=np.array([x['parameter'] for x in row['roots']]);lo=[];hi=[]
        for f,t in zip(ff,t0):
            a,b=next((a,b) for a,b in f['intervals'] if (a is None or a<t) and (b is None or t<b))
            lo.append(-np.inf if a is None else a+1e-10*max(1,abs(a)))
            hi.append(np.inf if b is None else b-1e-10*max(1,abs(b)))
        def roots(t):return np.array([point(f,v) for f,v in zip(ff,t)])
        R0=roots(t0);co0=-R0@M@R0.T/2;adj,oa,old,q=adjacency(R0,G)
        c151=-(R0@M@G['X'][151])/2
        threshold=np.sqrt((1-c151*c151)[:,None]*(1-c151*c151)[None,:])-c151[:,None]*c151[None,:]
        finite=[(i,j) for i in range(len(ff)) for j in range(i+1,len(ff)) if co0[i,j]<1-1e-7]
        dotted=[(i,j) for i in range(len(ff)) for j in range(i+1,len(ff)) if co0[i,j]>=1-1e-7]
        options=[]
        for p in finite:
            allowed=[k for k in range(len(labels)) if not adj[p] or cs[k]<threshold[p]-1e-8]
            assert allowed
            opts=sorted([(k,float((cs[k]-co0[p])**2)) for k in allowed],key=lambda x:x[1])
            # This is an explicit bounded choice, not all angle assignments.
            options.append(opts[:2])
        best=1e99
        for choice in best_products(options,limit):
            target=np.array([cs[k] for k in choice])
            def residual(t):
                R=roots(t);co=-R@M@R.T/2;ar=-R[:,1:]/R[:,0,None]
                return np.r_[[co[p]-v for p,v in zip(finite,target)],
                    [max(1+1e-8-co[p],0) for p in dotted],
                    100*np.maximum(1e-8-np.max(Z@ar.T-1,axis=1),0)]
            fit=least_squares(residual,t0,bounds=(lo,hi),max_nfev=nfev,
                ftol=1e-11,xtol=1e-11,gtol=1e-11,x_scale='jac')
            err=float(np.max(abs(residual(fit.x))));best=min(best,err)
            rec=dict(iteration=ii,target_pairs=[list(p) for p in finite],target_labels=[labels[k] for k in choice],
                parameters=fit.x.tolist(),maximum_residual=err,nfev=int(fit.nfev))
            records.append(rec)
            if err>1e-7:continue
            R=roots(fit.x);co=-R@M@R.T/2;np.fill_diagonal(co,0)
            if not gate_mask(co).all() or not gate_mask(-(R@M@G['X'].T)/np.sqrt(2*G['norms'])).all():continue
            if np.max(abs(np.einsum('ij,jk,ik->i',R,M,R)-2))>1e-4:continue
            newadj,newoa,_,_=adjacency(R,G)
            if not np.all(newadj[adj]) or np.min(newadj.sum(1)+newoa.sum(1))<3:continue
            ar=-R[:,1:]/R[:,0,None]
            face=poly(base,V,np.r_[AA,ar@V],np.r_[bb,1-ar@base],ids+list(range(158,158+len(R))))
            success.append(dict(**rec,family_ids=[f['id'] for f in ff],roots=R.tolist(),constraints=ar.tolist(),
                minimum_potential_degree=int(np.min(newadj.sum(1)+newoa.sum(1))),face=face))
        print('CLUSTER',ii,'caps',len(ff),'finite_pairs',len(finite),'best_residual',best,'success',len(success),flush=True)
    src=[f'scripts/{TAG}.py',BASE,ACUTE,'scripts/fresh_h20_two_anchor_families_20261016.py',
        'scripts/fresh_h20_quad_adjacent_h3_closure_20261016.py','scripts/fresh_h20_release144_cap_20261016_residual_scout.py']
    return dict(status='PASS_BOUNDED_INCIDENCE_PRESERVING_CLUSTER_RECOVERY',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in src},
        scan=scan,eligible_iterations=eligible,selected_iterations=selected,common_sample_count=len(Z),
        assignment_limit=limit,optimizer_evaluation_limit=nfev,label_domain=labels,records=records,successes=success,
        scope='Numerical necessary intrinsic-neighbor screening, then bounded local angle recovery. The degree condition is necessary only, not a closed-facet certificate. No failure implies all-parameter exclusion; success requires independent exactification and full face closure.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--limit',type=int,default=32);ap.add_argument('--nfev',type=int,default=500)
    a=ap.parse_args();d=run(a.limit,a.nfev);p=ROOT/f'results/{TAG}.json'
    if a.verify:assert close(d,json.loads(p.read_text()))
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['records']),len(d['successes']),flush=True)
