#!/usr/bin/env python3
"""Discover and exactly replay ball/linear dual certificates at the null link.

Optimization supplies candidate multipliers only. Every accepted certificate is
checked with rational outward intervals and an exact height-balance identity.
"""
import argparse,json,math
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from fresh_h20_s61_multicap_20261014_exact import I
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_null_link_escape_gate_20261016'
GEO=f'results/{TAG}_geometry.json';PRO='results/fresh_h20_null_link_profiles_20261016.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def mid(x):return float((x.l+x.h)/2)
def enc(x):
    den=1<<160
    return[str(F((x.l*den).__floor__(),den)),str(F((x.h*den).__ceil__(),den))]
def square(x):return I(0 if x.l<=0<=x.h else min(x.l*x.l,x.h*x.h),max(x.l*x.l,x.h*x.h))
def norm(v):
    q=sum((square(x)for x in v),I(0))
    return I(0)if q.h==0 else q.sqrt()
def inverse(A):
    n=len(A);B=[row[:]+[F(i==j)for j in range(n)]for i,row in enumerate(A)]
    for i in range(n):
        k=next(k for k in range(i,n)if B[k][i]);B[i],B[k]=B[k],B[i]
        t=B[i][i];B[i]=[x/t for x in B[i]]
        for j in range(n):
            if j!=i:
                t=B[j][i];B[j]=[x-t*y for x,y in zip(B[j],B[i])]
    return [row[n:]for row in B]
def problem(geo,profile,mask):
    white=geo['white'];J=inverse([[F(x)for x in row]for row in geo['positive_Gram']])
    tail='wall'in profile
    if tail:coeff=[I(0)]*13;R=I(2).sqrt()
    else:
        mags={3:I(1),4:I(2).sqrt(),5:(I(1)+I(5).sqrt())/2,6:I(3).sqrt()}
        d=[I(0)]*13
        for i,m in zip(profile['walls'],profile['labels']):d[white.index(i)]=mags[m]
        coeff=[-sum((d[i]*J[i][j]for i in range(13)),I(0))for j in range(13)]
        cost=sum((I(F(c))*I(rad).sqrt()for c,rad in zip(profile['cost'],[1,2,3,5,6,10,15,30])),I(0))
        if profile['cost']==['2']+['0']*7:R=I(0)
        else:R=(2-cost).sqrt()
    def row(v):
        a=[I(*x)for x in v['free']]
        if tail:a.append(-I(*v['white'][white.index(profile['wall'])])/I(2).sqrt())
        a.append(I(*v['null']))
        p=sum((c*I(*x)for c,x in zip(coeff,v['white'])),I(0))
        return a,p
    rows=[];bounds=[];labels=[]
    for i,v in enumerate(geo['normals']):
        if i in geo['incident']:continue
        a,p=row(v);rows.append(a);bounds.append(-p);labels.append(['acute',i])
    for j,i in enumerate([2,5,61,87]):
        a,p=row(geo['normals'][i])
        if mask&(1<<j):rows.append(a);bounds.append(-1-p);labels.append(['nonright',i])
        else:rows.append([-x for x in a]);bounds.append(p);labels.append(['right_reverse',i])
    if tail:
        n=len(rows[0]);v=[I(0)]*n;v[-2]=I(1)
        rows.extend([v,[-x for x in v]]);bounds.extend([I(2).sqrt(),-I(F(3,2)).sqrt()]);labels.extend([['tail_upper'],['tail_lower']])
    c,e=row(geo['escape'])
    # Positive rational row scaling changes neither the feasible set nor proof.
    scales=[]
    for a,b in zip(rows,bounds):
        largest=max([1.]+[abs(mid(x))for x in a]+[abs(mid(b))]);scales.append(F(1,2**max(0,math.ceil(math.log2(largest)))))
    rows=[[x*s for x in a]for a,s in zip(rows,scales)];bounds=[b*s for b,s in zip(bounds,scales)]
    return rows,bounds,c,e,R,labels
def certificate(data,witness):
    A,b,c,e,R,labels=data;n=len(A)
    mu=[I(0)]*n
    for i,q in witness['rational_multipliers']:mu[i]=I(F(q));assert mu[i].l>=0
    pivot=witness['height_balance_pivot'];target=c[-1]if witness['kind']=='upper_bound'else I(0)
    mu[pivot]=(target-sum((mu[i]*A[i][-1]for i in range(n)if i!=pivot),I(0)))/A[pivot][-1]
    if mu[pivot].l<0:return None
    targetz=c[:-1]if witness['kind']=='upper_bound'else [I(0)]*(len(c)-1)
    residual=[x-sum((mu[i]*A[i][j]for i in range(n)),I(0))for j,x in enumerate(targetz)]
    bound=sum((mu[i]*b[i]for i in range(n)),I(0))+R*norm(residual)
    if witness['kind']=='upper_bound':bound+=e
    if bound.h>=0:return None
    return dict(bound=enc(bound),height_pivot_multiplier=enc(mu[pivot]),
                used_rows=[labels[i]for i,x in enumerate(mu)if x.h>0])
def discover(data):
    import cvxpy as cp
    import numpy as np
    A,b,c,e,R,labels=data;AA=np.array([[mid(x)for x in row]for row in A]);bb=np.array([mid(x)for x in b]);cc=np.array([mid(x)for x in c])
    v=cp.Variable(len(c));ineq=AA@v<=bb
    ball=cp.norm(v[:-1])<=mid(R)
    prob=cp.Problem(cp.Maximize(cc@v+mid(e)),[ineq,ball])
    options=dict(tol_gap_abs=1e-10,tol_gap_rel=1e-10,tol_feas=1e-10,max_iter=300)
    try:prob.solve(solver='CLARABEL',**options)
    except cp.error.SolverError:pass
    if prob.status in ['optimal','optimal_inaccurate'] and prob.value is not None and prob.value<-1e-7:
        kinds=[('upper_bound',ineq.dual_value)]
    else:kinds=[]
    slack=cp.Variable();phase_row=AA@v-slack<=bb
    phase=cp.Problem(cp.Minimize(slack),[phase_row,ball])
    if not kinds:
        try:phase.solve(solver='CLARABEL',**options)
        except cp.error.SolverError:pass
        if phase.status in ['optimal','optimal_inaccurate'] and phase.value is not None and phase.value>1e-8:
            kinds.append(('infeasible',phase_row.dual_value))
    for kind,raw in kinds:
        for digits in [8,11,14]:
            den=10**digits;q=[max(F(0),F(round(float(x)*den),den))for x in raw]
            candidates=sorted((i for i,x in enumerate(q)if x>0 and abs(mid(A[i][-1]))>1e-9),key=lambda i:-float(q[i])*abs(mid(A[i][-1])))
            for pivot in candidates[:12]:
                w=dict(kind=kind,height_balance_pivot=pivot,rational_multipliers=[[i,str(x)]for i,x in enumerate(q)if x and i!=pivot])
                valid=certificate(data,w)
                if valid is not None:return dict(**w,verified=valid)
    return dict(kind='UNRESOLVED',optimizer_status=prob.status,optimizer_value=None if prob.value is None or not math.isfinite(prob.value)else float(prob.value),phase_status=phase.status,phase_value=None if phase.value is None or not math.isfinite(phase.value)else float(phase.value))
def run(verify=False):
    geo=json.loads((ROOT/GEO).read_text());pro=json.loads((ROOT/PRO).read_text())
    for d in [geo,pro]:
        for p,h in d['sources'].items():assert sha(p)==h,p
    profiles=pro['profiles']+pro['continuous_tails'];records=[];unresolved=[]
    old=json.loads((ROOT/f'results/{TAG}.json').read_text())if verify else None
    for k,profile in enumerate(profiles):
        for mask in range(16):
            if k==mask==0:
                w=dict(kind='zero_line_exclusion',dependency=GEO)
            else:
                data=problem(geo,profile,mask)
                if verify:
                    w=old['records'][16*k+mask]['certificate']
                    if w['kind']!='UNRESOLVED':assert certificate(data,w)==w['verified']
                else:w=discover(data)
            if w['kind']=='UNRESOLVED':unresolved.append([k,mask])
            records.append(dict(profile=k,mask=mask,certificate=w))
        print('verified profile',k,'unresolved',len(unresolved),flush=True)
    src=[GEO,PRO,f'scripts/{TAG}.py','scripts/fresh_h20_s61_multicap_20261014_exact.py']
    out=dict(status='PASS_EXACT_ALL_THROUGH_NULL_COXETER_ESCAPE_GATE'if not unresolved else 'PARTIAL_EXACT_NULL_LINK_DUAL_CERTIFICATES',
             sources={p:sha(p)for p in src},profile_domains=len(profiles),cases=len(records),records=records,unresolved=unresolved,
             retained_walls=135,extra_center_constraint=False,
             scope='All-through-nu acute norm2 Coxeter E1-cutter exclusion only if every case certified. Norm-ball relaxation and minimum-contact gates at2/5/61/87; all other nonincident contacts merely acute. No finite label cutoff, no old-center preservation hypothesis.')
    if verify:assert out==old
    return out
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run(args.verify)
    if not args.verify:(ROOT/f'results/{TAG}.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'unresolved',len(out['unresolved']),flush=True)
