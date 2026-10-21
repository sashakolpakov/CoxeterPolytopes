#!/usr/bin/env python3
"""Bounded numerical nonright-S61 profile and two-anchor experiments.

Old physical pairings are evaluated at100 decimal digits before conversion to
double precision. Exact support identities are imposed explicitly. No numeric
upper bound is promoted to a theorem.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import mpmath as mp
import numpy as np
import cvxpy as cp
import fresh_h20_s61_multicap_20261014_scout as C
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_nonright_s61_20261014_scout'
BASE='results/fresh_h20_s61_multicap_20261014_exact.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def load():
    N,Q,H,dot,w,frame,P,hh,old=C.load();base=json.loads((ROOT/BASE).read_text())
    for p,h in base['sources'].items():assert sha(p)==h,p
    def f(x):q=F(x);return mp.mpf(q.numerator)/q.denominator
    J=mp.matrix([[f(x)for x in row]for row in base['support_inverse']]);S=base['support']
    omega=[sum((J[j,k]*N[S[k]]for k in range(18)),mp.zeros(21,1))for j in range(18)]
    coeff=[mp.matrix([(f(x[0])+f(x[1]))/2 for x in row])for row in base['chain']['frame_coefficient_bounds']]
    N +=[sum((a*v for a,v in zip(co,frame)),mp.zeros(21,1))for co in coeff]
    PP=mp.matrix([[dot(v,n)for v in frame]for n in N])
    WP=mp.matrix([[dot(v,n)for n in N]for v in omega])
    for i in S:
        for j in range(3):assert abs(PP[i,j])<mp.mpf('1e-65');PP[i,j]=0
    for j in range(18):
        for k,i in enumerate(S):assert abs(WP[j,i]-(j==k))<mp.mpf('1e-65');WP[j,i]=int(j==k)
        for i in range(107,133):assert abs(WP[j,i])<mp.mpf('1e-65');WP[j,i]=0
    for i,co in enumerate(coeff,107):
        for j,sg in enumerate([-1,1,1]):assert abs(PP[i,j]-sg*co[j])<mp.mpf('1e-65');PP[i,j]=sg*co[j]
    E=mp.matrix([f(x)for x in base['strict_escape']['rational_point']]);ep=mp.matrix([dot(v,E)for v in frame]);AE=mp.matrix([-dot(v,E)for v in omega])
    return dict(N=N,Q=Q,H=H,dot=dot,frame=frame,P=PP,W=WP,J=J,S=S,E=E,ep=ep,AE=AE,base=base,f=f)
def classify(g,tol=2e-7):
    if abs(g)<tol:return 'right'
    if g>tol:return 'obtuse'
    if g<-2-tol:return 'dotted'
    if abs(g+2)<tol:return 'parallel'
    val=np.pi/np.arccos(-g/2);m=int(round(val))
    if m>=3 and abs(g+2*np.cos(np.pi/m))<tol:return f'pi/{m}'
    return 'forbidden'
def run():
    d=load();P=np.array(d['P'].tolist(),float);W=np.array(d['W'].tolist(),float);J=np.array(d['J'].tolist(),float);ep=np.array(list(d['ep']),float);AE=np.array(list(d['AE']),float);S=d['S']
    # Exact-contact coefficients are1,sqrt2,phi,sqrt3. Every profile is listed.
    labs={3:1.,4:np.sqrt(2.),5:(1+np.sqrt(5.))/2,6:np.sqrt(3.)}
    profiles=[([j],[m])for j in range(18)for m in labs]
    profiles +=[(list(s),[3,3])for s in combinations(range(18),2)]
    y=cp.Variable(3);contact=cp.Parameter(18,nonneg=True);mass=cp.Parameter(nonneg=True)
    scale=np.maximum(np.linalg.norm(P,axis=1),np.linalg.norm(W,axis=0));scale=np.maximum(scale,1)
    acute=(P@y-W.T@contact)/scale<=0
    causal=cp.norm(cp.hstack([y[1:],mass]),2)<=y[0]
    prob=cp.Problem(cp.Maximize(ep@y+AE@contact),[acute,causal])
    rows=[];skipped=[];proposals=[];trials=0;positives=[]
    for pid,(ss,ls)in enumerate(profiles):
        c=np.zeros(18);c[ss]=[labs[x]for x in ls];q=2-c@J@c
        if q>1e-10:skipped.append(dict(id=pid,support=ss,walls=[S[i]for i in ss],labels=ls,residual_norm=q,reason='spacelike_residual_not_in_this_SOC_model'));continue
        contact.value=c;mass.value=np.sqrt(max(0,-q))
        prob.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_gap_rel=1e-9,tol_feas=1e-9,max_iter=300)
        row=dict(id=pid,support=ss,walls=[S[i]for i in ss],labels=ls,residual_norm=float(q),status=prob.status,relaxed_escape_upper_bound=float(prob.value)if np.isfinite(prob.value)else None)
        if y.value is not None:
            yy=y.value;g=P@yy-W.T@c;row.update(residual_frame=yy.tolist(),norm_slack=float(-yy[0]**2+yy[1:]@yy[1:]-q),max_acute_pair=float(max(g)),dual_oldwalls=acute.dual_value.tolist())
        if prob.value is None or not np.isfinite(prob.value)or prob.value<=1e-6:rows.append(row);continue
        positives.append(pid)
        # Only the nearest eight nonsupport old constraints at this optimizer
        # are used as contact anchors. This is a bounded event scout.
        anchors=sorted((i for i in range(133)if i not in S),key=lambda i:abs(g[i]))[:8]
        row['anchor_pool']=anchors
        for ii,jj in combinations(anchors,2):
            L=P[[ii,jj],:];cross=np.cross(L[0],L[1]);cn=np.linalg.norm(cross)
            if cn<1e-10:continue
            v=cross/cn;qa=-v[0]**2+v[1:]@v[1:]
            for mi,mj in product([2,3,4,5,6],repeat=2):
                rhs=(W.T@c)[[ii,jj]]-2*np.array([np.cos(np.pi/mi),np.cos(np.pi/mj)])
                p=L.T@np.linalg.solve(L@L.T,rhs);qb=2*(-p[0]*v[0]+p[1:]@v[1:]);qc=-p[0]**2+p[1:]@p[1:]-q
                if abs(qa)<1e-12:ts=[]if abs(qb)<1e-12 else[-qc/qb]
                else:
                    disc=qb*qb-4*qa*qc
                    ts=[]if disc< -1e-8 else[(-qb+sg*np.sqrt(max(0,disc)))/(2*qa)for sg in[-1,1]]
                for t in ts:
                    trials+=1;z=p+t*v
                    if z[0]<=0:continue
                    cut=ep@z+AE@c
                    if cut<=1e-6:continue
                    gp=P@z-W.T@c
                    if max(gp)>2e-7:continue
                    kinds=[classify(x)for x in gp]
                    if any(x in['obtuse','forbidden']for x in kinds):continue
                    proposals.append(dict(profile=pid,contacts=c.tolist(),anchors=[ii,jj],anchor_labels=[mi,mj],residual_frame=z.tolist(),escape_pair=float(cut),all_pairs=gp.tolist(),relations=kinds))
        rows.append(row)
        print('positive',pid,ss,ls,float(prob.value),'proposals',len(proposals),flush=True)
    paths=[BASE,'scripts/fresh_h20_s61_multicap_20261014_scout.py','scripts/fresh_h20_coupled_roles_20261014.py',f'scripts/{TAG}.py']
    return dict(status='BOUNDED_NUMERICAL_NONRIGHT_S61_SOC_AND_TWO_ANCHOR_SCOUT',sources={p:sha(p)for p in paths},support=S,profiles_total=len(profiles),profiles=rows,skipped=skipped,positive_relaxed_profiles=positives,quadratic_branches_tested=trials,proposals=proposals,scope='Only72 single-node labels3/4/5/6 and153 two-node all3 profiles. Positive causal-residual SOC optima are proposals, not Coxeter certificates. Spacelike residuals are explicitly omitted. At each positive optimum only8 nearest anchors and labels2..6 are tested; no all-profile or all-angle exclusion. Input pairings at100 digits; event solves are numerical and require exactification.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();r=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert r==json.loads(p.read_text())
    else:p.write_text(json.dumps(r,indent=2)+'\n')
    print(r['status'],len(r['positive_relaxed_profiles']),r['quadratic_branches_tested'],len(r['proposals']),flush=True)
