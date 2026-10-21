#!/usr/bin/env python3
"""Numerical continuous all-right PD7 face-cap feasibility, no cap census claim."""
from pathlib import Path
import json,time,itertools,argparse
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog,minimize
from fresh_h9_makarov_joint_completion_20261013_profiles import setup,ROOT,P,CORE
OUT=ROOT/f'results/{P}_faces.json'
def cliques(adj,n):
    out=[]
    def visit(v,rest):
        if len(v)==n:out.append(v);return
        if len(v)+len(rest)<n:return
        for i in sorted(rest):visit(v+[i],{j for j in rest&adj[i]if j>i})
    visit([],set(range(len(adj))));return out
def build(latest=False):
    start=time.time();src,Q,W,c,E=setup();target=E
    if latest:
        target=np.array(json.loads((ROOT/f'results/{P}_all_relaxed_volume.json').read_text())['best_point'])
        target=.999*target+.001*c/(-c@Q@c)
    tags={tuple(row['pair']):row['type']for row in src['all_pair_gates']}
    adj=[set()for _ in W]
    for (i,j),tag in tags.items():
        if tag!='dotted':adj[i].add(j);adj[j].add(i)
    supports=cliques(adj,7);stats={'cliques':len(supports),'pd':0,'acute_can_cut':0,'dotted_can_cut':0};rows=[]
    scores=[]
    for supp in supports:
        G=W[supp]@Q@W[supp].T
        if min(np.linalg.eigvalsh(G))<1e-8:continue
        stats['pd']+=1;V=null_space(W[supp]@Q);H=V.T@Q@V
        outside=sorted(set(range(27))-set(supp));A=W[outside]@Q@V;cc=c@Q@V;ee=target@Q@V
        res=linprog(-ee,A_ub=A,b_ub=np.zeros(len(A)),A_eq=[cc],b_eq=[-1],bounds=[(None,None)]*3,method='highs')
        if res.status==3:
            res=linprog(np.zeros(3),A_ub=np.vstack([A,-ee]),b_ub=np.r_[np.zeros(len(A)),-1],A_eq=[cc],b_eq=[-1],bounds=[(None,None)]*3,method='highs')
        if not res.success or ee@res.x<1e-9:continue
        stats['acute_can_cut']+=1;scores.append((float(ee@res.x),supp,V,H,A,cc,ee,res.x))
    print('LP',stats,'seconds',time.time()-start,flush=True)
    # Most promising supports first; this is a bounded numerical optimization.
    for score,supp,V,H,A,cc,ee,x in sorted(scores,key=lambda r:-r[0])[:150]:
        starts=[x]
        for i in range(16):
            direction=np.random.default_rng(74100+i).normal(size=3)
            lp=linprog(-direction,A_ub=A,b_ub=np.zeros(len(A)),A_eq=[cc],b_eq=[-1],bounds=[(None,None)]*3,method='highs')
            if lp.success:starts.append(lp.x)
        best=None
        for xx in starts:
            q=xx@H@xx
            if q<=1e-10:continue
            xx=xx*np.sqrt(2/q)
            sol=minimize(lambda z:-ee@z,xx,jac=lambda z:-ee,method='SLSQP',constraints=[
                {'type':'eq','fun':lambda z:z@H@z-2,'jac':lambda z:2*H@z},
                {'type':'ineq','fun':lambda z:-2-A@z,'jac':lambda z:-A},
                {'type':'ineq','fun':lambda z:-1e-8-cc@z,'jac':lambda z:-cc}],options={'maxiter':250,'ftol':1e-11})
            z=sol.x
            if abs(z@H@z-2)>1e-5 or max(A@z)>-2+1e-6 or cc@z>=0:continue
            if best is None or ee@z>best[0]:best=(float(ee@z),z,sol.success)
        if best is not None:
            cut,z,success=best;r=V@z
            rows.append({'support':supp,'normal':r.tolist(),'target_pairing':cut,'all27_pairings':(W@Q@r).tolist(),'norm':float(r@Q@r),'optimizer_success':bool(success)})
            if cut>1e-8:stats['dotted_can_cut']+=1;print('CUT',supp,cut,flush=True)
    out={'status':'NUMERICAL_PD7_RIGHT_SUPPORT_FEASIBILITY','counts':stats,'optimized_support_limit':150,'caps':rows,'target':target.tolist(),'seconds':time.time()-start,'scope':'All finite-pair seven-cliques receive numerical PD and acute-cone LP gates. At most150 promising supports receive bounded nonlinear dotted feasibility optimization. Neither optimization failure nor this all-right profile is an exhaustive cap exclusion.'}
    outpath=OUT.with_name(OUT.stem+'_latest.json')if latest else OUT
    outpath.write_text(json.dumps(out,indent=2)+'\n');print('DONE',stats,len(rows),time.time()-start,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--latest',action='store_true');a=p.parse_args();build(a.latest)
