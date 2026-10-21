#!/usr/bin/env python3
"""Bounded numerical center-hull coverage diagnostic for one exact dotted cap.

No global optimizer or all-source exclusion is claimed. The exact geometric
coverage criterion is in the main report; this file tests a few points only.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import mpmath as mp
import numpy as np
from scipy.optimize import minimize
import fresh_h9_mixed_d4_coupling_20261014_review as A

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_mixed_d4_global_caps_20261014_coverage_scout'
BASE='results/fresh_h9_mixed_d4_coupling_20261014.json'
CAP='results/fresh_h9_mixed_d4_global_caps_20261014.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    b,d=[json.loads((ROOT/p).read_text())for p in[BASE,CAP]]
    for data in[b,d]:
        for p,h in data['sources'].items():assert sha(p)==h,p
    Q=A.dm(b['form']);N=A.dm(b['final_normals']);r=A.dv(d['raw_cap']);q=A.dec(d['raw_cap_norm'])
    gg=[A.pair(r,n,Q)for n in N]
    mp.mp.dps=80
    def f(x):return mp.mpf(x.numerator)/x.denominator
    def ev(x):return f(x.a.a)+f(x.a.b)*mp.sqrt(5)+(f(x.b.a)+f(x.b.b)*mp.sqrt(5))*mp.sqrt((16+4*mp.sqrt(5))/11)
    g=mp.matrix([ev(x)for x in gg]);qr=ev(q)
    G=mp.matrix([[ev(A.pair(n,m,Q))for m in N]for n in N])
    PP=G-g*g.T/qr;k=[mp.sqrt(-PP[i,i])for i in range(12)]
    # For X=sum lambda_i*(projected_n_i/k_i), an outside-all-balls witness
    # requires all lambda'C_j lambda>0. This is the exact inequality
    # q*(X.n_j)^2+g_j^2*X^2>0, divided by positive g_j^2*q.
    C=np.array([[[float(PP[i,j]/(k[i]*k[j]*qr)+PP[i,l]*PP[j,l]/(g[l]**2*k[i]*k[j]))for j in range(12)]for i in range(12)]for l in range(12)])
    def scores(w):return np.einsum('i,kij,j->k',w,C,w)
    samples=[]
    for size in[2,3,12]:
        for ids in combinations(range(12),size):
            w=np.zeros(12);w[list(ids)]=1/size
            v=scores(w);samples.append(dict(support=list(ids),minimum_scaled_ball_margin=float(min(v)),blocking_ball=int(np.argmin(v))))
    assert len(samples)==287
    opts=[]
    def inequalities(z):return scores(z[:-1])-z[-1]
    for i in range(13):
        w=np.ones(12)/12 if i==12 else np.eye(12)[i]
        start=np.r_[w,min(scores(w))]
        out=minimize(lambda z:-z[-1],start,method='SLSQP',bounds=[(0,1)]*12+[(None,None)],
            constraints=[{'type':'eq','fun':lambda z:sum(z[:-1])-1},{'type':'ineq','fun':inequalities}],
            options={'maxiter':300,'ftol':1e-12})
        v=scores(out.x[:-1])
        opts.append(dict(start=i,success=bool(out.success),message=str(out.message),iterations=int(out.nit),
            weights=out.x[:-1].tolist(),minimum_scaled_ball_margin=float(min(v))))
    maximum=max(r['minimum_scaled_ball_margin']for r in samples+opts)
    assert maximum<0
    paths=[BASE,CAP,'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',f'scripts/{TAG}.py']
    return dict(status='BOUNDED_NUMERICAL_CENTER_HULL_DIAGNOSTIC_NO_WITNESS',sources={p:sha(p)for p in paths},
        coefficient_precision_decimal=80,pair_sums=66,triple_sums=220,full_barycenters=1,
        samples=samples,local_optimizer_starts=13,local_optimizer_rows=opts,
        maximum_observed_minimum_margin=maximum,
        scope='Only287 fixed convex combinations and13 local optimization starts. Failure to find a point outside all balls is not a coverage proof or all-source obstruction. Local optimizer failures are retained explicitly.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:
        old=json.loads(p.read_text());assert d['sources']==old['sources']
        assert len(d['samples'])==len(old['samples'])==287 and len(d['local_optimizer_rows'])==13
        assert d['maximum_observed_minimum_margin']<0 and old['maximum_observed_minimum_margin']<0
        assert [x['support']for x in d['samples']]==[x['support']for x in old['samples']]
        assert max(abs(x['minimum_scaled_ball_margin']-y['minimum_scaled_ball_margin'])for x,y in zip(d['samples'],old['samples']))<1e-10
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],'best observed margin',d['maximum_observed_minimum_margin'],flush=True)
