#!/usr/bin/env python3
"""Small continuous placement relaxation; no Coxeter or coverage acceptance."""
import argparse,json,time
from pathlib import Path
from hashlib import sha256
import numpy as np
import mpmath as mp
from scipy.linalg import expm
from scipy.optimize import least_squares
from itertools import combinations
from fresh_h9_h8_facet_cover_20261014_geometry import data,SOURCE,BASE,CAP
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_h8_facet_cover_20261014_placement_scout'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def one(u):
 D=data(u=u);J=np.diag([-1.]+[1.]*8)
 # Source center at midpoint of diametral vertices; exact mp Gram--Schmidt.
 import mpmath as mp
 Q=D['Qs'];p,q=D['V'][32:34];h=(p+q)/mp.sqrt(-((p+q).T*Q*(p+q))[0]);basis=[h]
 for j in range(9):
  v=mp.eye(9)[:,j]
  v+=h*((v.T*Q*h)[0])
  for b in basis[1:]:v-=b*((v.T*Q*b)[0])
  n=(v.T*Q*v)[0]
  if n>mp.mpf('1e-50'):basis.append(v/mp.sqrt(n))
  assert len(basis)<=9
  if len(basis)==9:break
 B=mp.matrix([[v[i]for v in basis]for i in range(9)])
 Ns=np.array(D['Ns']*Q*B,dtype=float).reshape(11,9)@J
 Vs=np.array(mp.matrix([[v[i]for i in range(9)]for v in D['V']])*Q*B,dtype=float).reshape(41,9)@J
 CG=np.array([[-float(mp.cosh(D['DC'][i,j]))for j in range(12)]for i in range(12)])
 w,U=np.linalg.eigh(CG);ids=[int(np.argmin(w))]+list(np.argsort(w)[-8:]);C=U[:,ids]*np.sqrt(abs(w[ids]));C[:,0]*=np.sign(np.sum(C[:,0]));
 assert np.max(abs(C@J@C.T-CG))<1e-9
 # Move center hull barycenter to target origin.
 def boost(p):
  a=np.sqrt(1+np.dot(p,p));ret=np.eye(9);ret[0,0]=a;ret[0,1:]=p;ret[1:,0]=p;ret[1:,1:]+=np.outer(p,p)/(1+a);return ret
 c=np.sum(C,axis=0);c/=np.sqrt(-c@J@c);C=C@boost(-c[1:]).T
 rho=np.array([float(x)for x in D['rho']]);cr=np.cosh(rho)
 pairs=list(combinations(range(8),2))
 def rotate(v,z):
  # minimal proper rotation mapping v to z
  v=v/np.linalg.norm(v);z=z/np.linalg.norm(z);k=np.outer(z,v)-np.outer(v,z);dot=np.dot(v,z)
  if dot<-.999999:raise ValueError()
  return np.eye(8)+k+k@k/(1+dot)
 def matrix(par,base):
  K=np.zeros((8,8));
  for x,(i,j)in zip(par[8:],pairs):K[i,j]=x;K[j,i]=-x
  R=np.eye(9);R[1:,1:]=expm(K);return boost(par[:8])@R@base
 # squared/hyperbolic distances monotonic; violated normalized cosh coverage plus acute centers.
 def residual(par,base):
  M=matrix(par,base);N=Ns@M.T;V=Vs@M.T
  ct=C@J@N.T
  vv=-V@J@C.T/cr[None,:]-1
  return np.r_[np.maximum(ct,0).ravel(),np.maximum(np.min(vv,axis=1),0)]
 def diagnostics(par,base):
  M=matrix(par,base);N=Ns@M.T;V=Vs@M.T
  ct=C@J@N.T;vv=-V@J@C.T/cr[None,:]-1
  return dict(max_center_violation=float(ct.max()),max_vertex_uncovered=float(np.min(vv,axis=1).max()),cost=float(np.sum(residual(par,base)**2)),matrix=M.tolist(),vertex_ball_ids=np.argmin(vv,axis=1).tolist(),uncovered=[int(i)for i in np.where(np.min(vv,axis=1)>1e-7)[0]])
 # Start source diameter direction aligned with target center7 versus centerj;
 # shift midpoint along that geodesic to balance available cap radii.
 axis=(Vs[32]-Vs[33])[1:];axis/=np.linalg.norm(axis)
 rows=[];best=1e100;t0=time.time()
 for j in range(4):
  v,z=C[7],C[j];mid=v+z;mid/=np.sqrt(-mid@J@mid);M0=boost(mid[1:]);loc=(v-z)@boost(-mid[1:]).T
  rot=rotate(axis,loc[1:]);R=np.eye(9);R[1:,1:]=rot;base=M0@R
  # adjust source center toward larger ball7
  shift=(rho[7]-rho[j])/2;direction=loc[1:]/np.linalg.norm(loc[1:]);base=M0@boost(-np.sinh(shift)*direction)@R
  for seed in range(2):
   par=np.zeros(36)
   if seed:par[8:]=np.random.default_rng(100+j).normal(0,.08,28)
   out=least_squares(residual,par,args=(base,),max_nfev=250,ftol=1e-9,gtol=1e-9,xtol=1e-9,verbose=0)
   row=dict(ball=j,seed=seed,nfev=int(out.nfev),status=int(out.status),message=out.message,**diagnostics(out.x,base));rows.append(row)
   print('ROW',j,seed,out.nfev,row['cost'],row['max_center_violation'],row['max_vertex_uncovered'],len(row['uncovered']),'sec',time.time()-t0,flush=True)
   if row['cost']<best:best=row['cost']
   if best<1e-15:break
  if best<1e-15:break
 return dict(u=str(D['u']),t=str(D['t']),rows=rows,best=min(rows,key=lambda x:x['cost']),source_normals=Ns.tolist(),source_vertices=Vs.tolist(),centers=C.tolist(),radii=rho.tolist())

def run():
 cases=[one(u)for u in[None,'0.85','0.9','0.95','0.98']]
 paths=[SOURCE,BASE,CAP,'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py','scripts/fresh_h9_h8_facet_cover_20261014_geometry.py',f'scripts/{TAG}.py']
 return dict(status='BOUNDED_NUMERICAL_PLACEMENT_RELAXATION',sources={p:sha(p)for p in paths},cases=cases,
  scope='Five cap parameters, eight starts each,250 function evaluations perstart. All132 center/facet acuteness inequalities and41 vertex/ball-union membership inequalities are relaxed continuously; no discrete Coxeter-angle constraints. No feasible local solution is not a nonexistence proof. Even a feasible vertex placement would not certify whole-polytope coverage.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:
  old=json.loads(p.read_text());assert d['sources']==old['sources'] and len(d['cases'])==len(old['cases'])==5
  for x,y in zip(d['cases'],old['cases']):
   assert x['u']==y['u'] and len(x['rows'])==len(y['rows'])
   for a,b in zip(x['rows'],y['rows']):
    assert (a['ball'],a['seed'],a['status'])==(b['ball'],b['seed'],b['status'])
    for k in['cost','max_center_violation','max_vertex_uncovered']:assert abs(a[k]-b[k])<1e-7
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print('PASS',[(r['u'][:8],r['best']['cost'])for r in d['cases']],flush=True)
