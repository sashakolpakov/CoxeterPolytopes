#!/usr/bin/env python3
"""Bounded rotation-only numerical test of one exact-contact candidate direction."""
import argparse,json,time
from pathlib import Path
from hashlib import sha256
from itertools import combinations
import numpy as np
import mpmath as mp
from scipy.linalg import expm
from scipy.optimize import least_squares
from fresh_h9_h8_facet_cover_20261014_geometry import data,SOURCE,BASE
import fresh_h9_mixed_d4_coupling_20261014_review as A
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_h8_corner_anchors_20261014_single_scout'
LITERAL_V=['3.346178','4.2047','3.038572','1.952451','.909731','.150188','-1.149919','-.592165','-.912382','-.750588']
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 mp.mp.dps=90;D=data(90);b=json.loads((ROOT/BASE).read_text());J=np.diag([-1.]+[1.]*8)
 def f(x):return mp.mpf(x.numerator)/x.denominator
 def ev(x):return f(x.a.a)+f(x.a.b)*mp.sqrt(5)+(f(x.b.a)+f(x.b.b)*mp.sqrt(5))*mp.sqrt((16+4*mp.sqrt(5))/11)
 Q=mp.matrix([[ev(x)for x in row]for row in A.dm(b['form'])]);N=mp.matrix([[ev(x)for x in row]for row in A.dm(b['final_normals'])]);H=mp.matrix([ev(x)for x in A.dv(b['center'])]);V=mp.matrix([mp.mpf(x)for x in LITERAL_V]);delta=V-H;old8=N[8,:].T;a=(mp.sqrt(5)-1)/4
 def ip(x,y):return(x.T*Q*y)[0]
 g0,g1=ip(H,old8),ip(delta,old8);aa=g1*g1-(2+a)*ip(delta,delta);bb=2*g0*g1-2*(2+a)*ip(H,delta);cc=g0*g0-(2+a)*ip(H,H)
 roots=[(-bb+sg*mp.sqrt(bb*bb-4*aa*cc))/(2*aa)for sg in[-1,1]];t=min(roots,key=lambda x:abs(x-1));r=H+t*delta;q=ip(r,r);g=N*Q*r;G=N*Q*N.T;PP=G-g*g.T/q;k=[mp.sqrt(-PP[i,i])for i in range(12)];assert abs(k[8]**2-a)<mp.mpf('1e-65')
 CG=np.array([[-float(PP[i,j]/(k[i]*k[j]))for j in range(12)]for i in range(12)])
 # CG stores cosh distances; factor its negative Lorentz Gram.
 vals,U=np.linalg.eigh(-CG);ids=[int(np.argmin(vals))]+list(np.argsort(vals)[-8:]);C=U[:,ids]*np.sqrt(abs(vals[ids]));C[:,0]*=np.sign(sum(C[:,0]))
 def boost(p):
  v=np.sqrt(1+p@p);M=np.eye(9);M[0,0]=v;M[0,1:]=p;M[1:,0]=p;M[1:,1:]+=np.outer(p,p)/(1+v);return M
 C=C@boost(-C[8,1:]).T;assert np.max(abs(C[8]-np.r_[1,np.zeros(8)]))<1e-10
 rho=np.array([float(mp.asinh(mp.sqrt(G[i,i])/k[i]))for i in range(12)]);cr=np.cosh(rho)
 target=np.sum([v[1:]/np.linalg.norm(v[1:])for i,v in enumerate(C)if i!=8],axis=0);target/=np.linalg.norm(target)
 pairs=list(combinations(range(8),2));rows=[]
 def rotate(v,z):
  v=v/np.linalg.norm(v);z=z/np.linalg.norm(z);K=np.outer(z,v)-np.outer(v,z);return np.eye(8)+K+K@K/(1+v@z)
 for vi in[0,32,33]:
  Qs=D['Qs'];origin=D['V'][vi];basis=[origin]
  for j in range(9):
   v=mp.eye(9)[:,j];v+=origin*((v.T*Qs*origin)[0])
   for w in basis[1:]:v-=w*((v.T*Qs*w)[0])
   n=(v.T*Qs*v)[0]
   if n>mp.mpf('1e-55'):basis.append(v/mp.sqrt(n))
   if len(basis)==9:break
  B=mp.matrix([[v[i]for v in basis]for i in range(9)]);Ns=np.array(D['Ns']*Qs*B,dtype=float).reshape(11,9)@J;Vs=np.array(mp.matrix([[v[i]for i in range(9)]for v in D['V']])*Qs*B,dtype=float).reshape(41,9)@J
  support=D['source']['geometry']['vertices'][vi];du=-np.linalg.inv(Ns[support,1:])@np.ones(8);du/=np.linalg.norm(du);R0=rotate(du,target)
  def mat(par):
   K=np.zeros((8,8))
   for x,(i,j)in zip(par,pairs):K[i,j]=x;K[j,i]=-x
   M=np.eye(9);M[1:,1:]=expm(K)@R0;return M
  def residual(par):
   M=mat(par);NN=Ns@M.T;VV=Vs@M.T;ct=C@J@NN.T;bv=-VV@J@C.T/cr[None,:]-1
   return np.r_[np.maximum(ct,0).ravel(),np.maximum(np.min(bv,axis=1),0)]
  for seed in range(8):
   par=np.zeros(28)if seed==0 else np.random.default_rng(1000+vi*10+seed).normal(0,.35,28)
   out=least_squares(residual,par,max_nfev=400,ftol=1e-10,gtol=1e-10,xtol=1e-10);M=mat(out.x);NN=Ns@M.T;VV=Vs@M.T;ct=C@J@NN.T;bv=-VV@J@C.T/cr[None,:]-1
   row=dict(source_vertex=vi,start=seed,status=int(out.status),nfev=int(out.nfev),cost=float(residual(out.x)@residual(out.x)),max_center_violation=float(ct.max()),max_vertex_uncovered=float(np.min(bv,axis=1).max()),matrix=M.tolist(),vertex_ball_ids=np.argmin(bv,axis=1).tolist())
   rows.append(row);print(vi,seed,row['cost'],row['max_center_violation'],row['max_vertex_uncovered'],flush=True)
 paths=[SOURCE,BASE,'scripts/fresh_h9_h8_facet_cover_20261014_geometry.py','scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',f'scripts/{TAG}.py']
 return dict(status='BOUNDED_SINGLE_CORNER_ANCHOR_PLACEMENT_SCOUT',sources={p:sha(p)for p in paths},literal_V=LITERAL_V,t=str(t),raw_cap_norm=str(q),raw_cap=[str(x)for x in r],
  anchor_old_wall=8,anchor_normalized_projection_norm=str(k[8]**2),centers=C.tolist(),ball_radii=rho.tolist(),rows=rows,best=min(rows,key=lambda x:x['cost']),
  scope='One candidate direction,3 exceptional sourcecorner choices,8 deterministic local rotations each,400function evaluations. Anchored wall8 contacts are exact at the defining quadratic level; no other pair identities or full-base coverage are certified. No feasible relaxation found is only bounded numerical failure.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:
  o=json.loads(p.read_text());assert d['sources']==o['sources'] and d['t']==o['t'];assert len(d['rows'])==len(o['rows'])==24
  for x,y in zip(d['rows'],o['rows']):
   assert(x['source_vertex'],x['start'],x['status'])==(y['source_vertex'],y['start'],y['status'])
   for k in['cost','max_center_violation','max_vertex_uncovered']:assert abs(x[k]-y[k])<1e-7
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print('PASS best',d['best']['cost'],flush=True)
