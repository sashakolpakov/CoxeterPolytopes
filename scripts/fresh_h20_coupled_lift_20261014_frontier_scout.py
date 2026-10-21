import json,numpy as np
from pathlib import Path
from hashlib import sha256
ROOT=Path(__file__).resolve().parents[1]
import os
os.chdir(ROOT)
from fractions import Fraction as F
D=json.load(open('results/fresh_dimension_construction_routes_20261014.json'));sc=json.load(open('results/fresh_h20_coupled_lift_20261014_scout.json'));f=lambda t:float(F(t['a']))+float(F(t['b']))*2**.5
Q=np.array([[f(x)for x in r]for r in D['form']]);N=np.array([[f(x)for x in r]for r in D['normals']]);h=np.array([f(x)for x in D['strict_center']]);v=np.array([f(x)for x in D['strict_escape']['boundary_point']]);E=np.array([f(x)for x in D['strict_escape']['point']]);G=N@Q@N.T
labels=list(range(2,31));cs=np.cos(np.pi/np.array(labels));hits=[];trials=0
for si,rec in enumerate(sc['recipes']):
 S=rec['support'];B=G[np.ix_(S,S)];h0=h-N[S].T@np.linalg.solve(B,N[S]@Q@h)
 _,_,V=np.linalg.svd(np.vstack([N[S]@Q,v@Q,h0@Q]));w=V[-1];w/=np.sqrt(w@Q@w)
 qa=h0@Q@h0;qb=2*v@Q@h0;qc=v@Q@v;t0=(-qb-np.sqrt(qb*qb-4*qa*qc))/(2*qa);p=N@Q@v;a=N@Q@h0;b=N@Q@w;others=[i for i in range(60)if i not in S]
 for j in others:
  for sign in[-1,1]:
   for label,cos in [(0,1)]+list(zip(labels,cs)):
    z=-cos*np.sqrt(2*G[j,j])-sign*b[j]
    poly=[a[j]**2-z*z*qa,2*a[j]*p[j]-z*z*qb,p[j]**2-z*z*qc]
    for tr in np.roots(poly):
     if abs(tr.imag)>1e-8:continue
     t=float(tr.real)
     if not 0<t<t0*(1-1e-10):continue
     q=qa*t*t+qb*t+qc
     if q<=0 or abs((p[j]+t*a[j])/np.sqrt(q)-z)>1e-6:continue
     trials+=1;m=v+t*h0;caps=np.array([m+np.sqrt(q)*w,m-np.sqrt(q)*w]);rat=caps@Q@N.T/np.sqrt(2*q*np.diag(G));vals=rat[:,others]
     if np.any(vals>1e-7):continue
     finite=vals[(vals>-1+1e-7)&(vals<-1e-7)]
     if len(finite) and max(np.min(abs(finite[:,None]+cs[None,:]),axis=1))>1e-6:continue
     if np.min(caps@Q@E)<=0:continue
     hits.append(dict(support=si,t=t,anchor=j,sign=sign,label=int(label),q=float(q),finite_angles=finite.tolist(),caps=caps.tolist()))
print('trials',trials,'hits',len(hits));print([{k:v for k,v in x.items()if k!='caps'}for x in hits[:20]])
paths=['results/fresh_dimension_construction_routes_20261014.json','results/fresh_h20_coupled_lift_20261014_scout.json','scripts/fresh_h20_coupled_lift_20261014_frontier_scout.py']
json.dump(dict(status='NUMERICAL_BOUNDED_PAIR_FRONTIER',sources={p:sha256(Path(p).read_bytes()).hexdigest()for p in paths},labels=labels,trials=trials,hits=hits,scope='Numerical finite-anchor events for four supports using the original irrational source center; labels2..30 plus parallel only. No exact new cap or closure claim. This is a separate proposal branch from the frozen rational-center92.'),open('results/fresh_h20_coupled_lift_20261014_frontier_scout.json','w'),indent=2)
