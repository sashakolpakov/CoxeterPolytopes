#!/usr/bin/env python3
"""High-precision numerical geometry for a fixed H8-source attachment pilot."""
import json
from pathlib import Path
import mpmath as mp
import numpy as np
import fresh_h9_mixed_d4_coupling_20261014_review as A
ROOT=Path(__file__).resolve().parents[1]
SOURCE='data/known_examples/bugaenko_compact_h8.json'
BASE='results/fresh_h9_mixed_d4_coupling_20261014.json'
CAP='results/fresh_h9_mixed_d4_global_caps_20261014.json'
def data(dps=80,u=None):
 mp.mp.dps=dps
 src,b,d=[json.loads((ROOT/p).read_text())for p in[SOURCE,BASE,CAP]]
 assert src['dimension']==8 and len(src['simple_roots'])==11 and len(src['geometry']['vertices'])==41
 phi=(1+mp.sqrt(5))/2
 def sf(x):return mp.mpf(x['a'])+mp.mpf(x['b'])*phi
 Qs=mp.matrix([[sf(x)for x in row]for row in src['lattice_gram']]);Ns=mp.matrix([[sf(x)for x in row]for row in src['simple_roots']]);H=mp.matrix([sf(x)for x in src['geometry']['interior_point']])
 def pair(x,y,Q):return (x.T*Q*y)[0]
 GS=Ns*Qs*Ns.T
 assert max(abs(GS[i,j]-sf(src['wall_gram'][i][j]))for i in range(11)for j in range(11))<mp.mpf('1e-60')
 V=[]
 for ids in src['geometry']['vertices']:
  sub=mp.matrix([[Ns[i,j]for j in range(9)]for i in ids]);G=sub*Qs*sub.T;c=mp.lu_solve(G,sub*Qs*H);v=H-sub.T*c
  assert pair(v,v,Qs)<0 and max(Ns*Qs*v)<mp.mpf('1e-60')
  V.append(v/mp.sqrt(-pair(v,v,Qs)))
 D=mp.matrix([[mp.acosh(max(1,-pair(v,w,Qs)))for w in V]for v in V])
 Q=A.dm(b['form']);N=A.dm(b['final_normals']);r=A.dv(d['raw_cap']);q=A.dec(d['raw_cap_norm'])
 def f(x):return mp.mpf(x.numerator)/x.denominator
 def ev(x):return f(x.a.a)+f(x.a.b)*mp.sqrt(5)+(f(x.b.a)+f(x.b.b)*mp.sqrt(5))*mp.sqrt((16+4*mp.sqrt(5))/11)
 G=mp.matrix([[ev(A.pair(n,m,Q))for m in N]for n in N])
 H0=A.dv(b['center']);E0=[A.E(x)for x in d['literal_escape']]
 hh,he,ee=[ev(A.pair(x,y,Q))for x,y in[(H0,H0),(H0,E0),(E0,E0)]]
 h0=[ev(A.pair(H0,n,Q))for n in N];e0=[ev(A.pair(E0,n,Q))for n in N]
 def roots(a,b,c):
  disc=b*b-4*a*c
  if disc<0:return []
  return [(-b-mp.sqrt(disc))/(2*a),(-b+mp.sqrt(disc))/(2*a)]
 aa=hh-2*he+ee;bb=2*(he-hh);cc=hh
 tnull=min(t for t in roots(aa,bb,cc)if 0<t<1)
 pars=[]
 for i in range(12):
  z=e0[i]-h0[i]
  pars.extend((t,i)for t in roots(z*z-G[i,i]*aa,2*h0[i]*z-G[i,i]*bb,h0[i]*h0[i]-G[i,i]*cc)if tnull<t<1)
 tparallel,parallel_wall=min(pars)
 t=mp.mpf(d['interpolation_parameter'])if u is None else tnull+mp.mpf(u)*(tparallel-tnull)
 g=mp.matrix([h0[i]+t*(e0[i]-h0[i])for i in range(12)]);qr=aa*t*t+bb*t+cc
 if u is None:
  assert abs(qr-ev(q))<mp.mpf('1e-60')
  assert max(abs(g[i]-ev(A.pair(r,n,Q)))for i,n in enumerate(N))<mp.mpf('1e-60')
 assert qr>0 and max(g)<0
 PP=G-g*g.T/qr;k=[mp.sqrt(-PP[i,i])for i in range(12)]
 DC=mp.matrix([[mp.acosh(max(1,-PP[i,j]/(k[i]*k[j])))for j in range(12)]for i in range(12)])
 rho=[mp.asinh(mp.sqrt(G[i,i])/k[i])for i in range(12)]
 bound=mp.matrix([[DC[i,j]+rho[i]+rho[j]for j in range(12)]for i in range(12)])
 return dict(source=src,Qs=Qs,Ns=Ns,H=H,V=V,D=D,DC=DC,rho=rho,bound=bound,G=G,PP=PP,k=k,g=g,q=qr,t=t,u=(t-tnull)/(tparallel-tnull),tnull=tnull,tparallel=tparallel,parallel_wall=parallel_wall)
if __name__=='__main__':
 d=data();print('diameter',max(d['D']),'union diameter',max(d['bound']))
