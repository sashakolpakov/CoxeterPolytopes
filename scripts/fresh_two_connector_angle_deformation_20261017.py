#!/usr/bin/env python3
"""Continuous endpoint-angle relaxation of A-X-Y-R, numerical discovery only.

For fixed A and endpoint cosines, enumerate rectangle-boundary maxima and
all roots of the degree-eight interior stationary polynomial. Only three
outer variables (A and two endpoint cosines) use a bounded optimizer.
"""
import argparse,json,math
from collections import Counter
from hashlib import sha256
from pathlib import Path
import numpy as np
from numpy.polynomial import Polynomial as Poly
from scipy.optimize import differential_evolution
import fresh_two_connector_continuous_pair_20261017 as P
D,T,N,S,H,B,G,M,ROOT=P.D,P.T,P.N,P.S,P.H,P.B,P.G,P.M,P.ROOT
TAG='fresh_two_connector_angle_deformation_20261017'

def conics(anchor,theta,ns,h,old,extra):
 f=np.array([ns[8],anchor]);gram=f@M@f.T
 if np.max(np.abs(gram-2*np.eye(2)))>2e-6:return []
 inv=np.linalg.inv(gram)
 def proj(v):return v-(v@M@f.T)@inv@f
 part=-theta*anchor/2;delta=2-theta*theta/2
 if delta<=0:return []
 k0=proj(h);q=float(k0@M@k0)
 if q>=-1e-10:return []
 k0/=np.sqrt(-q);options=[]
 for j in [30,44,63,7,11]:
  v=proj(ns[j]);v+=float(v@M@k0)*k0;q=float(v@M@v)
  if q>1e-9:options.append((q,j,v))
 if not options:return []
 q,seed,k1=max(options,key=lambda z:z[0]);k1/=np.sqrt(q);out=[]
 for sg in [-1,1]:
  aa=np.sqrt(delta)/2*(k0+sg*k1);bb=np.sqrt(delta)/2*(-k0+sg*k1)
  ints=[(0.,float('inf'))]
  for i,v in old+[('extra',extra),('center',h)]:
   if i==8:continue
   a,b,c=[float(w@M@v)for w in [aa,bb,part]];bound=0 if i=='center' else 2
   if abs(a)+abs(b)<1e-9:
    if abs(c)<1e-8:continue
    if c+bound>1e-9:ints=[];break
   else:ints=D.meet(ints,D.qset(a,c+bound,b))
   if not ints:break
  if ints:out.append(dict(part=part,aa=aa,bb=bb,sign=sg,seed=seed,intervals=ints))
 return out

def maximum(x,y,ix,iy):
 a,b,c,d=[float(v@M@w)for v,w in [(x['aa'],y['aa']),(x['aa'],y['bb']),(x['bb'],y['aa']),(x['bb'],y['bb'])]]
 e,f,g,h=[float(v@M@w)for v,w in [(x['aa'],y['part']),(x['bb'],y['part']),(x['part'],y['aa']),(x['part'],y['bb'])]]
 i=float(x['part']@M@y['part']);sl,sh=ix;zl,zh=iy
 if min(sl,zl)<=0 or a>=-1e-10 or d>=-1e-10 or min(b,c)<=1e-10:return None
 if not np.isfinite(sh) and a*zl+b/zl+e>=-1e-10:return None
 if not np.isfinite(zh) and a*sl+c/sl+g>=-1e-10:return None
 def val(s,z):return a*s*z+b*s/z+c*z/s+d/(s*z)+e*s+f/s+g*z+h/z+i
 pts=[]
 for s in [sl]+([sh]if np.isfinite(sh)else []):
  A=a*s+c/s+g;B0=b*s+d/s+h
  for z in [zl]+([zh]if np.isfinite(zh)else []):pts.append((val(s,z),s,z,'corner'))
  if A<0 and B0<0:
   z=math.sqrt(B0/A)
   if zl<z<zh:pts.append((val(s,z),s,z,'boundary_z'))
 for z in [zl]+([zh]if np.isfinite(zh)else []):
  A=a*z+b/z+e;B0=c*z+d/z+f
  if A<0 and B0<0:
   s=math.sqrt(B0/A)
   if sl<s<sh:pts.append((val(s,z),s,z,'boundary_s'))
 # Let A(s)=P/s, B(s)=Q/s. Interior z-max is C-2sqrt(PQ)/s.
 # Its stationary equation implies (e*s^2-f)^2*PQ-(s*(PQ)'-2PQ)^2=0.
 L=Poly([c,g,a])*Poly([d,h,b]);ss=Poly([0,1]);der=ss*L.deriv()-2*L
 polynomial=Poly([-f,0,e])**2*L-der**2
 co=polynomial.coef;scale=max(abs(co))
 if scale<1e-20:return None
 roots=Poly(co/scale).roots();nreal=0
 for r in roots:
  if abs(r.imag)>1e-6*max(1,abs(r.real)):continue
  s=float(r.real)
  if not sl<s<sh:continue
  A=a*s+c/s+g;B0=b*s+d/s+h
  if A>=0 or B0>=0:continue
  z=math.sqrt(B0/A)
  if not zl<z<zh:continue
  nreal+=1;pts.append((val(s,z),s,z,'stationary_polynomial'))
 value,s,z,where=max(pts)
 return dict(maximum=value,s=s,z=z,where=where,coefficients=[a,b,c,d,e,f,g,h,i],stationary_candidates=nreal)

def run():
 z,labels,ns,h,profiles,primary=H.setup();old=[(i,ns[i])for i in labels if i not in S.DELETE+D.WHITE]
 parent=json.loads((ROOT/'results'/f'{S.TAG}.json').read_text());rb=[r for r in parent['accepted_low_caps']if r['recipe']['profile']==0 and r['recipe']['labels'][0]==2]
 k0,k1=T.a_curve(ns,h);counts=Counter();rows=[];potential=[]
 for rr in rb:
  r=np.array(rr['projection']);best=None;opts=[]
  for asign in [-1,1]:
   aa=np.sqrt(.5)*(k0+asign*k1);bb=np.sqrt(.5)*(-k0+asign*k1)
   for lo,hi in P.intervals(aa,bb,old,{7,8},h,[r]):
    ul=max(1/32,lo/hi if np.isfinite(hi)else 0);uh=1-1e-6
    if ul>=uh:continue
    def objective(v):
     nonlocal best
     u,tx,ty=v;t=lo/u;aroot=aa*t+bb/t;counts['outer_evaluations']+=1
     if abs(aroot@M@aroot-2)>3e-5:counts['unstable_A']+=1;return 1e8
     xs=conics(aroot,tx,ns,h,old,r);ys=conics(r,ty,ns,h,old,aroot);found=[]
     for x in xs:
      for y in ys:
       for ix in x['intervals']:
        for iy in y['intervals']:
         counts['conic_rectangles']+=1;ans=maximum(x,y,ix,iy)
         if ans is None:counts['unresolved_rectangles']+=1;continue
         found.append(ans['maximum'])
         if best is None or ans['maximum']>best['maximum']:
          best=dict(R_index=rr['index'],R_recipe=rr['recipe'],maximum=ans['maximum'],A_t=t,A_sign=asign,A_projection=aroot.tolist(),theta_AX=tx,theta_YR=ty,X_interval=P.interval_data([ix])[0],Y_interval=P.interval_data([iy])[0],pair=ans,
           X_conic=[x[k].tolist()for k in ['part','aa','bb']],Y_conic=[y[k].tolist()for k in ['part','aa','bb']])
     return -max(found)if found else 1e8
    fit=differential_evolution(objective,[(ul,uh),(0.,1.999),(0.,1.999)],seed=1700+rr['index'],maxiter=18,popsize=6,polish=False,tol=1e-8,atol=1e-9)
    opts.append(dict(A_sign=asign,u_domain=[ul,uh],success=bool(fit.success),nfev=int(fit.nfev),x=fit.x.tolist(),objective=float(fit.fun)))
  rows.append(dict(R_index=rr['index'],R_recipe=rr['recipe'],best=best,optimizations=opts))
  if best and best['maximum']>-2+1e-5:potential.append(best)
  print('ANGLE_DEFORMATION',rr['index'],None if best is None else best['maximum'],dict(counts),flush=True)
 paths=[Path(__file__),Path(P.__file__),Path(D.__file__),Path(T.__file__),Path(H.__file__),ROOT/'results'/f'{S.TAG}.json',H.PARENT,H.BOUNDARY]
 return dict(status='BOUNDED_CONTINUOUS_ENDPOINT_ANGLE_RELAXATION',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},counts=dict(counts,R_seeds=len(rb)),records=rows,potential_domains=potential,deleted_walls=S.DELETE,exact_old_deletion_cost=parent['exact_old_deletion_cost'],
  scope='Numerical discovery only. Endpoint products are relaxed continuously from0 to1.999; no Coxeter claim for relaxed angles. All X/Y native-dotted conic intervals are tested at each outer parameter; degree-eight roots give interior stationary candidates numerically. Outer A and endpoint-angle domains are sampled by one fixed bounded optimizer per R seed. No all-parameter exclusion, exact geometry, B/C completion, or finite-volume claim.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
 if args.verify:assert out==json.loads(path.read_text())
 else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
 print(out['status'],out['counts'],'potential',len(out['potential_domains']),flush=True)
