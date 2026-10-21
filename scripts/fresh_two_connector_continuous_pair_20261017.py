#!/usr/bin/env python3
"""Continuous right-angle A-X-Y-R discovery on old8; numerical only.

For fixed A,R, the full native-dotted X/Y conic intervals are intersected.
Their mutual product maximum is computed analytically on each rectangle.
Only the outer A parameter is searched by a bounded scalar routine.
"""
import argparse,json,math
from collections import Counter
from hashlib import sha256
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
import fresh_zero_arc_cap_two_connectors_20261017 as D

T,N,S,H,B,G,M,ROOT=D.T,D.N,D.S,D.H,D.B,D.G,D.M,D.ROOT
TAG='fresh_two_connector_continuous_pair_20261017'

def interval_data(ints):return [[a,None if not np.isfinite(b)else b]for a,b in ints]

def intervals(aa,bb,old,assigned,h,extra=()):
 out=[(0.,float('inf'))]
 for i,v in old:
  if i in assigned:continue
  a=float(aa@M@v);b=float(bb@M@v)
  if abs(a)+abs(b)<1e-9:continue
  out=D.meet(out,D.qset(a,2.,b))
  if not out:return []
 for v in extra:
  out=D.meet(out,D.qset(float(aa@M@v),2.,float(bb@M@v)))
  if not out:return []
 out=D.meet(out,D.qset(float(aa@M@h),0.,float(bb@M@h)))
 return out

def conics(anchor,ns,h,old,extra):
 f=np.array([ns[8],anchor]);gram=f@M@f.T
 assert np.max(np.abs(gram-2*np.eye(2)))<2e-6
 inv=np.linalg.inv(gram)
 def proj(v):return v-(v@M@f.T)@inv@f
 k0=proj(h);q=float(k0@M@k0)
 if q>=-1e-10:return []
 k0/=np.sqrt(-q);options=[]
 for j in [30,44,63,7,11]:
  v=proj(ns[j]);v+=float(v@M@k0)*k0;q=float(v@M@v)
  if q>1e-9:options.append((q,j,v))
 if not options:return []
 q,seed,k1=max(options,key=lambda z:z[0]);k1/=np.sqrt(q)
 out=[]
 for sg in[-1,1]:
  aa=np.sqrt(.5)*(k0+sg*k1);bb=np.sqrt(.5)*(-k0+sg*k1)
  ints=intervals(aa,bb,old,{8},h,[extra])
  if ints:out.append(dict(aa=aa,bb=bb,intervals=ints,sign=sg,seed=seed))
 return out

def pair_max(x,y,ix,iy):
 """All rectangle boundary extrema; any interior critical point is a saddle.

 f=a*s*z+b*s/z+c*z/s+d/(s*z). Future/past null orientation gives
 a,d<=0 and b,c>=0. The interior stationary Hessian has eigenvalues
 4*a*s*z and4*b*s/z, of opposite signs when nondegenerate.
 """
 a,b,c,d=[float(v@M@w)for v,w in[(x['aa'],y['aa']),(x['aa'],y['bb']),(x['bb'],y['aa']),(x['bb'],y['bb'])]]
 if a>=-1e-9 or d>=-1e-9 or b<=1e-9 or c<=1e-9:return None
 sl,sh=ix;zl,zh=iy;points=[]
 def val(s,z):return a*s*z+b*s/z+c*z/s+d/(s*z)
 for s in [sl]+([sh]if np.isfinite(sh)else[]):
  A=a*s+c/s;B0=b*s+d/s
  points.append((val(s,zl),s,zl))
  if np.isfinite(zh):points.append((val(s,zh),s,zh))
  elif A>=-1e-10:return dict(unbounded_or_degenerate=True,coefficients=[a,b,c,d])
  if A<0 and B0<0:
   z=math.sqrt(B0/A)
   if zl<z<zh:points.append((val(s,z),s,z))
 for z in [zl]+([zh]if np.isfinite(zh)else[]):
  A=a*z+b/z;B0=c*z+d/z
  points.append((val(sl,z),sl,z))
  if np.isfinite(sh):points.append((val(sh,z),sh,z))
  elif A>=-1e-10:return dict(unbounded_or_degenerate=True,coefficients=[a,b,c,d])
  if A<0 and B0<0:
   s=math.sqrt(B0/A)
   if sl<s<sh:points.append((val(s,z),s,z))
 value,s,z=max(points)
 return dict(maximum=value,s=s,z=z,coefficients=[a,b,c,d],boundary_candidates=len(points))

def run():
 z,labels,ns,h,profiles,primary=H.setup()
 old=[(i,ns[i])for i in labels if i not in S.DELETE+D.WHITE]
 actual=[i for i in z['actual_sides']if i not in S.DELETE]
 parent=json.loads((ROOT/'results'/f'{S.TAG}.json').read_text())
 rb=[r for r in parent['accepted_low_caps']if r['recipe']['profile']==0 and r['recipe']['labels'][0]==2]
 assert len(rb)==12
 k0,k1=T.a_curve(ns,h);counts=Counter();rows=[];accepted=[]
 for rr in rb:
  r=np.array(rr['projection']);local=[];best=None
  for asign in[-1,1]:
   aa=np.sqrt(.5)*(k0+asign*k1);bb=np.sqrt(.5)*(-k0+asign*k1)
   aints=intervals(aa,bb,old,{7,8},h,[r])
   for ai,(lo,hi)in enumerate(aints):
    assert lo>0
    # t=lo/u on an unbounded interval. The compact numerical search keeps
    # u>=1/64 and stays strictly inside the native-dotted endpoint.
    ul=max(1/64,lo/hi if np.isfinite(hi) else 0);uh=1-1e-7
    if ul>=uh:continue
    cache={}
    def evaluate(u):
     nonlocal best
     key=float(u)
     if key in cache:return cache[key]
     counts['A_evaluations']+=1;t=lo/u;aroot=aa*t+bb/t
     if abs(aroot@M@aroot-2)>3e-5:counts['unstable_A_norm']+=1;cache[key]=None;return None
     xs=conics(aroot,ns,h,old,r);ys=conics(r,ns,h,old,aroot);candidates=[]
     for x in xs:
      for y in ys:
       for ix in x['intervals']:
        for iy in y['intervals']:
         counts['complete_conic_rectangles']+=1;q=pair_max(x,y,ix,iy)
         if q is None or q.get('unbounded_or_degenerate'):
          counts['unresolved_rectangles']+=1;continue
         candidates.append((q['maximum'],q,x,y,ix,iy))
     ans=max(candidates,key=lambda v:v[0])if candidates else None
     if ans is not None and(best is None or ans[0]>best['maximum']):
      value,q,x,y,ix,iy=ans
      best=dict(maximum=value,R_index=rr['index'],R_recipe=rr['recipe'],A_sign=asign,A_t=t,A_u=u,A_interval=[lo,None if not np.isfinite(hi)else hi],
       X_sign=x['sign'],Y_sign=y['sign'],X_interval=interval_data([ix])[0],Y_interval=interval_data([iy])[0],pair=q,
       A_projection=aroot.tolist(),X_conic=[x['aa'].tolist(),x['bb'].tolist()],Y_conic=[y['aa'].tolist(),y['bb'].tolist()])
     cache[key]=ans;return ans
    def objective(u):
     a=evaluate(u);return -a[0] if a is not None else 1e9
    # Eight bounded outer intervals; no unbounded optimizer or label grid.
    cuts=np.linspace(ul,uh,9);opt=[]
    for left,right in zip(cuts,cuts[1:]):
     for u in[left,(left+right)/2,right]:evaluate(float(u))
     fit=minimize_scalar(objective,bounds=(left,right),method='bounded',options={'maxiter':45,'xatol':1e-9})
     evaluate(float(fit.x));opt.append(dict(interval=[float(left),float(right)],success=bool(fit.success),evaluations=int(fit.nfev),u=float(fit.x),objective=float(fit.fun)))
    local.append(dict(A_sign=asign,A_interval=[lo,None if not np.isfinite(hi)else hi],compact_u_domain=[ul,uh],optimization=opt))
  if best and best['maximum']>-2+1e-6:
   counts['potential_finite_pair_domains']+=1
   # Retain every promising domain for exact reconstruction; a numerical
   # maximum at a parallel boundary is not an accepted Coxeter chamber.
   accepted.append(best)
  rows.append(dict(R_index=rr['index'],R_recipe=rr['recipe'],A_domains=local,best=best))
  print('CONTINUOUS_PAIR',rr['index'],None if best is None else best['maximum'],dict(counts),flush=True)
 paths=[Path(__file__),Path(D.__file__),Path(T.__file__),Path(H.__file__),ROOT/'results'/f'{S.TAG}.json',H.PARENT,H.BOUNDARY]
 return dict(status='BOUNDED_OUTER_A_CONTINUOUS_TWO_CONNECTOR_PAIR_DIAGNOSTIC',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
  deleted_walls=S.DELETE,boundary_cycle=S.CYCLE,exact_old_deletion_cost=parent['exact_old_deletion_cost'],
  topology='A-X-Y-R alongold8; AX=YR=2, X/YrightS6+8, XYcontactleftfree. No A-R,A-Y,X-R adjacency.',
  counts=dict(counts,R_seeds=len(rb)),records=rows,potential_domains=accepted,
  scope='For each numerically tested A, all native-dotted conic parameter rectangles for X/Y are maximized analytically in floating arithmetic. The outer A parameter is bounded and numerically optimized, not exhaustively certified. R uses12 low native-admissible right8/profile0 seeds. Degenerate or unbounded pair cases remain unresolved. No interval theorem, exact wall, completed B/C patch, full boundary improvement or finite-volume claim.')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert out==json.loads(path.read_text())
 else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
 print(out['status'],out['counts'],'potential',len(out['potential_domains']),flush=True)
