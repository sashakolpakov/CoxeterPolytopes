#!/usr/bin/env python3
"""Exact continuous bridge gate over the literal low R endpoint domain.

The X conic parameter is never sampled. Both conic sheets and all isolated
high X-contact labels are covered. Frozen physical interval helpers are shared.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import fresh_new_cavity_construction_gate_20261017 as Q
import fresh_broken_cap_adjacency_20261017 as T

G,P,U,I,ROOT=Q.G,Q.P,Q.U,Q.I,Q.ROOT
TAG='fresh_zero_arc_cap_continuous_20261017'
TARGETS=[4,7,8,9,10,11,30,44,53,63]

def scalar(x):return x if isinstance(x,I)else I(x)
def vec(i):return {i:I(1)}
def scale(c,v):return {i:U.rnd(c*x)for i,x in v.items()}
def add(*vs):
 out={}
 for v in vs:
  for i,x in v.items():out[i]=U.rnd(out.get(i,I())+x)
 return out

def combine(cs,vs):return add(*(scale(c,v)for c,v in zip(cs,vs)))
def isnegative(x):return x.h<0

def quadratic_negative(a,b,c,left):
 """Sufficient exact coefficient certificate on t>=left by shifting t."""
 co=[a,U.rnd(2*a*left+b),U.rnd(a*left*left+b*left+c)]
 return all(x.h<0 for x in co),[x.enc()for x in co]

class Engine:
 def __init__(self):
  self.base=Q.Engine();self.counts=Counter();self.rbase=self.base;self.rpairs={};self.frame_cache={}
  self.a_frame=[vec(7),vec(8),vec(11)]
  ka=self.project(vec(53),self.a_frame);q=self.dot(ka,ka);assert q.h<0
  ka=scale(1/U.sqrt(-q),ka);sp=scale(1/U.sqrt(I(2)),vec(11))
  self.au=scale(1/U.sqrt(I(2)),add(ka,scale(-1,sp)))
  self.av=scale(1/U.sqrt(I(2)),add(scale(-1,ka),scale(-1,sp)))
  # The other spatial sheet pairs strictly positively with retained old11.
  bounds=[]
  for j in TARGETS:
   if j in[7,8,11]:continue
   a=self.dot(self.au,vec(j));b=self.dot(self.av,vec(j))
   if a.h<0 and b.l>0:
    root=U.sqrt(-b/a);bounds.append((root.l,j,root,a,b))
  assert bounds
  self.amin,self.awall,self.aroot,self.aa,self.ab=max(bounds,key=lambda x:x[0])
  self.a_certificate=dict(branch_exclusion='opposite spatial sheet has A.old11>0',
   selected_lower_bound_wall=self.awall,coefficient_t=self.aa.enc(),coefficient_inverse_t=self.ab.enc(),
   exact_root_interval=self.aroot.enc(),relaxed_lower_parameter=str(self.amin))
 def dot(self,a,b):
  def pair(i,j):
   if i==67 or j==67:return self.rpairs[j if i==67 else i]
   return self.base.pair(i,j)
  return U.rnd(sum((x*y*pair(i,j)for i,x in a.items()for j,y in b.items()),I()))
 def project(self,v,frame):
  gram=[[self.dot(a,b)for b in frame]for a in frame];iv,_=P.C53.inverse(gram)
  cs=U.imv(iv,[self.dot(v,f)for f in frame]);return add(v,scale(-1,combine(cs,frame)))
 def r_roots(self,rec):
  cp=G.profile(rec['profile']);rq=U.rnd(2-G.white(cp,cp))
  rc=[U.rnd(-Q.theta(m)-G.white(cp,self.base.sp.nodes[i]['cp']))for m,i in zip(rec['labels'],Q.RF)]
  coeff=U.imv(self.base.inv,rc);rem=U.rnd(rq-U.idot(rc,coeff));assert rem.l>0
  scl=U.sqrt(rem/self.base.kq);out=[]
  for sg in[-1,1]:
   pairs={67:rq}
   for j in TARGETS:
    tc=[self.base.pair(i,j)for i in Q.RF]
    kb=self.base.kq if j==self.base.seed else U.rnd(self.base.pair(self.base.seed,j)-U.idot(self.base.kc,tc))
    pairs[j]=U.rnd(U.idot(coeff,tc)+sg*scl*kb)
   for j,c in zip(Q.RF,rc):pairs[j]=c
   bad=[]
   for j in TARGETS:
    if j in Q.RF:continue
    val=U.rnd(pairs[j]+G.white(cp,self.base.sp.nodes[j]['cp']));ex=Q.exclusion(val)
    if ex:bad.append(dict(wall=j,**ex))
   out.append((sg,pairs,cp,rq,bad))
  return out
 def bridge(self,r,rcp,rq,rec,box):
  pi,m8,mr=rec;cp=G.profile(pi);xq=U.rnd(2-G.white(cp,cp))
  frame=[vec(8),r]
  if self.current not in self.frame_cache:
   gram=[[I(2),-Q.theta(self.r8)],[-Q.theta(self.r8),rq]];iv,_=P.C53.inverse(gram)
   k0=self.project(vec(53),frame);q0=self.dot(k0,k0);assert q0.h<0;k0=scale(1/U.sqrt(-q0),k0)
   seeds=[]
   for seed in[7,11,30,44,63]:
    kv=self.project(vec(seed),frame);kv=add(kv,scale(self.dot(kv,k0),k0));qs=self.dot(kv,kv)
    if qs.l>0:seeds.append((float((qs.h-qs.l)/qs.l),seed,kv,qs))
   assert seeds
   _,kernel_seed,k1,q1=min(seeds,key=lambda x:x[0]);k1=scale(1/U.sqrt(q1),k1)
   self.frame_cache[self.current]=(iv,k0,k1,kernel_seed)
  iv,k0,k1,kernel_seed=self.frame_cache[self.current]
  cx=[-box.get('X8',Q.theta(m8)),U.rnd(-box.get('XR',Q.theta(mr))-G.white(cp,rcp))]
  part=combine(U.imv(iv,cx),frame);d=U.rnd(xq-U.idot(cx,U.imv(iv,cx)))
  # All actual finite corner labels give d>0. Closed high boxes may meet d=0.
  if d.h<0:return dict(reason='negative_corner_residual',d=d.enc())
  d=I(max(F(0),d.l),d.h);branches=[]
  for sg in[-1,1]:
   xu=scale(1/U.sqrt(I(2)),add(k0,scale(sg,k1)))
   xvunit=scale(1/U.sqrt(I(2)),add(scale(-1,k0),scale(sg,k1)))
   xv=scale(d/2,xvunit)
   # X=part+s*xu+xv/s.  Its norm residual is exactly d, also as d tends to0.
   inequalities=[];impossible=[];lower=[]
   for j in TARGETS:
    if j==8:continue
    a=self.dot(xu,vec(j));b=U.rnd(d*self.dot(xvunit,vec(j))/2);c=U.rnd(self.dot(part,vec(j))+G.white(cp,self.base.sp.nodes[j]['cp'])+(U.sqrt(2*xq) if j in[30,63] else 0))
    inequalities.append((j,a,b,c))
    if a.l>0 and b.l>=0 and c.l>=0:impossible.append(dict(wall=j,coefficients=[a.enc(),c.enc(),b.enc()]))
    if a.h<0 and b.l>=0:
     disc=U.rnd(c*c-4*a*b)
     if disc.l>0:
      root=U.rnd((-c-U.sqrt(disc))/(2*a))
      if root.l>0:lower.append((root.l,j,root,a,b,c))
   if impossible:
    branches.append(dict(sign=sg,reason='strict_positive_native_quadratic',witness=impossible[0]));continue
   if not lower:return None
   s0,j,root,a,b,c=max(lower,key=lambda x:x[0]);t0=self.amin
   # f(t,s)=t*F(s)+G(s)/t, F=a*s+c+b/s, G=e*s+g+h/s.
   fa=self.dot(self.au,xu);fb=U.rnd(d*self.dot(self.au,xvunit)/2);fc=self.dot(self.au,part)
   ga=self.dot(self.av,xu);gb=U.rnd(d*self.dot(self.av,xvunit)/2);gc=self.dot(self.av,part)
   fneg,certF=quadratic_negative(fa,fc,fb,s0)
   ht=[U.rnd(t0*t0*x-y)for x,y in[(fa,ga),(fc,gc),(fb,gb)]]
   hneg,certH=quadratic_negative(*ht,s0)
   da=U.rnd(t0*fa+ga/t0);db=U.rnd(t0*fb+gb/t0)
   deriv=U.rnd(da*s0*s0-db)
   value=U.rnd(t0*(fa*s0+fc+fb/s0)+(ga*s0+gc+gb/s0)/t0)
   if not(fneg and hneg and da.h<0 and deriv.h<0 and value.h<-2):return None
   branches.append(dict(sign=sg,reason='AX_strictly_dotted_on_entire_parameter_quadrant',
    X_lower_bound=dict(wall=j,coefficients=[a.enc(),c.enc(),b.enc()],root_interval=root.enc(),relaxed_lower=str(s0)),
    F_shifted_coefficients=certF,t_derivative_shifted_coefficients=certH,
    s_derivative_leading=da.enc(),s_derivative_at_boundary=deriv.enc(),maximum_AX=value.enc()))
  return dict(d=d.enc(),X_projected_square=xq.enc(),intrinsic_nonintersection_threshold=(-U.sqrt(2*xq)).enc(),spacelike_seed=kernel_seed,branches=branches)
 def cover(self,r,cp,rq,rec,box,depth=0):
  self.counts['boxes']+=1
  for k,v in box.items():
   ex=Q.exclusion(-v)
   if ex:self.counts['angle_gap_leaves']+=1;return dict(box={k:v.enc()for k,v in box.items()},assigned_gap=dict(parameter=k,**ex))
  proof=self.bridge(r,cp,rq,rec,box)
  if proof is not None:
   self.counts['leaves']+=1;self.counts['max_depth']=max(self.counts['max_depth'],depth)
   return dict(box={k:v.enc()for k,v in box.items()},proof=proof)
  assert box and depth<20,('unresolved',self.current,rec,depth,{k:v.enc()for k,v in box.items()})
  k=max(box,key=lambda x:box[x].h-box[x].l);v=box[k];mid=(v.l+v.h)/2;children=[]
  for w in[I(v.l,mid),I(mid,v.h)]:
   nb=dict(box);nb[k]=w;children.append(self.cover(r,cp,rq,rec,nb,depth+1))
  return dict(split=k,mid=str(mid),children=children)

def run():
 e=Engine();dom=T.S.domain(Q.RF);rows=[];high=I(Q.theta(6).l,F(2))
 for idx,rec in enumerate(dom['low']):
  rr=[];e.r8=rec['labels'][0]
  for sign,pairs,cp,rq,bad in e.r_roots(rec):
   e.rpairs=pairs;r=vec(67)
   if bad:rr.append(dict(sign=sign,native_obstruction=bad));e.counts['excluded_R_branches']+=1;continue
   e.counts['R_branches']+=1;e.current=(idx,sign);cases=[]
   for pi,m8,mr in product(range(7),range(2,7),range(2,7)):
    if T.corner(pi,m8,mr,rec['profile'],e.r8)is None:continue
    if not any(T.corner(pi,m8,ma,0,2)is not None for ma in range(2,7)):continue
    e.counts['X_cases']+=1;box={k:high for k,m in[('X8',m8),('XR',mr)]if m==6}
    cases.append(dict(recipe=[pi,m8,mr],certificate=e.cover(r,cp,rq,(pi,m8,mr),box)))
   rr.append(dict(sign=sign,X_cases=cases))
  rows.append(dict(R_recipe=rec,branches=rr))
  print('ZERO_ARC_CONTINUOUS',idx+1,'of',len(dom['low']),dict(e.counts),flush=True)
 paths=[Path(__file__),Path(Q.__file__),Path(T.__file__),Path(T.S.__file__),Path(G.__file__),Path(G.E.X.__file__),G.CONTROL,T.H.PARENT,T.H.BOUNDARY]
 parent=json.loads((ROOT/'results'/f'{T.S.TAG}.json').read_text())
 return dict(status='PASS_EXACT_CONTINUOUS_ZERO_ARC_LOW_R_BRIDGE_GATE',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
  deleted_walls=T.S.DELETE,boundary_cycle=T.S.CYCLE,exact_old_deletion_cost=parent['exact_old_deletion_cost'],counts=dict(e.counts),A_conic=e.a_certificate,records=rows,
  scope='Every low R recipe (labels2..5) covering8-44 and44-63, every finite X-8-R corner and every continuous A parameter. All high X8/XR labels represented by full closed cosine intervals. X30 andX63 are required intrinsically nonintersecting: their projected products are at most -sqrt(2*qX), including parallel equality as a relaxation. Finite native contacts satisfying this threshold are included. The remaining tested old contacts are merely acute. A-X is proved strictly below-2, so no finite A-X corner. High labels at R, alternative old-wall releases, more connector caps and X intersecting old30 or63 intrinsically are not excluded. Shared exact physical interval arithmetic; no new chamber.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
 print(out['status'],out['counts'],flush=True)
