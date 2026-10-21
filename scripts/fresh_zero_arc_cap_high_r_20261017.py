#!/usr/bin/env python3
"""All ten symbolic high-R families for the continuous single-X bridge.

Shares the frozen physical interval arithmetic and finite-corner classifier.
The R8-high half uses a new conic chart with no division by 4-cosine**2.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import product
from pathlib import Path
import fresh_zero_arc_cap_continuous_20261017 as L

Q,T,G,P,U,I,ROOT=L.Q,L.T,L.G,L.P,L.U,L.I,L.ROOT
vec,scale,add=L.vec,L.scale,L.add
combine=L.combine
TARGETS=L.TARGETS
TAG='fresh_zero_arc_cap_high_r_20261017'

@lru_cache(None)
def cosine(m):return U.rnd(Q.theta(m))

def exclusion(value):
 if value.l>0:return dict(kind='strict_obtuse',product=value.enc())
 if value.h<=-2:return None
 for m in range(2,65):
  lower,upper=-cosine(m+1),-cosine(m)
  if lower.h<value.l and value.h<upper.l:
   return dict(kind='strict_consecutive_angle_gap',product=value.enc(),labels=[m,m+1],
     lower_threshold=lower.enc(),upper_threshold=upper.enc())
 return None

class Engine(L.Engine):
 def __init__(self):
  super().__init__();self.counts=Counter()
  # The old acute lower endpoint A.53=0 is not Coxeter with old9.
  # Hence every actual A has A.53<=-1, not merely A.53<=0.
  assert self.awall==53 and self.aa.h<0 and self.ab.l>0
  endpoint9=U.rnd(self.dot(self.au,vec(9))*self.aroot+self.dot(self.av,vec(9))/self.aroot)
  gap=exclusion(endpoint9);assert gap and gap['kind']=='strict_consecutive_angle_gap'
  tightened=U.rnd((-1-U.sqrt(I(1)-4*self.aa*self.ab))/(2*self.aa))
  assert tightened.l>self.amin
  self.a_certificate.update(right53_endpoint_old9_gap=gap,nonright53_root_interval=tightened.enc(),
    strengthened_relaxed_lower_parameter=str(tightened.l),
    reason='A.53 is strictly decreasing; its only allowed value above-1 is0, whose unique parameter is excluded by old9')
  self.amin=tightened.l
  for i in [8,30,63]:
   v=self.base.pair(44,i);assert v.l==v.h==0,(i,v.enc())
  self.a44=self.dot(self.au,vec(44));self.b44=self.dot(self.av,vec(44))
  self.A044=U.rnd(self.amin*self.a44+self.b44/self.amin)
  assert self.A044.h<0
 def r_roots_box(self,rec,box):
  assert rec['profile']==0
  rc=[-box['R'+str(j)] if j in rec['high_positions'] else -Q.theta(m)
      for j,m in enumerate(rec['labels'])]
  co=U.imv(self.base.inv,rc);rem=U.rnd(I(2)-U.idot(rc,co))
  if rem.h<0:return [],dict(reason='R_no_real_norm_branch',remainder=rem.enc())
  scl=Q.sqrt_closed(I(max(F(0),rem.l),rem.h)/self.base.kq);out=[]
  for sg in [-1,1]:
   pairs={67:I(2)}
   for j in L.TARGETS:
    tc=[self.base.pair(i,j) for i in Q.RF]
    kb=self.base.kq if j==self.base.seed else U.rnd(self.base.pair(self.base.seed,j)-U.idot(self.base.kc,tc))
    pairs[j]=U.rnd(U.idot(co,tc)+sg*scl*kb)
   for j,c in zip(Q.RF,rc):pairs[j]=c
   bad=[]
   for j in L.TARGETS:
    if j not in Q.RF:
     ex=exclusion(pairs[j])
     if ex:bad.append(dict(wall=j,**ex))
   out.append((sg,pairs,bad))
  return out,dict(remainder=rem.enc(),scale=scl.enc())
 def bridge(self,r,rcp,rq,rec,box):
  pi,m8,mr=rec;cp=G.profile(pi);xq=U.rnd(2-G.white(cp,cp))
  frame=[vec(8),r]
  if self.current not in self.frame_cache:
   gram=[[I(2),-Q.theta(self.r8)],[-Q.theta(self.r8),rq]];iv,_=P.C53.inverse(gram)
   if self.family['labels']==[2,6,2]:
    # R=-c*44/2+sg*sqrt(delta/(2*kq))*Kseed.  The regular
    # orthonormal complement below avoids subtracting large projections
    # and dividing by sqrt(delta) at the parallel endpoint.
    c=-self.rpairs[44];delta=U.rnd(4-c*c);delta=I(max(F(0),delta.l),delta.h)
    g=self.base.pair(8,63)
    k0=scale(1/U.sqrt(g*g/2-2),add(vec(63),scale(-g/2,vec(8))))
    kv=add(vec(self.base.seed),scale(-1,combine(self.base.kc,[vec(i)for i in Q.RF])))
    k1=add(scale(Q.sqrt_closed(delta)/(2*U.sqrt(I(2))),vec(44)),
           scale(self.rsign*c/(2*U.sqrt(self.base.kq)),kv))
    kernel_seed='regularized_R44_chart'
   else:
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
    a=self.dot(xu,vec(j));b=U.rnd(d*self.dot(xvunit,vec(j))/2);c=U.rnd(self.dot(part,vec(j))+G.white(cp,self.base.sp.nodes[j]['cp'])+(U.sqrt(2*xq) if j in[30,44,63] else 0))
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
   fneg,certF=L.quadratic_negative(fa,fc,fb,s0)
   ht=[U.rnd(t0*t0*x-y)for x,y in[(fa,ga),(fc,gc),(fb,gb)]]
   hneg,certH=L.quadratic_negative(*ht,s0)
   da=U.rnd(t0*fa+ga/t0);db=U.rnd(t0*fb+gb/t0)
   deriv=U.rnd(da*s0*s0-db)
   value=U.rnd(t0*(fa*s0+fc+fb/s0)+(ga*s0+gc+gb/s0)/t0)
   if not(fneg and hneg and da.h<0 and deriv.h<0 and value.h<-2):return None
   branches.append(dict(sign=sg,reason='AX_strictly_dotted_on_entire_parameter_quadrant',
    X_lower_bound=dict(wall=j,coefficients=[a.enc(),c.enc(),b.enc()],root_interval=root.enc(),relaxed_lower=str(s0)),
    F_shifted_coefficients=certF,t_derivative_shifted_coefficients=certH,
    s_derivative_leading=da.enc(),s_derivative_at_boundary=deriv.enc(),maximum_AX=value.enc()))
  return dict(d=d.enc(),X_projected_square=xq.enc(),intrinsic_nonintersection_threshold=(-U.sqrt(2*xq)).enc(),spacelike_seed=kernel_seed,branches=branches)
 def regular_high8(self,box):
  """All X profiles at once, using the projected finite-corner threshold."""
  c=box['R0'];delta=U.rnd(4-c*c)
  # The actual high cosine lies below2; include the parallel limit in the box.
  delta=I(max(F(0),delta.l),delta.h);assert delta.h>=0
  g=self.base.pair(8,63);h=self.rpairs[63]
  K=U.rnd(2*(g*g+c*g*h+h*h-delta))
  if K.l<=0:return None
  sq=U.sqrt(K)
  N=scale(1/sq,add(scale(delta,vec(63)),scale(-(2*g+c*h),vec(8)),scale(-(c*g+2*h),vec(67))))
  n30=self.dot(N,vec(30))
  if n30.l>=0:return dict(reason='X63_orientation_forces_X30_nonnegative',K=K.enc(),N30=n30.enc())
  bounds=[((1/sq).l,63,1/sq)]
  if n30.h<0:bounds.append(((1/(-n30)).l,30,1/(-n30)))
  w0,wall,root=max(bounds,key=lambda z:z[0]);assert w0>0
  r0=U.sqrt(I(F(1,4))/(w0*w0)+delta/2);r=I(0,r0.h)
  an=self.dot(self.au,N);bn=self.dot(self.av,N);t0=self.amin
  f=U.rnd(an-r*self.a44)
  hder=U.rnd(t0*t0*an-bn-r*(t0*t0*self.a44-self.b44))
  A0n=U.rnd(t0*an+bn/t0)
  wder=U.rnd(A0n-r*self.A044)
  value=U.rnd(w0*A0n-U.sqrt(I(F(1,4))+delta*w0*w0/2)*self.A044)
  if not(f.h<0 and hder.h<0 and wder.h<0 and value.h<-1):return None
  return dict(reason='all_X_profiles_fail_projected_finite_AX_corner',delta=delta.enc(),K=K.enc(),N30=n30.enc(),
    opposite_X_sheet='strictly positive product with old44',
    w_lower=dict(wall=wall,root_interval=root.enc(),relaxed_lower=str(w0)),
    ratio_s_over_w_enclosure=r.enc(),auN=an.enc(),avN=bn.enc(),A0_old44=self.A044.enc(),
    F_over_w=f.enc(),t_derivative_over_w=hder.enc(),w_derivative_upper=wder.enc(),
    maximum_AX_divided_by_sqrt_2q=value.enc())
 def evaluate(self,rec,xrec,box):
  for k,v in box.items():
   ex=exclusion(-v)
   if ex:return dict(reason='assigned_cosine_gap',parameter=k,**ex)
  roots,rproof=self.r_roots_box(rec,box)
  if not roots:return rproof
  branches=[]
  for sg,pairs,bad in roots:
   if bad:
    branches.append(dict(sign=sg,native_obstruction=bad[0]));continue
   self.rpairs=pairs;self.r8=rec['labels'][0];self.family=rec;self.rsign=sg
   self.current=(sg,tuple((k,tuple(v.enc())) for k,v in sorted(box.items())))
   self.frame_cache.clear()
   try:
    proof=(self.regular_high8(box) if self.r8==6 else
      self.bridge(vec(67),G.profile(0),I(2),xrec,{k:v for k,v in box.items() if k.startswith('X')}))
   except (AssertionError,ZeroDivisionError):return None
   if proof is None:return None
   branches.append(dict(sign=sg,proof=proof))
  return dict(R=rproof,branches=branches)
 def cover_high(self,rec,xrec,box,depth=0):
  self.counts['nodes']+=1
  if self.counts['nodes']%250==0:print('HIGH_R_BOXES',dict(self.counts),'depth',depth,flush=True)
  proof=self.evaluate(rec,xrec,box)
  if proof is not None:
   self.counts['leaves']+=1;self.counts['max_depth']=max(self.counts['max_depth'],depth)
   if proof.get('reason')=='assigned_cosine_gap':self.counts['assigned_gap_leaves']+=1
   return dict(box={k:v.enc() for k,v in box.items()},proof=proof)
  assert depth<32 and self.counts['nodes']<50000,('unresolved',rec,xrec,depth,{k:v.enc()for k,v in box.items()})
  k=max(box,key=lambda z:box[z].h-box[z].l);v=box[k];mid=(v.l+v.h)/2
  children=[]
  for w in [I(v.l,mid),I(mid,v.h)]:
   nb=dict(box);nb[k]=w;children.append(self.cover_high(rec,xrec,nb,depth+1))
  self.counts['splits']+=1
  return dict(split=k,mid=str(mid),children=children)

def run(selected=None):
 e=Engine();dom=T.S.domain(Q.RF);assert len(dom['tails'])==10
 high=I(Q.theta(6).l,F(2));rows=[]
 assert [(p,a,b)for p in range(7)for a in range(2,7)for b in range(2,7)if T.S.finite_corner(7,8,p,a,b)]==[(0,2,2)]
 for idx,rec in enumerate(dom['tails']):
  if selected is not None and idx not in selected:continue
  cases=[];rbox={'R'+str(j):high for j in rec['high_positions']}
  xc=[(pi,m8,mr)for pi,m8,mr in product(range(7),range(2,7),range(2,7))
      if T.corner(pi,m8,mr,rec['profile'],rec['labels'][0]) is not None
      and any(T.corner(pi,m8,ma,0,2) is not None for ma in range(2,7))]
  if rec['labels'][0]==6:
   assert xc==[(pi,2,2)for pi in range(7)],xc
   cases.append(dict(X_recipes=xc,certificate=e.cover_high(rec,None,rbox)))
  else:
   for xr in xc:
    box=dict(rbox)
    for k,m in [('X8',xr[1]),('XR',xr[2])]:
     if m==6:box[k]=high
    cases.append(dict(X_recipes=[xr],certificate=e.cover_high(rec,xr,box)))
  rows.append(dict(R_index=idx,R_recipe=rec,X_recipe_count=len(xc),cases=cases))
  print('HIGH_R_FAMILY_DONE',idx,rec['labels'],'Xcases',len(xc),dict(e.counts),flush=True)
 paths=[Path(__file__),Path(L.__file__),Path(Q.__file__),Path(T.__file__),Path(T.S.__file__),Path(G.__file__),Path(G.E.X.__file__),G.CONTROL,T.H.PARENT,T.H.BOUNDARY,
        ROOT/'results'/f'{L.TAG}.json',ROOT/'results'/f'{L.TAG}.md']
 parent=json.loads((ROOT/'results'/f'{L.TAG}.json').read_text())
 return dict(status='PASS_EXACT_CONTINUOUS_ZERO_ARC_HIGH_R_BRIDGE_GATE' if selected is None else 'PASS_SELECTED_HIGH_R_FAMILIES',
  sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
  counts=dict(e.counts,R_families=len(rows),X_recipes=sum(r['X_recipe_count']for r in rows)),
  A_conic=e.a_certificate,deleted_walls=parent['deleted_walls'],boundary_cycle=parent['boundary_cycle'],exact_old_deletion_cost=parent['exact_old_deletion_cost'],
  regular_chart=dict(high8_R_indices=[5,6,7,8,9],kernel_direction='N=(delta*63-(2g+c*h)*8-(c*g+2h)*R)/sqrt(K)',
    identities=['N^2=-delta','N.44=0','N.63=-sqrt(K)'],X_normalized='X/sqrt(2q)=w*N-sqrt(1/4+delta*w*w/2)*44',
    scope='Actual c<2; interval certificates include the parallel limit delta=0. No division by delta. All seven X A2 profiles are covered together.'),
  records=rows,scope='All ten high-R families, both norm branches, all endpoint labels and all continuous A and X parameters in the stated single zero-arc-X incidence. X must be intrinsically nonintersecting old30,44,63, with projected threshold -sqrt(2*qX). Combined with the frozen low-R package this closes that incidence, not other patches or a finite-volume construction.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--families');args=ap.parse_args()
 selected=None if args.families is None else set(map(int,args.families.split(',')))
 assert not(args.verify and selected is not None)
 out=json.loads(json.dumps(run(selected)))
 if selected is None:
  p=ROOT/'results'/f'{TAG}.json'
  if args.verify:assert out==json.loads(p.read_text())
  else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
 print(out['status'],out['counts'],flush=True)
