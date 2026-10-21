#!/usr/bin/env python3
"""Bounded finite-endpoint discovery on the new cap49/old33 pi/10 ridge."""
import argparse,json
from collections import Counter
from fractions import Fraction
from pathlib import Path
from hashlib import sha256
import mpmath as mp
import fresh_nonarithmetic_escape49_20261016_polar as P
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_ridge_pair_20261016_endpoints';CL='results/fresh_nonarithmetic_affine_filtered_20261016_cluster.json';CO='results/fresh_nonarithmetic_affine_filtered_20261016_cone.json'
mp.mp.dps=80
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def midpoint(a):
 q=(Fraction(a[0])+Fraction(a[1]))/2;return mp.mpf(q.numerator)/q.denominator

def run():
 j=P.B.read(P.BATCH);cl=P.B.read(CL);co=P.B.read(CO)
 Q=mp.matrix([[P.scalar(P.D.encfield(P.D.lift(P.A.B.db(x))))for x in row]for row in j['form']]);N=[mp.matrix([P.scalar(x)for x in row])for row in j['normals']]
 N=[n*mp.sqrt(2/(n.T*Q*n)[0])for n in N]+[mp.matrix([midpoint(x)for x in c['normal_interval']])for c in cl['caps']]
 H=mp.matrix([P.scalar(x)for x in j['center']]);E=mp.matrix([P.scalar(x)for x in co['strict_escape']['point']]);E=E/(-(E.T*Q*H)[0]);H=H/mp.sqrt(-(H.T*Q*H)[0])
 supports=[];candidates=[];stats=Counter();diagnostics=[]
 for cap in[49,50,51]:
  common=[i for i in range(49)if abs((N[i].T*Q*N[33])[0])<mp.mpf('1e-40')and abs((N[i].T*Q*N[cap])[0])<mp.mpf('1e-40')]
  supports.append(dict(cap=cap,common_old_orthogonal=common))
 # First cap's complete one-root spherical extensions of H4+A2.
 S=[0,1,2,3,5,6,33,49];W=mp.matrix([[N[i][j]for i in S]for j in range(10)]);G=W.T*Q*W;J=G**-1;assert min(mp.eigsy(G,eigvals_only=True))>0
 T=H-W*J*W.T*Q*H;T=T/mp.sqrt(-(T.T*Q*T)[0]);assert(T.T*Q*H)[0]<0
 d=N[7]-W*J*W.T*Q*N[7];d=d+(d.T*Q*T)[0]*T;d=d/mp.sqrt((d.T*Q*d)[0]);assert abs((d.T*Q*H)[0])<mp.mpf('1e-35')
 profiles=[(2,2)]+[(m,2)for m in[3,4,5]]+[(2,m)for m in[3,4,5]]
 for m5,m6 in profiles:
  c=mp.matrix([0,0,0,0,-2*mp.cos(mp.pi/m5)if m5!=2 else 0,-2*mp.cos(mp.pi/m6)if m6!=2 else 0,0,0]);p=W*J*c;rho=2-(p.T*Q*p)[0];assert rho>0
  counts=Counter()
  for anchor in range(52):
   if anchor in S:continue
   at=(T.T*Q*N[anchor])[0];ad=(d.T*Q*N[anchor])[0];pp=(p.T*Q*N[anchor])[0]
   for m in range(2,21):
    val=(-2*mp.cos(mp.pi/m)if m!=2 else 0)-pp
    # u=p+x*T+y*d; -x²+y²=rho and at*x+ad*y=val.
    if abs(ad)<mp.mpf('1e-40'):
     if abs(at)<mp.mpf('1e-40'):counts['degenerate_anchor']+=1;continue
     x=val/at;ys=[mp.sqrt(rho+x*x),-mp.sqrt(rho+x*x)];branches=[(x,y)for y in ys]
    else:
     aa=at*at-ad*ad;bb=-2*val*at;cc=val*val-ad*ad*rho
     if abs(aa)<mp.mpf('1e-40'):
      if abs(bb)<mp.mpf('1e-40'):counts['degenerate_equation']+=1;continue
      xs=[-cc/bb]
     else:
      disc=bb*bb-4*aa*cc
      if disc<0:continue
      xs=[(-bb+sgn*mp.sqrt(disc))/(2*aa)for sgn in[-1,1]]
     branches=[(x,(val-at*x)/ad)for x in xs]
    for branch,(x,y)in enumerate(branches):
     counts['real_branches']+=1;r=p+x*T+y*d;assert abs((r.T*Q*r)[0]-2)<mp.mpf('1e-30')
     hp=(r.T*Q*H)[0];ep=(r.T*Q*E)[0]
     if hp>=0:counts['wrong_center']+=1;continue
     vals=[(r.T*Q*n)[0]for n in N]
     if max(vals)>mp.mpf('1e-30'):counts['obtuse']+=1;continue
     counts['acute']+=1
     if ep<=0:counts['noncutting']+=1
     else:counts['acute_cutting']+=1
     bad=[];finite=[]
     for i,g in enumerate(vals):
      if abs(g)<mp.mpf('1e-30'):finite.append([i,2]);continue
      if g<-2-mp.mpf('1e-30'):continue
      if abs(g+2)<mp.mpf('1e-30'):bad.append([i,'near_parallel']);continue
      mm=mp.pi/mp.acos(-g/2);mi=int(mp.nint(mm))
      if abs(mm-mi)<mp.mpf('1e-25'):finite.append([i,mi])
      else:bad.append([i,mp.nstr(mm,30)])
     row=dict(profile=[m5,m6],anchor=anchor,anchor_label=m,branch=branch,escape_pairing=mp.nstr(ep,40),center_pairing=mp.nstr(hp,40),finite=finite,bad=bad,point=[mp.nstr(z,60)for z in r])
     if not bad:candidates.append(row);print('PROPOSAL',row,flush=True)
     elif len(diagnostics)<100:diagnostics.append(row)
  stats.update(counts);supports.append(dict(profile=[m5,m6],counts=dict(counts)));print('profile',m5,m6,dict(counts),flush=True)
 paths=[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py',P.BATCH,CL,CO]
 return dict(status='BOUNDED_NONCUTTING_PI10_RIDGE_ENDPOINT_DISCOVERY',sources={p:sha(p)for p in paths},ridge=[33,49],spherical_support=S,profile_labels=profiles,anchor_labels=list(range(2,21)),support_checks=supports,stats=dict(stats),candidates=candidates,first100_failed_angle_candidates=diagnostics,scope='Numerical80-digit discovery for the seven one-root finite extensions of the literal H4+A2 support, with a finite anchor-label pilot2..20. Not all normals rightto the pi/10 ridge and not a universal nonexistence result. Exactify any surviving candidate before geometric claims.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved bounded pi10 ridge scout')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['candidates']))
