#!/usr/bin/env python3
"""Bounded mixed-contact conic discovery on the renewed escape's active star."""
import argparse,json,gzip,itertools,math
from collections import Counter
from pathlib import Path
from hashlib import sha256
from fractions import Fraction
import mpmath as mp
import numpy as np
import fresh_nonarithmetic_escape49_20261016_polar as P
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_ridge_bridge_20261016_general_scout'
PROFILE_LABELS=[3,4,5,6,10];ANCHOR_LABELS=list(range(2,21))
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def smooth(x):return float(format(float(x),'.13g'))
def run():
 d=P.B.read(P.BATCH);cone=P.B.read('results/fresh_nonarithmetic_affine_filtered_20261016_cone.json');cl=P.B.read('results/fresh_nonarithmetic_affine_filtered_20261016_cluster.json');mp.mp.dps=80
 Q=mp.matrix([[P.scalar(P.D.encfield(P.D.lift(P.A.B.db(x))))for x in row]for row in d['form']])
 N=[mp.matrix([P.scalar(x)for x in row])for row in d['normals']]
 def mid(x):
  q=(Fraction(x[0])+Fraction(x[1]))/2;return mp.mpf(q.numerator)/q.denominator
 N += [mp.matrix([mid(x)for x in row['normal_interval']])for row in cl['caps']]
 N=[n*mp.sqrt(2/(n.T*Q*n)[0])for n in N]
 H=mp.matrix([P.scalar(x)for x in d['center']]);E=mp.matrix([P.scalar(x)for x in cone['strict_escape']['point']]);E/=-(E.T*Q*H)[0]
 G=mp.matrix([[(n.T*Q*m)[0]for m in N]for n in N]);h=mp.matrix([(H.T*Q*n)[0]for n in N]);e=mp.matrix([(E.T*Q*n)[0]for n in N])
 active=cone['strict_escape']['active'];assert len(active)==9
 stats=Counter();supports=[];accepted={};near=[]
 for omitted in active:
  S=[i for i in active if i!=omitted];GS=mp.matrix([[G[i,j]for j in S]for i in S])
  try:L=mp.cholesky(GS)
  except ValueError:supports.append(dict(support=S,positive_definite_numeric=False));continue
  J=GS**-1;cols=mp.matrix([[N[j][i]for j in S]for i in range(10)])
  T=H-cols*J*mp.matrix([h[i]for i in S]);tq=(T.T*Q*T)[0];assert tq<0;T/=mp.sqrt(-tq)
  best=None
  for i in range(52):
   w=N[i]-cols*J*mp.matrix([G[j,i]for j in S]);w+=((w.T*Q*T)[0])*T;ww=(w.T*Q*w)[0]
   if best is None or ww>best[0]:best=(ww,w)
  assert best[0]>0;D=best[1]/mp.sqrt(best[0]);assert abs((D.T*Q*T)[0])<mp.mpf('1e-50')
  cross=np.array([[float(G[i,j])for j in range(52)]for i in S]);JJ=np.array([[float(x)for x in row]for row in J.tolist()])
  tt=np.array([float((T.T*Q*n)[0])for n in N]);dd=np.array([float((D.T*Q*n)[0])for n in N]);tt[S]=0;dd[S]=0
  hv=np.array([float(h[i])for i in S]);ev=np.array([float(e[i])for i in S]);th=float((T.T*Q*H)[0]);dh=float((D.T*Q*H)[0]);te=float((T.T*Q*E)[0]);de=float((D.T*Q*E)[0])
  out=[49,50,51];anchors=np.repeat(out,len(ANCHOR_LABELS));labels=np.tile(ANCHOR_LABELS,len(out));target=-2*np.cos(np.pi/labels)
  qa=dd[anchors]**2-tt[anchors]**2;validq=abs(qa)>1e-12;stats['skipped_null_anchor_slots']+=int((~validq).sum())
  profiles=[[]]+[[(i,m)]for i in S for m in PROFILE_LABELS]+[[(i,m),(j,n)]for i,j in itertools.combinations(S,2)for m in PROFILE_LABELS for n in PROFILE_LABELS]
  counts=Counter()
  for profile in profiles:
   cv=np.zeros(8)
   for i,m in profile:cv[S.index(i)]=-2*math.cos(math.pi/m)
   coeff=cv@JJ;rho=2-coeff@cv
   pv=coeff@cross;w=target-pv[anchors];disc=w*w-rho*qa
   valid=validq&(disc>1e-10);stats['profiles']+=1;counts['profiles']+=1
   for sign in[-1,1]:
    take=np.flatnonzero(valid);s=np.sqrt(disc[take]);j=anchors[take];q=qa[take]
    aa=(-w[take]*tt[j]+sign*dd[j]*s)/q;bb=(w[take]*dd[j]-sign*tt[j]*s)/q
    vals=pv[None,:]+aa[:,None]*tt+bb[:,None]*dd
    hp=coeff@hv+aa*th+bb*dh;ep=coeff@ev+aa*te+bb*de
    stats['real_isolated_branches']+=len(take);counts['branches']+=len(take)
    keep=np.flatnonzero((np.max(vals,axis=1)<1e-6)&(hp<-1e-7)&(ep>1e-7))
    counts['acute_escape_cutters']+=len(keep)
    for k in keep:
     ps=vals[k];finite=(ps>-2+1e-7)&(ps<-1e-7);ml=np.pi/np.arccos(np.clip(-ps/2,-1+1e-15,1-1e-15));bad=np.flatnonzero(finite&(abs(ml-np.rint(ml))>1e-5)).tolist();parallel=np.flatnonzero(abs(ps+2)<1e-7).tolist()
     rec=dict(support=S,contacts=profile,anchor=int(anchors[take[k]]),anchor_label=int(labels[take[k]]),branch=sign,escape_pairing=smooth(ep[k]),center_pairing=smooth(hp[k]),normalized_wall_pairings=[smooth(x)for x in ps],bad_finite=[[i,smooth(ml[i])]for i in bad],near_parallel=parallel)
     if not bad and not parallel:
      key=tuple(np.round(ps,7));accepted.setdefault(key,rec)
     elif len(bad)+len(parallel)<=2:near.append(rec)
  supports.append(dict(support=S,positive_definite_numeric=True,**counts));print('support',S,dict(counts),'accepted',len(accepted),'near',len(near),flush=True)
 bank=sorted(accepted.values(),key=lambda r:-r['escape_pairing']);near=sorted(near,key=lambda r:(len(r['bad_finite'])+len(r['near_parallel']),-r['escape_pairing']))
 return dict(status='BOUNDED_GENERAL_BRIDGE_CUTTER_SCOUT',sources={p:sha(p)for p in[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py',P.BATCH,'results/fresh_nonarithmetic_affine_filtered_20261016_cone.json','results/fresh_nonarithmetic_affine_filtered_20261016_cluster.json']},escape_active=active,profile_labels=PROFILE_LABELS,anchor_labels=ANCHOR_LABELS,max_nonzero_support_contacts=2,counts=dict(stats),support_results=supports,candidate_count=len(bank),candidates=bank,near_miss_count=len(near),near_misses=near[:100],scope='Numerical only: each positive-definite8-subset of the nine renewed-escape active walls; at most two nonright support contacts at listed labels; one of the three newcaps as a finite anchor at labels2..20; both signs of Schur residual are retained; both isolated norm2 conic branches. Near-null anchor projections and discriminant<=1e-10 are omitted. All52wall gates are numerical; failed finite fits are unresolved, not all-label exclusions. No finite-volume or complete-family claim.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved target scout replay')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['candidate_count'],d['near_miss_count'])
