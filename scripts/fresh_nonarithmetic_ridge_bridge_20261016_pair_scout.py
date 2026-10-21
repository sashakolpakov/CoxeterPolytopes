#!/usr/bin/env python3
"""Two finite anchors and seven-root supports: bounded bridge-cutter discovery."""
import argparse,json,itertools,math
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import mpmath as mp
import numpy as np
import fresh_nonarithmetic_escape49_20261016_polar as P
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_ridge_bridge_20261016_pair_scout';CL='results/fresh_nonarithmetic_affine_filtered_20261016_cluster.json';CO='results/fresh_nonarithmetic_affine_filtered_20261016_cone.json'
PROFILE=[3,4,5,6,10];L33=[3,4,5,6,7,8,10,12,20];LCAP=[2,3,4,5,6,7,8,10,12,20]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def smooth(x):return float(format(float(x),'.13g'))
def run():
 old=P.B.read(P.BATCH);cl=P.B.read(CL);co=P.B.read(CO);mp.mp.dps=80
 Q=mp.matrix([[P.scalar(P.D.encfield(P.D.lift(P.A.B.db(x))))for x in row]for row in old['form']]);N=[mp.matrix([P.scalar(x)for x in row])for row in old['normals']]
 def mid(x):
  q=(Fraction(x[0])+Fraction(x[1]))/2;return mp.mpf(q.numerator)/q.denominator
 N += [mp.matrix([mid(x)for x in r['normal_interval']])for r in cl['caps']];N=[n*mp.sqrt(2/(n.T*Q*n)[0])for n in N]
 H=mp.matrix([P.scalar(x)for x in old['center']]);E=mp.matrix([P.scalar(x)for x in co['strict_escape']['point']]);E/=-(E.T*Q*H)[0]
 G=mp.matrix([[(n.T*Q*m)[0]for m in N]for n in N]);h=mp.matrix([(H.T*Q*n)[0]for n in N]);e=mp.matrix([(E.T*Q*n)[0]for n in N]);active=co['strict_escape']['active'];stats=Counter();records=[];accepted={};near=[]
 for S in itertools.combinations(active,7):
  GS=mp.matrix([[G[i,j]for j in S]for i in S])
  try:mp.cholesky(GS)
  except ValueError:stats['nonPD7supports']+=1;continue
  profiles=[[]]+[[(i,m)]for i in S for m in PROFILE]
  for cap in[49,50,51]:
   ids=list(S)+[33,cap];M=mp.matrix([[G[i,j]for j in ids]for i in ids]);det=mp.det(M)
   if abs(det)<mp.mpf('1e-35'):stats['unresolved_singular9']+=1;continue
   J=M**-1;W=mp.matrix([[N[i][j]for i in ids]for j in range(10)]);k=H-W*J*mp.matrix([h[i]for i in ids]);qk=(k.T*Q*k)[0]
   if abs(qk)<mp.mpf('1e-35'):stats['unresolved_null_kernel']+=1;continue
   k/=mp.sqrt(abs(qk));ks=1 if qk>0 else-1
   JJ=np.array([[float(x)for x in row]for row in J.tolist()]);cross=np.array([[float(G[i,j])for j in range(52)]for i in ids]);kp=np.array([float((k.T*Q*n)[0])for n in N]);kp[ids]=0;kh=float((k.T*Q*H)[0]);ke=float((k.T*Q*E)[0]);hv=np.array([float(h[i])for i in ids]);ev=np.array([float(e[i])for i in ids]);counts=Counter()
   for profile in profiles:
    cc=[];labels=[]
    for m,n in itertools.product(L33,LCAP):
     c=np.zeros(9)
     for i,l in profile:c[ids.index(i)]=-2*math.cos(math.pi/l)
     c[7]=-2*math.cos(math.pi/m);c[8]=-2*math.cos(math.pi/n)if n!=2 else 0
     cc.append(c);labels.append((m,n))
    cv=np.array(cc);coef=cv@JJ;rho=(2-np.sum(cv*coef,axis=1))/ks;take=np.flatnonzero(rho>1e-9);stats['contact_recipes']+=len(cv);counts['contact_recipes']+=len(cv)
    for sign in[-1,1]:
     lam=sign*np.sqrt(rho[take]);ps=coef[take]@cross+lam[:,None]*kp;hp=coef[take]@hv+lam*kh;ep=coef[take]@ev+lam*ke;stats['real_branches']+=len(take);counts['real_branches']+=len(take)
     keep=np.flatnonzero((hp<-1e-7)&(ep>1e-7)&(np.max(ps,axis=1)<1e-6));counts['acute_cutters']+=len(keep)
     for t in keep:
      vals=ps[t];finite=(vals>-2+1e-7)&(vals<-1e-7);ml=np.pi/np.arccos(np.clip(-vals/2,-1+1e-15,1-1e-15));bad=np.flatnonzero(finite&(abs(ml-np.rint(ml))>1e-5)).tolist();par=np.flatnonzero(abs(vals+2)<1e-7).tolist();m,n=labels[take[t]]
      row=dict(support=list(S),contacts=profile,anchors=[33,cap],anchor_labels=[m,n],branch=sign,escape_pairing=smooth(ep[t]),center_pairing=smooth(hp[t]),normalized_wall_pairings=[smooth(x)for x in vals],bad_finite=[[i,smooth(ml[i])]for i in bad],near_parallel=par)
      if not bad and not par:accepted.setdefault(tuple(np.round(vals,7)),row)
      elif len(bad)+len(par)<=2:near.append(row)
   records.append(dict(support=list(S),anchors=[33,cap],**counts));print('support',S,'cap',cap,dict(counts),'accepted',len(accepted),flush=True)
 paths=[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py',P.BATCH,CL,CO]
 return dict(status='BOUNDED_TWO_ANCHOR_BRIDGE_CUTTER_DISCOVERY',sources={p:sha(p)for p in paths},active=active,profile_labels=PROFILE,wall33_labels=L33,newcap_labels=LCAP,stats=dict(stats),support_records=records,candidates=sorted(accepted.values(),key=lambda r:-r['escape_pairing']),near_miss_count=len(near),near_misses=sorted(near,key=lambda r:(len(r['bad_finite'])+len(r['near_parallel']),-r['escape_pairing']))[:100],scope='Numerical discovery only: PD7subsets of the renewed active9, at most one nonright support contact, two finite anchors33 and one of49..51 at listed labels, both isolated norm branches. Singular9frames and tiny discriminants remain untested. No all-label exclusion or exact cap acceptance.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved two-anchor scout')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['candidates']),d['near_miss_count'])
