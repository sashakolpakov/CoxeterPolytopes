#!/usr/bin/env python3
"""Bounded three-anchor repairs of all six replacement H3 exterior rays."""
import argparse,itertools,json,math
from collections import Counter
from hashlib import sha256
import mpmath as mp
import numpy as np
import fresh_nonarithmetic_h3_closure_20261016_scout as H3
P=H3.P;ROOT=P.ROOT;TAG='fresh_nonarithmetic_h3_closure_20261016_repair_scout'
TRIAL='results/fresh_nonarithmetic_h3_closure_20261016_trial.json'
LABELS=list(range(2,21))
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def smooth(x):return float(format(float(x),'.13g'))
def run():
 Q,N,H,E=P.geometry();cap=P.P.B.read(H3.CAP);N.append(mp.matrix([P.mid(x)for x in cap['normal_interval']]));tr=P.P.B.read(TRIAL)
 N.extend(mp.matrix([mp.mpf(x)for x in row])for row in tr['normal_intervals_or_proposals'])
 V=[mp.matrix([mp.mpf(x)for x in row['point']])for row in tr['rays']if row['kind']=='hyperideal_rays'];targets=[row for row in tr['rays']if row['kind']=='hyperideal_rays']
 G=mp.matrix([[(n.T*Q*m)[0]for m in N]for n in N]);hh=mp.matrix([(H.T*Q*n)[0]for n in N]);stats=Counter();records=[];accepted={};near=[]
 for target,row in enumerate(targets):
  anchors=[i for i in row['active']if i not in H3.S];assert len(anchors)==3;frame=H3.S+anchors;M=mp.matrix([[G[i,j]for j in frame]for i in frame]);J=M**-1;W=mp.matrix([[N[i][j]for i in frame]for j in range(10)]);k=H-W*J*mp.matrix([hh[i]for i in frame]);qk=(k.T*Q*k)[0];assert abs(qk)>mp.mpf('1e-30');ks=1 if qk>0 else-1;k/=mp.sqrt(abs(qk))
  JJ=np.array(J.tolist(),dtype=float);cross=np.array([[float(G[i,j])for j in range(len(N))]for i in frame]);kp=np.array([float((k.T*Q*n)[0])for n in N]);kp[frame]=0;hv=np.array([float(hh[i])for i in frame]);kh=float((k.T*Q*H)[0]);rv=np.array([[float((N[i].T*Q*v)[0])for v in V]for i in frame]);kr=np.array([float((k.T*Q*v)[0])for v in V]);counts=Counter()
  labels=list(itertools.product(LABELS,repeat=3));cv=np.zeros((len(labels),9));cv[:,6:]=np.array([[-2*math.cos(math.pi/m)if m!=2 else 0 for m in labs]for labs in labels]);coef=cv@JJ;rho=(2-np.sum(cv*coef,axis=1))/ks;take=np.flatnonzero(rho>1e-10);counts['recipes']=len(labels)
  for branch in[-1,1]:
   lam=branch*np.sqrt(rho[take]);vals=coef[take]@cross+lam[:,None]*kp;hp=coef[take]@hv+lam*kh;cover=coef[take]@rv+lam[:,None]*kr;counts['real_branches']+=len(take)
   keep=np.flatnonzero((np.max(vals,axis=1)<1e-6)&(hp<-1e-7)&(cover[:,target]>1e-8));counts['acute_target_cutters']+=len(keep)
   for j in keep:
    ps=vals[j];finite=(ps>-2+1e-7)&(ps<-1e-7);ml=np.pi/np.arccos(np.clip(-ps/2,-1+1e-15,1-1e-15));bad=np.flatnonzero(finite&(abs(ml-np.rint(ml))>1e-5)).tolist();par=np.flatnonzero(abs(ps+2)<1e-7).tolist();r=dict(target=target,anchors=anchors,labels=list(labels[take[j]]),branch=branch,pairings=[smooth(x)for x in ps],cut_face_rays=np.flatnonzero(cover[j]>1e-8).tolist(),face_pairings=[smooth(x)for x in cover[j]],center_pairing=smooth(hp[j]),bad=[[i,smooth(ml[i])]for i in bad],near_parallel=par)
    if not bad and not par:accepted.setdefault(tuple(np.round(ps,7)),r)
    elif len(bad)+len(par)<=1:near.append(r)
  stats.update(counts);records.append(dict(target=target,anchors=anchors,**counts));print('target',target,anchors,dict(counts),'accepted',len(accepted),'near',len(near),flush=True)
 paths=[f'scripts/{TAG}.py',f'scripts/{H3.TAG}.py',TRIAL,H3.CAP,H3.CONE,P.CL,P.CO,P.P.BATCH]
 return dict(status='BOUNDED_RESIDUAL_H3_REPLACEMENT_RAY_REPAIR_SCOUT',sources={p:sha(p)for p in paths},right_white_support=H3.S,anchor_labels=LABELS,stats=dict(stats),targets=records,candidates=list(accepted.values()),near_miss_count=len(near),near_misses=near,scope='Numerical discovery only. Every triple of labels2..20 on the three incident anchors of each of six replacement H3 exterior rays, both isolated norm2 branches, and all59 ambient wall gates are tested. No finite-volume, exact-angle, no-affine, or all-label claim is made.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved replacement-ray repair scout')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['candidates']))
