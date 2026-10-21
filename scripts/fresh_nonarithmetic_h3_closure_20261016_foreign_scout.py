#!/usr/bin/env python3
"""Residual-H3 three-anchor conics, keeping every ambient wall gate."""
import argparse,itertools,json,math
from collections import Counter
from hashlib import sha256
from pathlib import Path
import mpmath as mp
import numpy as np
import fresh_nonarithmetic_h3_closure_20261016_scout as H3
P=H3.P;ROOT=P.ROOT;TAG='fresh_nonarithmetic_h3_closure_20261016_foreign_scout';INPUT=f'results/{H3.TAG}.json';LABELS=[2,3,4,5,6,7,8,10,12]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def smooth(v):return float(format(float(v),'.13g'))
def run():
 Q,N,H,E=P.geometry();cap=P.P.B.read(H3.CAP);sc=P.P.B.read(INPUT);N.append(mp.matrix([P.mid(x)for x in cap['normal_interval']]))
 polars=[1,3,6]
 N.extend(mp.matrix([mp.mpf(x)for x in sc['proposals'][i]['point']])for i in polars)
 V=[mp.matrix([mp.mpf(x)for x in row['point']])for row in sc['proposals']]
 G=mp.matrix([[(n.T*Q*m)[0]for m in N]for n in N]);hh=mp.matrix([(H.T*Q*n)[0]for n in N]);stats=Counter();rows=[];accepted={};near=[]
 target=5;incident=[30,47,52];facewalls=sorted({i for rows in sc['face_census'].values()for ids in rows for i in ids if i not in H3.S}|{53,54,55})
 triples=sorted({tuple(sorted(pair+(j,)))for pair in itertools.combinations(incident,2)for j in facewalls if j not in incident})
 for anchors in triples:
  frame=H3.S+list(anchors);M=mp.matrix([[G[i,j]for j in frame]for i in frame]);
  if abs(mp.det(M))<mp.mpf('1e-30'):stats['unresolved_singular_frames']+=1;continue
  J=M**-1;W=mp.matrix([[N[i][j]for i in frame]for j in range(10)]);k=H-W*J*mp.matrix([hh[i]for i in frame]);qk=(k.T*Q*k)[0];assert abs(qk)>mp.mpf('1e-30');ks=1 if qk>0 else-1;k/=mp.sqrt(abs(qk))
  JJ=np.array(J.tolist(),dtype=float);cross=np.array([[float(G[i,j])for j in range(len(N))]for i in frame]);kp=np.array([float((k.T*Q*n)[0])for n in N]);kp[frame]=0;hv=np.array([float(hh[i])for i in frame]);kh=float((k.T*Q*H)[0]);rv=np.array([[float((N[i].T*Q*v)[0])for v in V]for i in frame]);kr=np.array([float((k.T*Q*v)[0])for v in V]);counts=Counter()
  labels=list(itertools.product(LABELS,repeat=3));cv=np.zeros((len(labels),9));cv[:,6:]=np.array([[-2*math.cos(math.pi/m)if m!=2 else 0 for m in row]for row in labels]);coef=cv@JJ;rho=(2-np.sum(cv*coef,axis=1))/ks;take=np.flatnonzero(rho>1e-10);counts['recipes']=len(labels)
  for sign in[-1,1]:
   lam=sign*np.sqrt(rho[take]);vals=coef[take]@cross+lam[:,None]*kp;hp=coef[take]@hv+lam*kh;cover=coef[take]@rv+lam[:,None]*kr;counts['real_branches']+=len(take)
   keep=np.flatnonzero((np.max(vals,axis=1)<1e-6)&(hp<-1e-7)&(cover[:,target]>1e-7));counts['acute_target_cutters']+=len(keep)
   for j in keep:
    ps=vals[j];finite=(ps>-2+1e-7)&(ps<-1e-7);ml=np.pi/np.arccos(np.clip(-ps/2,-1+1e-15,1-1e-15));bad=np.flatnonzero(finite&(abs(ml-np.rint(ml))>1e-5)).tolist();par=np.flatnonzero(abs(ps+2)<1e-7).tolist();r=dict(target=target,anchors=anchors,labels=list(labels[take[j]]),branch=sign,pairings=[smooth(x)for x in ps],cut_old_face_rays=np.flatnonzero(cover[j]>1e-7).tolist(),old_face_pairings=[smooth(x)for x in cover[j]],center_pairing=smooth(hp[j]),bad=[[i,smooth(ml[i])]for i in bad],near_parallel=par)
    if not bad and not par:accepted.setdefault(tuple(np.round(ps,7)),r)
    elif len(bad)+len(par)<=2:near.append(r)
  stats.update(counts);rows.append(dict(target=target,anchors=anchors,**counts));print('target',target,dict(counts),'accepted',len(accepted),'near',len(near),flush=True)
 paths=[f'scripts/{TAG}.py',f'scripts/{H3.TAG}.py',INPUT,H3.CAP,H3.CONE,P.CL,P.CO,P.P.BATCH]
 return dict(status='BOUNDED_RESIDUAL_H3_FOREIGN_ANCHOR_CUTTER_SCOUT',sources={p:sha(p)for p in paths},right_white_support=H3.S,assumed_extra_polars=polars,anchor_labels=LABELS,stats=dict(stats),targets=rows,candidates=list(accepted.values()),near_miss_count=len(near),near_misses=near[:100],scope='Numerical discovery only. The remaining exterior ray5 supplies two incident anchors, with a third anchor chosen from all current face walls and the three proposed polars; every listed label and both isolated norm2 branches are tested. All56 ambient wall gates include three proposed right polars. These polars and surviving candidates need exactification; higher labels, other anchors and changes in white contacts are outside this bounded search.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved residualH3 conic scout')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['candidates']))
