#!/usr/bin/env python3
"""Two-anchor polygon-side scout on the common S6+30 two-face."""
import argparse,itertools,json,math
from collections import Counter
from hashlib import sha256
import mpmath as mp
import numpy as np
import fresh_nonarithmetic_h3_closure_20261016_scout as H3
P=H3.P;ROOT=P.ROOT;TAG='fresh_nonarithmetic_h3_closure_20261016_h2_scout';WHITE=H3.S+[30];LABELS=list(range(2,21))
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def smooth(x):return float(format(float(x),'.13g'))
def geometry():
 Q,N,H,E=P.geometry();cap=P.P.B.read(H3.CAP);N.append(mp.matrix([P.mid(x)for x in cap['normal_interval']]));W=mp.matrix([[N[i][j]for i in WHITE]for j in range(10)]);J=(W.T*Q*W)**-1;proj=lambda n:n-W*J*W.T*Q*n
 F=mp.matrix([[proj(N[i])[j]for i in[7,33,44]]for j in range(10)]);G=F.T*Q*F;Gi=G**-1;coords=lambda n:Gi*F.T*Q*n;Gn=np.array(G.tolist(),dtype=float);val,vec=np.linalg.eigh(Gn);assert sum(val<0)==1;B=vec@np.diag(1/np.sqrt(abs(val)));Bi=np.linalg.inv(B);g=np.diag(np.sign(val));to=lambda n:Bi@np.array(list(coords(n)),dtype=float)
 sc=P.P.B.read(f'results/{H3.TAG}.json');V=[mp.matrix([mp.mpf(x)for x in r['point']])for r in sc['proposals']]
 return np.array([to(n)for n in N]),g,to(H),np.array([to(v)for v in V]),sc
def run():
 N,g,H,V,sc=geometry();anchors=[i for i in range(len(N))if i not in WHITE and np.linalg.norm(N[i])>1e-10];labels=list(itertools.product(LABELS,repeat=2));cv=np.array([[-2*math.cos(math.pi/m)if m!=2 else 0 for m in row]for row in labels]);stats=Counter();out={};near=[];domains=[]
 for ai,(a,b)in enumerate(itertools.combinations(anchors,2)):
  W=N[[a,b]].T;G=W.T@g@W
  if abs(np.linalg.det(G))<1e-8*max(1,np.max(abs(G))**2):stats['degenerate_anchor_pairs']+=1;continue
  J=np.linalg.inv(G);k=np.cross(g@W[:,0],g@W[:,1]);q=k@g@k
  if abs(q)<1e-11*(k@k):stats['near_null_kernel']+=1;continue
  k/=math.sqrt(abs(q));ks=1 if q>0 else-1;coef=cv@J;p=coef@W.T;rho=(2-np.einsum('ij,jk,ik->i',p,g,p))/ks;ix=np.flatnonzero(rho>1e-9);stats['anchor_pairs']+=1;stats['recipes']+=len(labels);hit=0
  for branch in[-1,1]:
   R=p[ix]+branch*np.sqrt(rho[ix])[:,None]*k;ps=R@g@N.T;hp=R@g@H;cover=R@g@V.T;stats['real_branches']+=len(ix);keep=np.flatnonzero((np.max(ps,axis=1)<1e-6)&(hp<-1e-7)&(np.max(cover,axis=1)>1e-6));stats['acute_cutters']+=len(keep)
   for j in keep:
    row=ps[j];finite=(row>-2+1e-7)&(row<-1e-7);ml=np.pi/np.arccos(np.clip(-row/2,-1+1e-15,1-1e-15));bad=np.flatnonzero(finite&(abs(ml-np.rint(ml))>1e-5)).tolist();par=np.flatnonzero(abs(row+2)<1e-7).tolist();r=dict(anchors=[a,b],labels=list(labels[ix[j]]),branch=branch,pairings=[smooth(x)for x in row],cut_original_H3_rays=np.flatnonzero(cover[j]>1e-6).tolist(),H3_pairings=[smooth(x)for x in cover[j]],bad=[[i,smooth(ml[i])]for i in bad],near_parallel=par)
    if not bad and not par:out.setdefault(tuple(np.round(row,6)),r);hit+=1
    elif len(bad)+len(par)<=1:near.append(r)
  if hit:domains.append(dict(anchors=[a,b],proposals=hit))
  if ai%100==0:print('anchor pair',ai,dict(stats),'proposals',len(out),'near',len(near),flush=True)
 paths=[f'scripts/{TAG}.py',H3.CAP,H3.CONE,f'results/{H3.TAG}.json',f'scripts/{H3.TAG}.py',P.CL,P.CO,P.P.BATCH]
 return dict(status='BOUNDED_TWO_ANCHOR_S6_PLUS30_POLYGON_SIDE_SCOUT',sources={p:sha(p)for p in paths},white_support=WHITE,anchor_labels=LABELS,stats=dict(stats),accepted_domains=domains,candidates=list(out.values()),near_misses=near,scope='Numerical discovery only. Every independent pair of projected old53 normals and every assigned pair of labels2..20 is tested, keeping both real norm2 branches, old acuteness, center orientation and all numerical Coxeter gates. The right seven-root support is fixed. Higher labels, one-contact continuous families and nonright white profiles are outside this search; no exact or finite-volume claim is made.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved H2 side scout')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['candidates']))
