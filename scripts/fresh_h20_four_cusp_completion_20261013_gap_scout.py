#!/usr/bin/env python3
"""Three-anchor proposals in the residualH3 at the saved remote gap.

Numerical scouting only. Seventeen fixed short walls give a PD17 block,
leaving a4-dimensional Lorentzian complement; three finite anchors and
the norm equation propose both remaining branches.
"""
import argparse,json,itertools
import numpy as np
from scipy.linalg import null_space
from fresh_h20_four_cusp_completion_20261013 import state,OUT as SOURCE,ROOT,read,vec,B
from fresh_h20_coupled_completion_20261013 import num
SHORT=[0,1,2,3,4,5,8,9,13,14,15,20,21,22,24,25,26]
OUT='data/enumerations/fresh_h20_four_cusp_completion_20261013_gap_scout.json'

def run(wide=False,nextgap=False,drop=None):
 Q,h,hc,roots=state();roots.extend((B,vec(c['normal']))for c in read(SOURCE)['new_caps']);q=np.array([[float(x.a)for x in row]for row in Q]);R=np.array([[num(f,x)for x in r]for f,r in roots]);ns=np.einsum('ij,jk,ik->i',R,q,R)
 hn=np.array([float(x.a)for x in h]);target=np.array([float(x.a)for x in vec(read(SOURCE)['strict_escape']['vector'])])
 short=SHORT
 G=R[short]@q@R[short].T;assert np.linalg.eigvalsh(G)[0]>1e-6
 U=null_space(R[short]@q);ev,T=np.linalg.eigh(U.T@q@U);assert np.sum(ev<0)==1
 V=U@T/np.sqrt(abs(ev));J=np.diag(np.sign(ev));P=R@q@V;C=P/np.sqrt(ns[:,None]);ct=target@q@V;ch=hn@q@V
 anchors=[i for i in range(27)if i not in short]+[33,49,51,55]
 if wide:
  seen=set();anchors=[]
  for i in range(71):
   c=C[i];nn=np.linalg.norm(c)
   if nn<1e-8:continue
   key=tuple(np.round(c/nn,7))
   if key not in seen:seen.add(key);anchors.append(i)
  anchors+=[173,174,175,176,177,178,214,215,216,217,218,219,220,221,222,223]
 nanchor=20-len(short)
 labels=[2,3,4,5,6,8,10,12,0];cs=np.array([0 if m==2 else 1 if m==0 else np.cos(np.pi/m)for m in labels]);profiles=np.array(list(itertools.product(range(len(labels)),repeat=nanchor)));rhs=-cs[profiles]
 allowed=np.array([0]+[np.cos(np.pi/m)for m in range(3,121)])
 found=[];seen=set();tested=0;realbranches=0;cutbranches=0;acute=0
 for abc in itertools.combinations(anchors,nanchor):
  A=C[list(abc)]
  if np.linalg.matrix_rank(A,tol=1e-8)!=nanchor:continue
  K=A.T@np.linalg.inv(A@A.T);a=rhs@K.T;v=null_space(A)[:,0];qa=v@J@v
  if abs(qa)<1e-9:continue
  qb=2*a@J@v;qc=np.einsum('ij,jk,ik->i',a,J,a)-1;dis=qb*qb-4*qa*qc;ok=np.flatnonzero(dis>1e-9);tested+=len(profiles);realbranches+=2*len(ok)
  if not len(ok):continue
  for branch in(-1,1):
   t=(-qb[ok]+branch*np.sqrt(dis[ok]))/(2*qa);Y=a[ok]+t[:,None]*v
   mask=(Y@ct>1e-7)&(Y@ch<0);Y=Y[mask];pp=ok[mask];cutbranches+=len(Y)
   if not len(Y):continue
   z=-(Y@C.T);mask=np.max(-z,axis=1)<1e-7;Y=Y[mask];z=z[mask];pp=pp[mask];acute+=len(Y)
   if not len(Y):continue
   good=np.ones(len(Y),dtype=bool)
   for col in z.T:
    good&=(col>=1-1e-7)|(np.min(abs(col[:,None]-allowed),axis=1)<1e-7)
   for yy,profile in zip(Y[good],pp[good],strict=True):
    r=V@yy;key=tuple(np.round(r,7))
    if key in seen:continue
    seen.add(key);row={'anchors':list(abc),'labels':[labels[i]for i in profiles[profile]],'branch':branch,'normal':r.tolist(),'norm':float(r@q@r),'target_pairing':float(r@q@target),'center_pairing':float(r@q@hn),'finite_neighbors':np.flatnonzero(-R@q@r/np.sqrt(ns)<1-1e-7).tolist()};found.append(row)
    if len(found)<=10 or len(found)%100==0:print('FOUND',len(found),abc,row['labels'],row['target_pairing'],flush=True)
 out={'status':'NUMERICAL_SHORT_BLOCK_FINITE_ANCHOR_GAP_PROPOSALS','source':SOURCE,'short_block':short,'anchor_indices':anchors,'labels':labels,'profiles_tested':tested,'real_branches':realbranches,'target_cut_branches':cutbranches,'acute_cut_branches':acute,'new_candidates':found,'scope':'Bounded floating proposals at the224-wall escape, outside the nine fixed highest-root families; no exact compatibility or completeness claim.'}
 path=OUT.replace('.json','_next'+(''if drop is None else'_drop'+str(drop))+'.json')if nextgap else OUT.replace('.json','_wide.json')if wide else OUT;(ROOT/path).write_text(json.dumps(out,indent=2)+'\n');print('DONE',tested,realbranches,cutbranches,acute,len(found),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--wide',action='store_true');p.add_argument('--nextgap',action='store_true');p.add_argument('--drop',type=int);a=p.parse_args();run(a.wide,a.nextgap,a.drop)
