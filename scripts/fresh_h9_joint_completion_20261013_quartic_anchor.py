#!/usr/bin/env python3
"""Four free seed distances, a quartic anchor and a Q(phi) anchor.

Splitting the quartic anchor into its rational and radical coefficient rows
gives two Q(phi)-linear conditions. A third anchor condition leaves a line,
so the norm equation remains quadratic over Q(phi).
"""
import argparse,heapq,json
from itertools import combinations,product
import numpy as np
from scripts.fresh_h9_joint_completion_20261013 import setup,ROOT,read
from scripts.fresh_h8_cap_completion_20261013 import gate
from scripts.fresh_h8_alternate_caps_20261013 import *
from scripts.fresh_h8_release_h9_free_caps_20261013 import esign,emul,esub,eencode
OUTPUT=ROOT/'data/enumerations/fresh_h9_joint_completion_20261013_quartic_anchor.json'

def run(min_retained=18,both=True,dotted=False,quartic_only=False):
 G,I,center,nodes,polys,cn,In,hn,original=setup();quartic=1039;chosen=sorted([i for i in original if i!=1038]+[quartic]);anchors=[i for i in chosen if nodes[i]['d']==Z];a0=nodes[quartic]['a'];a1=nodes[quartic]['b'];assert nodes[quartic]['d']!=Z;mags=(Z,O,P,scale(P,2),add(O,scale(P,2)),add(scale(O,2),P))if dotted else(Z,O,P);cases=[(None,m,0)for m in range(len(mags))]if quartic_only else[(k,m,n)for k in anchors for m,n in product(range(len(mags)),range(3))];rows=np.array([[[real(x)for x in row]for row in((a0,a1)if quartic_only else(a0,a1,nodes[k]['a']))]for k,m,n in cases]);mag=tuple(real(x)for x in mags);rhs=np.array([[mag[m],0.]if quartic_only else[mag[m],0.,mag[n]]for k,m,n in cases]);target=read('data/enumerations/fresh_h9_joint_completion_20261013_targets.json');ys=[tuple(map(parse,r['certificate']['y']))for r in target['targets']];Y=np.array([[real(x)for x in y]for y in ys]);IY=In@Y.T;IC=In@np.array([real(x)for x in center]);selected_cn=cn[chosen];selected_norms=np.array([real(nodes[i]['q'])for i in chosen]);allowed=np.array([float((lo+hi)/8)for _,_,lo,hi in polys]);patterns=[(supp,labs)for k in range(4)for supp in combinations(range(7 if quartic_only else 6),k)for labs in product((0,1),repeat=k)];heaps={};serial=0;stats={'free_sets':0,'free_distances':3 if quartic_only else 4,'anchor_line_candidates':0,'joint_cut_candidates':0,'retention_candidates':0,'best_retention':0,'targets_required':[5,22]if both else'either5or22'}
 for free in combinations(range(10),3 if quartic_only else 4):
  stats['free_sets']+=1;remaining=[i for i in range(10)if i not in free];assigned=np.zeros(len(cases),dtype=bool)
  for u in free:
   solved=tuple(i for i in free if i!=u);matrices=rows[:,:,solved];det=np.linalg.det(matrices);active=(abs(det)>1e-9)&~assigned;assigned|=active
   if not active.any():continue
   ids=np.flatnonzero(active);r=rows[ids];inv=np.linalg.inv(matrices[ids]);w1=np.zeros((len(ids),10));w1[:,u]=1;w1[:,solved]=-np.einsum('bij,bj->bi',inv,r[:,:,u]);A=np.einsum('bi,ij,bj->b',w1,In,w1)
   for positions,labels in patterns:
    f=np.zeros(10)
    for j,l in zip(positions,labels):f[remaining[j]]=(1.,real(P))[l]
    rr=rhs[ids]-np.einsum('bij,j->bi',r,f);sol=np.einsum('bij,bj->bi',inv,rr);w0=np.broadcast_to(f,(len(ids),10)).copy();w0[:,solved]=sol;B=2*np.einsum('bi,ij,bj->b',w1,In,w0);C=np.einsum('bi,ij,bj->b',w0,In,w0)-2;D=B*B-4*A*C;good=(D>1e-9)&(abs(A)>1e-9)
    if not good.any():continue
    for branch in(-1,1):
     x=np.full(len(ids),np.nan);x[good]=(-B[good]+branch*np.sqrt(D[good]))/(2*A[good]);ws=w0+x[:,None]*w1;cuts=ws@IY;nc=(cuts>1e-7).sum(axis=1);ok=good&np.all(ws[:,list(free)]>2+1e-7,axis=1)&(ws@IC< -1e-7);stats['anchor_line_candidates']+=int(ok.sum());focus=(cuts[:,5]>1e-7)&(cuts[:,22]>1e-7)if both else(cuts[:,5]>1e-7)|(cuts[:,22]>1e-7);ok&=focus&(nc>=2);stats['joint_cut_candidates']+=int(ok.sum());idx=np.flatnonzero(ok)
     if not len(idx):continue
     ps=-ws[idx]@selected_cn.T;rat=ps*ps/(2*selected_norms);near=np.min(np.abs(rat[:,:,None]-allowed),axis=2)<1e-6;compatible=(ps<=1e-7)&((rat>1+1e-7)|near);nr=compatible.sum(axis=1);stats['best_retention']=max(stats['best_retention'],int(nr.max()))
     for jj in np.flatnonzero(nr>=min_retained):
      j=int(idx[jj]);stats['retention_candidates']+=1;k,m,n=cases[int(ids[j])];mask=tuple(chosen[z]for z in range(len(chosen))if compatible[jj,z]);length=np.linalg.norm(ws[j]@In);score=(int(nc[j]),int(nr[jj]),float((min(cuts[j,5],cuts[j,22])if both else max(cuts[j,5],cuts[j,22]))/length));key=k,mask;serial+=1;item=(score,serial,(u,solved,tuple(remaining[z]for z in positions),labels,k,m,n,branch));h=heaps.setdefault(key,[])
      if len(h)<3:heapq.heappush(h,item)
      elif score>h[0][0]:heapq.heapreplace(h,item)
  if stats['free_sets']%30==0:print(stats,flush=True)
 caps=[];seen=set()
 for score,_,spec in sorted([t for h in heaps.values()for t in h],reverse=True):
  u,solved,support,labels,k,m,n,branch=spec;f=[Z]*10
  for j,l in zip(support,labels):f[j]=(O,P)[l]
  f=tuple(f);ar=(a0,a1)if quartic_only else(a0,a1,nodes[k]['a']);M,N=mags[m],(Z,O,P)[n];targetrhs=(M,Z)if quartic_only else(M,Z,N);mat=tuple(tuple(row[j]for j in solved)for row in ar);inv=inverse(mat);constants=mv(inv,tuple(sub(t,dot(row,f))for t,row in zip(targetrhs,ar)));directions=mv(inv,tuple(neg(row[u])for row in ar));w0=list(f);w1=[Z]*10;w1[u]=O
  for j,v,t in zip(solved,constants,directions):w0[j]=v;w1[j]=t
  w0,w1=tuple(w0),tuple(w1);assert all(dot(row,w0)==t and dot(row,w1)==Z for row,t in zip(ar,targetrhs));A=pair(w1,I,w1);B=scale(pair(w0,I,w1),2);C=sub(pair(w0,I,w0),scale(O,2));D=sub(mul(B,B),scale(mul(A,C),4));assert sign(D)>0;x=(div(neg(B),scale(A,2)),div(scale(O,branch),scale(A,2)));pa=tuple(neg(add(a,mul(b,x[0])))for a,b in zip(w0,w1));pb=tuple(neg(mul(b,x[1]))for b in w1);ca,cb=mv(I,pa),mv(I,pb);key=(D,ca,cb)
  if key in seen:continue
  seen.add(key);nd={'q':scale(O,2),'d':D,'a':ca,'b':cb,'pa':pa,'pb':pb};norm=(add(dot(pa,ca),mul(D,dot(pb,cb))),scale(dot(pa,cb),2));assert esign(esub(norm,(scale(O,2),Z)),D)==0
  for j in range(10):
   assert esign((pa[j],pb[j]),D)<=0;sq=emul((pa[j],pb[j]),(pa[j],pb[j]),D)
   if j in(u,*solved):assert esign(esub(sq,(scale(O,4),Z)),D)>0
   else:assert any(esign(esub(sq,(z,Z)),D)==0 for z in(Z,O,add(O,P)))
  assert esign((dot(ca,center),dot(cb,center)),D)>0;checks=[{'old_node':i,**gate(nd,nodes[i],polys)}for i in chosen];kept=[r['old_node']for r in checks if r['compatible']];dropped=[r['old_node']for r in checks if not r['compatible']];cut=[j for j,y in enumerate(ys)if esign((dot(ca,y),dot(cb,y)),D)<0]
  if len(kept)<min_retained or quartic not in kept or(not quartic_only and k not in kept):continue
  if both:assert 5 in cut and 22 in cut
  else:assert 5 in cut or 22 in cut
  caps.append({'free_seed_indices':[u,*solved],'fixed_pairing_magnitudes':vec(f),'quartic_anchor_node':quartic,'rational_anchor_node':k,'anchor_pairing_magnitudes':vec([M]if quartic_only else[M,N]),'line_constant':vec(w0),'line_direction':vec(w1),'norm_polynomial_descending':vec([A,B,C]),'discriminant':enc(D),'parameter':eencode(x),'pairings_rational':vec(pa),'pairings_radical':vec(pb),'coefficients_rational':vec(ca),'coefficients_radical':vec(cb),'exact_retained_nodes':kept,'exact_dropped_nodes':dropped,'exact_pair_checks':checks,'cut_target_indices':cut,'numeric_score':list(score)})
 out={'schema':'quartic-anchor-four-free-distances-H9-v1','target_source':'data/enumerations/fresh_h9_joint_completion_20261013_targets.json','selected_configuration':chosen,'quartic_anchor':quartic,'quartic_anchor_pairing_magnitudes':vec(mags),'rational_anchors':[]if quartic_only else anchors,'minimum_retained':min_retained,'proposal_counts':stats,'caps':caps,'quartic_only':quartic_only,'scope':'The quartic anchor equation is split into its two Q(phi) coefficient equations. Three free distances, or four with a third rational-anchor equation, leave a line. Only accepted proposals and their exact multi-gap/retention checks are certified; the parameter family is not exhaustively excluded.'};OUTPUT.write_text(json.dumps(out,indent=2)+'\n');print({'accepted_caps':len(caps),'stats':stats},flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--min-retained',type=int,default=18);p.add_argument('--either',action='store_true');p.add_argument('--dotted',action='store_true');p.add_argument('--quartic-only',action='store_true');p.add_argument('--output');a=p.parse_args()
 if a.output:OUTPUT=ROOT/a.output
 run(a.min_retained,not a.either,a.dotted,a.quartic_only)
