#!/usr/bin/env python3
"""Joint escape cutting with two anchors in one coherent H9 configuration.

Numerical proposals enforce multiple gaps and preserve most existing walls.
Accepted roots, retention masks, angles, and cuts are then checked exactly.
"""
import argparse,heapq,json
from itertools import combinations,product
import numpy as np
from scipy.spatial import HalfspaceIntersection
from scripts.fresh_h8_cap_completion_20261013 import load,certify,node,gate,ROOT,read
from scripts.fresh_h8_alternate_caps_20261013 import *
from scripts.fresh_h8_release_h9_free_caps_20261013 import eadd,esub,emul,esign,eencode

OUTPUT=ROOT/'data/enumerations/fresh_h9_joint_completion_20261013.json'
TARGETS=ROOT/'data/enumerations/fresh_h9_joint_completion_20261013_targets.json'

def setup():
 G,I,center,nodes,polys,cn,qn,hn,_=load();old=read('data/enumerations/fresh_h8_cap_completion_20261013_bottleneck_caps.json');nodes+=list(node(c,G)for c in old['caps']);cn=np.array([[real(a)+np.sqrt(real(n['d']))*real(b)for a,b in zip(n['a'],n['b'])]for n in nodes]);configuration=read('data/enumerations/fresh_h8_cap_completion_20261013_extensions.json')['cases'][0]['selected_node_indices']
 return G,I,center,nodes,polys,cn,qn,hn,configuration

def targets():
 G,I,center,nodes,polys,cn,qn,hn,chosen=setup();A=np.vstack((np.eye(10),cn[chosen]));half=np.c_[A[:,9,None]-A[:,:9],-A[:,9]];half/=np.linalg.norm(half[:,:9],axis=1)[:,None];xx=HalfspaceIntersection(half,hn[:9],qhull_options='Qx').intersections;ys=np.c_[xx,1-xx.sum(axis=1)];qs=np.einsum('bi,ij,bj->b',ys,qn,ys);good=np.flatnonzero(qs>1e-7);picked=[];rows=[];distance=np.full(len(ys),np.inf)
 for iteration in range(24):
  j=max(good,key=lambda j:qs[j])if not picked else max(good,key=lambda j:(distance[j],qs[j]))
  if j in picked:break
  picked.append(j);distance=np.minimum(distance,np.linalg.norm(ys-ys[j],axis=1));row=certify(ys[j],chosen,I,center,nodes,hn)
  if row:rows.append({'proposal_vertex_index':int(j),'numeric_vertex_norm':float(qs[j]),'certificate':row})
 out={'source_configuration':'data/enumerations/fresh_h8_cap_completion_20261013_extensions.json','configuration_index':0,'selected_node_indices':chosen,'numeric_vertex_count':len(ys),'numerically_spacelike_vertices':len(good),'targets':rows,'scope':'Farthest-point sampling of numerical vertices proposes geometrically separated gaps. Every saved target is an exact strict future spacelike witness of the full32-facet configuration.'};TARGETS.write_text(json.dumps(out,indent=2)+'\n');print({'targets':len(rows),'vertices':len(ys)},flush=True)

def run(max_support=3,min_retained=18,keep=2):
 if not TARGETS.exists():targets()
 G,I,center,nodes,polys,cn,In,Cn,chosen=setup();data=json.loads(TARGETS.read_text());witnesses=[tuple(map(parse,r['certificate']['y']))for r in data['targets']];anchors=[i for i in chosen if nodes[i]['d']==Z];anchor_pairs=list(combinations(anchors,2));cases=[(a,b,l,m)for a,b in anchor_pairs for l,m in product(range(3),repeat=2)];aa=np.array([[real(x)for x in nodes[a]['a']]for a,b,l,m in cases]);bb=np.array([[real(x)for x in nodes[b]['a']]for a,b,l,m in cases]);mag=np.array([0.,1.,real(P)]);mm=np.array([[mag[l],mag[m]]for a,b,l,m in cases]);Y=np.array([[real(x)for x in y]for y in witnesses]);IY=In@Y.T;IC=In@np.array([real(x)for x in center]);selected_cn=cn[chosen];selected_norms=np.array([real(nodes[i]['q'])for i in chosen]);allowed=np.array([float((lo+hi)/8)for _,_,lo,hi in polys]);heaps={};serial=0;stats={'free_triples':0,'multi_gap_center_candidates':0,'retention_candidates':0,'best_numeric_retention':0,'best_numeric_target_count':0};patterns=[(supp,labs)for k in range(max_support+1)for supp in combinations(range(7),k)for labs in product((0,1),repeat=k)]
 for free in combinations(range(10),3):
  stats['free_triples']+=1;remaining=[i for i in range(10)if i not in free];assigned=np.zeros(len(cases),dtype=bool)
  for u,v,w in (free,(free[1],free[0],free[2]),(free[2],free[0],free[1])):
   det=aa[:,v]*bb[:,w]-aa[:,w]*bb[:,v];active=(abs(det)>1e-9)&~assigned;assigned|=active
   if not active.any():continue
   ids=np.flatnonzero(active);a,b=aa[ids],bb[ids];d=det[ids];m=mm[ids];tv=(-a[:,u]*b[:,w]+a[:,w]*b[:,u])/d;tw=(-a[:,v]*b[:,u]+a[:,u]*b[:,v])/d;w1=np.zeros((len(ids),10));w1[:,u]=1;w1[:,v]=tv;w1[:,w]=tw;A=np.einsum('bi,ij,bj->b',w1,In,w1)
   for positions,labels in patterns:
    f=np.zeros(10)
    for j,l in zip(positions,labels):f[remaining[j]]=(1.,real(P))[l]
    rhs0=m[:,0]-a@f;rhs1=m[:,1]-b@f;Lv=(rhs0*b[:,w]-a[:,w]*rhs1)/d;Lw=(a[:,v]*rhs1-rhs0*b[:,v])/d;w0=np.broadcast_to(f,(len(ids),10)).copy();w0[:,v]=Lv;w0[:,w]=Lw;B=2*np.einsum('bi,ij,bj->b',w1,In,w0);C=np.einsum('bi,ij,bj->b',w0,In,w0)-2;D=B*B-4*A*C;good=(D>1e-9)&(abs(A)>1e-9)
    if not good.any():continue
    for branch in(-1,1):
     x=np.full(len(ids),np.nan);x[good]=(-B[good]+branch*np.sqrt(D[good]))/(2*A[good]);ws=w0+x[:,None]*w1;cuts=ws@IY;nc=(cuts>1e-7).sum(axis=1);ok=good&np.all(ws[:,list(free)]>2+1e-7,axis=1)&(nc>=2)&(ws@IC< -1e-7);stats['multi_gap_center_candidates']+=int(ok.sum())
     idx=np.flatnonzero(ok)
     if not len(idx):continue
     ps=-ws[idx]@selected_cn.T;rat=ps*ps/(2*selected_norms);near=np.min(np.abs(rat[:,:,None]-allowed),axis=2)<1e-6;compat=(ps<=1e-7)&((rat>1+1e-7)|near);nr=compat.sum(axis=1);stats['best_numeric_retention']=max(stats['best_numeric_retention'],int(nr.max()));stats['best_numeric_target_count']=max(stats['best_numeric_target_count'],int(nc[idx].max()))
     for jj in np.flatnonzero(nr>=min_retained):
      j=int(idx[jj]);stats['retention_candidates']+=1;anchor1,anchor2,ml,mm2=cases[int(ids[j])];mask=tuple(chosen[k]for k in range(len(chosen))if compat[jj,k]);cnum=int(nc[j]);length=np.linalg.norm(ws[j]@In);positive=sorted(cuts[j][cuts[j]>1e-7],reverse=True);score=(cnum,int(nr[jj]),float(positive[1]/length));key=(anchor1,anchor2,mask);serial+=1;item=(score,serial,(u,v,w,tuple(remaining[z]for z in positions),labels,anchor1,anchor2,ml,mm2,branch));h=heaps.setdefault(key,[])
      if len(h)<keep:heapq.heappush(h,item)
      elif score>h[0][0]:heapq.heapreplace(h,item)
  if stats['free_triples']%20==0:print(stats,flush=True)
 caps=[];seen=set()
 for score,_,spec in sorted([r for h in heaps.values()for r in h],reverse=True):
  u,v,w,support,labels,anchor1,anchor2,mi,mj,branch=spec;f=[Z]*10
  for j,l in zip(support,labels):f[j]=(O,P)[l]
  f=tuple(f);ar=nodes[anchor1]['a'];br=nodes[anchor2]['a'];M,N=(Z,O,P)[mi],(Z,O,P)[mj];det=sub(mul(ar[v],br[w]),mul(ar[w],br[v]));assert det!=Z;rhs0=sub(M,dot(ar,f));rhs1=sub(N,dot(br,f));Lv=div(sub(mul(rhs0,br[w]),mul(ar[w],rhs1)),det);Lw=div(sub(mul(ar[v],rhs1),mul(rhs0,br[v])),det);Tv=div(sub(mul(ar[w],br[u]),mul(ar[u],br[w])),det);Tw=div(sub(mul(ar[u],br[v]),mul(ar[v],br[u])),det);w0=list(f);w0[v]=Lv;w0[w]=Lw;w0=tuple(w0);w1=tuple(O if j==u else Tv if j==v else Tw if j==w else Z for j in range(10));A=pair(w1,I,w1);B=scale(pair(w0,I,w1),2);C=sub(pair(w0,I,w0),scale(O,2));D=sub(mul(B,B),scale(mul(A,C),4));assert sign(D)>0;x=(div(neg(B),scale(A,2)),div(scale(O,branch),scale(A,2)));pa=tuple(neg(add(a,mul(b,x[0])))for a,b in zip(w0,w1));pb=tuple(neg(mul(b,x[1]))for b in w1);ca,cb=mv(I,pa),mv(I,pb);key=(D,ca,cb)
  if key in seen:continue
  seen.add(key);n={'q':scale(O,2),'d':D,'a':ca,'b':cb,'pa':pa,'pb':pb};norm=(add(dot(pa,ca),mul(D,dot(pb,cb))),scale(dot(pa,cb),2));assert esign(esub(norm,(scale(O,2),Z)),D)==0
  for j in range(10):
   assert esign((pa[j],pb[j]),D)<=0;sq=emul((pa[j],pb[j]),(pa[j],pb[j]),D)
   if j in(u,v,w):assert esign(esub(sq,(scale(O,4),Z)),D)>0
   else:assert any(esign(esub(sq,(t,Z)),D)==0 for t in(Z,O,add(O,P)))
  assert esign((dot(ca,center),dot(cb,center)),D)>0
  for ai,t in((anchor1,M),(anchor2,N)):assert esign(eadd((dot(pa,nodes[ai]['a']),dot(pb,nodes[ai]['a'])),(t,Z)),D)==0
  checks=[{'old_node':i,**gate(n,nodes[i],polys)}for i in chosen];retained=[r['old_node']for r in checks if r['compatible']];dropped=[r['old_node']for r in checks if not r['compatible']];cut=[j for j,y in enumerate(witnesses)if esign((dot(ca,y),dot(cb,y)),D)<0]
  if len(retained)<min_retained or len(cut)<2 or anchor1 not in retained or anchor2 not in retained:continue
  caps.append({'free_seed_indices':[u,v,w],'fixed_pairing_magnitudes':vec(f),'anchor_node_indices':[anchor1,anchor2],'anchor_pairing_magnitudes':vec([M,N]),'line_constant':vec(w0),'line_direction':vec(w1),'norm_polynomial_descending':vec([A,B,C]),'discriminant':enc(D),'parameter':eencode(x),'pairings_rational':vec(pa),'pairings_radical':vec(pb),'coefficients_rational':vec(ca),'coefficients_radical':vec(cb),'exact_retained_nodes':retained,'exact_dropped_nodes':dropped,'exact_pair_checks':checks,'cut_target_indices':cut,'numeric_score':list(score)})
 out={'schema':'joint-anchor-multigap-H9-caps-v1','base_source':'data/enumerations/fresh_h8_alternate_3_6_20261013.json','target_source':str(TARGETS.relative_to(ROOT)),'selected_configuration':chosen,'anchor_node_indices':anchors,'max_fixed_support':max_support,'minimum_retained_caps':min_retained,'proposal_counts':stats,'caps':caps,'scope':'Bounded numerical proposals with two anchors inside one actual configuration; exact accepted roots, seed angles, retention/replacement masks and multi-witness cuts. No exhaustive cap or volume claim.'};OUTPUT.write_text(json.dumps(out,indent=2)+'\n');print({'accepted_caps':len(caps),'stats':stats},flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--targets-only',action='store_true');p.add_argument('--support',type=int,default=3);p.add_argument('--min-retained',type=int,default=18);p.add_argument('--keep',type=int,default=2);a=p.parse_args()
 if a.targets_only:targets()
 else:run(a.support,a.min_retained,a.keep)
