#!/usr/bin/env python3
"""Three dotted seed pairings, two anchor-angle relations, exact quadratic norm.

Numerical bounded proposals only. All accepted caps have exact norm and angle
checks, contain the strict center, and cut one of the designated witnesses.
"""
import argparse,heapq,json
from itertools import combinations,product
import numpy as np
from scripts.fresh_h8_alternate_selection_20261013 import load,test
from scripts.fresh_h8_alternate_caps_20261013 import *
from scripts.fresh_h8_release_h9_free_caps_20261013 import eadd,esub,emul,esign,eencode,bsign
OUTPUT=ROOT/'data/enumerations/fresh_h8_threefree_caps_20261013.json'

def run(max_support=3,per_anchor=3,all_targets=False,target_file=None,output=None):
 base,pol,raw,G,I,center,nodes=load();graph=json.loads((ROOT/'data/enumerations/fresh_h8_alternate_selection_20261013.json').read_text());proof=json.loads((ROOT/'results/fresh_h8_alternate_core_cover_independent_20261013.json').read_text());witnesses=[tuple(map(parse,r['escape_y']))for r in proof['five_complete_escape_clauses']]
 if target_file:
  proof=json.loads((ROOT/target_file).read_text());witnesses=[tuple(map(parse,r['certificate']['y']))for r in proof['clauses']]
 polar=348;anchors=[i for i in range(348)if f'{i},{polar}'in graph['edges']];assert len(anchors)==60
 In=np.array([[real(x)for x in row]for row in I]);Cn=np.array([real(x)for x in center]);Y=np.array([[real(x)for x in y]for y in witnesses]);IY=In@Y.T;IC=In@Cn
 targets_to_keep=list(range(len(witnesses)))if all_targets else[0]
 cases=[(a,l,m)for a in anchors for l,m in product(range(3),repeat=2)];aa=np.array([[real(x)for x in nodes[polar]['a']]for _ in cases]);bb=np.array([[real(x)for x in nodes[a]['a']]for a,l,m in cases]);mag=np.array([0.,1.,real(P)]);mm=np.array([[mag[l],mag[m]]for a,l,m in cases]);heaps={(a,t):[]for a in anchors for t in targets_to_keep};serial=0;stats={'quadratic_positive_center_cutters':0,'multi_core_cutters':0,'free_triples':0,'cutters_per_core_witness':[0]*len(witnesses)};patterns=[(supp,labs)for k in range(max_support+1)for supp in combinations(range(7),k)for labs in product((0,1),repeat=k)]
 for free in combinations(range(10),3):
  stats['free_triples']+=1;remaining=[i for i in range(10)if i not in free];assigned=np.zeros(len(cases),dtype=bool)
  for u,v,w in (free,(free[1],free[0],free[2]),(free[2],free[0],free[1])):
   det=aa[:,v]*bb[:,w]-aa[:,w]*bb[:,v];active=(abs(det)>1e-9)&~assigned;assigned|=active
   if not active.any():continue
   ids=np.flatnonzero(active);a,b=aa[ids],bb[ids];d=det[ids];targets=mm[ids]
   tv=(-a[:,u]*b[:,w]+a[:,w]*b[:,u])/d;tw=(-a[:,v]*b[:,u]+a[:,u]*b[:,v])/d
   w1=np.zeros((len(ids),10));w1[:,u]=1;w1[:,v]=tv;w1[:,w]=tw;A=np.einsum('bi,ij,bj->b',w1,In,w1)
   for positions,labels in patterns:
    f=np.zeros(10)
    for j,l in zip(positions,labels):f[remaining[j]]=(1.,real(P))[l]
    rhs0=targets[:,0]-a@f;rhs1=targets[:,1]-b@f;Lv=(rhs0*b[:,w]-a[:,w]*rhs1)/d;Lw=(a[:,v]*rhs1-rhs0*b[:,v])/d
    w0=np.broadcast_to(f,(len(ids),10)).copy();w0[:,v]=Lv;w0[:,w]=Lw
    B=2*np.einsum('bi,ij,bj->b',w1,In,w0);C=np.einsum('bi,ij,bj->b',w0,In,w0)-2;D=B*B-4*A*C;good=(D>1e-9)&(abs(A)>1e-9)
    if not good.any():continue
    for branch in (-1,1):
     x=np.full(len(ids),np.nan);x[good]=(-B[good]+branch*np.sqrt(D[good]))/(2*A[good]);ws=w0+x[:,None]*w1
     cuts=ws@IY;orient=ws@IC;ok=good&np.all(ws[:,list(free)]>2+1e-7,axis=1)&np.any(cuts[:,targets_to_keep]>1e-7,axis=1)&(orient< -1e-7)
     stats['quadratic_positive_center_cutters']+=int(ok.sum());stats['multi_core_cutters']+=int((ok&(np.sum(cuts>1e-7,axis=1)>1)).sum())
     for target in range(len(witnesses)):stats['cutters_per_core_witness'][target]+=int((ok&(cuts[:,target]>1e-7)).sum())
     for j in np.flatnonzero(ok):
      anchor,ml,mm2=cases[int(ids[j])];ncut=int(np.sum(cuts[j]>1e-7));length=np.linalg.norm(ws[j]@In)
      for target in targets_to_keep:
       if cuts[j,target]<=1e-7:continue
       score=(ncut,float(cuts[j,target]/length));serial+=1
       item=(score,serial,(u,v,w,tuple(remaining[z]for z in positions),labels,anchor,ml,mm2,branch));h=heaps[anchor,target]
       if len(h)<per_anchor:heapq.heappush(h,item)
       elif item[0]>h[0][0]:heapq.heapreplace(h,item)
  if stats['free_triples']%20==0:print({'free_triples':stats['free_triples'],'cutters':stats['quadratic_positive_center_cutters']},flush=True)
 old_dotted=[]
 for node in nodes:
  ids=[]
  for j in range(10):
   sq=emul((node['pa'][j],node['pb'][j]),(node['pa'][j],node['pb'][j]),node['d'])
   if esign(esub(emul(sq,(div(scale(O,2),node['q']),Z),node['d']),(scale(O,4),Z)),node['d'])>0:ids.append(j)
  old_dotted.append(tuple(ids))
 caps=[];seen=set()
 for score,_,spec in sorted([r for h in heaps.values()for r in h],reverse=True):
  u,v,w,support,labels,anchor,mi,mj,branch=spec;f=[Z]*10
  for j,l in zip(support,labels):f[j]=(O,P)[l]
  f=tuple(f);ar=nodes[polar]['a'];br=nodes[anchor]['a'];M,N=(Z,O,P)[mi],(Z,O,P)[mj];det=sub(mul(ar[v],br[w]),mul(ar[w],br[v]));assert det!=Z
  rhs0=sub(M,dot(ar,f));rhs1=sub(N,dot(br,f));Lv=div(sub(mul(rhs0,br[w]),mul(ar[w],rhs1)),det);Lw=div(sub(mul(ar[v],rhs1),mul(rhs0,br[v])),det);Tv=div(sub(mul(ar[w],br[u]),mul(ar[u],br[w])),det);Tw=div(sub(mul(ar[u],br[v]),mul(ar[v],br[u])),det)
  w0=list(f);w0[v]=Lv;w0[w]=Lw;w0=tuple(w0);w1=tuple(O if j==u else Tv if j==v else Tw if j==w else Z for j in range(10))
  A=pair(w1,I,w1);B=scale(pair(w0,I,w1),2);C=sub(pair(w0,I,w0),scale(O,2));D=sub(mul(B,B),scale(mul(A,C),4));assert sign(D)>0 and A!=Z
  x=(div(neg(B),scale(A,2)),div(scale(O,branch),scale(A,2)));pa=tuple(neg(add(a,mul(b,x[0])))for a,b in zip(w0,w1));pb=tuple(neg(mul(b,x[1]))for b in w1);ca,cb=mv(I,pa),mv(I,pb);key=(D,ca,cb)
  if key in seen:continue
  seen.add(key);norm=(add(dot(pa,ca),mul(D,dot(pb,cb))),scale(dot(pa,cb),2));assert esign(esub(norm,(scale(O,2),Z)),D)==0
  dotted=[]
  for j,(a,b)in enumerate(zip(pa,pb)):
   assert esign((a,b),D)<=0;squared=emul((a,b),(a,b),D)
   if esign(esub(squared,(scale(O,4),Z)),D)>0:dotted.append(j)
   else:assert any(esign(esub(squared,(q,Z)),D)==0 for q in (Z,O,scale(O,2),scale(O,3),add(O,P),add(scale(O,2),P)))
  assert sorted(dotted)==sorted((u,v,w))
  node={'d':D,'q':scale(O,2),'a':ca,'b':cb,'pa':pa,'pb':pb}
  for index,target in ((polar,M),(anchor,N)):
   assert esign(eadd((dot(pa,nodes[index]['a']),dot(pb,nodes[index]['a'])),(target,Z)),D)==0 and test(node,nodes[index])[0]
  assert esign((dot(ca,center),dot(cb,center)),D)>0
  vals=[(dot(ca,y),dot(cb,y))for y in witnesses];cut=[j for j,p in enumerate(vals)if esign(p,D)<0];assert set(cut)&set(targets_to_keep)
  same_old=[]
  for oi,nodeold in enumerate(nodes):
   if old_dotted[oi]!=tuple(sorted(dotted)):continue
   factor=div(scale(O,2),nodeold['q']);equal=True
   for j in range(10):
    newq=emul((pa[j],pb[j]),(pa[j],pb[j]),D);oldq=emul((nodeold['pa'][j],nodeold['pb'][j]),(nodeold['pa'][j],nodeold['pb'][j]),nodeold['d'])
    diff=(sub(newq[0],mul(factor,oldq[0])),newq[1],neg(mul(factor,oldq[1])),Z)
    if bsign(diff,D,nodeold['d'])!=0:equal=False;break
   if equal:same_old.append(oi)
  if same_old:continue
  caps.append({'free_seed_indices':[u,v,w],'free_original_walls':[base['retained_original_walls'][i]for i in(u,v,w)],'fixed_pairing_magnitudes':vec(f),'anchor_node_indices':[polar,anchor],'anchor_pairing_magnitudes':vec([M,N]),'line_constant':vec(w0),'line_direction':vec(w1),'norm_polynomial_descending':vec([A,B,C]),'discriminant':enc(D),'parameter':eencode(x),'pairings_rational':vec(pa),'pairings_radical':vec(pb),'coefficients_rational':vec(ca),'coefficients_radical':vec(cb),'cut_core_witness_indices':cut,'core_constraint_values':[eencode(p)for p in vals],'numeric_score':[score[0],score[1]]})
 result={'schema':'H8-three-free-two-anchor-H9-caps-v1','base_source':'data/enumerations/fresh_h8_alternate_3_6_20261013.json','old_dictionary_source':'data/enumerations/fresh_h8_alternate_selection_20261013.json','witness_source':'results/fresh_h8_alternate_core_cover_independent_20261013.json','fixed_pairing_maximum_support':max_support,'kept_proposals_per_anchor_and_target':per_anchor,'target_core_witnesses':targets_to_keep,'second_anchor_node_indices':anchors,'proposal_counts':stats,'caps':caps,'outside_old_dictionary_reason':'Exact normalized seed-pairing-square comparison, after dotted-support filtering, proves every accepted cap is a different mirror from all493 old entries.','scope':'Bounded numerical proposals, with exact accepted root/seed-angle/two-anchor/core-cut checks. No volume conclusion and no exhaustive family claim.'}
 if target_file:result['witness_source']=target_file
 (ROOT/output if output else OUTPUT).write_text(json.dumps(result,indent=2)+'\n');print({'exact_caps':len(caps),'stats':stats},flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--support',type=int,default=3);p.add_argument('--per-anchor',type=int,default=3);p.add_argument('--all-targets',action='store_true');p.add_argument('--target-file');p.add_argument('--output');a=p.parse_args();run(a.support,a.per_anchor,a.all_targets,a.target_file,a.output)
