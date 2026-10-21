#!/usr/bin/env python3
"""Two vertex-free bridge H2 faces: numerical residual and bounded anchor scout."""
import json,math,itertools
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from collections import Counter
import numpy as np
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e7_bridge_face_20261014_scout';BASE='results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json'
LABELS=[2,3,4,5,6,8,10,12,16,24,30,60,'parallel']
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 x=json.loads((ROOT/BASE).read_text())
 for p,h in x['sources'].items():assert sha(p)==h,p
 rad=np.sqrt(x['field']['base_square_root_basis']);dd=sum(float(F(a))*b for a,b in zip(x['field']['discriminant'],rad))**.5
 def val(a):return sum(float(F(q))*r for q,r in zip(a[0],rad))+dd*sum(float(F(q))*r for q,r in zip(a[1],rad))
 def vec(v):return np.array([val(a)for a in v])
 N=np.array([vec(v)for v in x['normals']]);E=vec(x['strict_escape']['point']);H=vec(x['strict_center']);Q=np.array([[float(F(a))for a in r]for r in x['form']]);norms=np.array(list(map(float,x['norms'])));G=N@Q@N.T;NQ=N@Q
 finite=[i for i in range(95)if i!=87 and -G[i,87]<math.sqrt(norms[i]*norms[87])-1e-7]
 pd18=[];pd19=[]
 for rank,out in[(19,pd19),(18,pd18)]:
  for ss in itertools.combinations(finite,rank-1):
   ids=sorted(list(ss)+[87]);M=G[np.ix_(ids,ids)]
   if np.min(np.linalg.eigvalsh(M))>1e-7:out.append(ids)
 assert len(pd19)==0 and len(pd18)==2
 reports=[];hits=[];stats=Counter()
 for S in pd18:
  U=np.linalg.svd(N[S]@Q)[2][-3:].T;M=U.T@Q@U;A=NQ@U;e=E@Q@U;gi=np.linalg.inv(G[np.ix_(S,S)])
  V=E-(E@Q@N[S].T)@gi@N[S];P=H-(H@Q@N[S].T)@gi@N[S]
  projection_norms=norms-np.einsum('ij,jk,ik->i',G[:,S],gi,G[:,S]);assert max(projection_norms)<1e-7
  assert P@Q@P<0 and max((NQ@P)[[i for i in range(95)if i not in S]])< -1e-7
  lp=linprog(np.ones(95),A_eq=A.T,b_eq=e,bounds=(0,None),method='highs')
  other=[i for i in range(95)if i not in S];rays=[];keys=set();anchors=set()
  def label_gate(pair):
   labels=[];bad=[]
   for i,c in enumerate(-pair/np.sqrt(2*norms)):
    if abs(c)<2e-7:labels.append(2)
    elif c>=1-2e-7:labels.append('parallel_or_dotted')
    elif c<.5-2e-7:bad.append([i,float(c)]);labels.append(None)
    else:
     m=math.pi/math.acos(min(c,1));mi=round(m)
     if abs(m-mi)>2e-5:bad.append([i,float(c)]);labels.append(None)
     else:labels.append(mi)
   return labels,bad
  for i,j in itertools.combinations(other,2):
   z=np.cross(A[i],A[j]);q=z@M@z
   if q<1e-18:continue
   z*=math.sqrt(2/q)
   for z1 in[z,-z]:
    pair=A@z1
    if max(pair)>1e-6:continue
    key=tuple(np.round(z1,7))
    if key in keys:continue
    keys.add(key);cut=e@z1;labels,bad=label_gate(pair)
    active=[k for k in other if abs(pair[k])<1e-6];anchors.update(active)
    rays.append(dict(right_anchors=active,normal=(U@z1).tolist(),escape_pairing=float(cut),bad_finite_pairs=bad))
  # The finite anchor pool is chosen from this exact mathematical type of
  # polar-cone ray, but its enumeration and all gates here remain numeric.
  anchors=sorted(anchors);localhits=[];hitkeys=set()
  for i,j in itertools.combinations(anchors,2):
   AA=A[[i,j]];cross=np.cross(*AA)
   if np.linalg.norm(cross)<1e-10:continue
   cross/=np.linalg.norm(cross);pinv=AA.T@np.linalg.inv(AA@AA.T);qa=cross@M@cross
   for mi,mj in itertools.product(LABELS,repeat=2):
    rhs=-np.array([math.sqrt(2*norms[k])*(1 if m=='parallel'else math.cos(math.pi/m))for k,m in[(i,mi),(j,mj)]])
    part=pinv@rhs;qb=2*part@M@cross;qc=part@M@part-2;disc=qb*qb-4*qa*qc
    if disc< -1e-8:continue
    roots=[-qc/qb]if abs(qa)<1e-12 and abs(qb)>1e-12 else[]if abs(qa)<1e-12 else[(-qb+sign*math.sqrt(max(0,disc)))/(2*qa)for sign in[-1,1]]
    for t in roots:
     stats['normalized_proposals']+=1;z=part+t*cross
     if e@z>1e-7:stats['E_cutters']+=1
     pair=A@z
     if max(pair)>2e-6:continue
     stats['acute_proposals']+=1
     if e@z>1e-7:stats['acute_E_cutters']+=1
     labels,bad=label_gate(pair)
     if bad:continue
     r=U@z;key=tuple(np.round(r,7))
     if key in hitkeys:continue
     hitkeys.add(key);assert abs(r@Q@r-2)<1e-5
     rec=dict(support=S,anchors=[i,j],anchor_labels=[mi,mj],normal=r.tolist(),escape_pairing=float(r@Q@E),all_labels=labels)
     localhits.append(rec);hits.append(rec);
     if rec['escape_pairing']>0:print('NUMERIC_COXETER_CUTTER',i,j,mi,mj,rec['escape_pairing'],flush=True)
  reports.append(dict(support=S,minimum_Gram_eigenvalue=float(min(np.linalg.eigvalsh(G[np.ix_(S,S)]))),strict_face_point=P.tolist(),strict_face_point_norm=float(P@Q@P),maximum_old_projection_norm=float(max(projection_norms)),escape_projection_norm=float(V@Q@V),positive_dual_LP_success=bool(lp.success),acute_polar_extreme_rays=rays,finite_anchor_pool=anchors,accepted_normal_count=len(localhits)))
 coherent=[]
 for i,j in itertools.combinations(range(len(hits)),2):
  a=np.array(hits[i]['normal']);b=np.array(hits[j]['normal']);g=a@Q@b
  if g>1e-7:continue
  if g< -2+1e-7:kind='parallel_or_dotted'
  elif abs(g)<1e-7:kind='right'
  else:
   c=-g/2
   if c<.5-1e-7:continue
   m=math.pi/math.acos(c)
   if abs(m-round(m))>2e-5:continue
   kind='pi'+str(round(m))
  coherent.append(dict(candidates=[i,j],pairing=float(g),relation=kind))
 paths=[BASE,f'scripts/{TAG}.py']
 return dict(status='NUMERICAL_BRIDGE_FACE_FEASIBILITY_AND_BOUNDED_ANCHORS',sources={p:sha(p)for p in paths},bridge_finite_neighbors=finite,PD19_supports=pd19,PD18_supports=pd18,faces=reports,anchor_labels=LABELS,counts=dict(stats,accepted_Coxeter_normals=len(hits),accepted_escape_cutters=sum(r['escape_pairing']>1e-7 for r in hits),compatible_candidate_pairs=len(coherent)),candidates=hits,compatible_pairs=coherent,scope='All arithmetic/eigenvalue/LP and angle gates here are numerical. No exact construction or all-real/label exclusion. PD18-support completeness is a bounded190/1140subset numerical check among bridge finite neighbors. Two old finite anchors sampled only from the polar-ray active-wall pool. Primary95 and other frozen sources unchanged.')
if __name__=='__main__':
 d=run();(ROOT/f'results/{TAG}.json').write_text(json.dumps(d,indent=2)+'\n');print(d['counts']);print([(r['support'][-2:],r['finite_anchor_pool'])for r in d['faces']])
