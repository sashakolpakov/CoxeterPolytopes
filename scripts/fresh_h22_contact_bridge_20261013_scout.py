#!/usr/bin/env python3
"""Bounded numerical all-right edge/finite-anchor bridges of the pi/6 core."""
import json,sys,math
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import fresh_h22_multiquadratic_20261013 as P
V=P.V
TAG='fresh_h22_contact_bridge_20261013'
OUT=V.P/'results'/f'{TAG}_scout.json'
LABELS=(2,3,4,5,6,8)
def load():
 d=V.read(P.OUT)
 for p,s in d['sources'].items():assert V.sha(p)==s
 M,R,h=V.setup();recipes=V.read(V.REC)['recipes'];qold=np.array([float(V.pair(M,w,w))for w in R]);caps={};qq={}
 for k in d['selected_original_types']:
  rec=recipes[k];u=np.array(list(map(Q,rec['u'])),float);v=np.array(list(map(Q,rec['v'])),float);q=Q(rec['norm']);qu=V.pair(M,list(map(Q,rec['u'])),list(map(Q,rec['u'])));qv=V.pair(M,list(map(Q,rec['v'])),list(map(Q,rec['v'])));rad=q/qv-3*q*q*qu/(4*qv)
  r=-float(q)*math.sqrt(3)/2*u+rec['radical_branch']*math.sqrt(float(rad))*v
  if r[-1]<0:r[-1]*=-1
  caps[k]=r;qq[k]=float(q)
 c=V.read(P.COUPLED);u,v=np.array([[float(Q(a))for a in rr]for rr in c['rational_duals']]);a,b=[np.array([float(Q(x))for x in recipes[17][s]])for s in('u','v')]
 Z=(math.sqrt(31)/3*a-1.5*math.sqrt(3)*b)*d['normalized_frame_orientation']
 caps[23]=(209*Z-241*u)/120;caps[24]=(241*Z-418*v)/math.sqrt(29281);qq[23]=qq[24]=2.
 types=d['selected_original_types']+[23,24];W=list(map(lambda x:np.array(x,float),R));norms=list(qold);meta=[{'parent':i}for i in range(812)]
 for k in types:
  r=caps[k];f=r.copy();f[-1]*=-1;W.extend((r,f));norms.extend((qq[k],qq[k]));meta.extend(({'type':k,'height':1},{'type':k,'height':-1}))
 W=np.array(W);M=np.array(M,float);h=np.array(h,float);n=np.array(norms);G=W@M@W.T;np.fill_diagonal(G,n)
 # Incident Gram entries use the already certified finite labels, avoiding
 # cancellation in the large-coordinate cap recipes.
 perm=d['parent_reflection_permutation']
 for rec in d['cap_records']:
  k=rec['type'];idx=812+2*types.index(k)
  for j0,m in rec['finite_parent_labels'].items():
   j=int(j0);p=-math.sqrt(qq[k]*n[j])*math.cos(math.pi/m)if m!=2 else 0.
   for a,b in((idx,j),(idx+1,perm[j])):G[a,b]=G[b,a]=p
 for row in d['all_signed_cap_pairs']:
  if row['kind']=='right_by_exact_frame_identity':i=812+row['i'];j=812+row['j'];G[i,j]=G[j,i]=0.
 vertices=[A*Z+B*u+C*v for A,B,C in[x['frame_coefficients']for x in d['compact_H2_vertices']]]
 return M,W,h,G,n,meta,d,vertices
def solve(A,B,d):
 # -x^2+y^2=2, A*x+B*y=d. No one-dimensional family completeness is
 # inferred when a numerically degenerate coefficient is skipped.
 if abs(B)<1e-11:
  if abs(A)<1e-11:return []
  x=d/A;return [(x,sg*math.sqrt(2+x*x))for sg in(-1,1)]
 aa=A*A-B*B;bb=-2*d*A;cc=d*d-2*B*B
 if abs(aa)<1e-11:
  if abs(bb)<1e-11:return []
  xx=[-cc/bb]
 else:
  D=bb*bb-4*aa*cc
  if D< -1e-9:return []
  D=max(D,0.);xx=[(-bb+sg*math.sqrt(D))/(2*aa)for sg in(-1,1)]
 return [(x,(d-A*x)/B)for x in xx]
def run():
 M,W,h,G,n,meta,doc,vertices=load();WM=W@M;hc=M@h;target=np.array(list(map(Q,doc['strict_escape']['vector'])),float);tc=M@target;bounds=np.sqrt(2*n);special={i for i,x in enumerate(meta)if x.get('type')in(17,23,24)and x.get('height')==1};rows=[];accepted=[];seen=set();counts=Counter();best=[]
 for vi,source in enumerate(doc['compact_H2_vertices']):
  inc=source['incident'];A=G[np.ix_(inc,inc)];eig=np.linalg.eigvalsh(A);assert eig[0]>1e-9
  AI=np.linalg.inv(A);point=vertices[vi];v=point/math.sqrt(-source['norm']);assert abs(v@M@v+1)<2e-7
  for k,released in enumerate(inc):
   edge=[i for i in inc if i!=released];du=AI[k]@W[inc];e=du/math.sqrt(AI[k,k]);assert abs(e@M@e-1)<2e-6 and abs(e@M@v)<2e-6
   null=v-e;ps=WM@null;others=[i for i in range(860)if i not in edge];maxp=max(ps[others]);row={'vertex':vi,'released':released,'edge':edge,'special_contacts':sorted(special.intersection(edge)),'maximum_nonedge_null_pairing':float(maxp),'null_active_extra':[i for i in others if abs(ps[i])<1e-6],'trials':0,'accepted':0};rows.append(row)
   if not row['special_contacts']:counts['no_special_edge_contact']+=1;continue
   va=WM@v;ea=WM@e
   for anchor in range(210,511):
    if anchor in edge:continue
    for label in LABELS:
     d=-math.sqrt(2*n[anchor])*math.cos(math.pi/label)if label!=2 else 0.
     for x,y in solve(va[anchor],ea[anchor],d):
      row['trials']+=1;counts['roots']+=1
      if x<=1e-9 or y>=-1e-9:counts['not_outgoing']+=1;continue
      r=x*v+y*e
      if r[-1]<=1e-8 or r@hc>=-1e-8:counts['height_or_center']+=1;continue
      pair=x*va+y*ea;pair[edge]=0.;pair[anchor]=d
      if max(pair)>2e-6:counts['positive_pair']+=1;continue
      cos=-pair/bounds;intersect=np.where(cos<1-2e-6)[0];mm=np.rint(np.pi/np.arccos(np.clip(cos[intersect],-1,1))).astype(int);err=abs(cos[intersect]-np.cos(np.pi/np.maximum(mm,2)))
      if len(err)and max(err)>2e-6:counts['forbidden_finite_angle']+=1;continue
      if np.any(abs(cos-1)<2e-6):counts['parallel_or_near_parallel']+=1
      key=tuple(np.round(r,8))
      if key in seen:counts['duplicate']+=1;continue
      seen.add(key);row['accepted']+=1;counts['accepted']+=1
      contact={str(i):int(m)for i,m in zip(intersect,mm,strict=True)};rec={'index':len(accepted),'vertex':vi,'released':released,'edge':edge,'anchor':anchor,'label':label,'x':float(x),'y':float(y),'normal':r.tolist(),'computed_norm':float(r@M@r),'center_pair':float(r@hc),'height':float(r[-1]),'escape_pair':float(r@tc),'finite_contacts':contact,'minimum_dotted_cosine':float(min(cos[cos>=1]))if np.any(cos>=1)else None};accepted.append(rec)
   print('edge',vi,released,'null',maxp,'accepted',row['accepted'],flush=True)
 sources={p:V.sha(p)for p in(P.OUT,V.REC,P.COUPLED,'scripts/'+Path(__file__).name)}
 out={'status':'BOUNDED_NUMERICAL_ALL_RIGHT_EDGE_ANCHOR_SCOUT','sources':sources,'facets':860,'vertices':4,'edge_rows':rows,'counts':dict(counts),'anchors':list(range(210,511)),'prescribed_anchor_labels':list(LABELS),'accepted':accepted,'scope':'Only all-right profiles on the21-wall edges of four saved vertices, with one of301 positive-height rational parent anchors and six selected finite labels. Numerical proposals require exact certification. Degenerate equations and other labels/profiles are not excluded.'}
 OUT.write_text(json.dumps(out,indent=2)+'\n');print(dict(counts),'cutters',sum(r['escape_pair']>0 for r in accepted),OUT,flush=True)
if __name__=='__main__':run()
