#!/usr/bin/env python3
"""Bounded spherical endpoint profiles on the current46 open H5 edges."""
import json,itertools,math
from pathlib import Path
from collections import Counter
import numpy as np
import fresh_h9_h5_closure_20261013_scout as S
P=S.P;TAG='fresh_h9_h5_closure_next_20261013';OLD='data/enumerations/fresh_h9_h5_closure_20261013_verify.json';LABELS=(2,3,4,5,6,8,10,12,16)
def load():
 R,G,c,names,previous,F,p,C,_=S.load();doc=json.loads((P/OLD).read_text())
 for src,sha in doc['sources'].items():assert S.B.V.sha(src)==sha
 for row in doc['new_caps']:
  names.append(row['specification']['name']);R=np.vstack([R,np.array([[(int(a)+int(b))/(2*10**80)for a,b in row['normal']]])])
 assert len(R)==123;face=doc['H5_face']['facets'];J=np.diag([-1.,1,1,1,1,1]);RF=np.c_[-R@G@p,R@G@C];phys=np.column_stack([p,C]);y=np.array([float(S.B.K.parse(x)[0])+float(S.B.K.parse(x)[1])*(1+math.sqrt(5))/2 for x in doc['strict_escape']['point']]);target=-np.linalg.solve(G,y)
 return R,G,c,names,doc,p,C,J,RF,phys,target
def affine_free(five,b,gram):
 A=np.block([[gram,b[:,None]],[b[None,:],np.array([[2.]])]])
 for mask in range(1,32):
  ss=[5]+[i for i in range(5)if mask>>i&1];B=A[np.ix_(ss,ss)]
  if np.any(B< -2-1e-7):continue
  ev=np.linalg.eigvalsh(B)
  if abs(ev[0])<2e-7:return False
 return True
def run():
 R,G,c,names,doc,p,C,J,RF,phys,target=load();face=doc['H5_face']['facets'];fids=[names.index(n)for n in face];fixed=doc['H5_face']['fixed_walls'];fixedids=[names.index(n)for n in fixed];edges=[e['walls']for e in doc['H5_face']['edges']if len(e['endpoints'])==1];edgekeys={tuple(sorted(e)):i for i,e in enumerate(edges)};records=[];seen=set();counts=Counter();stats=[];allb=-2*np.cos(np.pi*np.reciprocal(np.array(LABELS,float)));allb[0]=0.;center=np.array([1.,0,0,0,0,0])
 for ei,edge in enumerate(edges):
  ids=[names.index(n)for n in edge];E=RF[ids];M=E@J@E.T;MI=np.linalg.inv(M);assert np.linalg.eigvalsh(M)[0]>0
  options=[]
  for i in range(4):options.append([k for k,b in enumerate(allb)if b*b*MI[i,i]<2-1e-7])
  good=[]
  for profile in itertools.product(*options):
   b4=allb[list(profile)];q=2-b4@MI@b4
   if q>1e-7:good.append((tuple(LABELS[k]for k in profile),b4))
  count=0;counts['spherical_profiles']+=len(good)
  for anchor in face:
   if anchor in edge:continue
   ai=names.index(anchor);five=edge+[anchor];jj=ids+[ai];W=RF[jj];H=W@J@W.T
   if abs(np.linalg.det(H))<1e-9:counts['dependent_anchor_spans']+=1;continue
   HI=np.linalg.inv(H);z=center-(center@J@W.T)@HI@W;qz=z@J@z
   if abs(qz)<1e-9:counts['null_complement']+=1;continue
   for profile,b4 in good:
    for m,banchor in zip(LABELS,allb):
     counts['metric_profiles']+=1;b=np.r_[b4,banchor];part=b@HI@W;rad=(2-part@J@part)/qz
     if rad<=1e-10:continue
     for branch in(-1,1):
      root=part+branch*math.sqrt(rad)*z
      if root[0]<=1e-7:continue
      ps=RF@J@root;non=np.ones(len(R),bool);non[jj+fixedids]=False
      if max(ps[non])>=-2-2e-6:continue
      if not affine_free(five,b,H):counts['affine_rejected']+=1;continue
      r=phys@root;key=tuple(np.round(r,8))
      if key in seen:counts['duplicate']+=1;continue
      seen.add(key);covers=[];newvertices=[]
      A=np.block([[H,b[:,None]],[b[None,:],np.array([[2.]])]])
      for subset in itertools.combinations(range(5),4):
       if np.linalg.eigvalsh(A[np.ix_(list(subset)+[5],list(subset)+[5])])[0]>1e-7:
        ss=tuple(sorted(five[i]for i in subset));newvertices.append(list(ss))
        if ss in edgekeys:covers.append(edgekeys[ss])
      assert ei in covers;count+=1;cut=float(r@G@target);records.append({'index':len(records),'edge_index':ei,'edge':edge,'profile':list(profile),'anchor':anchor,'label':m,'branch':branch,'normal':r.tolist(),'target_cut':cut,'closed_old_edges':sorted(covers),'new_vertex_fourwalls':newvertices,'nonincident_dotted_margin':float(-2-max(ps[non])),'complement_norm':float(qz)})
  stats.append({'edge_index':ei,'spherical_profiles':len(good),'new_candidates':count});print('edge',ei,'profiles',len(good),'accepted',count,'total',len(records),flush=True)
 counts['accepted']=len(records);counts['target_cutters']=sum(r['target_cut']>1e-8 for r in records);counts['multiple_old_edges']=sum(len(r['closed_old_edges'])>1 for r in records)
 out={'status':'BOUNDED_NUMERICAL_NONRIGHT_ENDPOINT_PROFILES_ON46_OPEN_H5_EDGES','sources':{p:S.B.V.sha(p)for p in(OLD,'scripts/'+TAG+'_scout.py')},'retained_facets':123,'labels':list(LABELS),'edges':edges,'statistics':stats,'counts':dict(counts),'candidates':records,'scope':'Bounded finite labels on spherical fourwall edge profiles, one additional retained face anchor, strict dotted other walls and numerical affine screen. No all-label completeness or exact candidate assertion.'};(P/'data/enumerations'/f'{TAG}_scout.json').write_text(json.dumps(out,indent=2)+'\n');print(dict(counts),flush=True)
if __name__=='__main__':run()
