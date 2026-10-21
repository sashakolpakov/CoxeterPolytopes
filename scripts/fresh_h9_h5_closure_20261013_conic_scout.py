#!/usr/bin/env python3
"""Four-contact continuous circle screen, then a fifth finite anchor."""
import argparse,json,math
from itertools import product
import numpy as np
from scipy.linalg import null_space
import fresh_h9_h5_closure_20261013_scout as S
P=S.P;PRE=S.PRE;TAU=2*np.pi
LABELS=(2,3,4,5,6,8);ANCHORS=(2,3,4,5,6,8,10,12,16)
def intersect(a,b):return[(max(x,u),min(y,v))for x,y in a for u,v in b if max(x,u)<=min(y,v)+1e-12]
def circle_constraint(a,b,c):
 r=np.hypot(a,b)
 if r<1e-10:return[(0.,TAU)]if c>=-1e-9 else[]
 if c>=r-1e-12:return[(0.,TAU)]
 if c<=-r-1e-12:return[]
 t=math.acos(np.clip(c/r,-1,1));phase=math.atan2(b,a);l=(phase+t)%TAU;u=(phase+TAU-t)%TAU
 return[(l,u)]if l<=u else[(0.,u),(l,TAU)]
def run(stage=0):
 R,G,c,names,final,F,p,C,target=S.load();J=np.diag([-1.,1,1,1,1,1]);RF=np.c_[-R@G@p,R@G@C];tf=np.r_[-target@G@p,target@G@C];phys=np.column_stack([p,C]);fixed=final['H5_face']['fixed_walls'];face=final['H5_face']['facets'];support=['seed6','nextH4','facePrism','roof3','adaptive0'];records=[];interval_count=0;tested=0;seen=set();stats=[]
 if stage:
  state=json.loads((P/'data/enumerations'/f'{PRE}_state{stage}.json').read_text());R=np.vstack([R,np.array([r['normal']for r in state['selected']])]);names+= [r['name']for r in state['selected']];face=face+[r['name']for r in state['selected']];target=np.array(state['target_point']);support=state['support'];RF=np.c_[-R@G@p,R@G@C];tf=np.r_[-target@G@p,target@G@C]
 for omitted in support:
  four=[n for n in support if n!=omitted];ids=[names.index(n)for n in four];E=RF[ids];M=E@J@E.T;MI=np.linalg.inv(M);U=null_space(E@J);UU=U.T@J@U;ev,vec=np.linalg.eigh(UU)
  if min(ev)<=0:stats.append({'omitted':omitted,'skipped':'nonpositive residual plane'});continue
  U=U@vec@np.diag(1/np.sqrt(ev));proj=MI@E;ab=RF@J@U;abt=tf@J@U;count=0;localint=0
  for labels in product(LABELS,repeat=4):
   tested+=1;bs=-2*np.cos(np.pi/np.array(labels));part=bs@proj;q=2-part@J@part
   if q<=1e-10:continue
   amp=np.sqrt(q);pc=RF@J@part;pt=tf@J@part;AB=amp*ab;BT=amp*abt
   feasible=circle_constraint(-BT[0],-BT[1],pt-1e-8)
   if not feasible:continue
   for i in range(len(R)):
    if names[i]in fixed or i in ids:continue
    feasible=intersect(feasible,circle_constraint(AB[i,0],AB[i,1],-pc[i]))
    if not feasible:break
   if not feasible:continue
   interval_count+=1;localint+=1
   for anchor in face:
    if anchor in four:continue
    ai=names.index(anchor);a,b=AB[ai];rad=np.hypot(a,b)
    if rad<1e-9:continue
    for m in ANCHORS:
     want=-2*np.cos(np.pi/m);z=(want-pc[ai])/rad
     if abs(z)>1+1e-10:continue
     alpha=math.acos(np.clip(z,-1,1));phase=math.atan2(b,a)
     for th in((phase+alpha)%TAU,(phase-alpha)%TAU):
      if not any(l-1e-9<=th<=u+1e-9 for l,u in feasible):continue
      vv=np.array([np.cos(th),np.sin(th)]);ps=pc+AB@vv;non=np.ones(len(R),bool);non[ids+[ai]+[names.index(n)for n in fixed]]=False
      if np.max(ps[non])>=-2-1e-7:continue
      r=phys@(part+amp*U@vv)
      if r@G@c>=-1e-8 or r@G@target<=1e-8:continue
      key=tuple(np.round(r,8))
      if key in seen:continue
      seen.add(key);records.append({'four':four,'profile':list(labels),'anchor':anchor,'label':m,'angle_parameter':float(th),'normal':r.tolist(),'target_cut':float(r@G@target),'nonincident_dotted_margin':float(np.max(ps[non])+2),'continuous_acute_intervals':[list(x)for x in feasible]});count+=1
  stats.append({'omitted':omitted,'continuous_profiles':localint,'accepted':count});print(omitted,'continuous',localint,'accepted',count,'total',len(records),flush=True)
 out={'status':'BOUNDED_NUMERICAL_FOUR_CONTACT_CONIC_SCOUT','profile_labels':list(LABELS),'anchor_labels':list(ANCHORS),'profiles_tested':tested,'continuous_profiles':interval_count,'statistics':stats,'candidates':records,'source_sha256':S.B.V.sha(S.OLD+'_final_verify.json'),'scope':'Continuous acuteness/target-cut screen on five chosen four-contact circles, followed by a bounded finite-label anchor set. No exhaustive cap or finite-volume claim.'}
 out['stage']=stage;out['retained_wall_count']=len(R);out['target_support']=support
 suffix=''if stage==0 else f'_stage{stage}'
 (P/'data/enumerations'/f'{PRE}_conic_scout{suffix}.json').write_text(json.dumps(out,indent=2)+'\n');print('DONE',len(records),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--stage',type=int,default=0);a=ap.parse_args();run(a.stage)
