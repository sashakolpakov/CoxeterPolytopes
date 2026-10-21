#!/usr/bin/env python3
"""Extract literal angle recipes and inspect both norm-quadratic branches."""
import json,sys
import numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import fresh_h9_global_h5_20261013_scout as S
P=S.P;TAG=S.TAG;FIXED=['seed2','compact1','compact2','compact3'];COMMON=['seed1']+FIXED
T=['seed6','crosspair0','nextH4','compact0'];BASES=['wall17','facePrism']
def edge(ci):return[n for i,n in enumerate(T)if i!=ci//2]+[BASES[ci%2]]+FIXED
def run():
 R,G,c,names,doc,p,C,NR,refs,sources,counts=S.load();lookup=dict(zip(names,R));extra=[]
 oldbank=json.loads((P/'data/enumerations/fresh_h9_h5_closure_next_20261013_scout.json').read_text())['candidates']
 for j,i in enumerate((7315,4818,9678)):
  row=oldbank[i];name='alternate_endpoint'+str(j);lookup[name]=np.array(row['normal']);extra.append({'name':name,'incident':FIXED+row['edge']+[row['anchor']],'profile':[2]*4+row['profile']+[row['label']],'branch':row['branch'],'source_row':i})
 def extract(tag,j,row):
  if tag=='fresh_h9_higher_face_20261013_scout':
   return ['crosspair0','nextH4','compact0','facePrism']+[n for n in COMMON if n!=row['free']]+['wall17'],[2]*8+[row['label']]
  if tag=='fresh_h9_higher_face_20261013_branch_scout':
   if j==0:return edge(1)+['wall17'],[2]*8+[16]
   if j<6:
    anchor='roof1'if j==4 else'roof0';return edge(row['corner'])+[anchor],row.get('profile',[2]*8)+[row['anchor_label']]
   return row['edge']+['roof'+str(row['anchor_existing'])],row['profile']+[row['anchor_label']]
  if tag=='fresh_h9_higher_face_20261013_profile_scout':return row['edge']+['roof'+str(row['anchor_existing'])],row['profile']+[row['anchor_label']]
  if tag=='fresh_h9_higher_face_20261013_finite_truncation_scout':return row['right_fixed']+row['incident'],[2]*4+row['profile']
  if tag=='fresh_h9_higher_face_20261013_truncation_scout':return row['right']+[row['anchor']],[2]*8+[row['label']]
  ss=FIXED+row.get('four',row.get('edge'))+[row['anchor']];ll=[2]*4+row['profile']+[row['label']]
  if tag.endswith('_stage1_scout')and'_closure_next_'in tag:ss=[n.replace('endpoint','alternate_endpoint')if n.startswith('endpoint')else n for n in ss]
  return ss,ll
 cache={};recipes=[];bad=[];floating_mismatch=[];escape=json.loads((P/'data/enumerations'/f'{TAG}_scout.json').read_text())['best']['strict_point'];x=np.array(escape)
 for tag,key in S.BANKS:
  rows=json.loads((P/'data/enumerations'/f'{tag}.json').read_text())[key]
  for j,row in enumerate(rows):
   ss,ll=extract(tag,j,row);E=np.array([lookup[n]for n in ss]);b=-2*np.cos(np.pi/np.array(ll));r=np.array(row['normal']);res=max(abs(E@G@r-b));assert res<2e-4,(tag,j,res)
   kk=tuple(ss)
   if kk not in cache:
    M=E@G@E.T;MI=np.linalg.inv(M);cp=E@G@c;cx=c@G@x;ex=E@G@x;qz=c@G@c-cp@MI@cp;zx=cx-cp@MI@ex;cache[kk]=(MI,cp,ex,qz,zx)
   MI,cp,ex,qz,zx=cache[kk];qp=b@MI@b;rad=(2-qp)/qz;assert rad>0,(tag,j,rad);rt=np.sqrt(rad);branches=[]
   part=b@MI@E;z=c-cp@MI@E;distance={sg:float(np.max(np.abs(part+sg*rt*z-r)))for sg in(-1,1)};chosen=min(distance,key=distance.get);assert distance[chosen]<2e-4,(tag,j,distance)
   for sg in(-1,1):
    center=b@MI@cp+sg*rt*qz;target=b@MI@ex+sg*rt*zx
    if center<-1e-7:
     branches.append({'branch':sg,'center':float(center),'target':float(target)})
     if target>=-1e-7:bad.append((tag,j,sg,float(target)))
   assert branches
   recipes.append({'source':tag,'row':j,'incident':ss,'profile':ll,'branch':chosen,'numerical_formula_error':distance[chosen],'center_oriented_branches':branches})
  print(tag,len(rows),'bad cumulative',len(bad),flush=True)
 cc=c/np.sqrt(-c@G@c);xx=(1-1e-6)*x+1e-6*cc;y=-G@xx;target=[{'a':str(int(round(v*10**12)))+'/1000000000000','b':'0'}for v in y]
 out={'status':'NUMERICAL_BOTH_BRANCH_RECIPE_SCREEN','sources':sources,'extra_anchors':extra,'recipes':recipes,'groups':len(cache),'bad_branches':bad,'target_point':escape,'exact_escape_proposal':target,'scope':'Recipe extraction and numerical diagnostic only. Exact verification must reconstruct the selected algebraic branches and all interval signs.'}
 (P/'data/enumerations'/f'{TAG}_recipes.json').write_text(json.dumps(out,indent=2)+'\n');print('recipes',len(recipes),'groups',len(cache),'bad',len(bad),bad[:10],flush=True)
if __name__=='__main__':run()
