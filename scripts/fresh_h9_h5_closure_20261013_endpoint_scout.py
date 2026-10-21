#!/usr/bin/env python3
"""Couple a spherical H3+A1 edge endpoint to the newly retained closure cap."""
import json,itertools
import numpy as np
from scipy.linalg import null_space
import fresh_h9_h5_closure_20261013_scout as S
P=S.P;PRE=S.PRE;LABELS=(2,3,4,5,6,8,10,12,16)
def run():
 R,G,c,names,final,F,p,C,oldtarget=S.load();state=json.loads((P/'data/enumerations'/f'{PRE}_state3.json').read_text());R=np.vstack([R,[r['normal']for r in state['selected']]]);names+=[r['name']for r in state['selected']];target=np.array(state['target_point']);fixed=final['H5_face']['fixed_walls'];four=['seed6','nextH4','facePrism','roof5'];ids=[names.index(n)for n in four];E=R[ids];M=E@G@E.T;assert np.linalg.eigvalsh(M)[0]>0
 good=[]
 for labs in itertools.product(LABELS,repeat=4):
  b=-2*np.cos(np.pi/np.array(labs));A=np.block([[M,b[:,None]],[b[None,:],np.array([[2.]])]])
  if np.linalg.eigvalsh(A)[0]>1e-7:good.append((labs,b))
 results=[]
 for anchor in ['closure0','closure1','closure2','closure3']:
  base=fixed+four+[anchor];ii=[names.index(n)for n in base];W=R[ii];H=W@G@W.T;HI=np.linalg.inv(H);cc=(c@G@W.T)@HI@W;z=c-cc;qz=z@G@z
  if qz<=0:continue
  for labs,b4 in good:
   for m in LABELS:
    b=np.r_[np.zeros(4),b4,-2*np.cos(np.pi/m)];part=b@HI@W;rad=(2-part@G@part)/qz
    if rad<=0:continue
    for sign in(-1,1):
     r=part+sign*np.sqrt(rad)*z;ps=R@G@r;non=np.ones(len(R),bool);non[ii]=False
     if r@G@c>=-1e-8 or max(ps[non])>=-2-1e-7:continue
     vs=S.geometry(np.vstack([R,r]),G,p,C);results.append({'four':four,'profile':list(labs),'anchor':anchor,'label':m,'normal':r.tolist(),'branch':sign,'old_target_cut':float(r@G@oldtarget),'current_target_cut':float(r@G@target),'maximum_norm':max(v['norm']for v in vs),'finite_vertices':sum(v['norm']<-1e-7 for v in vs),'projective_vertices':vs});print('FOUND',anchor,labs,m,results[-1]['finite_vertices'],results[-1]['maximum_norm'],flush=True)
 out={'status':'BOUNDED_NUMERICAL_SPHERICAL_ENDPOINT_COUPLING','old_spherical_block':four,'finite_profiles':len(good),'profiles':[list(x[0])for x in good],'anchors':['closure0','closure1','closure2','closure3'],'labels':list(LABELS),'accepted':results,'scope':'Spherical profile gate and all122 retained-wall numeric compatibility. No exact or all-label completeness claim.'};(P/'data/enumerations'/f'{PRE}_endpoint_scout.json').write_text(json.dumps(out,indent=2)+'\n');print('DONE',len(good),len(results),flush=True)
if __name__=='__main__':run()
