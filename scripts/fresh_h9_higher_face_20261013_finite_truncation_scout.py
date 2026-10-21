#!/usr/bin/env python3
"""Finite contact profiles on the entire five-wall remote face, not a fixed all-right plane."""
import json,sys
from itertools import product
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scripts import fresh_h9_higher_face_20261013_linked_scout as L
S=L.S

def run():
 _,_,_,_,_,names,_,_,R,G,c=S.load();d=S.V.read(S.PREFIX+'_selected_scout.json');R=np.vstack([R,np.array([x['normal']for x in d['recipes']])]);names+=['roof'+str(i)for i in range(7)];fixed=[n for n in S.COMMON if n!='seed1'];pd=S.V.read(S.PREFIX+'_polar_scout.json');labels=np.array(list(product((2,3,4,5,6,8,10,12,16),repeat=5)));bb=-2*np.cos(np.pi/labels);out=[];stats=[]
 for vi,row in enumerate(sorted(pd['rows'],key=lambda x:-x['norm'])):
  other=[n for n in row['incident']if n not in fixed]
  if len(other)!=5:continue
  E=R[[names.index(n)for n in fixed+other]];M=E@G@E.T;Mi=np.linalg.inv(M);bs=np.c_[np.zeros((len(bb),4)),bb];co=bs@Mi;part=co@E;qp=np.einsum('bi,bi->b',bs,co);ok=qp<2-1e-8;idx=np.where(ok)[0];qz=row['norm'];z=np.array(row['point'])/np.sqrt(qz);rr=part[ok]+np.sqrt(2-qp[ok])[:,None]*z
  ps=rr@G@R.T;ms=np.maximum(2,np.rint(np.pi/np.arccos(np.clip(-ps/2,-1+1e-15,1-1e-15))));finite=(ps<=0)&(ps>=-2)&(np.abs(ps+2*np.cos(np.pi/ms))<1e-7);valid=np.all((ps<-2+1e-7)|(np.abs(ps)<1e-7)|finite,axis=1)&(rr@G@c<0)
  keep=np.where(valid)[0]
  keep=sorted(keep,key=lambda k:(sum(labels[idx[k]]!=2),max(labels[idx[k]]),sum(labels[idx[k]])))[:100]
  for k in keep:out.append({'vertex':vi,'right_fixed':fixed,'incident':other,'profile':labels[idx[k]].tolist(),'normal':rr[k].tolist(),'target':row['point'],'target_cut':float(rr[k]@G@np.array(row['point']))})
  stats.append({'vertex':vi,'profiles':len(labels),'norm_equation_roots':len(idx),'accepted':int(valid.sum()),'archived_simplest':len(keep)});print(stats[-1],flush=True)
 (S.V.ROOT/(S.PREFIX+'_finite_truncation_scout.json')).write_text(json.dumps({'status':'NUMERICAL_FULL_REMOTE_FACE_CONTACT_PROFILES','labels':[2,3,4,5,6,8,10,12,16],'stats':stats,'candidates':out},indent=2)+'\n');print('total',len(out),flush=True)
if __name__=='__main__':run()
