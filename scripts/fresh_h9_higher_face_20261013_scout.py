#!/usr/bin/env python3
"""Numerical all-right H4-corner endpoint families for five H5 face lifts."""
import json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_remote_frontier_20261013_second_verify as S
V=S.V;K=S.K;N=S.N
PREFIX='data/enumerations/fresh_h9_higher_face_20261013'
COMMON=['seed1','seed2','compact1','compact2','compact3']
H4=['crosspair0','nextH4','compact0','facePrism']
def load():
 d=V.read(S.PREFIX+'_caps.json');G,I,h,g,old,names,c,pairs,_=S.old_state();r,_=S.make(d['branch'],g,old,names,c,pairs)
 for n in names:pairs['|'.join(sorted((S.NEW,n))) ]={'kind':'finite','label':S.PROFILE[n]}if n in S.PROFILE else{'kind':'dotted'}
 old=old+[r];names=names+[S.NEW];ff=lambda z:sum(z)/(2*V.SCALE)
 return G,I,h,g,old,names,c,pairs,np.array([[ff(z)for z in row]for row in old]),np.array([[ff(V.golden(z))for z in row]for row in G]),np.array([ff(z)for z in c])
def run():
 G,I,h,g,old,names,c,pairs,R,Gn,cn=load();b=R[names.index('wall17')];A=R[[names.index(n)for n in H4]];M=A@Gn@A.T
 bp=b-(b@Gn@A.T)@np.linalg.inv(M)@A;qb=bp@Gn@bp;assert qb<0
 v=bp/np.sqrt(-qb);AA=R@Gn@v;records=[];found=[]
 for free in COMMON:
  f=R[names.index(free)];uu=f/np.sqrt(2);BB=R@Gn@uu;edge=H4+[n for n in COMMON if n!=free]
  low=0.;upper=1.;bounds=[]
  for j,(a,bb)in enumerate(zip(AA,BB)):
   if names[j]in edge:continue
   if abs(a)<1e-8:
    if bb<-1e-8:upper=-1;bounds.append({'wall':names[j],'impossible':True})
   elif a<0:
    lo=bb/a
    if lo>low:low=lo
    if lo>0:bounds.append({'wall':names[j],'lower_tanh':float(lo)})
   else:upper=min(upper,bb/a)
  row={'free_common_wall':free,'edge':edge,'tanh_necessary_interval':[low,upper],'strongest_bounds':sorted(bounds,key=lambda z:z.get('lower_tanh',1e9),reverse=True)[:6]}
  contacts=[]
  for m in range(16,257):
   x=2*np.cos(np.pi/m);rr=(-x/qb)*bp-np.sqrt(1-x*x/(2*qb))*f;ps=R@Gn@rr
   if max(ps)>1e-7:continue
   bad=[]
   for j,p in enumerate(ps):
    if names[j]in edge or names[j]=='wall17':continue
    if p>-2+1e-8 and abs(p)>1e-7:
     ang=np.pi/np.arccos(np.clip(-p/2,-1,1));mi=int(round(ang))
     if mi<2 or abs(p+2*np.cos(np.pi/mi))>1e-7:bad.append([names[j],float(p)])
   contacts.append({'label':m,'bad_pairs':bad})
   if not bad:found.append({'free':free,'label':m,'normal':rr.tolist()})
  row['finite_wall17_angle_trials']=contacts;records.append(row)
 out={'status':'NUMERICAL_HIGHER_FACE_ENDPOINT_SCREEN','source':S.PREFIX+'_caps.json','source_sha256':V.sha(S.PREFIX+'_caps.json'),'H4':H4,'common_A1':COMMON,'wall17_projection_norm':float(qb),'seed6_to_wall17_pair_ratio':float((R[names.index('seed6')]@Gn@bp)/qb),'continuous_families':records,'proposals':found,'scope':'Numerical acuteness intervals for all-right H4+4A1 endpoint caps. Finite wall17 labels16..256 are a bounded scout only; no exhaustive finite-label claim.'}
 (V.ROOT/(PREFIX+'_scout.json')).write_text(json.dumps(out,indent=2)+'\n')
 print('qb',qb,'proposals',len(found));print([(r['free_common_wall'],r['tanh_necessary_interval'],r['strongest_bounds'][:2],len(r['finite_wall17_angle_trials']))for r in records])
if __name__=='__main__':run()
