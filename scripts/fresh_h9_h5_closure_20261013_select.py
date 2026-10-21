#!/usr/bin/env python3
"""Bounded coherent cap selection driven by full H5 projective geometry."""
import argparse,json
from itertools import combinations
import numpy as np
from scipy.optimize import linprog
import fresh_h9_h5_closure_20261013_scout as S
P=S.P;PRE=S.PRE
def run(stage=0):
 R,G,c,names,final,F,p,C,target=S.load();suffix=''if stage==0 else f'_stage{stage}';data=json.loads((P/'data/enumerations'/f'{PRE}_conic_scout{suffix}.json').read_text());caps=data['candidates'];N=np.array([r['normal']for r in caps]);bad=[]
 if stage:
  state=json.loads((P/'data/enumerations'/f'{PRE}_state{stage}.json').read_text());R=np.vstack([R,np.array([r['normal']for r in state['selected']])]);names+=[r['name']for r in state['selected']]
 oldG=R@G@R.T
 for ci,rec in enumerate(caps):
  ids=[names.index(n)for n in rec['four']+[rec['anchor']]];B=oldG[np.ix_(ids,ids)];bb=-2*np.cos(np.pi/np.array(rec['profile']+[rec['label']]));A=np.block([[np.array([[2.]]),bb[None,:]],[bb[:,None],B]])
  affine=False
  for mask in range(1,32):
   ix=[0]+[j+1 for j in range(5)if mask>>j&1];M=A[np.ix_(ix,ix)]
   if np.any(M<-2-1e-6):continue
   eig=np.linalg.eigvalsh(M)
   if abs(eig[0])<1e-7:affine=True;break
  if affine:bad.append(ci)
 good=[i for i in range(len(N))if i not in set(bad)];unfiltered_good=len(good)
 if len(good)>1400:
  vv=S.geometry(R,G,p,C);pt=np.array([v['point']for v in vv if v['norm']>1e-7]);cs=N[good]@G@pt.T;score=np.sum(cs>1e-7,axis=1);ranked=sorted(range(len(good)),key=lambda i:(-score[i],-caps[good[i]]['target_cut']));keep=ranked[:1000];seen=set()
  for i in ranked:
   r=caps[good[i]];key=(tuple(r['four']),tuple(r['profile']),r['anchor'])
   if key in seen:continue
   seen.add(key);keep.append(i)
   if len(keep)>=1400:break
  good=[good[i]for i in dict.fromkeys(keep)]
 NN=N[good];Pnew=NN@G@NN.T;compat=Pnew<-2-1e-6
 for m in(2,3,4,5,6,8,10,12,16):compat|=abs(Pnew+2*np.cos(np.pi/m))<1e-7
 np.fill_diagonal(compat,False);print('AFFINEFREE',len(good),'edges',int(compat.sum()//2),flush=True)
 def geom(ids):
  rr=np.vstack([R,NN[ids]])if ids else R;vs=S.geometry(rr,G,p,C);vs.sort(key=lambda r:-r['norm']);return vs
 initial=geom([]);pts=np.array([v['point']for v in initial if v['norm']>1e-7]);cuts=NN@G@pts.T
 starts=sorted(range(len(good)),key=lambda i:(-np.sum(cuts[i]>1e-7),-caps[good[i]]['target_cut']))[:30]
 # Diversify the first step across literal finite-contact profiles.
 seen=set();extra=[]
 for i in sorted(range(len(good)),key=lambda i:-caps[good[i]]['target_cut']):
  rec=caps[good[i]];key=(tuple(rec['four']),tuple(rec['profile']),rec['anchor'],rec['label'])
  if key in seen:continue
  seen.add(key);extra.append(i)
  if len(extra)==20:break
 starts=list(dict.fromkeys(starts+extra));runs=[]
 for si,start in enumerate(starts):
  chosen=[start];available=set(np.flatnonzero(compat[start]));vs=geom(chosen);trace=[{'added':good[start],'maximum_norm':vs[0]['norm'],'projective_vertices':len(vs)}]
  for step in range(9):
   if vs[0]['norm']<=-1e-7 or not available:break
   worst=np.array(vs[0]['point']);cutters=[i for i in available if NN[i]@G@worst>1e-7]
   if not cutters:break
   points=np.array([v['point']for v in vs if v['norm']>1e-7]);scores=NN[cutters]@G@points.T;ordered=sorted(range(len(cutters)),key=lambda j:(-np.sum(scores[j]>1e-7),-scores[j,0]))[:12];trials=[]
   for j in ordered:
    i=cutters[j];nv=geom(chosen+[i]);trials.append((nv[0]['norm'],sum(v['norm']>1e-7 for v in nv),i,nv))
   _,_,best,nv=min(trials,key=lambda z:z[:3]);chosen.append(best);available&=set(np.flatnonzero(compat[best]));vs=nv;trace.append({'added':good[best],'maximum_norm':vs[0]['norm'],'projective_vertices':len(vs)})
  runs.append({'selected':[good[i]for i in chosen],'trace':trace,'maximum_norm':vs[0]['norm'],'spacelike_vertices':sum(v['norm']>1e-7 for v in vs),'finite_vertices':sum(v['norm']<-1e-7 for v in vs),'leading_vertices':vs[:12]});print('RUN',si,'caps',len(chosen),'maxq',vs[0]['norm'],'positive',runs[-1]['spacelike_vertices'],flush=True)
  if vs[0]['norm']<-1e-7:break
 runs.sort(key=lambda r:(r['maximum_norm'],r['spacelike_vertices']));out={'status':'BOUNDED_NUMERICAL_COHERENT_H5_CAP_SELECTION','source_sha256':S.B.V.sha('data/enumerations/'+PRE+'_conic_scout'+suffix+'.json'),'individually_affine_profiles_rejected':bad,'all_numerically_affine_free_candidates':unfiltered_good,'remaining_candidates':good,'compatibility_edges':int(compat.sum()//2),'runs':runs,'scope':'Numerical angle/affine and complete projective-vertex geometry proposals. Above1400 candidates, heuristic coverage/profile selection reduces the graph. No exact candidate or compactness claim until replay.'}
 out['stage']=stage
 (P/'data/enumerations'/f'{PRE}_select{suffix}.json').write_text(json.dumps(out,indent=2)+'\n');print('BEST',runs[0]['selected'],runs[0]['maximum_norm'],flush=True)
 if stage:
  best=runs[0];chosen=[]
  for j,i in enumerate(best['selected']):
   rec=dict(caps[i]);rec['name']='closure'+str(len(state['selected'])+j);chosen.append(rec)
  fixed=final['H5_face']['fixed_walls'];newnames=names+[r['name']for r in chosen];v=best['leading_vertices'][0];ns={'selected':state['selected']+chosen,'target_point':v['point'],'support':[newnames[i]for i in v['active']if newnames[i]not in fixed],'prior_geometry':best};(P/'data/enumerations'/f'{PRE}_state{stage+1}.json').write_text(json.dumps(ns,indent=2)+'\n')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--stage',type=int,default=0);a=ap.parse_args();run(a.stage)
