#!/usr/bin/env python3
"""Coherent metric endpoints selected by full current H5 geometry."""
import json
from collections import Counter
import numpy as np
import fresh_h9_h5_closure_next_20261013_scout as S
P=S.P;TAG=S.TAG
def run():
 R,G,c,names,doc,p,C,J,RF,phys,target=S.load();data=json.loads((P/'data/enumerations'/f'{TAG}_scout.json').read_text());rec=data['candidates'];NR=np.array([x['normal']for x in rec]);Ncov=NR@G;npair=Ncov@p;assert max(npair)<0;scaled=Ncov/(-npair[:,None]);initial=S.S.geometry(R,G,p,C);pool=np.arange(len(rec));starts=sorted(range(len(rec)),key=lambda i:(len(rec[i]['closed_old_edges'])+(rec[i]['target_cut']>1e-8)*2,rec[i]['target_cut']),reverse=True);seeds=[];seenprofiles=set()
 for i in starts:
  key=(rec[i]['edge_index'],tuple(rec[i]['profile']),rec[i]['anchor'])
  if key in seenprofiles:continue
  seenprofiles.add(key);seeds.append(i)
  if len(seeds)==8:break
 # Also begin with the strongest gap cutter irrespective of old-edge count.
 seeds.insert(0,max(range(len(rec)),key=lambda i:rec[i]['target_cut']))
 cases=[]
 for ci,seed in enumerate(seeds):
  chosen=[];allowed=np.ones(len(rec),bool);verts=initial;trace=[]
  for step in range(12):
   qq=np.array([x['norm']for x in verts]);order=np.argsort(qq)[::-1];worst=verts[int(order[0])];points=np.array([verts[int(i)]['point']for i in order[:60]]);eligible=np.flatnonzero(allowed)
   if step==0:proposals=[seed]
   else:
    cuts=scaled[eligible]@points.T;good=cuts[:,0]>1e-8;eligible=eligible[good];cuts=cuts[good]
    if not len(eligible):break
    weights=np.maximum(qq[order[:60]],0);score=np.maximum(cuts,0)@weights;inds=np.argsort(score)[-16:][::-1];proposals=eligible[inds].tolist()
   options=[]
   for i in proposals:
    vs=S.S.geometry(np.vstack([R,NR[chosen+[i]]]),G,p,C);mx=max(v['norm']for v in vs);f=sum(v['norm']<-1e-7 for v in vs);options.append((mx,-f,i,vs))
   mx,negf,i,verts=min(options,key=lambda x:(x[0],x[1]));chosen.append(i);allowed&=(Ncov@NR[i]<-2-2e-6);allowed[i]=False;trace.append({'candidate':i,'maximum_norm':mx,'finite_vertices':-negf,'spacelike_vertices':sum(v['norm']>1e-7 for v in verts),'remaining_compatible':int(sum(allowed))});print('case',ci,'step',step,'cap',i,'q',mx,'finite',-negf,'available',sum(allowed),flush=True)
   if mx< -1e-7:break
  covers=sorted(set(j for i in chosen for j in rec[i]['closed_old_edges']));leading=sorted(verts,key=lambda v:v['norm'],reverse=True)[:12];case={'selected_indices':chosen,'closed_original_edges':covers,'trace':trace,'leading_vertices':leading,'maximum_norm':trace[-1]['maximum_norm'],'finite_vertices':trace[-1]['finite_vertices']};cases.append(case)
  out={'status':'NUMERICAL_COHERENT_ENDPOINT_GEOMETRY_LOOKAHEAD','sources':{p:S.S.B.V.sha(p)for p in(S.OLD,'data/enumerations/'+TAG+'_scout.json','scripts/'+TAG+'_select.py')},'cases':cases,'scope':'Numerical compatible selections. Every selected cap retains123; new-new relations are required strictly dotted. Full H5 projective geometry guides each step. Exact checks and finite-volume certification remain required.'};(P/'data/enumerations'/f'{TAG}_selection.json').write_text(json.dumps(out,indent=2)+'\n')
 print('BEST',min(enumerate(cases),key=lambda x:x[1]['maximum_norm'])[0],min(x['maximum_norm']for x in cases),flush=True)
if __name__=='__main__':run()
