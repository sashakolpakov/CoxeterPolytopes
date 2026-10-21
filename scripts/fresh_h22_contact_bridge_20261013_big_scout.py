#!/usr/bin/env python3
"""Separate four-anchor extension of the frozen parent-bank interval scout."""
import json,math
from collections import Counter
import numpy as np
import fresh_h22_contact_bridge_20261013_scout as S
import fresh_h22_contact_bridge_20261013_continuous as C
TAG=S.TAG+'_big_scout';OUT=S.V.P/'results'/f'{TAG}.json'
def run():
 M,W,h,G,n,meta,doc,vertices=S.load();WM=W@M;primary=S.V.read('results/'+S.TAG+'_verify.json');anchors=[i for i,r in enumerate(meta)if r.get('height')==1 and r.get('type')in(0,2,6,22)];assert len(anchors)==4;rows=[];counts=Counter()
 for row in primary['rows']:
  if row.get('conclusion')=='EMPTY_ACUTE_UPPER_CAP_INTERVAL':continue
  vi,released=row['vertex'],row['released'];src=doc['compact_H2_vertices'][vi];inc=src['incident'];AI=np.linalg.inv(G[np.ix_(inc,inc)]);k=inc.index(released);v=vertices[vi]/math.sqrt(-src['norm']);e=(AI[k]@W[inc])/math.sqrt(AI[k,k]);a=WM@(v-e);b=WM@(v+e);lo=max(1.,int(row['lower_cutoff'][0])/int(primary['interval_denominator']));proofs=[]
  for anchor in anchors:
   counts['anchor_tests']+=1;ii=C.finite_interval(a[anchor],b[anchor],n[anchor]);assert ii is not None;al,ah=max(lo,ii[0]),ii[1];rec={'anchor':anchor,'type':meta[anchor]['type'],'T_interval':[al,ah]}
   if al>=ah:rec['kind']='EMPTY';counts['empty']+=1
   else:
    mn,mx=C.ranges(a,b,n,al,ah);bad=[]
    for m in range(2,100):
     inds=np.where((mn>math.cos(math.pi/m)+1e-7)&(mx<math.cos(math.pi/(m+1))-1e-7))[0]
     for j in inds:bad.append({'wall':int(j),'gap_labels':[m,m+1],'cosine_range':[float(mn[j]),float(mx[j])]})
    rec['kind']='FORBIDDEN_RANGE'if bad else'UNRESOLVED';counts[rec['kind']]+=1;rec['witness']=bad[0]if bad else None
   proofs.append(rec)
  rows.append({'vertex':vi,'released':released,'anchors':proofs})
 out={'status':'NUMERICAL_FOUR_REMAINING_BIG_COMPONENT_ANCHORS','sources':{p:S.V.sha(p)for p in('results/'+S.TAG+'_verify.json',S.P.OUT,'scripts/'+TAG+'.py')},'anchors':anchors,'counts':dict(counts),'rows':rows,'scope':'Separate four-anchor continuous diagnostic; does not modify the frozen301-parent certificate. Numerical only until its inequalities are replayed exactly.'}
 OUT.write_text(json.dumps(out,indent=2)+'\n');print(dict(counts),anchors)
if __name__=='__main__':run()
