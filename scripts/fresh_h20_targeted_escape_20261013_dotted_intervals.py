#!/usr/bin/env python3
"""Numerical continuous-height intervals with onlyright/dotted cap relations.

Unlike finite-angle events, these intervals do not force an extra anchor
equality. All possible sign changes are endpoints of explicit quadratics.
"""
import json
from pathlib import Path
import numpy as np
from fresh_h20_targeted_escape_20261013_second_volume import extended_setup
from fresh_h20_targeted_escape_20261013_vertex_bank import F,VERTICES
ROOT=Path(__file__).resolve().parents[1]
OUT='data/enumerations/fresh_h20_targeted_escape_20261013_dotted_intervals.json'
GEOM='data/enumerations/fresh_h20_targeted_escape_20261013_bank_geometry.json'

def run():
 data=json.loads((ROOT/VERTICES).read_text());geometry=json.loads((ROOT/GEOM).read_text());w=np.array(geometry['whole_bank_geometry']['point'])
 Q,R,h,_,_,_,_=extended_setup();cov=R@Q;ns=np.einsum('ij,jk,ik->i',R,Q,R);z=cov[:,-1]
 found=[];intervals=[];best=None
 for i,row in enumerate(data['vertices']):
  v=np.array([float(F(x))for x in row['projection']]+[0.]);scale=np.sqrt(-2/float(F(row['norm'])));v*=scale
  a=cov@v;ident=(abs(a)<1e-8)&(abs(z)<1e-8)
  aa=z*z-3*ns;bb=2*a*z;cc=a*a+2*ns
  critical=[-np.sqrt(2/3),np.sqrt(2/3)]
  for j in range(len(R)):
   if ident[j]:continue
   if abs(z[j])>1e-9:critical.append(-a[j]/z[j])
   if abs(aa[j])<1e-9:
    if abs(bb[j])>1e-9:critical.append(-cc[j]/bb[j])
   else:
    d=bb[j]**2-4*aa[j]*cc[j]
    if d>=-1e-8:
     critical.extend([(-bb[j]+s*np.sqrt(max(0,d)))/(2*aa[j])for s in(-1,1)])
  critical=sorted(set(round(float(x),12)for x in critical if np.isfinite(x)))
  critical=[critical[0]-max(1,abs(critical[0]))]+critical+[critical[-1]+max(1,abs(critical[-1]))]
  for lo,hi in zip(critical,critical[1:]):
   if hi-lo<1e-10:continue
   t=(lo+hi)/2;q=-2+3*t*t
   if q<=1e-9:continue
   ps=a+z*t
   if not np.all(ps[~ident]<-1e-8):continue
   if not np.all(ps[~ident]**2-ns[~ident]*q>1e-8):continue
   # The cut inequality is affine in t. Prefer a well-separated point
   # toward the endpoint that maximizes it, without claiming exactendpoints.
   wa=v@Q@w;wb=3*w[-1]
   choices=[lo+(hi-lo)*f for f in(.1,.25,.5,.75,.9)]
   t=max(choices,key=lambda t:wa+wb*t);q=-2+3*t*t;r=v.copy();r[-1]=t;r*=np.sqrt(2/q)
   cut=r@Q@w;record={'vertex_index':i,'interval':[lo,hi],'height':t,'normal':r.tolist(),
    'target_pairing':float(cut),'right_walls':np.flatnonzero(ident).tolist()}
   intervals.append(record)
   if best is None or cut>best['target_pairing']:best=record
   if cut>1e-7:found.append(record)
 out={'status':'NUMERICAL_CONTINUOUS_RIGHT_DOTTED_HEIGHT_INTERVAL_SCOUT','vertex_source':VERTICES,'target_source':GEOM,
  'intervals':intervals,'cutters':found,'best':best,
  'scope':'All2046shortvertexdirections, numeric quadratic criticalpoints partitionheightline forright/dottedrelationsagainst87walls. Eachsaved openintervalpasses strictmiddletest; finiteprecision endpoints andnodiscovery/volumeclaim. This explorescontinuousfamilies omittedbyfiniteangleevents.'}
 (ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n');print({'intervals':len(intervals),'cutters':len(found),'best':None if best is None else{ k:best[k]for k in('vertex_index','interval','target_pairing')}})

if __name__=='__main__':run()
