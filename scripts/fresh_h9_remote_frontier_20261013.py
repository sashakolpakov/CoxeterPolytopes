#!/usr/bin/env python3
"""Freeze the new compact face and test the whole105-wall chamber."""
import json,sys
from pathlib import Path
from fractions import Fraction as Q
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.spatial import HalfspaceIntersection
from scripts import fresh_h9_remote_frontier_20261013_verify as F
V=F.V;K=F.K

def run():
 d=V.read(F.PREFIX+'_caps.json');out=F.check(d);print('EXACT FACE PASS',flush=True)
 G,I,h,g,old,names,c,pairs,_=F.old_state();r,_=F.make(d['branch'],g,old,names,c,pairs);allroots=old+[r]
 ff=lambda x:sum(x)/(2*V.SCALE)
 R=np.array([[ff(x)for x in row]for row in allroots]);In=np.array([[ff(V.golden(x))for x in row]for row in I]);hn=np.array([ff(V.golden(x))for x in h]);hn/=hn.sum()
 H=np.c_[R[:,9,None]-R[:,:9],-R[:,9]];H/=np.linalg.norm(H[:,:9],axis=1)[:,None]
 vv=HalfspaceIntersection(H,hn[:9],qhull_options='Qx').intersections;ys=np.c_[vv,1-vv.sum(axis=1)];qs=np.einsum('bi,ij,bj->b',ys,In,ys);found=None
 for j in np.argsort(-qs):
  if qs[j]<=0:break
  for eps in(1e-4,1e-6,1e-8):
   y=tuple((Q(round(float(x)*10**12),10**12),Q(0))for x in ys[j]+eps*hn)
   if K.sign(K.pair(y,I,y))<=0 or K.sign(K.pair(y,I,h))>=0:continue
   if min(V.dot(a,list(map(V.golden,y)))[0]for a in allroots)>0:found=(y,int(j));break
  if found:break
 assert found,'No escape found; a separate full-volume proof is required.'
 y,j=found;d['escape']=list(map(K.enc,y));d['numerical_geometry']={'projective_vertices':len(ys),'spacelike_vertices':int((qs>1e-8).sum()),'maximum_norm':float(qs.max()),'escape_active_walls':[(names+[F.NEW])[k]for k in np.flatnonzero(abs(R@ys[j])<1e-7)],'scope':'Numerical census only proposes the exact replayed strict escape.'}
 (V.ROOT/(F.PREFIX+'_caps.json')).write_text(json.dumps(d,indent=2)+'\n');out=F.check(d);(V.ROOT/(F.PREFIX+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
 print('FINAL',out['status'],d['numerical_geometry'],flush=True)

if __name__=='__main__':run()
