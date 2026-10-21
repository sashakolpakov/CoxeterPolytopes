#!/usr/bin/env python3
"""Actual H5 remote projective vertices, and their polar truncation test."""
import json,sys
from pathlib import Path
from itertools import combinations
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.linalg import null_space
from scipy.spatial import HalfspaceIntersection
from scripts import fresh_h9_higher_face_20261013_linked_scout as L
S=L.S

def run():
 _,_,_,_,_,names,_,_,R,G,c=S.load();d=S.V.read(S.PREFIX+'_selected_scout.json');rr=np.array([x['normal']for x in d['recipes']]);R=np.vstack([R,rr]);names=names+[f'roof{i}'for i in range(len(rr))];fixed=[names.index(n)for n in S.COMMON if n!='seed1'];F=R[fixed];p=c-(c@G@F.T)@np.linalg.inv(F@G@F.T)@F;p/=np.sqrt(-p@G@p);B=null_space(np.vstack([F@G,p@G]));M=B.T@G@B;B=B@np.linalg.inv(np.linalg.cholesky(M).T);A=R@G@B;b=R@G@p;keep=np.linalg.norm(A,axis=1)>1e-7;H=np.c_[A[keep],b[keep]];vv=HalfspaceIntersection(H,np.zeros(5),qhull_options='Qx').intersections;rows=[];valid=[]
 for v in vv:
  x=p+B@v;q=x@G@x
  if q<1e-7:continue
  r=np.sqrt(2/q)*x;pair=R@G@r;bad=[{'wall':names[i],'pair':float(y)}for i,y in enumerate(pair)if not L.gate(y)];row={'point':x.tolist(),'norm':float(q),'incident':[names[i]for i,y in enumerate(R@G@x)if abs(y)<1e-6],'bad':bad,'normal':r.tolist()};rows.append(row)
  if not bad:valid.append(r)
 rows.sort(key=lambda x:(len(x['bad']),-x['norm']));out={'status':'NUMERICAL_POLAR_SCREEN','source':S.PREFIX+'_selected_scout.json','rows':rows,'accepted':len(valid)}
 if valid:out['geometry_all_accepted']=L.geom(np.vstack([R,valid]),G,c,fixed)
 (S.V.ROOT/(S.PREFIX+'_polar_scout.json')).write_text(json.dumps(out,indent=2)+'\n');print('accepted',len(valid),'best',[(len(x['bad']),x['norm'],x['incident'],x['bad'])for x in rows[:5]])
if __name__=='__main__':run()
