#!/usr/bin/env python3
"""Numerical coherent-face scout from the frozen118-wall checkpoint."""
import json,sys,hashlib
from pathlib import Path
from itertools import combinations,product
import numpy as np
from scipy.linalg import null_space
from scipy.spatial import HalfspaceIntersection
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_higher_face_20261013_scout as B
P=B.V.ROOT;PRE='fresh_h9_h5_closure_20261013';OLD='data/enumerations/fresh_h9_higher_face_20261013'
def load():
 _,_,_,_,_,names,_,_,R,G,c=B.load();names=list(names)
 for suffix,extra in[('_verify.json',[f'roof{i}'for i in range(7)]),('_truncation_verify.json',['adaptive0','adaptive1','secondAdaptive0']),('_final_verify.json',['roofPolar','roofMetric'])]:
  d=B.V.read(OLD+suffix)
  for f,h in d['sources'].items():assert B.V.sha(f)==h
  nr=np.array([[(int(a)+int(b))/(2*B.V.SCALE)for a,b in row]for row in d['normals']]);R=np.vstack([R,nr]);names+=extra
 assert len(R)==118
 final=d;fixed=final['H5_face']['fixed_walls'];F=R[[names.index(x)for x in fixed]];p=c-(c@G@F.T)@np.linalg.inv(F@G@F.T)@F;p/=np.sqrt(-p@G@p)
 C=null_space(np.vstack([F@G,p@G]));C=C@np.linalg.inv(np.linalg.cholesky(C.T@G@C).T)
 target=np.array([(int(a)+int(b))/(2*B.V.SCALE)for a,b in final['H5_face']['strict_spacelike_escape']['point']])
 return R,G,c,names,final,F,p,C,target
def gate(z):
 if z<-2+1e-7:return True
 return min(abs(z+2*np.cos(np.pi/m))for m in(2,3,4,5,6,7,8,10,12,16))<2e-7
def geometry(R,G,p,C):
 A=R@G@C;b=R@G@p;keep=np.linalg.norm(A,axis=1)>1e-7;H=np.c_[A[keep],b[keep]];vv=HalfspaceIntersection(H,np.zeros(5),qhull_options='Qx').intersections
 return [{'point':(p+C@v).tolist(),'norm':float(v@v-1),'active':np.flatnonzero(abs(R@G@(p+C@v))<1e-6).tolist()}for v in vv]
def run():
 R,G,c,names,final,F,p,C,target=load();verts=geometry(R,G,p,C);polars=[]
 for v in verts:
  if v['norm']<=1e-7:continue
  x=np.array(v['point']);r=np.sqrt(2/v['norm'])*x;ps=R@G@r;bad=[i for i,z in enumerate(ps)if not gate(z)]
  polars.append({'normal':r.tolist(),'norm':v['norm'],'active':[names[i]for i in v['active']],'bad':[{'wall':names[i],'pair':float(ps[i])}for i in bad],'target_cut':float(r@G@target)})
 polars.sort(key=lambda r:(len(r['bad']),-r['norm']))
 out={'status':'NUMERICAL_CURRENT_H5_GAP_DIAGNOSTIC','source_sha256':B.V.sha(OLD+'_final_verify.json'),'finite_vertices':sum(v['norm']<-1e-7 for v in verts),'spacelike_vertices':sum(v['norm']>1e-7 for v in verts),'maximum_norm':max(v['norm']for v in verts),'projective_vertices':verts,'polars':polars,'scope':'Numerical finite checkpoint diagnostic, not a compactness or family-exclusion theorem.'}
 (P/'data/enumerations'/f'{PRE}_scout.json').write_text(json.dumps(out,indent=2)+'\n');print('GEOM',len(verts),out['finite_vertices'],out['spacelike_vertices'],out['maximum_norm'],'bestpolars',[(r['active'],len(r['bad']),r['target_cut'])for r in polars[:8]],flush=True)
if __name__=='__main__':run()
