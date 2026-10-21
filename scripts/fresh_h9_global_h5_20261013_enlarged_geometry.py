#!/usr/bin/env python3
"""Numerical full H5 relaxation, deliberately ignoring all conflicts."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.spatial import HalfspaceIntersection
sys.path.insert(0,str(Path(__file__).resolve().parent))
import fresh_h9_global_h5_20261013_scout as S
def run():
 R,G,c,names,d,p,C,NR,refs,sources,counts=S.load();path='data/enumerations/'+S.TAG+'_opposite_scout.json';rows=json.loads((S.P/path).read_text())['candidates'];W=np.vstack([R,NR,np.array([r['normal']for r in rows])]);A=W@G@C;b=W@G@p;keep=np.linalg.norm(A,axis=1)>1e-7;assert max(b[keep])<0;H=np.c_[A[keep]/(-b[keep,None]),-np.ones(sum(keep))]
 hs=HalfspaceIntersection(H,np.zeros(5),qhull_options='Qx');v=hs.intersections;q=np.sum(v*v,axis=1)-1;i=int(np.argmax(q));offset=float(max(hs.dual_equations[:,-1]));assert offset<-1e-8
 out={'status':'NUMERICAL_COMPLETE_ENLARGED_BANK_HALFSPACE_GEOMETRY','halfspaces':len(H),'vertices':len(v),'finite':int(sum(q<-1e-7)),'spacelike':int(sum(q>1e-7)),'max_norm':float(q[i]),'max_point':(p+C@v[i]).tolist(),'active_entries':np.flatnonzero(abs(W@G@(p+C@v[i]))<1e-6).tolist(),'affine_slice_bounded':True,'dual_origin_maximum_facet_offset':offset,'scope':'Floating point complete projective vertex census and dual-hull boundedness test only. Conflicts are ignored; this is not a Coxeter polytope or an exact compactness certificate.'}
 (S.P/'data/enumerations'/f'{S.TAG}_enlarged_geometry.json').write_text(json.dumps(out,indent=2)+'\n');print(out,flush=True)
if __name__=='__main__':run()
