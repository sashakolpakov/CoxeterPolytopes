#!/usr/bin/env python3
"""Whole-cone gate for all six-cap images under the full seed action."""
import json,sys,itertools
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.optimize import linprog
from scripts import fresh_h9_prism_completion_20261013_scout as S
V=S.V

def run():
    d=S.setup();G=np.array(d['gram']);R=np.array(d['normals']);basis=list(range(8))+[8,12];allimages=[]
    for p,q,bs,uv in itertools.product(itertools.permutations(range(4)),itertools.permutations(range(4)),range(2),range(2)):
        perm=[((i//8)^uv)*8+(((i%8)//4)^bs)*4+(p if (i%8)//4==0 else q)[i%4]for i in range(16)]
        M=R[[perm[i]for i in basis]];allimages.extend(R[90:]@M)
    images=np.array(allimages);_,inds=np.unique(np.round(images,6),axis=0,return_index=True);images=images[np.sort(inds)]
    U=[0,2,1,3,4,6,5,7,8,9];A=np.vstack([R,R[90:,U],images]);I=np.linalg.inv(G);rng=np.random.default_rng(20261013);best=None;history=[]
    for start in range(40):
        objective=rng.normal(size=10)
        for step in range(15):
            res=linprog(-objective,A_ub=-A,b_ub=np.zeros(len(A)),A_eq=np.ones((1,10)),b_eq=[1],bounds=[(None,None)]*10,method='highs')
            assert res.success,res.message
            y=res.x;value=float(y@I@y)
            if best is None or value>best['norm']:
                best={'norm':value,'point':y.tolist(),'start':start,'step':step,'active_old':np.flatnonzero(abs(A[:102]@y)<1e-7).tolist(),'all_normal_minimum':float((A@y).min())};print('BEST',best,flush=True)
            nxt=I@y
            if np.linalg.norm(nxt-objective)<1e-10:break
            objective=nxt
        history.append({'start':start,'norm':value})
        if best['norm']>1e-6:break
    out={'status':'NUMERICAL_FULL_SEED_ORBIT_CONE_GATE','source':S.PREFIX+'_caps.json','source_sha256':V.sha(S.PREFIX+'_caps.json'),'group_elements':2304,'raw_images':13824,'numerically_distinct_images':len(images),'best':best,'history':history,
         'scope':'Every six-cap image under the full16-seed action is included before compatibility filtering. Numerical LP only proposes a possible exact common escape; failure to find one would not prove compactness.'}
    (V.ROOT/(S.PREFIX+'_full_orbit_geometry.json')).write_text(json.dumps(out,indent=2)+'\n');print('DONE',out['numerically_distinct_images'],best['norm'])

if __name__=='__main__':run()
