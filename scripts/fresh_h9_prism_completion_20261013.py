#!/usr/bin/env python3
"""Freeze the symmetry-coupled remote insertion and test whole-chamber volume."""
import json,sys
from pathlib import Path
from fractions import Fraction as Q
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.spatial import HalfspaceIntersection
from scripts import fresh_h9_prism_completion_20261013_verify as T
V=T.V;K=T.K

def run():
    source=T.B.PREFIX+'_caps.json';prior=V.read(source)
    d={'source':source,'source_sha256':V.sha(source),'permutation':T.PERM,'source_normals':T.SOURCE_NAMES,'target':prior['escape']}
    path=T.PREFIX+'_caps.json';(V.ROOT/path).write_text(json.dumps(d,indent=2)+'\n')
    out=T.check(d);print('INITIAL',out['new_pair_counts'],out['new_component_frontier'],flush=True)
    G,I,h,g,old,names,roots,center,_,_=T.reconstruct(d);ff=lambda x:sum(x)/(2*V.SCALE)
    R=np.array([[ff(x)for x in r]for r in old+roots]);In=np.array([[ff(V.golden(x))for x in row]for row in I]);hn=np.array([ff(V.golden(x))for x in h]);hn/=hn.sum()
    H=np.c_[R[:,9,None]-R[:,:9],-R[:,9]];H/=np.linalg.norm(H[:,:9],axis=1)[:,None]
    vv=HalfspaceIntersection(H,hn[:9],qhull_options='Qx').intersections;ys=np.c_[vv,1-vv.sum(axis=1)];qs=np.einsum('bi,ij,bj->b',ys,In,ys)
    found=None
    for j in np.argsort(-qs):
        if qs[j]<=0:break
        for eps in(1e-4,1e-6,1e-8):
            y=tuple((Q(round(float(x)*10**12),10**12),Q(0))for x in ys[j]+eps*hn)
            if K.sign(K.pair(y,I,y))<=0 or K.sign(K.pair(y,I,h))>=0:continue
            if min(V.dot(r,list(map(V.golden,y)))[0]for r in old+roots)>0:found=(y,int(j));break
        if found:break
    assert found,'No escape found: actual complete volume verification required.'
    y,j=found;d['escape']=list(map(K.enc,y));d['numerical_geometry']={'vertices':len(ys),'spacelike_vertices':int((qs>1e-8).sum()),'maximum_norm':float(qs.max()),'escape_proposal_support':[int(i)for i in np.flatnonzero(abs(ys[j])>1e-7)],'escape_active_walls':[(names+T.NEW)[i]for i in np.flatnonzero(abs(R@ys[j])<1e-7)],'scope':'Numerical vertices only propose an exactly verified strict escape.'}
    (V.ROOT/path).write_text(json.dumps(d,indent=2)+'\n');out=T.check(d)
    (V.ROOT/(T.PREFIX+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
    print('FINAL',out['status'],d['numerical_geometry'],flush=True)

if __name__=='__main__':run()
