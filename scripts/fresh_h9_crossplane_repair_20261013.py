#!/usr/bin/env python3
"""Produce the exact compact-compatible metric repair and test its chamber."""
import json
from fractions import Fraction as Q
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.spatial import HalfspaceIntersection
from scripts import fresh_h9_crossplane_repair_20261013_verify as T
V=T.V;K=V.K

def run():
    source=T.B.PREFIX+'_caps.json';data={'status':'EXACT_CROSSPLANE_METRIC_REPAIR','source':source,'source_sha256':V.sha(source),
        'constraints':T.CONSTRAINTS,'branch':-1,'target_point':V.read(source)['escape'],
        'definition':{'A':'Nine retained normals in the stated constraint order.','M':'Gram(A).','p':'The prescribed pairings0,-1,-phi.',
            'u':'sum_i (M^-1 p)_i A_i','v':'center minus its orthogonal projection onto span(A)',
            'normal':'u-sqrt((2-u²)/v²)*v'},
        'scope':'Adaptive replacement of finite-angle equations yields one exact cap. All88 old walls retained; no exhaustive cap-family claim.'}
    result=T.check(data);print({'exact_cap':True,'affines':result['affine_and_vertex_census']['new_affine_circuits'],
        'new_vertices':len(result['affine_and_vertex_census']['new_finite_vertices'])},flush=True)
    G,I,h,g,old,names,r,center,_=T.reconstruct(data)
    f=lambda p:sum(p)/(2*V.SCALE)
    C=np.array([[f(p)for p in row]for row in old+[r]]);In=np.array([[f(V.golden(x))for x in row]for row in I]);hn=np.array([f(V.golden(x))for x in h]);hn/=hn.sum()
    half=np.c_[C[:,9,None]-C[:,:9],-C[:,9]];half/=np.linalg.norm(half[:,:9],axis=1)[:,None]
    vv=HalfspaceIntersection(half,hn[:9],qhull_options='Qx').intersections;ys=np.c_[vv,1-vv.sum(axis=1)];qs=np.einsum('bi,ij,bj->b',ys,In,ys);accepted=None
    for j in np.argsort(-qs):
        if qs[j]<=0:break
        for eps in(1e-4,1e-6,1e-8):
            p=ys[j]+eps*hn;y=tuple((Q(round(float(x)*10**12),10**12),Q(0))for x in p)
            if K.sign(K.pair(y,I,y))<=0 or K.sign(K.pair(y,I,h))>=0:continue
            values=[V.dot(a,list(map(V.golden,y)))for a in old+[r]]
            if min(p[0]for p in values)>0:accepted=(y,j);break
        if accepted:break
    assert accepted,'No exact escape found; investigate complete volume.'
    y,j=accepted;data['escape']=list(map(K.enc,y));data['numerical_geometry']={'vertices':len(ys),'spacelike_vertices':int((qs>1e-8).sum()),
        'maximum_vertex_norm':float(qs.max()),'escape_proposal_support':[int(i)for i in np.flatnonzero(ys[j]>1e-7)],
        'scope':'Numerical geometry only proposes the exact strict escape; no compactness claim from this census.'}
    result=T.check(data);path=T.PREFIX+'_caps.json';(V.ROOT/path).write_text(json.dumps(data,indent=2)+'\n');result['source_sha256']=V.sha(path)
    (V.ROOT/(T.PREFIX+'_verify.json')).write_text(json.dumps(result,indent=2)+'\n');print({'facets':89,'status':'EXACT_ESCAPE','geometry':data['numerical_geometry']},flush=True)

if __name__=='__main__':run()
