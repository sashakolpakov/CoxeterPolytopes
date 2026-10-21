#!/usr/bin/env python3
"""A coupled four-cap remote vertex insertion, followed by an exact escape gate."""
import json
import numpy as np
from scipy.spatial import HalfspaceIntersection
from scripts.fresh_h9_remote_caps_20261013_base import *
from scripts import fresh_h9_remote_caps_20261013_verify as V

def run():
    G,I,h,old,geometry,target,frame,norms,chosen,wp=source_frame();k=64
    ds=[div(scale(O,Q(1,2)),K.neg(norms[0])if i==0 else q)for i,q in enumerate(norms)]
    data={'status':'EXACT_COUPLED_REMOTE_VERTEX_INSERTION','source':V.OLD,'source_sha256':sha(V.OLD),
          'target_source':REMOTE,'target_source_sha256':sha(REMOTE),'target_point':list(map(enc,target)),
          'frame':[[enc(x)for x in row]for row in frame],'frame_norms':list(map(enc,norms)),
          'orthogonal_old_walls':list(S_NAMES),'gram_schmidt_sources':chosen,
          'projected_target_norm':enc(pairing(K,wp,G,wp)),'boost':str(k),'radicands':list(map(enc,ds)),
          'signs':list(map(list,SIGNS)),
          'scope':'Four coupled caps in the rank-five spherical-face complement, retaining all79 existing walls. This is one explicit parameter choice, not an exhaustive family search.'}
    cert=V.check(data);print({'exact_new_caps':4,'actual_finite_vertex':True},flush=True)
    _,_,_,_,_,caps,point,old_roots=V.reconstruct(data)
    C=np.array([[(x[0]+x[1])/(2*V.SCALE)for x in r]for r in old_roots+caps]);In=np.array([[float(value(K,x))for x in r]for r in I])
    hn=np.array([float(value(K,x))for x in h]);hn/=hn.sum()
    half=np.c_[C[:,9,None]-C[:,:9],-C[:,9]];half/=np.linalg.norm(half[:,:9],axis=1)[:,None]
    vertices=HalfspaceIntersection(half,hn[:9],qhull_options='Qx').intersections
    ys=np.c_[vertices,1-vertices.sum(axis=1)];qs=np.einsum('bi,ij,bj->b',ys,In,ys)
    accepted=None
    for index in np.argsort(-qs):
        if qs[index]<=0:break
        for epsilon in (1e-4,1e-6,1e-8):
            yp=ys[index]+epsilon*hn;y=tuple(scale(O,Q(round(float(x)*10**12),10**12))for x in yp)
            if K.sign(pairing(K,y,I,y))<=0 or K.sign(pairing(K,y,I,h))>=0:continue
            yi=list(map(V.golden,y));vals=[V.dot(r,yi)for r in old_roots+caps]
            if min(x[0]for x in vals)>0:accepted=(y,index);break
        if accepted:break
    assert accepted,'No exact escape found; investigate full volume independently.'
    y,index=accepted;data['escape']=list(map(enc,y));data['numerical_geometry']={
        'vertices':len(ys),'spacelike_vertices':int((qs>1e-8).sum()),'maximum_vertex_norm':float(qs.max()),
        'escape_proposal_support':[int(i)for i in np.flatnonzero(ys[index]>1e-7)],
        'scope':'Numerical vertex data propose an exact strict escape; absence of a numerical escape would not prove finite volume.'}
    cert=V.check(data)
    (ROOT/(PREFIX+'_caps.json')).write_text(json.dumps(data,indent=2)+'\n')
    cert['source_sha256']=sha(PREFIX+'_caps.json');(ROOT/(PREFIX+'_verify.json')).write_text(json.dumps(cert,indent=2)+'\n')
    print({'facets':83,'new_actual_vertex':'9A1','geometry':'EXACT_STRICT_ESCAPE','numerical':data['numerical_geometry']},flush=True)

if __name__=='__main__':run()
