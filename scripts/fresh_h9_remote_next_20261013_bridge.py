#!/usr/bin/env python3
"""Exact circle-coupled bridge proposal and volume gate."""
import json
from fractions import Fraction as Q
import numpy as np
from scipy.spatial import HalfspaceIntersection
from scripts import fresh_h9_remote_next_20261013_bridge_verify as B
V=B.V;K=V.K

def run():
    source=B.N.PREFIX+'_caps.json';data={'status':'EXACT_REMOTE_INDEFINITE_FACE_BRIDGE','source':source,'source_sha256':V.sha(source),
        'orthogonal_old_walls':B.N.S_NAMES+['wall22','nextH4'],'finite_anchor':'wall38','finite_label':3,'branch':-1,
        'target_point':V.read(source)['escape'],
        'definition':{'plane':'Orthogonal complement of S6+wall22+nextH4, whose Gram has signature(7,1).',
                      'u':'Orthogonal projection of wall38 onto the positive normal plane.',
                      'v':'Projected strict center minus its component along u.',
                      'd':'u^2','e':'v^2','q':'2-1/d',
                      'bridge0':'-u/d-sqrt(q/e)*v','bridge1':'-sqrt(q/d)*u+v/sqrt(d*e)'},
        'scope':'One exact orthogonal cap pair fixed by a pi3 anchor relation, preserving the86-wall chamber.'}
    result=B.check(data);print({'pair':True,'two_actual_vertices':True,'all86_retained':True},flush=True)
    G,I,h,g,roots,names,caps,center,params=B.reconstruct(data)
    C=np.array([[(p[0]+p[1])/(2*V.SCALE)for p in row]for row in roots+caps]);In=np.array([[sum(V.golden(x))/(2*V.SCALE)for x in row]for row in I])
    hn=np.array([sum(V.golden(x))/(2*V.SCALE)for x in h]);hn/=hn.sum();half=np.c_[C[:,9,None]-C[:,:9],-C[:,9]]
    half/=np.linalg.norm(half[:,:9],axis=1)[:,None];vv=HalfspaceIntersection(half,hn[:9],qhull_options='Qx').intersections
    ys=np.c_[vv,1-vv.sum(axis=1)];qs=np.einsum('bi,ij,bj->b',ys,In,ys);accepted=None
    for index in np.argsort(-qs):
        if qs[index]<=0:break
        for eps in(1e-4,1e-6,1e-8):
            yp=ys[index]+eps*hn;y=tuple((Q(round(float(x)*10**12),10**12),Q(0))for x in yp)
            if K.sign(K.pair(y,I,y))<=0 or K.sign(K.pair(y,I,h))>=0:continue
            vals=[V.dot(row,list(map(V.golden,y)))for row in roots+caps]
            if min(p[0]for p in vals)>0:accepted=(y,index);break
        if accepted:break
    assert accepted,'No exact escape found: investigate complete volume.'
    y,index=accepted;data['escape']=list(map(K.enc,y));data['numerical_geometry']={'vertices':len(ys),
        'spacelike_vertices':int((qs>1e-8).sum()),'maximum_vertex_norm':float(qs.max()),
        'escape_proposal_support':[int(i)for i in np.flatnonzero(ys[index]>1e-7)],
        'scope':'Numerical proposals only; the strict exact witness proves infinite volume.'}
    result=B.check(data);path=B.PREFIX+'_caps.json';(V.ROOT/path).write_text(json.dumps(data,indent=2)+'\n')
    result['source_sha256']=V.sha(path);(V.ROOT/(B.PREFIX+'_verify.json')).write_text(json.dumps(result,indent=2)+'\n')
    print({'facets':88,'status':'EXACT_STRICT_ESCAPE','numerical':data['numerical_geometry']},flush=True)

if __name__=='__main__':run()
