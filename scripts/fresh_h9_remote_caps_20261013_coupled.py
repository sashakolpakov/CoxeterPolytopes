#!/usr/bin/env python3
"""Tune a remote two-cap insertion to an actual H3 endpoint on an old edge."""
import json
import numpy as np
from scipy.spatial import HalfspaceIntersection
from scripts.fresh_h9_remote_caps_20261013_base import *
from scripts import fresh_h9_remote_caps_20261013_verify as V
from scripts import fresh_h9_remote_caps_20261013_coupled_verify as W

def run():
    G,I,h,old,geometry=selected();previous=read(PREFIX+'_caps.json');target=tuple(map(parse,previous['escape']))
    center=tuple(K.neg(x)for x in matvec(K,I,h));w=tuple(K.neg(x)for x in matvec(K,I,target))
    def project(v):
        for n in W.S_NAMES:
            p=div(dot(K,v,old[n]['p']),scale(O,2));v=tuple(sub(x,mul(p,r))for x,r in zip(v,old[n]['c'],strict=True))
        return v
    t=project(center);wp=project(w);a=K.neg(pairing(K,t,G,t));q=pairing(K,wp,G,wp);pt=pairing(K,wp,G,t)
    assert K.sign(a)>0 and K.sign(q)>0
    v=tuple(add(x,mul(div(pt,a),y))for x,y in zip(wp,t,strict=True));b=pairing(K,v,G,v)
    e=project(old['seed8']['c'])
    for u,nu in((t,K.neg(a)),(v,b)):
        p=div(pairing(K,e,G,u),nu);e=tuple(sub(x,mul(p,y))for x,y in zip(e,u,strict=True))
    c=pairing(K,e,G,e);assert K.sign(c)>0
    norms=[K.neg(a),b,c];frame=[t,v,e];ds=[div(O,z)for z in(a,b,c)]
    data={'status':'EXACT_FINITE_ANGLE_REMOTE_PAIR','source':PREFIX+'_caps.json','source_sha256':sha(PREFIX+'_caps.json'),
          'target_point':list(map(enc,target)),'frame':[[enc(x)for x in row]for row in frame],
          'frame_norms':list(map(enc,norms)),'orthogonal_old_walls':W.S_NAMES,'radicands':list(map(enc,ds)),
          'projected_target_norm':enc(q),'gram_schmidt_source':'seed8','finite_anchor':'wall17',
          'finite_cap_sign':1,'finite_label':3,'boost_branch':-1,
          'scope':'One coupled boost chosen by an exact angle equation, retaining all83 previous walls; not a bounded-bank exhaustive search.'}
    cert=W.check(data);print({'exact_caps':2,'actual_vertices':['H3+6A1','9A1'],'finite_edge':True},flush=True)
    _,_,_,_,roots,names,caps,point,poly=W.reconstruct(data)
    C=np.array([[(x[0]+x[1])/(2*V.SCALE)for x in row]for row in roots+caps]);In=np.array([[float(value(K,x))for x in row]for row in I])
    hn=np.array([float(value(K,x))for x in h]);hn/=hn.sum()
    half=np.c_[C[:,9,None]-C[:,:9],-C[:,9]];half/=np.linalg.norm(half[:,:9],axis=1)[:,None]
    vertices=HalfspaceIntersection(half,hn[:9],qhull_options='Qx').intersections;ys=np.c_[vertices,1-vertices.sum(axis=1)]
    qs=np.einsum('bi,ij,bj->b',ys,In,ys);accepted=None
    for index in np.argsort(-qs):
        if qs[index]<=0:break
        for epsilon in(1e-4,1e-6,1e-8):
            yn=ys[index]+epsilon*hn;y=tuple(scale(O,Q(round(float(x)*10**12),10**12))for x in yn)
            if K.sign(pairing(K,y,I,y))<=0 or K.sign(pairing(K,y,I,h))>=0:continue
            yi=list(map(V.golden,y));vals=[V.dot(r,yi)for r in roots+caps]
            if min(x[0]for x in vals)>0:accepted=(y,index);break
        if accepted:break
    assert accepted,'No exact escape: investigate volume independently.'
    y,index=accepted;data['escape']=list(map(enc,y));data['numerical_geometry']={'vertices':len(ys),
        'spacelike_vertices':int((qs>1e-8).sum()),'maximum_vertex_norm':float(qs.max()),
        'escape_proposal_support':[int(i)for i in np.flatnonzero(ys[index]>1e-7)],
        'scope':'Exploratory numerical vertex census; only the strict exact escape certifies infinite volume.'}
    cert=W.check(data);path=W.PREFIX+'_caps.json';(ROOT/path).write_text(json.dumps(data,indent=2)+'\n')
    cert['source_sha256']=sha(path);(ROOT/(W.PREFIX+'_verify.json')).write_text(json.dumps(cert,indent=2)+'\n')
    print({'facets':85,'geometry':'EXACT_STRICT_ESCAPE','numerical':data['numerical_geometry']},flush=True)

if __name__=='__main__':run()
