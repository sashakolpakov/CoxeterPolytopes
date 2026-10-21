#!/usr/bin/env python3
"""Metric-coupled H4 endpoint, then a certified volume gate."""
import json
from fractions import Fraction as Q
import numpy as np
from scipy.spatial import HalfspaceIntersection
from scripts import fresh_h9_remote_next_20261013_verify as N
V=N.V;K=V.K

def run():
    source=N.W.PREFIX+'_caps.json';old=V.read(source)
    data={'status':'EXACT_REMOTE_H4_ENDPOINT_CONSTRUCTION','source':source,'source_sha256':V.sha(source),
          'target_point':old['escape'],'orthogonal_walls':N.S_NAMES+['wall17'],'pi3_wall':'coupled0',
          'pi10_wall':'coupled1','negative_unit_seed':'seed0','norm':2,
          'definition':{'orthogonal_frame':'Previous coupled stage has P^2=-1 and R^2=C^2=1; R=(coupled0+coupled1)/2 and C=(coupled0-coupled1)/2.',
                        'angle_constant':'L=sqrt(2+phi)>0','b':'-(1+L)/2','c':'(L-1)/2',
                        'a':'-[b(R,wall17)+c(C,wall17)]/(P,wall17)',
                        'u':'-sqrt((2+a^2-b^2-c^2)/2)','normal':'aP+bR+cC+u*seed0'},
          'scope':'Exact cap with two coupled finite angles closing a specific old edge and cutting the current target; not an exhaustive cap search.'}
    certificate=N.check(data);print({'new_cap':True,'actual_vertex':'H4+5A1','all85_retained':True},flush=True)
    G,I,h,g,roots,names,cap,_,parameters,coords=N.reconstruct(data)
    C=np.array([[(x[0]+x[1])/(2*V.SCALE)for x in row]for row in roots+[cap]])
    In=np.array([[sum(V.golden(x))/(2*V.SCALE)for x in row]for row in I]);hn=np.array([sum(V.golden(x))/(2*V.SCALE)for x in h]);hn/=hn.sum()
    half=np.c_[C[:,9,None]-C[:,:9],-C[:,9]];half/=np.linalg.norm(half[:,:9],axis=1)[:,None]
    vertices=HalfspaceIntersection(half,hn[:9],qhull_options='Qx').intersections;ys=np.c_[vertices,1-vertices.sum(axis=1)]
    qs=np.einsum('bi,ij,bj->b',ys,In,ys);accepted=None
    for index in np.argsort(-qs):
        if qs[index]<=0:break
        for eps in(1e-4,1e-6,1e-8):
            yp=ys[index]+eps*hn;y=tuple((Q(round(float(x)*10**12),10**12),Q(0))for x in yp)
            if K.sign(K.pair(y,I,y))<=0 or K.sign(K.pair(y,I,h))>=0:continue
            vals=[V.dot(row,list(map(V.golden,y)))for row in roots+[cap]]
            if min(x[0]for x in vals)>0:accepted=(y,index);break
        if accepted:break
    assert accepted,'No exact escape: investigate volume rather than declaring success.'
    y,index=accepted;data['escape']=list(map(K.enc,y));data['numerical_geometry']={
        'vertices':len(ys),'spacelike_vertices':int((qs>1e-8).sum()),'maximum_vertex_norm':float(qs.max()),
        'escape_proposal_support':[int(i)for i in np.flatnonzero(ys[index]>1e-7)],
        'scope':'Numerical proposals only; exact strict escape is the volume certificate.'}
    certificate=N.check(data);path=N.PREFIX+'_caps.json';(V.ROOT/path).write_text(json.dumps(data,indent=2)+'\n')
    certificate['artifact_sha256']=V.sha(path);(V.ROOT/(N.PREFIX+'_verify.json')).write_text(json.dumps(certificate,indent=2)+'\n')
    print({'facets':86,'status':'EXACT_STRICT_ESCAPE','numerical':data['numerical_geometry']},flush=True)

if __name__=='__main__':run()
