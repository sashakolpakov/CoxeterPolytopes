#!/usr/bin/env python3
"""Freeze the finite tetrahedral face and run a whole-chamber volume gate."""
import json
from fractions import Fraction as Q
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.spatial import HalfspaceIntersection
from scripts import fresh_h9_compact_next_20261013_verify as T
V=T.V;K=T.K

def run():
    source=T.P.PREFIX+'_caps.json'
    d={'status':'EXACT_FOUR_CAP_COMPACT_TETRAHEDRAL_FACE','source':source,'source_sha256':V.sha(source),
       'orthogonal_old_walls':T.S,'anchor':T.ANCHOR,'anchor_label':3,'frame_sources':[0,3],'boost':'13/4','target':V.read(source)['escape'],
       'definition':{'a':'Orthogonal projection of nextH4 to S-perp; a²=q=(4-2phi)/3.',
         'frame':'Orthonormalize the projected center, projected target, and coordinate vectors0,3 in S-perp intersect a-perp; signature(-1,1,1,1).',
         'P':'cosh(t)*e0+sinh(t)*e1, with exp(t)=13/4',
         'B':'sinh(t)*e0+cosh(t)*e1',
         'compact0':'-a/q+sqrt(1/q-2)*P',
         'compact1':'sqrt(2/3)*B+2/sqrt(3)*e2',
         'compact2':'sqrt(2/3)*B-1/sqrt(3)*e2+e3',
         'compact3':'sqrt(2/3)*B-1/sqrt(3)*e2-e3'},
       'scope':'One explicit four-cap block. The finite3-5-3 face is certified; wholeH9 volume requires a separate gate.'}
    out=T.check(d);print({'exact_caps':4,'actual_finite_vertices':len(out['actual_finite_vertices']),'new_affine':out['affine_census']['new_affine_circuits']},flush=True)
    G,I,h,g,old,names,caps,center,_=T.reconstruct(d);ff=lambda x:sum(x)/(2*V.SCALE)
    C=np.array([[ff(x)for x in r]for r in old+caps]);In=np.array([[ff(V.golden(x))for x in row]for row in I]);hn=np.array([ff(V.golden(x))for x in h]);hn/=hn.sum()
    H=np.c_[C[:,9,None]-C[:,:9],-C[:,9]];H/=np.linalg.norm(H[:,:9],axis=1)[:,None]
    vv=HalfspaceIntersection(H,hn[:9],qhull_options='Qx').intersections;ys=np.c_[vv,1-vv.sum(axis=1)];qs=np.einsum('bi,ij,bj->b',ys,In,ys)
    accepted=None
    for j in np.argsort(-qs):
        if qs[j]<=0:break
        for eps in(1e-4,1e-6,1e-8):
            y=tuple((Q(round(float(x)*10**12),10**12),Q(0))for x in ys[j]+eps*hn)
            if K.sign(K.pair(y,I,y))<=0 or K.sign(K.pair(y,I,h))>=0:continue
            pp=[V.dot(r,list(map(V.golden,y)))for r in old+caps]
            if min(p[0]for p in pp)>0:accepted=(y,j);break
        if accepted:break
    assert accepted,'No exact escape found; perform a complete volume certificate.'
    y,j=accepted;d['escape']=list(map(K.enc,y));d['numerical_geometry']={'vertices':len(ys),'spacelike_vertices':int((qs>1e-8).sum()),'maximum_norm':float(qs.max()),
        'escape_proposal_support':[int(i)for i in np.flatnonzero(abs(ys[j])>1e-7)],'escape_active_walls':[names[i]if i<len(names)else T.NEW[i-len(names)]for i in np.flatnonzero(abs(C@ys[j])<1e-7)],
        'scope':'Floating vertex census proposes the exactly checked escape. These counts do not certify volume.'}
    out=T.check(d);path=T.PREFIX+'_caps.json';(V.ROOT/path).write_text(json.dumps(d,indent=2)+'\n');out['source_sha256']=V.sha(path)
    (V.ROOT/(T.PREFIX+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n');print({'status':'EXACT_ESCAPE','geometry':d['numerical_geometry']},flush=True)

if __name__=='__main__':run()
