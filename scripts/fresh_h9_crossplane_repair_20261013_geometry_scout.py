#!/usr/bin/env python3
"""Compare coherent pair/triple repairs against all original numeric vertices."""
from itertools import combinations
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.spatial import HalfspaceIntersection
from scripts import fresh_h9_crossplane_repair_20261013_verify as T
V=T.V

def run():
    data=V.read(T.PREFIX+'_caps.json');G,I,h,g,old,names,r,center,_=T.reconstruct(data);f=lambda p:sum(p)/(2*V.SCALE)
    C=np.array([[f(x)for x in a]for a in old]);GG=np.array([[f(V.golden(x))for x in a]for a in G]);In=np.linalg.inv(GG)
    hn=np.array([f(V.golden(x))for x in h]);hn/=hn.sum();half=np.c_[C[:,9,None]-C[:,:9],-C[:,9]];half/=np.linalg.norm(half[:,:9],axis=1)[:,None]
    vv=HalfspaceIntersection(half,hn[:9],qhull_options='Qx').intersections;ys=np.c_[vv,1-vv.sum(axis=1)];q=np.einsum('bi,ij,bj->b',ys,In,ys)
    proposals=V.read(T.PREFIX+'_pair2_scout.json')['candidates'];cases=[];unique=[]
    for i,record in enumerate(proposals):
        pair=np.array(record['normals'])
        if any(np.max(abs(pair-oldpair))<1e-6 for oldpair in unique):continue
        unique.append(pair);sets=[pair]
        # A third seed permutation is tested only as a separate coherent set.
        D=C[5]-C[7];third=pair[0]-((pair[0]@GG@D)/2)*D
        vals=np.array([2*np.cos(np.pi/m)for m in(2,3,4,5,6,8,10,12)]);ps=C@GG@third;er=np.abs(ps[:,None]+vals).min(axis=1);er=np.where(ps<-2+1e-7,0,er)
        if er.max()<1e-7 and np.max(abs(pair@GG@third))<1e-7:sets.append(np.vstack([pair,third]))
        for caps in sets:
            survive=np.flatnonzero(np.all(caps@ys.T>=-1e-7,axis=0));j=survive[np.argmax(q[survive])]
            A=np.vstack([C,caps]);hh=np.c_[A[:,9,None]-A[:,:9],-A[:,9]];hh/=np.linalg.norm(hh[:,:9],axis=1)[:,None]
            vv2=HalfspaceIntersection(hh,hn[:9],qhull_options='Qx').intersections;zz=np.c_[vv2,1-vv2.sum(axis=1)];qq=np.einsum('bi,ij,bj->b',zz,In,zz);k=int(np.argmax(qq))
            row={'proposal':i,'caps':len(caps),'normals':caps.tolist(),'surviving_old_spacelike_vertices':int(sum(q[survive]>1e-8)),
                 'largest_surviving_old_norm':float(q[j]),'old_escape':ys[j].tolist(),'vertices':len(zz),'spacelike_vertices':int(sum(qq>1e-8)),
                 'maximum_norm':float(qq[k]),'new_escape':zz[k].tolist(),'new_escape_active':[int(x)for x in np.flatnonzero(abs(A@zz[k])<1e-7)],
                 'new_finite_vertices_numeric':int(sum(qq< -1e-8))-22}
            cases.append(row);print({k:v for k,v in row.items()if k not in('normals','old_escape','new_escape','new_escape_active')},flush=True)
    out={'status':'NUMERICAL_COHERENT_PAIR_TRIPLE_GEOMETRY_COMPARISON','source':T.PREFIX+'_pair2_scout.json','source_sha256':V.sha(T.PREFIX+'_pair2_scout.json'),
         'old88_vertices':len(ys),'cases':cases,'scope':'Numeric complete projectivevertex scouts and directoldvertexsurvival. Everycasehasaspacelikeproposal; exactfinitevolumeorcompactness conclusions require exactverification.'}
    (V.ROOT/(T.PREFIX+'_geometry_scout.json')).write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':run()
