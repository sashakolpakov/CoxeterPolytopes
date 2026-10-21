#!/usr/bin/env python3
"""Bounded three-anchor search in the H2 face below the saved H22 gap."""
import itertools,json
from collections import Counter
import numpy as np
import sympy as sp
from fresh_h22_gap_polar_20261013 import ROOT,read,setup

OUT=ROOT/'data/enumerations/fresh_h22_gap_face_20261013.json'

def valid(c,tol=2e-6):
    allowed=np.cos(np.pi/np.arange(2,121))
    j=np.searchsorted(allowed,c).clip(1,len(allowed)-1)
    err=np.minimum(abs(c-allowed[j]),abs(c-allowed[j-1]))
    return (c>=1-tol)|(err<tol)

def run():
    Q,W,h,U,D,ht=setup();source=read('data/enumerations/fresh_h22_gap_polar_20261013.json')
    x=np.array(source['point']);S=[i for name,i in source['active']if name=='old']
    B=sp.Matrix.hstack(*sp.Matrix((W[S]@Q).tolist()).nullspace())
    basis=np.array(B,dtype=float);C=W@Q@basis;F=basis.T@Q@basis
    nn=np.sqrt(np.diag(W@Q@W.T));coords=x[:22]@Q@basis
    labels=[3,4,6]
    grids=np.array(list(itertools.product(labels,repeat=3)))
    values=np.cos(np.pi/grids)
    anchors=[5,12,26]
    counts=Counter();seen=set();oldcaps=[];survivors=[]
    for serial,ids in enumerate(itertools.combinations(anchors,3)):
        M=C[list(ids)]/nn[list(ids),None]
        if abs(np.linalg.det(M))<1e-9:counts['dependent_triples']+=1;continue
        v=-values@np.linalg.inv(M).T;counts['anchor_profiles']+=len(v)
        q=np.einsum('ij,jk,ik->i',v,F,v);ok=q<1-1e-8
        z=-np.sqrt(np.maximum(0,1-q));ok&=(v@coords+z*x[22]>1e-7)
        rows=np.flatnonzero(ok)
        if not len(rows):continue
        cosine=-(v[rows]@C.T)/nn;good=np.all(valid(cosine),axis=1)
        for j in rows[good]:
            key=tuple(np.round(v[j],8))
            if key in seen:continue
            seen.add(key);p=basis@v[j];height=z[j]
            row={'anchors':ids,'labels':grids[j].tolist(),'old_coordinates':v[j].tolist(),
                 'normal':np.r_[p,height].tolist(),'old_norm':float(q[j]),
                 'target_pairing':float(v[j]@coords+height*x[22])}
            oldcaps.append(row);counts['distinct_old_compatible_cutters']+=1
            center=np.array([h@Q@p+4*height,h@Q@p-4*height])/np.sqrt(2)
            if not np.all(valid(-center)):counts['center_failure']+=1;continue
            ps=U@Q@p/2
            if not np.all(valid(-(ps+np.sqrt(2.5)*abs(height)))) or not np.all(valid(-(ps-np.sqrt(2.5)*abs(height)))):
                counts['timelike_orbit_failure']+=1;continue
            ps=D@Q@p
            if not np.all(valid(-(ps+ht*abs(height))/ht)) or not np.all(valid(-(ps-ht*abs(height))/ht)):
                counts['cusp_orbit_failure']+=1;continue
            survivors.append(row);counts['full_survivors']+=1
            print('survivor',len(survivors),ids,grids[j],row['target_pairing'],flush=True)
        if serial%100==0:print('triples',serial,dict(counts),flush=True)
    return {'status':'BOUNDED_NUMERICAL_H22_GAP_FACE_SEARCH','old_orthogonal_walls':S,
            'basis':[[str(a)for a in B.row(i)]for i in range(22)],
            'anchors':anchors,'labels':labels,'counts':dict(counts),
            'old_compatible_cutters':oldcaps,'survivors':survivors,
            'scope':'Exactly27 profiles on the selected three anchors5,12,26, then all210 oldwalls and both signs of all retained symmetricorbits, centerlabel3, numerical gate labels2..120. A preliminary larger scout was interrupted after locating this triple and supplies no completeness assertion. Neither all-parameter nor all-label completeness.'}

if __name__=='__main__':
    out=run();OUT.write_text(json.dumps(out,indent=2)+'\n');print(out['counts'])
