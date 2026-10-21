#!/usr/bin/env python3
"""Discovery: unequal two-block H9 seeds and symmetry-reduced cap domains.

The full prescribed finite-pairing bank is used first, even when its caps
are mutually incompatible. An escape for that intersection rejects every
selection from that bank. Floats only; no exhaustive algebraic conclusion.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import linprog
from scipy.spatial import HalfspaceIntersection,QhullError

import fresh_h9_multiblock_seed_20261013 as base

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/enumerations/fresh_h9_asymmetric_blocks_20261013_scout.json'


def reduced_geometry(ns,G,R):
    vals,Q=np.linalg.eigh(G)
    idx=list(np.flatnonzero(vals < -1e-7))+list(np.flatnonzero(vals > 1e-7))
    assert len(idx)==10 and sum(vals < -1e-7)==1
    N=Q[:,idx]*np.sqrt(abs(vals[idx]))
    N[:,0]*=np.sign(sum(N[:,0]))
    W=np.vstack([N,R@N])
    assert min(W[:,0])>1e-8
    # Simultaneous permutations of U_i,V_i within each group are seed
    # isometries and preserve this complete orbit bank. Every point has
    # a representative in this ordered sector, with the same Lorentz norm.
    orders=[];start=0
    for n in ns:
        orders.extend(N[i]-N[i+1] for i in range(start,start+n-1))
        start+=n
    A=np.vstack([W,orders]);hs=np.c_[A[:,1:],-A[:,0]]
    lengths=np.linalg.norm(hs[:,:9],axis=1)
    lp=linprog(np.r_[np.zeros(9),-1.],A_ub=np.c_[hs[:,:9],lengths],
               b_ub=-hs[:,9],bounds=[(None,None)]*9+[(0,None)],method='highs')
    assert lp.success and lp.x[-1]>1e-7
    option='default'
    try:V=HalfspaceIntersection(hs,lp.x[:9]).intersections
    except QhullError:
        option='Qs'
        try:V=HalfspaceIntersection(hs,lp.x[:9],qhull_options='Qs').intersections
        except QhullError:
            option='QJ'
            V=HalfspaceIntersection(hs,lp.x[:9],qhull_options='QJ').intersections
    q=np.sum(V*V,axis=1);j=int(np.argmax(q));v=V[j]
    residual=W[:,1:]@v-W[:,0]
    assert max(residual)<1e-6
    y=-(N*np.r_[-1.,np.ones(9)])@np.r_[1.,v]
    return {'qhull_option':option,'ordered_sector_vertices':len(V),'outside':int(sum(q>1+1e-7)),
            'ideal':int(sum(abs(q-1)<1e-7)),'max_klein_square':float(q[j]),
            'worst_point_seed_covector':y.tolist(),'worst_klein_point':v.tolist(),
            'max_wall_residual':float(max(residual)),
            'active_walls':np.flatnonzero(abs(residual)<1e-7).tolist(),
            'N':N.tolist()}


def run(output):
    started=time.monotonic();records=[];census=[]
    # Fixed rational transverse parameters. All admissible roots of each
    # prescribed norm equation are retained; no top-cap-count truncation.
    bs=[1/16,1/8,1/4,3/8,1/2,5/8,3/4]
    for ns in ([5,3],[6,2],[7,1]):
        P=base.profiles(ns)
        for b in bs:
            if b*max(ns)>=4:continue
            points=base.parameters(ns,b,P)
            census.append({'groups':ns,'b':b,'parameter_points':len(points)})
            for point in points:
                a=point['a'];G,R,refs=base.bank(ns,a,b,P,point['profiles'])
                geo=reduced_geometry(ns,G,R)
                rec={'groups':ns,'a':a,'b':b,'cap_count':len(R),
                     'profile_indices':point['profiles'],
                     'profile_types':[P[1][i].tolist()for i in point['profiles']],
                     'geometry_all_caps':geo}
                if geo['outside']==0:
                    H=R@G@R.T;C=base.compat(H)
                    rec['compatible_pair_count']=int(C.sum()//2)
                    rec['cap_covectors']=refs
                records.append(rec)
            result={'status':'NUMERICAL_ASYMMETRIC_BLOCK_DOMAIN_SCOUT',
                    'source_sha256':hashlib.sha256(Path(base.__file__).read_bytes()).hexdigest(),
                    'b_values':bs,'census':census,'records':records,
                    'seconds':round(time.monotonic()-started,3),
                    'scope':'Unequal two-block sizes [5,3],[6,2],[7,1], listed rational b, all roots a in (b+2,20) returned by every prescribed finite pairing norm equation, labels2,3,5 on all16 seed walls. Every orbit member included. Numerical geometry of the whole cap-bank intersection in a permutation sector; cap pairs need not be compatible. No exact exclusion or finite-volume claim.'}
            output.write_text(json.dumps(result,indent=2)+'\n')
            print({'groups':ns,'b':b,'points':len(points),'total':len(records),
                   'best_max_klein_square':min(r['geometry_all_caps']['max_klein_square']for r in records),
                   'seconds':round(time.monotonic()-started,2)},flush=True)
    print({'done':len(records),'no_outside':[i for i,r in enumerate(records)if r['geometry_all_caps']['outside']==0]},flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=OUT)
    args=parser.parse_args();run(args.output)
