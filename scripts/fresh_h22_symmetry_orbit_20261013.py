#!/usr/bin/env python3
"""Exact integer symmetry orbit; numerical H22 compatibility and volume scout."""
import argparse,json,hashlib
from collections import Counter
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.optimize import linprog
from fresh_h22_borcherds_extension_20261013 import prepare
ROOT=Path(__file__).resolve().parents[1]
SOURCE='data/enumerations/fresh_h22_two_anchor_20261013.json'
SYMS='data/enumerations/fresh_h20_global_interior_sections_20261013.json'
OUT='data/enumerations/fresh_h22_symmetry_orbit_20261013.json'
ARRAY='data/enumerations/fresh_h22_symmetry_orbit_20261013.npz'

def run():
    Q,W,G,h,X,A,b,Qn,R=prepare();src=json.loads((ROOT/SOURCE).read_text());row=src['caps'][1548]
    assert row['center_label']==5 and row['anchors']==[14,16] and row['labels']==[4,2]
    ids=row['orthogonal_support']+row['anchors'];v=sp.Matrix((W[ids]@Q).tolist()).inv()*sp.Matrix([0]*20+[-1,0]);u=np.array([int(2*x)for x in v],dtype=np.int64)
    assert all(sp.Rational(int(x),2)==y for x,y in zip(u,v))and u@Q@u==-6 and u@Q@h==-33
    sd=json.loads((ROOT/SYMS).read_text());generators=[np.array(r['matrix'],dtype=np.int64) for r in sd['ambient_symmetry_generators']]
    for g,r in zip(generators,sd['ambient_symmetry_generators']):
        assert np.array_equal(g.T@Q@g,Q)and np.array_equal(g@h,h)
        assert np.array_equal(W@g.T,W[r['wall_permutation']])
    seen={tuple(u):0};orbit=[u];cursor=0
    while cursor<len(orbit):
        old=orbit[cursor];cursor+=1
        for g in generators:
            v=g@old;key=tuple(v)
            if key not in seen:seen[key]=len(orbit);orbit.append(v)
        if cursor%20000==0:print('orbit',cursor,len(orbit),flush=True)
        assert len(orbit)<=241920
    U=np.array(orbit,dtype=np.int64);assert np.all(np.einsum('ij,jk,ik->i',U,Q,U)==-6)and np.all(U@Q@h==-33)
    np.savez_compressed(ROOT/ARRAY,twice_old_projections=U)
    # For r=(u/2,sqrt(5/2)), r²=1 and r.r'=u.Q.u'/4+5/2.
    first_pairs=U@(Q@U[0])/4+2.5;counts=Counter(map(str,first_pairs));allowed=np.cos(np.pi/np.arange(2,121))
    compatible=(first_pairs<=-1+1e-9)|(np.min(abs(-first_pairs[:,None]-allowed),axis=1)<1e-9)
    compatible[0]=False
    print('complete orbit',len(U),'compatible with representative',int(sum(compatible)),flush=True)
    oldproj=U/2;B=oldproj@Q@X;bb=-oldproj@Q@h/np.sqrt(14);height=np.sqrt(2.5)
    # Both height signs are in the bank; pair checks between entries are
    # deliberately relaxed for this necessary volume screen.
    bankA=np.vstack([np.column_stack([B,np.full(len(U),height)]),np.column_stack([B,np.full(len(U),-height)])]);bankb=np.r_[bb,bb]
    rng=np.random.default_rng(221541);runs=[]
    for m in(3,4,5,6,8,12):
        T=np.sqrt(14+1/(2*np.cos(np.pi/m)**2));AA=np.vstack([np.column_stack([A,np.zeros(210)]),[[0]*21+[T],[0]*21+[-T]],bankA]);bb0=np.r_[b,np.sqrt(14),np.sqrt(14),bankb]
        escape=None;samplemax=-1
        for k in range(100):
            lp=linprog(rng.normal(size=22),A_ub=AA,b_ub=bb0,bounds=[(None,None)]*22,method='highs')
            if not lp.success:continue
            q=float(lp.x@lp.x-1);samplemax=max(samplemax,q)
            if q>1e-7:
                point=np.r_[h/np.sqrt(14)+X@lp.x[:21],lp.x[-1]]
                escape={'norm':q,'point':point.tolist(),'active_original_walls':np.flatnonzero(abs(A@lp.x[:21]-b)<1e-6).tolist(),'active_orbit_entries':np.flatnonzero(abs(bankA@lp.x-bankb)<1e-6).tolist()};break
        runs.append({'center_label':m,'escape':escape,'largest_sample_norm':samplemax});print('m',m,'escape',None if escape is None else escape['norm'],flush=True)
    out={'status':'EXACT_SINGLE_ORBIT_WITH_NUMERICAL_VOLUME_SCREEN','source':SOURCE,'source_cap':1548,'generators':SYMS,'orbit_size':len(U),'array':ARRAY,'array_sha256':hashlib.sha256((ROOT/ARRAY).read_bytes()).hexdigest(),
        'twice_old_representative':u.tolist(),'old_projection_norm':'-3/2','new_height_squared':'5/2','new_norm':'1','center_pairing':'-33/2','same_height_pairing_counts_with_representative':dict(sorted(counts.items())),
        'same_height_compatible_with_representative':int(sum(compatible)),'runs':runs,'scope':'Complete orbit under six exactly checked generators, finite closure by integer vectors. Every orbit cap inherits its representative old-wall relations, not yet promoted to a standalone exact cap certificate here. Whole-orbit bank may contain incompatible pairs. Volume screens numerical.'}
    (ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':run()
