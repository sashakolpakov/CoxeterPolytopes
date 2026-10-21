#!/usr/bin/env python3
"""Numerical placement of thirteen local cusp-completion normals.

All nonincident constraints are required dotted. This supplies rational-height
proposals, not an exact construction or an optimality claim.
"""
import json,numpy as np
from fractions import Fraction as F
from fresh_h22_new_cusp_completion_20261013_verify import setup,P,read
from fresh_h22_new_cusp_completion_20261013_tight import shifted,DEN
PRE='fresh_h22_transported_cusp_completion_20261013'
def run():
    _,Q,W,h,meta,gens,D,V,B,A,ids,sizes=setup()
    proposal=read('data/enumerations/'+PRE+'_local_proposal.json')
    C=np.array([[float(F(x))for x in r]for r in proposal['bases']])
    d=np.array([float(F(x))for x in proposal['target']['delta']]);q=np.array([float(F(x))for x in proposal['base_norms']])
    cov=Q@C[:,:22].T;dc=Q@d[:22];threshold=np.zeros(13);anchors=[None]*13;counts={}
    def block(name,rows,z,n,signed=True):
        size=len(rows);counts[name]=size*(2 if signed else 1)
        for start in range(0,size,32768):
            rr=rows[start:start+32768];zz=np.broadcast_to(z,size)[start:start+len(rr)];nn=np.broadcast_to(n,size)[start:start+len(rr)]
            dd=rr@dc;bb=rr@cov
            for sign in ((1,-1)if signed else(1,)):
                slope=dd+sign*zz*d[22]*np.sqrt(11);pair=bb+sign*zz[:,None]*C[:,22]*np.sqrt(11)
                assert max(slope)<1e-6,(name,start,max(slope))
                incident=abs(slope)<1e-6
                if np.any(incident):
                    cs=-pair[incident]/np.sqrt(nn[incident,None]*q[None,:]);allowed=np.cos(np.pi/np.arange(2,121))
                    assert np.all((cs>=1-1e-6)|(np.min(abs(cs[:,:,None]-allowed),axis=2)<1e-6)),name
                values=np.full_like(pair,-np.inf)
                values[~incident]=(pair[~incident]+np.sqrt(nn[~incident,None]*q[None,:]))/-slope[~incident,None]
                where=np.argmax(values,axis=0);mx=values[where,np.arange(13)]
                for j in np.flatnonzero(mx>threshold):
                    threshold[j]=mx[j];anchors[j]={'family':name,'index':start+int(where[j]),'sign':sign,'threshold':float(mx[j])}
    wn=np.einsum('ij,ij->i',W,W@Q);block('old210',W.astype(float),0.,wn,False)
    core=read('results/fresh_h22_incidence_core_20261013.json');binary=read('results/fresh_h22_incidence_binary_20261013.json')
    block('core4',np.array(core['old_projections'],dtype=float),2.,2.)
    block('binary560',np.array(binary['accepted_old_projections'],dtype=float),2.,2.)
    z=np.load(P/'data/enumerations/fresh_h22_symmetry_orbit_20261013.npz');block('U',z['twice_old_projections']/2,np.sqrt(2.5),1.)
    z=np.load(P/'data/enumerations/fresh_h22_cusp_orbits_20261013.npz');heights=z['four_heights']/4;block('K',z['old_null_projections'].astype(float),heights,heights**2)
    z=np.load(P/'data/enumerations/fresh_h22_gap_face_orbit_20261013.npz');block('R',z['six_old_projections']/6,np.sqrt(11),2/3)
    an=np.einsum('ij,ij->i',A,A@Q)
    for name,rows,t,ix,index in shifted(D,V,B,A,ids,sizes,Q,True):
        nq=an if index==0 else [279/14,279/14,55/62496,55/62496][index-1]
        block(name,rows/DEN,t*np.sqrt(11)/10,nq)
    old=read('results/fresh_h22_incidence_binary_20261013_escape.json');x=np.array(old['numerator'])/old['denominator']
    ds=float(d[:22]@Q@x[:22]+d[22]*np.sqrt(11)*x[22]);basepair=C[:,:22]@Q@x[:22]+C[:,22]*np.sqrt(11)*x[22];cut=basepair/-ds
    shifts=np.ceil(threshold*100+1)/100
    out={'status':'NUMERICAL_THIRTEEN_LOCAL_CUSP_WALL_PLACEMENTS','counts':counts,'all_signed_wall_count':sum(counts.values()),'thresholds':threshold.tolist(),'anchors':anchors,'proposed_shifts':[str(F(float(t)).limit_denominator(100))for t in shifts],'escape_cut_thresholds':cut.tolist(),'proposed_escape_pairings':(basepair+shifts*ds).tolist(),'scope':'All nonincident relations are required dotted; floating-point bounds only. Exact field/radical pair gate and local bounded-link proof are required.'}
    (P/('data/enumerations/'+PRE+'_placements.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],flush=True);print('thresholds',threshold.tolist(),flush=True);print('cut shifts',cut.tolist(),flush=True);print('cutters',np.flatnonzero(basepair+shifts*ds>0).tolist(),flush=True)
if __name__=='__main__':run()
