#!/usr/bin/env python3
"""Numerical new-cap proposals on the five exactly certified vertex-free edges.

Finite ninth anchors and continuous all-dotted nonedge profiles are separate
bounded scouts.  The exact PD8/profile certificate is an input; every new
normal/angle/interval in this output still requires exactification.
"""
from pathlib import Path
from fractions import Fraction
import json,time,hashlib
import numpy as np
from scipy.linalg import null_space
from fresh_h9_makarov_joint_completion_20261013_profiles import setup,ROOT,P,CORE
from fresh_h9_makarov_joint_completion_20261013_continuous import intervals,intersect
SOURCE=ROOT/'results/fresh_h9_joint_frontier_20261013_all_edges.json'
OUT=ROOT/f'results/{P}_remote_edges.json'
LABELS=[2,3,4,5,6,8,10,12,16,24]
def build():
    start=time.time();src,Q,W,c,E=setup();front=json.loads(SOURCE.read_text())
    target=np.array(json.loads((ROOT/f'results/{P}_all_relaxed_volume.json').read_text())['best_point'])
    target=.999*target+.001*c/(-c@Q@c)
    targets=np.array([target,E,np.array(src['escape']['point'],dtype=float)])
    targets/=-(targets@Q@c)[:,None]
    finite=np.array([-2*np.cos(np.pi/m)for m in[2,3,4,5,6,7,8,9,10,12,15,16,17,20,24,30]])
    caps=[];keys=set();regions=[];stats={'finite_profiles':0,'roots':0,'angle_pass':0,'continuous_regions':0}
    def save(r,recipe):
        ps=r@Q@W.T
        if max(ps)>3e-7 or r@Q@c>=0:return
        ok=ps<-2+3e-7
        for val in finite:ok|=abs(ps-val)<3e-7
        if not ok.all():return
        stats['angle_pass']+=1
        key=tuple(np.rint(r*1e7).astype('int64'))
        if key in keys:return
        keys.add(key);caps.append(dict(id=len(caps),normal=r.tolist(),all27_pairings=ps.tolist(),target_pairings=(r@Q@targets.T).tolist(),**recipe))
    for edge in front['new_edges_without_finite_endpoint']:
        eid=edge['id'];supp=edge['support'];packet=edge['spherical_profiles']
        profiles=[r['contacts']for r in packet['finite_profiles']]
        profiles += [[[r['wall'],m]]for r in packet['dihedral_tail_profiles']for m in LABELS if m>=3]
        stats['finite_profiles']+=len(profiles);outside=sorted(set(range(27))-set(supp))
        G=W[supp]@Q@W[supp].T;Gi=np.linalg.inv(G)
        rhs=[]
        for contacts in profiles:
            p=np.zeros(8)
            for j,m in contacts:p[supp.index(j)]=-2*np.cos(np.pi/m)
            rhs.append(p)
        rhs=np.array(rhs)
        for anchor in outside:
            M=W[supp+[anchor]]@Q
            if np.linalg.matrix_rank(M,tol=1e-9)!=9:continue
            z=null_space(M)[:,0];aa=z@Q@z
            if abs(aa)<1e-9:continue
            T=M.T@np.linalg.inv(M@M.T)
            rr=[];meta=[]
            for pi,row in enumerate(rhs):
                for m in LABELS:rr.append(np.r_[row,-2*np.cos(np.pi/m)if m>2 else 0]);meta.append((pi,m))
            u=np.array(rr)@T.T;bb=2*u@Q@z;cc=np.einsum('ij,jk,ik->i',u,Q,u)-2;D=bb*bb-4*aa*cc
            good=np.flatnonzero(D>1e-9)
            for branch in[-1,1]:
                R=u[good]+((-bb[good]+branch*np.sqrt(D[good]))/(2*aa))[:,None]*z
                stats['roots']+=len(R)
                for r,idx in zip(R,good):
                    pi,m=meta[idx];save(r,{'kind':'finite_anchor','edge':eid,'support':supp,'contacts':profiles[pi],'anchor':anchor,'anchor_label':m,'quadratic_sign':branch})
        # A timelike center in the entire geodesic, no finite endpoint assumed.
        p=c-(W[supp]@Q@c)@Gi@W[supp];p/=np.sqrt(-(p@Q@p))
        plane=null_space(W[supp]@Q);b=plane[:,0]+(plane[:,0]@Q@p)*p
        if b@Q@b<1e-10:b=plane[:,1]+(plane[:,1]@Q@p)*p
        b/=np.sqrt(b@Q@b)
        assert abs(p@Q@b)<1e-6
        for pi,row in enumerate(rhs):
            u=row@Gi@W[supp];R=2-row@Gi@row
            if R<=1e-8:continue
            sq=np.sqrt(R);d=u@Q@W.T;A=p@Q@W.T;B=b@Q@W.T
            for branch in[-1,1]:
                allowed=[(0.,np.inf)]
                for j in outside:
                    allowed=intersect(allowed,intervals(sq*(A[j]+branch*B[j])/2,d[j]+2,sq*(-A[j]+branch*B[j])/2))
                    if not allowed:break
                for lo,hi in allowed:
                    stats['continuous_regions']+=1;rid=len(regions)
                    regions.append({'edge':eid,'support':supp,'contacts':profiles[pi],'branch':branch,'interval':[lo,None if not np.isfinite(hi)else hi],'timelike_unit':p.tolist(),'spacelike_unit':b.tolist()})
                    proposals=[lo+max(1e-5,abs(lo)*1e-5)]
                    if np.isfinite(hi):proposals +=[(lo+hi)/2,hi-max(1e-5,abs(hi)*1e-5)]
                    else:proposals +=[max(lo+1,lo*1.2),max(lo+2,lo*2)]
                    for target in targets:
                        a=p@Q@target;bb=b@Q@target
                        if abs(a+branch*bb)>1e-12 and (a-branch*bb)/(a+branch*bb)>0:proposals.append(np.sqrt((a-branch*bb)/(a+branch*bb)))
                    for x in proposals:
                        if not lo<x<hi or x<1e-6 or x>1e6:continue
                        xx=Fraction(float(x)).limit_denominator(1000000);x=float(xx)
                        if not lo<x<hi:continue
                        r=u+sq*((x-1/x)/2*p+branch*(x+1/x)/2*b)
                        save(r,{'kind':'continuous_all_dotted','region':rid,'edge':eid,'support':supp,'contacts':profiles[pi],'boost_exp_rational':str(xx),'branch':branch})
        print('edge',eid,'caps',len(caps),'stats',stats,flush=True)
    out={'status':'NUMERICAL_VERTEX_FREE_EDGE_CAP_SCOUT','source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in[CORE,SOURCE]},'labels':LABELS,'targets':targets.tolist(),'counts':stats,'regions':regions,'caps':caps,'target_cut_counts':(np.array([r['target_pairings']for r in caps])>1e-8).sum(axis=0).tolist()if caps else[0]*3,'seconds':time.time()-start,'scope':'Exact input supports and finite spherical profiles; sampled tail/ninth-anchor labels, numerical continuous intervals. New normals and all angle/global geometry predicates are numerical proposals only.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print('DONE',len(caps),out['target_cut_counts'],out['seconds'],flush=True)
if __name__=='__main__':build()
