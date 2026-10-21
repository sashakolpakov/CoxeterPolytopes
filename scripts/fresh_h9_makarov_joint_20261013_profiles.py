#!/usr/bin/env python3
"""Numerical endpoint-profile cap proposals; every recipe is exact-definable.

The spherical support census is exact input. Solving the ninth finite contact,
norm, all-wall angle gate and volume scouting here are numerical only.
"""
from pathlib import Path
from itertools import product
import argparse,hashlib,json,time
import numpy as np
from scipy.linalg import null_space
from fresh_h9_makarov_joint_20261013 import read,volume_scout
import fresh_h9_makarov_joint_20261013 as J

ROOT=Path(__file__).resolve().parents[1]
FRONT=ROOT/'results/fresh_h9_makarov_frontier_20261013.json'
OUT=ROOT/'results/fresh_h9_makarov_joint_20261013_profiles.json'

def build():
    start=time.time();data,Q,W,c=read();F=json.loads(FRONT.read_text())
    ns=np.einsum('ij,jk,ik->i',W,Q,W);labels=[2,3,4,5,6,8,10,12,16,24]
    gate=[2,3,4,5,6,7,8,9,10,12,15,16,17,20,24,30]
    finite=np.array([-2*np.cos(np.pi/m) for m in gate])
    caps=[];keys=set();stats={'edge_profiles':0,'quadratic_roots':0,'acute':0,'angle_pass':0}
    for edgeid,packet in F['H4_released_edge_profiles'].items():
        e=F['outgoing_edges'][int(edgeid)];support=e['support']
        profiles=[row['contacts'] for row in packet['finite_profiles']]
        profiles += [[[row['wall'],m]] for row in packet['dihedral_tail_profiles'] for m in labels if m>=3]
        stats['edge_profiles']+=len(profiles)
        for anchor in sorted(set(range(11))-set(support)):
            ids=support+[anchor];M=W[ids]@Q
            if np.linalg.matrix_rank(M,tol=1e-9)<9:continue
            pin=M.T@np.linalg.inv(M@M.T);z=null_space(M)[:,0];aa=z@Q@z
            if abs(aa)<1e-9:continue
            for pi,contacts in enumerate(profiles):
                rhs=np.zeros((len(labels),9));lookup=dict(contacts)
                for k,i in enumerate(support):
                    if i in lookup:rhs[:,k]=-2*np.cos(np.pi/lookup[i])*np.sqrt(ns[i]/2)
                rhs[:,8]=[-2*np.cos(np.pi/m)*np.sqrt(ns[anchor]/2) if m>2 else 0 for m in labels]
                r0=rhs@pin.T;bb=2*r0@Q@z;cc=np.einsum('ij,jk,ik->i',r0,Q,r0)-2
                dd=bb*bb-4*aa*cc
                for mi in np.flatnonzero(dd>1e-9):
                    for sign in [-1,1]:
                        r=r0[mi]+(-bb[mi]+sign*np.sqrt(dd[mi]))/(2*aa)*z
                        stats['quadratic_roots']+=1;ps=W@Q@r
                        if max(ps)>2e-7 or r@Q@c>=0:continue
                        stats['acute']+=1;p=ps/np.sqrt(ns/2)
                        if not np.all((p<-2+2e-7)|(np.min(abs(p[:,None]-finite),axis=1)<2e-7)):continue
                        stats['angle_pass']+=1
                        key=tuple(np.rint(r*1e7).astype('int64'))
                        if key in keys:continue
                        keys.add(key)
                        caps.append({'id':len(caps),'edge':int(edgeid),'support':support,'contacts':contacts,
                                     'anchor':anchor,'anchor_label':labels[mi],'quadratic_sign':sign,
                                     'normal':r.tolist(),'seed_pairings':ps.tolist()})
    out={'status':'NUMERICAL_ENDPOINT_PROFILE_PROPOSALS','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [FRONT,J.SOURCE]},
         'finite_anchor_labels':labels,'tail_labels_sampled':labels[1:],
         'counts':stats,'entries':len(caps),'caps':caps,'elapsed_seconds':time.time()-start,
         'scope':'Finite profile census input is complete, but anchor and isolated-tail labels here are bounded; numeric acceptance is not an exact certificate.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print('PROFILES',stats,len(caps),time.time()-start,flush=True)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--objectives',type=int,default=100);p.add_argument('--volume-only',action='store_true');a=p.parse_args()
    bank=json.loads(OUT.read_text())if a.volume_only else build()
    # Combine the earlier free-center bank, without forcing compatibility.
    initial=json.loads(J.OUT.read_text());combined={'caps':initial['caps']+bank['caps']}
    J.OUT=OUT
    volume_scout(combined,a.objectives)
