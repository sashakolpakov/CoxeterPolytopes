#!/usr/bin/env python3
"""New finite-endpoint profiles of the exact27-wall chamber: numerical scout."""
from pathlib import Path
from fractions import Fraction
import argparse,hashlib,json,time
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[1];P='fresh_h9_makarov_joint_completion_20261013'
CORE=ROOT/'results/fresh_h9_makarov_joint_20261013_exact.json'
FRONT=ROOT/'results/fresh_h9_joint_frontier_20261013.json'
OUT=ROOT/f'results/{P}_profiles.json'
def setup():
    src=json.loads(CORE.read_text());f=lambda x:float(Fraction(x[0]))+np.sqrt(5)*float(Fraction(x[1]))
    Q=np.array([[f(x)for x in row]for row in src['form']]);c=np.array([f(x)for x in src['strict_center']])
    W=[[f(x)for x in row]for row in src['seed_normals']]
    W += [[f(x[0])+np.sqrt(f(row['radicand']))*f(x[1])for x in row['normal']]for row in src['caps']]
    W=np.array(W);W*=np.sqrt(2/np.einsum('ij,jk,ik->i',W,Q,W))[:,None]
    oldbank=json.loads((ROOT/f'results/{P}_old_bank_gate.json').read_text())
    E=np.array(oldbank['escape']['point'],dtype=float);E/=-(E@Q@c)
    return src,Q,W,c,E
def build(all_targets=False):
    start=time.time();src,Q,W,c,E=setup();front=json.loads(FRONT.read_text())
    labels=[2,3,4,5,6,8,10,12,16,24];finite=np.array([-2*np.cos(np.pi/m)for m in[2,3,4,5,6,7,8,9,10,12,15,16,17,20,24,30]])
    caps=[];keys=set();counts={'systems':0,'profiles':0,'roots':0,'target_cut':0,'acute':0,'angle_pass':0}
    for eid,packet in front['open_edge_spherical_profiles'].items():
        edge=front['edges'][int(eid)];supp=edge['support']
        profiles=[r['contacts']for r in packet['finite_profiles']]
        profiles += [[[r['wall'],m]]for r in packet['dihedral_tail_profiles']for m in labels if m>=3]
        counts['profiles']+=len(profiles)
        for anchor in sorted(set(range(27))-set(supp)):
            ids=supp+[anchor];M=W[ids]@Q
            if np.linalg.matrix_rank(M,tol=1e-9)!=9:continue
            z=null_space(M)[:,0];aa=z@Q@z
            if abs(aa)<1e-9:continue
            T=M.T@np.linalg.inv(M@M.T);counts['systems']+=1
            rhs=[];meta=[]
            for pi,contacts in enumerate(profiles):
                p=np.zeros(9);lookup=dict(contacts)
                for j,w in enumerate(supp):
                    if w in lookup:p[j]=-2*np.cos(np.pi/lookup[w])
                for m in labels:
                    r=p.copy();r[-1]=-2*np.cos(np.pi/m)if m>2 else 0
                    rhs.append(r);meta.append((pi,m))
            u=np.array(rhs)@T.T;bb=2*u@Q@z;cc=np.einsum('ij,jk,ik->i',u,Q,u)-2;D=bb*bb-4*aa*cc
            good=np.flatnonzero(D>1e-9)
            for branch in[-1,1]:
                R=u[good]+((-bb[good]+branch*np.sqrt(D[good]))/(2*aa))[:,None]*z
                counts['roots']+=len(R);cut=R@Q@E>1e-8
                counts['target_cut']+=int(cut.sum())
                if all_targets:cut=np.ones(len(R),dtype=bool)
                R=R[cut];indices=good[cut]
                if not len(R):continue
                ps=R@Q@W.T;ok=(ps.max(axis=1)<3e-7)&(R@Q@c<0)
                counts['acute']+=int(ok.sum());R=R[ok];ps=ps[ok];indices=indices[ok]
                if not len(R):continue
                ok=(ps<-2+3e-7)
                for val in finite:ok |=abs(ps-val)<3e-7
                mask=ok.all(axis=1);counts['angle_pass']+=int(mask.sum())
                for r,p,idx in zip(R[mask],ps[mask],indices[mask]):
                    key=tuple(np.rint(r*1e6).astype('int64'))
                    if key in keys:continue
                    keys.add(key);pi,m=meta[idx]
                    caps.append({'id':len(caps),'edge':int(eid),'support':supp,'contacts':profiles[pi],
                                 'anchor':anchor,'anchor_label':m,'quadratic_sign':branch,'normal':r.tolist(),
                                 'all27_pairings':p.tolist(),'target_pairing':float(r@Q@E)})
        if int(eid)%20==0:print('edge',eid,'new',len(caps),'counts',counts,'secs',round(time.time()-start,1),flush=True)
    out={'status':'BOUNDED_NUMERICAL_NEW_CORNER_PROFILE_CAPS','source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in[CORE,FRONT]},
         'labels':labels,'counts':counts,'caps':caps,'seconds':time.time()-start,'requires_primary_target_cut':not all_targets,
         'scope':'All listed finite spherical profiles are tried, but tail and ninth-anchor labels are sampled. Every root numerically passes all27 pair gates; the primary-target filter is explicitly recorded. No exact geometry or exhaustive cap-family claim.'}
    outpath=OUT.with_name(OUT.stem+'_all.json')if all_targets else OUT
    outpath.write_text(json.dumps(out,indent=2)+'\n');print('DONE',len(caps),counts,out['seconds'],flush=True);return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--all-targets',action='store_true');a=p.parse_args();build(a.all_targets)
