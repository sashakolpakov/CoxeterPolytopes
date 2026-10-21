#!/usr/bin/env python3
"""Bounded numerical new-endpoint search at the new PD7+2A1 corner."""
import json,time,hashlib
from fractions import Fraction as F
import numpy as np
from scipy.linalg import null_space
from fresh_h9_makarov_joint_completion_20261013_profiles import setup,ROOT,P,CORE
from fresh_h9_joint_frontier_20261013 import spherical_profiles,K,Z,PHI,I
SOURCE=ROOT/f'results/{P}_coupled_exact.json'
OUT=ROOT/f'results/{P}_coupled_frontier.json'
def build():
    start=time.time();src,Q,W,c,E=setup();source=json.loads(SOURCE.read_text());f=lambda x:float(F(x[0]))+np.sqrt(5)*float(F(x[1]))
    W=np.vstack([W,[[f(x)*np.sqrt(f(row['radicand']))for x in row['projective_normal']]for row in source['new_caps']]])
    p=np.array([f(x)for x in source['finite_vertex']['point']]);incident=source['finite_vertex']['incident']
    targets=np.array([source['escape']['point'],json.loads((ROOT/f'results/{P}_old_bank_gate.json').read_text())['escape']['point']],dtype=float);targets/=-(targets@Q@c)[:,None]
    labels=[2,3,4,5,6,8,10,12,16,24];finite=np.array([-2*np.cos(np.pi/m)for m in[2,3,4,5,6,7,8,9,10,12,15,16,17,20,24,30]])
    G=[[None]*29 for _ in range(29)]
    for i in range(29):G[i][i]=K(2)
    for row in src['all_pair_gates']:
        if row['type']=='dotted':continue
        i,j=row['pair'];G[i][j]=G[j][i]={'m2':Z,'m3':K(-1),'m5':I.neg(PHI)}[row['type']]
    for row in source['new_pairs']:
        if row['type']=='dotted':continue
        i,j=row['walls'];G[i][j]=G[j][i]=Z
    caps=[];keys=set();edges=[];stats={'profile_instances':0,'roots':0,'angle_pass':0}
    for released in incident:
        supp=[i for i in incident if i!=released];packet=spherical_profiles([[G[i][j]for j in supp]for i in supp],supp)
        profiles=[r['contacts']for r in packet['finite_profiles']]+[[[r['wall'],m]]for r in packet['dihedral_tail_profiles']for m in labels if m>=3]
        stats['profile_instances']+=len(profiles);edges.append({'released':released,'support':supp,'spherical_profiles':packet})
        for anchor in sorted(set(range(29))-set(supp)):
            M=W[supp+[anchor]]@Q
            if np.linalg.matrix_rank(M,tol=1e-9)!=9:continue
            z=null_space(M)[:,0];aa=z@Q@z
            if abs(aa)<1e-9:continue
            T=M.T@np.linalg.inv(M@M.T);rhs=[];meta=[]
            for pi,contacts in enumerate(profiles):
                r=np.zeros(9)
                for wid,m in contacts:r[supp.index(wid)]=-2*np.cos(np.pi/m)
                for m in labels:
                    rr=r.copy();rr[-1]=-2*np.cos(np.pi/m)if m>2 else 0;rhs.append(rr);meta.append((pi,m))
            u=np.array(rhs)@T.T;bb=2*u@Q@z;cc=np.einsum('ij,jk,ik->i',u,Q,u)-2;D=bb*bb-4*aa*cc;good=np.flatnonzero(D>1e-9)
            for branch in[-1,1]:
                R=u[good]+((-bb[good]+branch*np.sqrt(D[good]))/(2*aa))[:,None]*z;stats['roots']+=len(R)
                ps=R@Q@W.T;acute=(ps.max(axis=1)<1e-6)&(R@Q@c<0);R=R[acute];ps=ps[acute];idxs=good[acute]
                if not len(R):continue
                ok=ps<-2+1e-6
                for x in finite:ok|=abs(ps-x)<1e-6
                for r,rr,idx in zip(R[ok.all(axis=1)],ps[ok.all(axis=1)],idxs[ok.all(axis=1)]):
                    stats['angle_pass']+=1;key=tuple(np.rint(r*1e6).astype('int64'))
                    if key in keys:continue
                    keys.add(key);pi,m=meta[idx];caps.append({'id':len(caps),'support':supp,'released':released,'contacts':profiles[pi],'anchor':anchor,'anchor_label':m,'quadratic_sign':branch,'normal':r.tolist(),'all29_pairings':rr.tolist(),'target_pairings':(r@Q@targets.T).tolist()})
        print('released',released,'caps',len(caps),flush=True)
    out={'status':'NUMERICAL_NINE_NEW_EDGE_PROFILE_SCOUT','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in[CORE,SOURCE]},'edges':edges,'labels':labels,'counts':stats,'caps':caps,'targets':targets.tolist(),'target_cut_counts':(np.array([r['target_pairings']for r in caps])>1e-8).sum(axis=0).tolist()if caps else[0,0],'seconds':time.time()-start,'scope':'Finite spherical profiles at the one new corner, sampled I2-tail and ninth-anchor labels. Root identities, compatibility and target cuts remain numerical proposals; no volume conclusion.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print('DONE',len(caps),stats,out['target_cut_counts'],time.time()-start,flush=True)
if __name__=='__main__':build()
