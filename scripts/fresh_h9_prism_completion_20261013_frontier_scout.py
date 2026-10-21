#!/usr/bin/env python3
"""Targeted finite endpoints on the newly inserted three-vertex component."""
import json,sys,itertools
from pathlib import Path
from fractions import Fraction
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scripts import fresh_h9_prism_completion_20261013_scout as S
T=S.B;V=S.V

def run():
    d=S.setup();G=np.array(d['gram']);old=np.array(d['normals']);perm=[0,2,1,3,4,6,5,7,8,9]
    R=np.vstack([old,old[90:,perm]]);names=d['names']+['mirror'+n[0].upper()+n[1:]for n in d['names'][90:]]
    source=S.PREFIX+'_caps.json';prior=V.read(S.PREFIX+'_verify.json');recipe=V.read(source)
    h=np.array(d['h']);y=np.array([float(Fraction(x['a']))for x in recipe['escape']]);cov=R@G
    labels=(2,3,4,5,6,8,10,12);rows=[];stats=[];seen=set()
    edges=[e for e in prior['new_component_frontier']['edges']if e['status']=='open_to_light_cone']
    for ei,edge in enumerate(edges):
        ids=[names.index(n)for n in edge['common_facets']];M=R[ids]@G@R[ids].T;Mi=np.linalg.inv(M)
        profiles=[np.zeros(8)]
        for i in range(8):
            for m in labels[1:]:
                p=np.zeros(8);p[i]=-2*np.cos(np.pi/m)
                if p@Mi@p<2-1e-8:profiles.append(p)
        for i,j in itertools.combinations(range(8),2):
            for m,n in itertools.product((3,4,5),repeat=2):
                p=np.zeros(8);p[i]=-2*np.cos(np.pi/m);p[j]=-2*np.cos(np.pi/n)
                if p@Mi@p<2-1e-8:profiles.append(p)
        pp=np.array(profiles);stat={'edge':ei,'profiles':len(pp),'equations':0,'target_cut_branches':0,'retains_all':0,'at_most_two_optional_conflicts':0}
        for ai in range(len(R)):
            if ai in ids:continue
            A=cov[ids+[ai]];_,s,v=np.linalg.svd(A)
            if s[-1]<1e-8:continue
            inv=A.T@np.linalg.inv(A@A.T);direction=v[-1];qa=direction@G@direction
            if abs(qa)<1e-10:continue
            for label in labels:
                rhs=np.column_stack([pp,np.full(len(pp),0 if label==2 else -2*np.cos(np.pi/label))]);base=rhs@inv.T
                qb=2*base@G@direction;qc=np.einsum('ij,jk,ik->i',base,G,base)-2;disc=qb*qb-4*qa*qc;stat['equations']+=len(pp)
                for branch in(-1,1):
                    inds=np.flatnonzero(disc>1e-10)
                    rr=base[inds]+(((-qb[inds]+branch*np.sqrt(disc[inds]))/(2*qa))[:,None])*direction
                    filt=(rr@h>1e-7)&(rr@y<-1e-7);rr=rr[filt];inds=inds[filt];stat['target_cut_branches']+=len(rr)
                    for r,pi in zip(rr,inds):
                        pair=cov@r
                        valid=(pair<-2-1e-7)|np.logical_or.reduce([abs(pair+2*np.cos(np.pi/m))<1e-6 for m in labels])
                        bad=np.flatnonzero(~valid)
                        if len(bad)>2 or any(j<16 for j in bad):continue
                        key=tuple(np.round(r,7))
                        if key in seen:continue
                        seen.add(key);stat['at_most_two_optional_conflicts']+=1;stat['retains_all']+=int(not len(bad))
                        rows.append({'edge_index':ei,'edge':edge,'anchor':names[ai],'anchor_label':label,'profile_labels':[2 if abs(x)<1e-8 else min(labels[1:],key=lambda m:abs(x+2*np.cos(np.pi/m)))for x in pp[pi]],'normal':r.tolist(),'target_cut':float(-r@y),'conflicts':[names[j]for j in bad],
                                     'finite_relations':{names[j]:min(labels,key=lambda m:abs(pair[j]+2*np.cos(np.pi/m)))for j in range(len(R))if valid[j]and pair[j]>-2+1e-6}})
        stats.append(stat);print(stat,flush=True)
        out={'status':'NUMERICAL_NEW_COMPONENT_TARGETED_ENDPOINT_SCOUT','source':source,'source_sha256':V.sha(source),'stats':stats,'candidates':rows,'scope':'New23 open edges only; profiles with at most two finite nonright pairings and one extra finite old anchor. Targeted to the new102-wall escape, preserving all16 seed walls and allowing at most two optional conflicts. Not exhaustive over all caps.'}
        (V.ROOT/(S.PREFIX+'_frontier_scout.json')).write_text(json.dumps(out,indent=2)+'\n')
    print('DONE',len(rows),'ALL_RETAINED',sum(not x['conflicts']for x in rows),flush=True)

if __name__=='__main__':run()
