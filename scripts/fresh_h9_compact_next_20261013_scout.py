#!/usr/bin/env python3
"""Finite-corner insertion with the boost fixed by an old-wall angle."""
import json
from itertools import permutations
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scripts import fresh_h9_crossplane_repair_20261013_pair_verify as P
V=P.V;K=P.K
PREFIX='data/enumerations/fresh_h9_compact_next_20261013'
BLOCKS=[['crosspair0','seed1','seed2','seed6','wall17'],['crosspair0','nextH4','seed1','seed2','wall17']]
SIGNS=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])

def run():
    d=V.read(P.PREFIX+'_caps.json');G,I,h,g,old,names,caps,center,_=P.reconstruct(d)
    ff=lambda x:sum(x)/(2*V.SCALE)
    C=np.array([[ff(x)for x in r]for r in old+caps]);names+=P.NEW;Q=np.array([[ff(V.golden(x))for x in row]for row in G])
    hh=np.array([ff(V.golden(x))for x in h]);cc=np.array([ff(x)for x in center]);yy=np.array([ff(V.golden(K.parse(x)))for x in d['escape']]);w=-np.linalg.solve(Q,yy)
    labels=[2,3,4,5,6];val=np.array([2*np.cos(np.pi/m)for m in labels]);counts={'frames':0,'real_boosts':0,'target_cutters':0,'compatible':0};accepted=[];seen=set();best=None
    for ss in BLOCKS:
        ids=[names.index(n)for n in ss];A=C[ids];M=A@Q@A.T;Mi=np.linalg.inv(M)
        def project(v):return v-(v@Q@A.T)@Mi@A
        t=project(cc);qn=t@Q@t;assert qn<0;t/=np.sqrt(-2*qn)
        v=project(w);v-=((v@Q@t)/(t@Q@t))*t;qn=v@Q@v;assert qn>0;v/=np.sqrt(2*qn)
        rest=[i for i in range(90)if i not in ids]
        for src in permutations(range(10),3):
            frame=[t,v]
            for j in src:
                z=project(np.eye(10)[j])
                for b in frame:z-=((z@Q@b)/(b@Q@b))*b
                q=z@Q@z
                if q<1e-8:break
                frame.append(z/np.sqrt(2*q))
            if len(frame)!=5:continue
            for sign in(-1,1):
                F=np.array(frame);F[-1]*=sign
                key=tuple(np.round(F[2:].ravel(),9))
                if key in seen:continue
                seen.add(key);counts['frames']+=1
                spatial=SIGNS@F[2:];aa=(t@Q@C.T+v@Q@C.T)/2;bb=(v@Q@C.T-t@Q@C.T)/2;dd=spatial@Q@C.T
                assert np.max(aa[rest])<0
                # A k²+(d+2cos(pi/m))k+B=0 forces one exact angle.
                for ci in range(4):
                    a=np.repeat(aa[rest],len(labels));b=np.repeat(bb[rest],len(labels));e=np.repeat(dd[ci,rest],len(labels))+np.tile(val,len(rest))
                    delta=e*e-4*a*b;mask=delta>1e-10
                    for branch in(-1,1):
                        ks=(-e[mask]+branch*np.sqrt(delta[mask]))/(2*a[mask]);js=np.flatnonzero(mask)
                        for k,idx in zip(ks,js,strict=True):
                            if k<=1+1e-8:continue
                            counts['real_boosts']+=1;ct=(k+1/k)/2;st=(k-1/k)/2;R=st*t+ct*v+spatial
                            cuts=-R@yy
                            if cuts.min()<=1e-6 or(R@hh).min()<=1e-6:continue
                            counts['target_cutters']+=1;pp=R@Q@C.T
                            error=np.abs(pp[:,:,None]+val).min(axis=2);error=np.where(pp<-2-1e-7,0,error);error[:,ids]=0
                            bad=np.argwhere(error>2e-7);score=(len(bad),float(np.sum(error*error)))
                            rec={'orthogonal_old_walls':ss,'frame_sources':list(src),'last_frame_sign':sign,'anchor_cap':ci,'anchor_old':names[rest[idx//len(labels)]],
                                 'anchor_label':labels[idx%len(labels)],'boost':float(k),'boost_branch':branch,'score':list(score),'target_cuts':cuts.tolist(),
                                 'normals':R.tolist(),'bad':[[int(i),names[int(j)]]for i,j in bad]}
                            if best is None or score<tuple(best['score']):best=rec
                            if len(bad):continue
                            # Exclude accidental finite contacts beyond the prescribed one.
                            finite=[(int(i),names[int(j)])for i,j in np.argwhere(pp>-2+1e-6)if j not in ids]
                            if finite!=[(ci,rec['anchor_old'])]:continue
                            counts['compatible']+=1;accepted.append(rec)
                if counts['frames']%20==0:print(counts,'accepted',len(accepted),flush=True)
    out={'status':'NUMERICAL_METRIC_FIXED_FOUR_CAP_CORNER_SCOUT','source':P.PREFIX+'_caps.json','source_sha256':V.sha(P.PREFIX+'_caps.json'),
         'blocks':BLOCKS,'labels':labels,'counts':counts,'candidates':accepted,'best_failed':best,
         'scope':'Two rank-five spherical blocks meeting crosspair0, coordinate Gram-Schmidt frames, two tetrahedral orientations, one prescribed finite anchor angle, and both norm branches. Numerical proposals only.'}
    (V.ROOT/(PREFIX+'_scout.json')).write_text(json.dumps(out,indent=2)+'\n');print('DONE',counts,flush=True)

if __name__=='__main__':run()
