#!/usr/bin/env python3
"""Coupled mirror pair fixed by the seed6/7 reflection, with metric exchanges."""
import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scripts import fresh_h9_crossplane_repair_20261013_verify as T
V=T.V;PREFIX=T.PREFIX

def run(rounds=10,width=80,mutual=2):
    d=V.read(PREFIX+'_caps.json');G,I,h,g,old,names,r,center,_=T.reconstruct(d)
    f=lambda p:sum(p)/(2*V.SCALE)
    C=np.array([[f(p)for p in row]for row in old]);GG=np.array([[f(V.golden(p))for p in row]for row in G])
    hh=np.array([f(V.golden(p))for p in h]);cc=np.array([f(p)for p in center]);D=C[6]-C[7];assert abs(D@GG@D-4)<1e-8
    mirrored=C-((C@GG@D)/2)[:,None]*D;bank=np.vstack([C,mirrored]);cov=bank@GG
    targets=np.array([[f(V.golden(V.K.parse(p)))for p in V.read(d['source'])['escape']],[f(V.golden(V.K.parse(p)))for p in d['escape']]])
    labels=(2,3,4,5,6,8,10,12);vals=np.array([2*np.cos(np.pi/m)for m in labels]);value=dict(zip(labels,vals,strict=True))
    coupling=np.sqrt(4+4*np.cos(np.pi/mutual));seen=set();counts=Counter();accepted=[];history=[]
    seed=tuple(sorted((names.index(n),m)for n,m in d['constraints']))
    def affine_pair(pair):
        P=pair@GG@C.T;roots=np.vstack([C,pair]);H=roots@GG@roots.T
        for new in(88,89):
            neigh=[i for i in range(90)if i!=new and H[new,i]>-2+1e-7]
            if any(abs(H[new,i]+2)<1e-7 for i in range(90)if i!=new):return['parallel']
            for size in range(2,len(neigh)+1):
                for S in combinations(neigh,size):
                    M=H[np.ix_(S+(new,),S+(new,))]
                    if M.min()<-2+1e-7:continue
                    ev,U=np.linalg.eigh(M)
                    if abs(ev[0])<1e-7 and ev[1]>1e-6 and(np.all(U[:,0]>1e-6)or np.all(U[:,0]<-1e-6)):return list(S)+[new]
        return[]
    def evaluate(state):
        if state in seen:return[]
        seen.add(state);A=np.vstack([bank[[i for i,m in state]],D]);p=np.r_[[-value[m]for i,m in state],coupling];M=A@GG@A.T
        try:Mi=np.linalg.inv(M)
        except np.linalg.LinAlgError:return[]
        if np.linalg.norm(M@Mi-np.eye(9))>1e-5:return[]
        u=p@Mi@A;v=cc-(cc@GG@A.T)@Mi@A;aa=v@GG@v;bb=2*u@GG@v;c=u@GG@u-2
        if abs(aa)<1e-10:return[]
        disc=bb*bb-4*aa*c
        if disc<1e-9:return[]
        rows=[]
        for branch in(-1,1):
            counts['norm_branches']+=1;t=(-bb+branch*np.sqrt(disc))/(2*aa);r=u+t*v;s=r-(r@GG@D)/2*D;pair=np.array([r,s])
            if np.max(abs(A@GG@r-p))>1e-5 or abs(r@GG@r-2)>1e-5:continue
            if np.min(pair@hh)<1e-6 or np.min(np.max(-pair@targets.T,axis=0))<1e-6:continue
            pp=bank@GG@r;error=np.abs(pp[:,None]+vals).min(axis=1);error=np.where(pp<-2+1e-8,0,error);bad=np.flatnonzero(error>1e-7)
            aff=affine_pair(pair)if not len(bad)else[];score=(len(bad)+3*bool(aff),float(error@error)+.01*bool(aff),float(error.max()))
            row={'state':state,'branch':branch,'score':score,'bad':bad,'error':error,'pairings':pp,'pair':pair,'affine':aff,'target_cuts':-pair@targets.T}
            rows.append(row)
            if not len(bad)and not aff:accepted.append(row)
        return rows
    beam=[]
    for j in range(9):beam+=evaluate(seed[:j]+seed[j+1:])
    assert beam;beam=sorted(beam,key=lambda r:r['score'])[:width];best=beam[0]
    for depth in range(rounds):
        proposals=[]
        for row in beam:
            state=row['state'];ids={i for i,m in state}
            if row['affine']and row['affine']!=['parallel']:
                for i in row['affine']:
                    if i>=88:continue
                    for j,(idx,m)in enumerate(state):
                        if idx%88!=i or m==2:continue
                        k=labels.index(m);proposals+=evaluate(tuple(sorted(state[:j]+state[j+1:]+((idx,labels[k-1]),))))
            bad=sorted((int(i)for i in row['bad']if i not in ids),key=lambda i:-row['error'][i])[:6]
            for i in bad:
                for k in np.argsort(abs(row['pairings'][i]+vals))[:2]:
                    for j in range(8):proposals+=evaluate(tuple(sorted(state[:j]+state[j+1:]+((i,labels[int(k)]),))))
        if not proposals:break
        beam=sorted(proposals,key=lambda r:r['score'])[:width]
        if beam[0]['score']<best['score']:best=beam[0]
        history.append({'round':depth+1,'states':len(seen),'best_score':list(best['score']),'bad':[(names[i%88],i>=88)for i in best['bad']],'accepted':len(accepted)})
        print(history[-1],flush=True)
        if accepted:break
    def save(row):
        return{'constraints':[{'old':names[i%88],'reflected':i>=88,'label':m}for i,m in row['state']],'branch':row['branch'],
               'score':list(row['score']),'bad':[{'old':names[i%88],'reflected':bool(i>=88)}for i in row['bad']],
               'normals':row['pair'].tolist(),'pairings':row['pairings'].tolist(),'affine':row['affine'],'target_cuts':row['target_cuts'].tolist()}
    out={'status':'NUMERICAL_REFLECTION_COUPLED_CAP_REPAIR','source':PREFIX+'_caps.json','source_sha256':V.sha(PREFIX+'_caps.json'),
        'retained_source':d['source'],'reflection_walls':['seed6','seed7'],'mutual_label':mutual,'round_limit':rounds,'beam_width':width,
        'counts':dict(counts),'states':len(seen),'history':history,'candidates':[save(r)for r in accepted],'best_failed':save(best),
        'scope':'Each proposal replaces the89th cap by a pair retaining theoriginal88walls, cutsbothsavedgaps, and has prescribedmutualangle viaonefixedlinearreflectioncoupling. Boundednumericalbeam; exactcheckingrequired.'}
    (V.ROOT/(PREFIX+'_pair'+str(mutual)+'_scout.json')).write_text(json.dumps(out,indent=2)+'\n');print('DONE',len(accepted),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--rounds',type=int,default=10);p.add_argument('--width',type=int,default=80);p.add_argument('--mutual',type=int,default=2);a=p.parse_args();run(a.rounds,a.width,a.mutual)
