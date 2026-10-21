#!/usr/bin/env python3
"""Beam search exchanging finite-angle constraints around a target cutter."""
import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scripts import fresh_h9_remote_next_20261013_bridge_verify as B
from scripts.fresh_h9_crossplane_repair_20261013_scout import INCIDENT,PREFIX
V=B.V

def run(rounds=8,width=48,compact=False):
    source=B.PREFIX+'_caps.json';d=V.read(source)
    G,I,h,g,old,names,caps,center,_=B.reconstruct(d);old+=caps;names+=['bridge0','bridge1']
    f=lambda p:sum(p)/(2*V.SCALE)
    C=np.array([[f(p)for p in r]for r in old]);GG=np.array([[f(V.golden(p))for p in r]for r in G]);cov=C@GG
    hh=np.array([f(V.golden(p))for p in h]);y=np.array([f(V.golden(V.K.parse(p)))for p in d['escape']]);cc=np.array([f(p)for p in center])
    labels=(2,3,4,5,6,8,10,12);values=np.array([2*np.cos(np.pi/m)for m in labels]);value=dict(zip(labels,values,strict=True))
    if compact:
        initial=V.read(PREFIX+'_exchange_scout.json')['candidates'][0]
        seed=tuple(sorted((names.index(n),m)for n,m in initial['constraints']))
    else:
        initial=V.read(PREFIX+'_face_scout.json')['best_failed'];seedlabels=dict(zip(initial['nonright_support'],initial['labels'],strict=True))
        seed=tuple(sorted((names.index(n),seedlabels.get(n,2))for n in INCIDENT))
    seen=set();counts=Counter();accepted=[];history=[];oldgram=C@GG@C.T
    def affine(P):
        parallel=np.flatnonzero(abs(P+2)<1e-7)
        if len(parallel):return[int(parallel[0])]
        neigh=[i for i in range(88)if P[i]>-2+1e-7]
        for size in range(2,len(neigh)+1):
            for S in combinations(neigh,size):
                H=oldgram[np.ix_(S,S)]
                if H.min()<-2+1e-7:continue
                M=np.block([[H,P[list(S),None]],[P[None,list(S)],np.array([[2.]])]])
                ev,U=np.linalg.eigh(M)
                if abs(ev[0])<1e-7 and ev[1]>1e-6 and np.min(abs(U[:,0]))>1e-6 and (np.all(U[:,0]>0)or np.all(U[:,0]<0)):
                    return list(S)
        return[]
    def evaluate(state):
        if state in seen:return[]
        seen.add(state);ids=[i for i,m in state];p=-np.array([value[m]for i,m in state]);A=C[ids];M=A@GG@A.T
        try:Mi=np.linalg.inv(M)
        except np.linalg.LinAlgError:return[]
        if np.linalg.norm(M@Mi-np.eye(9))>1e-5:return[]
        u=p@Mi@A;v=cc-(cc@GG@A.T)@Mi@A;aa=v@GG@v;bb=2*u@GG@v;cc0=u@GG@u-2
        if abs(aa)<1e-10:return[]
        disc=bb*bb-4*aa*cc0
        if disc<1e-9:return[]
        out=[]
        for branch in(-1,1):
            counts['norm_branches']+=1;t=(-bb+branch*np.sqrt(disc))/(2*aa);r=u+t*v
            if r@hh<=1e-6 or -r@y<=1e-6:continue
            P=cov@r
            if max(abs(P[ids]-p))>2e-5 or abs(r@GG@r-2)>2e-5:continue
            error=np.abs(P[:,None]+values).min(axis=1);error=np.where(P<-2+1e-8,0,error)
            bad=np.flatnonzero(error>1e-7);hard=int(sum(i<16 for i in bad));aff=affine(P)if compact and not len(bad)else[]
            score=(len(bad)+2*bool(aff),float(error@error)+.01*bool(aff),float(error.max()),hard)
            row={'state':state,'branch':branch,'score':score,'bad':bad,'affine':aff,'error':error,'pairings':P,'coefficients':r,'target_cut':float(-r@y)}
            out.append(row)
            if not len(bad)and not aff:accepted.append(row)
        return out
    beam=evaluate(seed);best=min(beam,key=lambda r:r['score'])
    for depth in range(rounds):
        proposals=[]
        for row in beam:
            state=row['state'];ids={i for i,m in state}
            if row['affine']:
                for i in row['affine']:
                    for j,(idx,m)in enumerate(state):
                        if idx!=i or m==2:continue
                        k=labels.index(m);proposals.extend(evaluate(tuple(sorted(state[:j]+state[j+1:]+((i,labels[k-1]),)))))
            bad=sorted((int(i)for i in row['bad']if i not in ids),key=lambda i:(i>=16,-row['error'][i]))[:5]
            for i in bad:
                nearest=np.argsort(abs(row['pairings'][i]+values))[:2]
                for k in nearest:
                    for j in range(9):
                        nextstate=tuple(sorted(state[:j]+state[j+1:]+((i,labels[int(k)]),)))
                        proposals.extend(evaluate(nextstate))
        if not proposals:break
        beam=sorted(proposals,key=lambda r:r['score'])[:width]
        if beam[0]['score']<best['score']:best=beam[0]
        history.append({'round':depth+1,'states':len(seen),'best_score':list(best['score']),'best_bad_walls':[names[i]for i in best['bad']],'accepted':len(accepted)})
        print(history[-1],flush=True)
        if accepted:break
    def save(row):
        return {'constraints':[[names[i],m]for i,m in row['state']],'branch':row['branch'],'score':list(row['score']),
                'bad_walls':[names[i]for i in row['bad']],'affine_walls':[names[i]for i in row['affine']],'pairings':list(map(float,row['pairings'])),
                'coefficients':list(map(float,row['coefficients'])),'target_cut':row['target_cut']}
    out={'status':'NUMERICAL_FINITE_CONSTRAINT_EXCHANGE_SCOUT','source':source,'source_sha256':V.sha(source),'round_limit':rounds,'beam_width':width,
         'states':len(seen),'counts':dict(counts),'history':history,'candidates':[save(r)for r in accepted],'best_failed':save(best),
         'scope':'Adaptive finite-angle constraint exchanges; bounded beam and label set. All proposed roots cut the fixed88-wall escape and preserve its strict center numerically. Exact cap and volume verification is required.'}
    (V.ROOT/(PREFIX+('_compact_exchange_scout.json'if compact else'_exchange_scout.json'))).write_text(json.dumps(out,indent=2)+'\n');print('DONE',len(accepted),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--rounds',type=int,default=8);p.add_argument('--width',type=int,default=48);p.add_argument('--compact',action='store_true');a=p.parse_args();run(a.rounds,a.width,a.compact)
