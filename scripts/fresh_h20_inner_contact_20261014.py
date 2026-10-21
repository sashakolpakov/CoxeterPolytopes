#!/usr/bin/env python3
"""Bounded central-contact inner-cap construction and joint probe-cover scout."""
import argparse
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import numpy as np
from fresh_h20_balanced_20261014 import data

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_inner_contact_20261014'


def prepare(central_labels=(2,),max_deviations=2):
    d=data();iv=np.linalg.inv(d['spatial_Gram']);L=d['cholesky']
    profiles=[np.full(19,3,dtype=np.int8)]
    for k in range(1,max_deviations+1):
        for inds in combinations(range(19),k):
            for labels in product((4,5,6),repeat=k):
                p=np.full(19,3,dtype=np.int8)
                for i,m in zip(inds,labels):p[i]=m
                profiles.append(p)
    labels=np.array(profiles);q=np.cos(np.pi/labels)
    aa=np.einsum('ij,ij->i',q@iv,q)-1;bb=2*q@(iv@np.ones(19))
    sums=np.column_stack((q[:,:6].sum(1),q[:,6:12].sum(1),q[:,12:].sum(1)))
    finite=np.cos(np.pi/np.arange(2,101))
    blocks=[]
    for m in central_labels:
        cc=np.cos(np.pi/m)if m else 1.
        linear=bb+4*cc;disc=linear*linear-768*aa
        rr=(linear-np.sqrt(np.maximum(0,disc)))/(2*aa)
        omitted=np.array([7,7,8])[None,:]/rr[:,None]-sums
        gate=(disc>=0)&(rr>0)&(rr<2)&np.all((omitted>=1-1e-8)|
                  (np.min(abs(omitted[:,:,None]-finite),axis=2)<1e-8),axis=1)
        blocks.append((np.where(gate)[0],labels[gate],q[gate],rr[gate],omitted[gate],
                       np.full(sum(gate),m,dtype=int)))
    indices,labels,q,rr,omitted,central=np.concatenate([b[0]for b in blocks]),\
        np.concatenate([b[1]for b in blocks]),np.concatenate([b[2]for b in blocks]),\
        np.concatenate([b[3]for b in blocks]),np.concatenate([b[4]for b in blocks]),\
        np.concatenate([b[5]for b in blocks])
    gaps=np.column_stack((omitted[:,0],q[:,:6],omitted[:,1],q[:,6:12],omitted[:,2],q[:,12:]))*rr[:,None]
    base=(1-rr[:,None]*q)@iv@L
    vertices=np.array([v['point']for v in d['vertices']]);target=.2*vertices
    powers=np.sum(base*base,axis=1)[:,None]+7.84-2*base@target.T-rr[:,None]**2
    cover=[np.where(row<-1e-8)[0] for row in powers]
    width=max(map(len,cover));local=np.array(list(product(range(7),range(7),range(8))))
    centers=[];radii=[];targets=[]
    for shift in local:
        rotated=np.concatenate([np.roll(gaps[:,lo:hi],int(s),axis=1)
                                for (lo,hi),s in zip(((0,7),(7,14),(14,22)),shift)],axis=1)
        c=(1-rotated[:,[j for j in range(22)if j not in(0,7,14)]])@iv@L
        centers.append(c);radii.append(rr)
        block=np.full((len(rr),width),392,dtype=np.int16)
        for i,ids in enumerate(cover):
            xyz=(local[ids]+shift)%np.array([7,7,8])
            block[i,:len(ids)]=xyz[:,0]*56+xyz[:,1]*8+xyz[:,2]
        targets.append(block)
    C=np.concatenate(centers);R=np.concatenate(radii);T=np.concatenate(targets)
    co=np.array([np.cos(np.pi/m)if m else 1. for m in central])
    assert np.max(abs(np.sum(C*C,axis=1)-R*R-4-4*R*np.tile(co,392)))<1e-7
    return C,R,T,dict(profiles=labels.tolist(),original_indices=indices.tolist(),
                       central_labels=central.tolist(),radii=rr.tolist(),seed_gaps=gaps.tolist(),base_probes=[x.tolist()for x in cover])


def run(args):
    C,R,T,base=prepare(args.central_labels,args.max_deviations);N=np.sum(C*C,axis=1)
    if args.keep_middle:
        from fresh_h20_balanced_integer_joint_20261014 import prepare as middle_bank
        d,MC,_=middle_bank()
        old=json.loads((ROOT/'results/fresh_h20_balanced_integer_joint_20261014_motif.json').read_text())
        chosen=np.array([MC[i,j]for i,j in old['best']['selected']])
    else:chosen=np.empty((0,19))
    finite=np.cos(np.pi/np.arange(2,61));initial=np.ones(len(R),dtype=bool)
    for c in chosen:
        cosine=(N+c@c-2*C@c-R*R-4)/(4*R)
        allowed=cosine>=1-1e-8
        for cc in finite:allowed|=abs(cosine-cc)<1e-8
        initial&=allowed
    rng=np.random.default_rng(20261014);attempts=[]
    incidence=[np.where(np.any(T==i,axis=1))[0] for i in range(392)]
    best=None
    for attempt in range(args.restarts):
        live=initial.copy();uncovered=np.ones(393,dtype=bool);uncovered[392]=False;selection=[]
        while uncovered.any():
            counts=np.bincount(T[live].ravel(),minlength=393)[:392]
            missing=np.where(uncovered[:392])[0]
            m=counts[missing].min()
            if not m:break
            worst=missing[counts[missing]==m]
            probe=int(worst[0]if attempt==0 else rng.choice(worst))
            ids=incidence[probe];ids=ids[live[ids]]
            gain=uncovered[T[ids]].sum(axis=1);best_ids=ids[gain==gain.max()]
            if attempt==0:
                # Prefer the smallest central-right radius at equal gain.
                k=int(best_ids[np.argmin(R[best_ids])])
            else:k=int(rng.choice(best_ids))
            selection.append(k);uncovered[T[k]]=False
            cosine=(N+N[k]-2*C@C[k]-R*R-R[k]**2)/(2*R*R[k])
            compatible=cosine>=1-1e-8
            for v in finite:compatible|=abs(cosine-v)<1e-8
            # Same literal cap is harmless; exclude it from further selection.
            live&=compatible
        row=dict(attempt=attempt,covered=int(392-uncovered[:392].sum()),
                 chosen=len(selection),selected=selection,
                 uncovered=np.where(uncovered[:392])[0].tolist(),
                 remaining_candidates=int(live.sum()))
        attempts.append(row)
        if best is None or row['covered']>best['covered']:best=row
        print('attempt',attempt,'covered',row['covered'],'caps',row['chosen'],flush=True)
        if row['covered']==392:break
    paths=[f'scripts/{TAG}.py','scripts/fresh_h20_balanced_20261014.py']
    if args.keep_middle:paths+=['results/fresh_h20_balanced_integer_joint_20261014_motif.json',
                               'scripts/fresh_h20_balanced_integer_joint_20261014.py']
    out=dict(status='BOUNDED_CENTRAL_CONTACT_INNER_CAP_JOINT_SCOUT',
             sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
             base_recipe_count=len(base['radii']),candidate_count=len(R),
             initial_candidates=int(initial.sum()),keep_middle=args.keep_middle,
             central_labels=args.central_labels,max_deviations=args.max_deviations,
             target='all392 z/5 probes',max_individual_probe_coverage=int(np.max((T<392).sum(1))),
             base=base,attempts=attempts,best=best,
             scope='Numerical scout with19 prescribed seed contacts, all3 except the stated maximum changed to4/5/6; stated central labels (0 means parallel), omitted labels2..100 or dotted. All392 cyclictransports. Greedy failure is not infeasibility; probe coverage is not finite volume. Exact final geometry must be independently certified.')
    suffix='middle'if args.keep_middle else'free'
    if args.central_labels!=[2]or args.max_deviations!=2:
        suffix+='_m'+'_'.join(map(str,args.central_labels))+'_d'+str(args.max_deviations)
    (ROOT/f'results/{TAG}_{suffix}.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--restarts',type=int,default=8)
    ap.add_argument('--central-labels',type=int,nargs='+',default=[2])
    ap.add_argument('--max-deviations',type=int,default=2)
    ap.add_argument('--keep-middle',action='store_true');run(ap.parse_args())
