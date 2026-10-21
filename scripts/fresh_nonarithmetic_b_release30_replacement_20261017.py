#!/usr/bin/env python3
"""Change the7/30/49 hinge, with all seven A2 profiles and whole-face scoring."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_nonarithmetic_b_release52_patch_20261017 as P
ROOT=P.ROOT;B=P.B;M=B.METRIC;TAG='fresh_nonarithmetic_b_release30_replacement_20261017'
ALPHABET=[2,3,4,5,6,8,10,12]
FRAMES=[[7,49,k]for k in[33,47,52,55,59,63,64,44]]+[[7,47,52],[49,47,52]]

def incidence(ns,actual):
    old=np.array([ns[i]for i in actual]);g=old@M@old.T;di=np.diag(g);eps=1e-8
    epairs=[(a,b)for a,b in combinations(range(len(actual)),2)if di[a]*di[b]-g[a,b]**2>eps*di[a]*di[b]]
    etri=[t for t in combinations(range(len(actual)),3)if np.linalg.eigvalsh(g[np.ix_(t,t)])[0]>eps]
    return old,g,epairs,etri

def compute():
    z,labels,ns,h,profiles=P.setup();retained=[i for i in labels if i!=30];actual=[i for i in z['actual_sides']if i!=30]
    old=np.array([ns[i]for i in retained]);poly,g,epairs,etri=incidence(ns,actual);d=np.diag(g)
    assert len(etri)==18
    baseends=2*len(epairs)-3*len(etri)
    ms=np.array(list(product(ALPHABET,repeat=3)));desired=-np.array([[B.theta(m)for m in row]for row in ms])
    rows=[];counts=Counter();scores=Counter();seen=set();best=[]
    for ids in FRAMES:
        f=np.array([ns[i]for i in ids]);fg=f@M@f.T
        if abs(np.linalg.det(fg))<1e-10:continue
        gi=np.linalg.inv(fg);k=np.linalg.svd(f@M)[2][-1];kq=k@M@k
        for pi,pr in enumerate(profiles):
            off=np.array(pr['offsets']);cc=desired-off[[labels.index(i)for i in ids]];pp=cc@gi@f
            rem=(pr['square']-np.einsum('ij,jk,ik->i',pp,M,pp))/kq;good=np.where(rem>1e-9)[0]
            for sign in[-1,1]:
                yy=pp[good]+sign*np.sqrt(rem[good])[:,None]*k;prod=yy@M@old.T+off[[labels.index(i)for i in retained]]
                for p,i in enumerate(ids):prod[:,retained.index(i)]=desired[good,p]
                center=yy@M@h+pr['center_offset'];counts['real_norm_branches']+=len(good)
                for local in np.where((center<-1e-8)&(prod.max(axis=1)<=2e-7))[0]:
                    rel=[B.relation(x)for x in prod[local]]
                    if any(r is None or r=='parallel_or_unresolved'for r in rel):counts['gap']+=1;continue
                    y=yy[local];key=(pi,)+tuple(np.round(y,7))
                    if key in seen:continue
                    seen.add(key);counts['accepted']+=1
                    products=poly@M@y;q=pr['square'];det=d*q-products**2
                    if np.any(np.abs(det)<1e-8*d*q):counts['intrinsic_parallel_unresolved']+=1;continue
                    ne=[j for j in range(len(actual))if det[j]>1e-8*d[j]*q];nv=[]
                    for a,b in epairs:
                        gg=np.array([[g[a,a],g[a,b],products[a]],[g[a,b],g[b,b],products[b]],[products[a],products[b],q]])
                        if np.linalg.eigvalsh(gg)[0]>1e-8:nv.append([actual[a],actual[b]])
                    ends=baseends+2*len(ne)-3*len(nv);scores[str(ends)]+=1
                    row=dict(index=len(rows),frame=ids,labels=ms[good[local]].tolist(),branch=sign,profile=pi,projection=y.tolist(),
                        finite_contacts=[[i,r]for i,r in zip(retained,rel)if isinstance(r,int)],
                        projected_new_edges=[actual[j]for j in ne],new_finite_corner_pairs=nv,
                        open_ends_from_complete_finite_incidence=ends,finite_count=18+len(nv))
                    rows.append(row)
                    if ends<15:
                        row['complete_numeric_boundary']=B.ledger(np.array([ns[i]for i in actual]+[y]),actual+[30])
                        assert row['complete_numeric_boundary']['open_ends']==ends
                        print('IMPROVEMENT',row,flush=True)
        print('FRAME',ids,'accepted',len(rows),'scores',dict(scores),flush=True)
    rows.sort(key=lambda r:(r['open_ends_from_complete_finite_incidence'],-r['finite_count']))
    return dict(status='BOUNDED_NUMERICAL_RELEASE30_HINGE_REPLACEMENT',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in[Path(__file__),Path(P.__file__),P.C.PARENT,
            ROOT/'results/fresh_nonarithmetic_replace47_close12_20261017.json']},
        deleted_wall=30,alphabet=ALPHABET,frames=FRAMES,profiles=profiles,counts=dict(counts),score_counts=dict(scores),
        after_deletion=dict(finite_vertices=18,geodesic_edges=len(epairs),open_ends=baseends),records=rows,
        incidence_formula='For an acute pointed H3 chamber with no ideal vertices, all positive-definite two-normal blocks give actual geodesic edges, all positive-definite triples give actual finite vertices, and openends=2E−3V. Values here are numerical diagnostics; interesting candidates also get full projective-cone enumeration.',
        scope='Finite three-contact event search for a genuinely changed7/30/49 hinge, all seven A2 profiles, all63 retained native gates. No all-label or continuous exclusion. No new chamber promoted without exact native and complete-boundary certificates.')

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=json.loads(json.dumps(compute()))
    path=ROOT/'results'/(TAG+'.json')
    if a.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],out['score_counts'],flush=True)
if __name__=='__main__':main()
