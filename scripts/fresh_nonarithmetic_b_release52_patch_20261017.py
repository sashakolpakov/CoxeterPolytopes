#!/usr/bin/env python3
"""Bounded coordinated discovery after deleting52, all seven A2 profiles.

Three-old-contact events seed a finite compatibility graph. Full native
products include the white A2 offsets. Only complete boundary improvement
can nominate a cluster; this floating-point scout proves no exclusion.
"""
import argparse,json
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_nonarithmetic_replace47_b_pair_20261017 as C

ROOT=C.ROOT;B=C.B;TAG='fresh_nonarithmetic_b_release52_patch_20261017'
ALPHABET=[2,3,4,5,6,8,10,12]
ANCHORS=[7,8,12,30,47,64,65]
J=np.array([[2.,1.],[1.,2.]])/3

def setup():
    import fresh_nonarithmetic_replace47_close12_20261017 as G
    z,labels,ns,h=C.setup();exact=G.geometry();pair=exact['pair'];H=exact['common_center'];N=exact['normals']
    oldwhite=np.array([[B.mid(G.literal(exact,j,i).enc())for j in[5,6]]for i in labels])
    whiteH=np.array([B.mid(pair(N[i],H).enc())for i in[5,6]])
    whiteS=[[B.mid(G.literal(exact,i,j).enc())for j in G.S]for i in G.S]
    productsH=np.array([B.mid(pair(N[i],H).enc())for i in G.S])
    tau=np.sqrt(-(B.mid(pair(H,H).enc())-productsH@np.linalg.inv(whiteS)@productsH))
    profiles=[dict(node=None,label=2,contacts=[0.,0.])]
    profiles +=[dict(node=i,label=m,contacts=[-B.theta(m)if i==j else 0. for j in[5,6]])for i,m in product([5,6],[3,4,5])]
    for pr in profiles:
        c=np.array(pr['contacts']);pr['square']=float(2-c@J@c);pr['offsets']=(oldwhite@J@c).tolist();pr['center_offset']=float(whiteH@J@c/tau)
    return z,labels,ns,h,profiles

def compute():
    z,labels,ns,h,profiles=setup();retained=[i for i in labels if i!=52];ix=[labels.index(i)for i in retained]
    old=np.array([ns[i]for i in retained]); actual=[i for i in z['actual_sides']if i!=52]
    ms=np.array(list(product(ALPHABET,repeat=3)));desired=-np.array([[B.theta(m)for m in row]for row in ms])
    counts=Counter();byprofile=Counter();rows=[];seen=set()
    for ids in combinations(ANCHORS,3):
        f=np.array([ns[i]for i in ids]);gram=f@B.METRIC@f.T
        if abs(np.linalg.det(gram))<1e-10:counts['singular_frame']+=1;continue
        inv=np.linalg.inv(gram);k=np.linalg.svd(f@B.METRIC)[2][-1];kq=k@B.METRIC@k
        if abs(kq)<1e-10:counts['null_kernel']+=1;continue
        for pi,pr in enumerate(profiles):
            off=np.array(pr['offsets']);cc=desired-off[[labels.index(i)for i in ids]];pp=cc@inv@f
            rem=(pr['square']-np.einsum('ij,jk,ik->i',pp,B.METRIC,pp))/kq
            finite=np.zeros(len(ms),dtype=bool)
            for a,b in combinations(range(3),2):
                g2=gram[np.ix_([a,b],[a,b])]
                if np.linalg.eigvalsh(g2)[0]<=1e-9:continue
                finite|=pr['square']-np.einsum('ij,jk,ik->i',cc[:,[a,b]],np.linalg.inv(g2),cc[:,[a,b]])>1e-8
            good=np.where(rem>1e-9)[0];counts['real_norm_branches']+=2*len(good)
            for sign in[-1,1]:
                yy=pp[good]+sign*np.sqrt(rem[good])[:,None]*k
                prod=yy@B.METRIC@old.T+off[ix]
                for pos,i in enumerate(ids):prod[:,retained.index(i)]=desired[good,pos]
                center=yy@B.METRIC@h+pr['center_offset']
                for local in np.where((center<-1e-8)&(prod.max(axis=1)<=2e-7))[0]:
                    rel=[B.relation(x)for x in prod[local]]
                    if any(r is None or r=='parallel_or_unresolved'for r in rel):counts['finite_angle_gap']+=1;continue
                    y=yy[local];key=(pi,)+tuple(np.round(y,7))
                    if key in seen:counts['duplicate']+=1;continue
                    seen.add(key);byprofile[str(pi)]+=1
                    row=dict(index=len(rows),profile=pi,frame=list(ids),labels=ms[good[local]].tolist(),branch=sign,
                        projection=y.tolist(),finite_contacts=[[i,r]for i,r in zip(retained,rel)if isinstance(r,int)],
                        has_old_old_new_finite_corner=bool(finite[good[local]]),center_product=float(center[local]))
                    rows.append(row)
        print('FRAME',ids,'pool',len(rows),dict(byprofile),flush=True)
    pairs=[];compatible=defaultdict(set);finite_edges=[];relation_counts=Counter()
    ys=np.array([r['projection']for r in rows]);cs=np.array([profiles[r['profile']]['contacts']for r in rows])
    if len(rows):
        products=ys@B.METRIC@ys.T+cs@J@cs.T
        for i,j in combinations(range(len(rows)),2):
            q=float(products[i,j]);rel=B.relation(q);relation_counts[str(rel)]+=1
            if rel is None or rel=='parallel_or_unresolved':continue
            compatible[i].add(j);compatible[j].add(i)
            if isinstance(rel,int):finite_edges.append((i,j));pairs.append(dict(indices=[i,j],product=q,label=rel))
    # Connected compatible triples must contain a finite two-edge tree.
    fn=defaultdict(set)
    for a,b in finite_edges:fn[a].add(b);fn[b].add(a)
    triples=set()
    for center,near in fn.items():
        for a,b in combinations(sorted(near),2):
            if b in compatible[a]:triples.add(tuple(sorted([a,b,center])))
    batches=[]
    for selection in sorted(triples):
        labs=actual+[1000+i for i in selection]
        norms=np.array([ns[i]for i in actual]+[rows[i]['projection']for i in selection])
        ledger=B.ledger(norms,labs)
        batches.append(dict(candidates=list(selection),boundary=ledger))
    batches.sort(key=lambda r:(r['boundary']['open_ends'],r['boundary']['exterior'],-r['boundary']['finite']))
    return dict(status='BOUNDED_NUMERICAL_RELEASE52_SEVEN_PROFILE_PATCH_SCOUT',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in[Path(__file__),Path(C.__file__),C.PARENT,
            ROOT/'scripts/fresh_nonarithmetic_replace47_close12_20261017.py',ROOT/'results/fresh_nonarithmetic_replace47_close12_20261017.json']},
        deleted_old_wall=52,alphabet=ALPHABET,anchor_pool=ANCHORS,profiles=profiles,retained_labels=retained,
        old_projected_normals={str(i):ns[i].tolist()for i in retained},common_projected_center=h.tolist(),
        baseline=B.ledger(np.array([ns[i]for i in z['actual_sides']]),z['actual_sides']),
        after_deletion=B.ledger(np.array([ns[i]for i in actual]),actual),
        counts=dict(counts),candidate_counts_by_profile=dict(byprofile),candidates=rows,
        pair_relation_counts=dict(relation_counts),finite_pairs=pairs,compatible_triples=len(triples),batches=batches,
        scope='All seven A2 facet profiles; finite angle samples2,3,4,5,6,8,10,12 on three distinct old anchors. No old-old-new finite corner is required: bridging facets with only new finite corners are retained. All63 retained native products and common-center signs checked numerically, including white offsets. Connected compatible triples scored on the entire numerical S6 projective cone. No exhaustive continuous-domain or all-label exclusion and no exact promoted geometry.')

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=json.loads(json.dumps(compute()))
    path=ROOT/'results'/(TAG+'.json')
    if a.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],len(out['candidates']),len(out['finite_pairs']),out['compatible_triples'],flush=True)
    if out['batches']:print('BEST',out['batches'][0],flush=True)
if __name__=='__main__':main()
