#!/usr/bin/env python3
"""Change the rigid12--52 pi/6 bond, then solve the changed cavity.

Bounded numerical construction diagnostic; all frozen predecessors read-only.
"""
import argparse
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations
from fractions import Fraction as F
import json
from pathlib import Path

import numpy as np
import fresh_nonarithmetic_release47_cycle_patch_20261017 as A

ROOT=A.ROOT;B=A.B;M=A.M;J=A.J
TAG='fresh_nonarithmetic_release47_replace12_frontier_20261017'
REPLACEMENT_LABELS=[3,4,5]


def cost():
    data=json.loads((ROOT/'results/fresh_nonarithmetic_b_release_costs_20261017.json').read_text())
    vs=[tuple(v) for v in data['finite_vertices']];edges=defaultdict(list);faces=defaultdict(list)
    for i,v in enumerate(vs):
        for e in combinations(v,8):edges[e].append(i)
        for k in range(8):
            for w in combinations(v,k):faces[w].append(i)
    assert all(53 not in v for v in vs)
    lost={i for i,v in enumerate(vs) if 12 in v or 47 in v};affected=[]
    for w,inds in faces.items():
        if not lost.intersection(inds) or len(inds)<10-len(w):continue
        if all(all(len(edges[tuple(j for j in vs[i] if j!=k)])==2 for k in vs[i] if k not in w) for i in inds):
            affected.append(dict(support=list(w),dimension=9-len(w),finite_vertices=inds))
    return dict(changed_or_deleted=[12,47,53],old53_finite_vertices=0,
        lost_old_finite_supports=[list(vs[i]) for i in sorted(lost)],
        old_compact_faces_touched=affected,counts_by_dimension=dict(Counter(x['dimension'] for x in affected)),
        scope='Exact loss of old certificates before checking a new geometry; an altered8-face can be rebuilt.')


def replacement_roots(data):
    z,labels,ns,h,profiles=data;frame=np.array([ns[i] for i in [7,8,52]])
    retained=[i for i in labels if i not in [12,47]];old=np.array([ns[i] for i in retained])
    trials=[];accepted=[]
    for m in REPLACEMENT_LABELS:
        for branch,r in A.C.R.roots(frame,[0.,0.,-B.theta(m)]):
            vals=old@M@r
            for i,val in [(7,0.),(8,0.),(52,-B.theta(m))]:vals[retained.index(i)]=val
            rel=[B.relation(v) for v in vals]
            blockers=[dict(wall=i,value=float(v),kind='obtuse' if v>0 else 'gap_or_nearparallel')
                for i,v,t in zip(retained,vals,rel) if t is None or t=='parallel_or_unresolved']
            row=dict(index=len(trials),label_to52=m,branch=branch,projection=r.tolist(),
                center_product=float(r@M@h),blockers=blockers,
                finite_old_contacts=[[i,t] for i,t in zip(retained,rel) if isinstance(t,int)])
            if all(b['wall']==53 for b in blockers) and r@M@h<-1e-8:
                row['status']='native_coherent_after53_release'
                row['additional_deleted_wall']=53
                newns={**ns,12:r};actual=[i for i in z['actual_sides'] if i!=47]
                row['boundary']=B.ledger(np.array([newns[i] for i in actual]),actual)
                accepted.append(row['index'])
            else:row['status']='rejected'
            trials.append(row)
    return trials,accepted


def three_cap_patches(seed,data):
    z,labels,ns,h,profiles=data
    rot=seed['rotation'];cyc=A.CYCLE[rot:]+A.CYCLE[:rot]
    L=np.array(seed['projection']);lp=seed['profile'];actual=[i for i in z['actual_sides'] if i!=47]
    pp,pc=A.solve_frame([cyc[2],cyc[3],-1],np.array([ns[cyc[2]],ns[cyc[3]],L]),
        [0,0,lp],[(0,1),(0,2)],data)
    qq,qc=A.solve_frame([cyc[3],cyc[0],-1],np.array([ns[cyc[3]],ns[cyc[0]],L]),
        [0,0,lp],[(0,1),(1,2)],data)
    counts=Counter();scores=[];improved=[]
    for p in pp:
        pv=np.array(p['projection'])
        for q in qq:
            qv=np.array(q['projection']);value=A.native(pv,p['profile'],qv,q['profile'],profiles)
            rel=B.relation(value);counts[str(rel)]+=1
            if not isinstance(rel,int) and rel!='dotted':continue
            boundary=B.ledger(np.array([ns[i] for i in actual]+[L,pv,qv]),actual+[66,67,68])
            row=dict(P=p['index'],Q=q['index'],mutual_relation=rel,product=value,
                finite=boundary['finite'],exterior=boundary['exterior'],open_ends=boundary['open_ends'])
            scores.append(row)
            if boundary['open_ends']<15:
                row={**row,'boundary':boundary};improved.append(row)
                print('CHANGED_PATCH_IMPROVEMENT',seed['index'],p['index'],q['index'],boundary['open_ends'],flush=True)
    print('CHANGED_PATCH_DONE',seed['index'],'P',len(pp),'Q',len(qq),'coherent',len(scores),'improved',len(improved),flush=True)
    return dict(seed=seed['index'],P=pp,Q=qq,P_counts=dict(pc),Q_counts=dict(qc),
        pair_counts=dict(counts),complete_boundary_scores=scores,improved=improved)


def compute():
    data=A.setup();z,labels,ns,h,profiles=data
    assert all(abs(p['offsets'][labels.index(12)])<1e-12 for p in profiles)
    redundant=next(r for r in z['projections'] if r['wall']==53)
    assert redundant['kind']=='future_timelike_redundant_on_H3'
    assert F(redundant['square'][1])<0 and F(redundant['projected_center_product'][1])<0
    trials,accepted=replacement_roots(data);records=[]
    for ind in accepted:
        replacement=trials[ind];newns={**ns,12:np.array(replacement['projection'])}
        newlabels=[i for i in labels if i!=53]
        newprofiles=[{**p,'offsets':[p['offsets'][labels.index(i)] for i in newlabels]} for p in profiles]
        newdata=(z,newlabels,newns,h,newprofiles);seeds=[];counts=Counter()
        for rot in range(4):
            cyc=A.CYCLE[rot:]+A.CYCLE[:rot];ids=cyc[:3]
            out,cc=A.solve_frame(ids,np.array([newns[i] for i in ids]),[0]*3,[(0,1),(1,2)],newdata)
            counts.update(cc)
            for row in out:
                row['index']=len(seeds);row['rotation']=rot;row['endpoint_duties']=A.endpoint_mask(row,newdata)
                seeds.append(row)
            print('CHANGED_SEED_FRAME',replacement['label_to52'],ids,len(out),flush=True)
        patches=[three_cap_patches(seed,newdata) for seed in seeds]
        records.append(dict(replacement=ind,seeds=seeds,seed_counts=dict(counts),patches=patches))
    paths=[Path(__file__),Path(A.__file__),ROOT/'results'/f'{A.TAG}.json',A.INPUT,A.C.PARENT,
           ROOT/'results/fresh_nonarithmetic_b_release_costs_20261017.json']
    return dict(status='BOUNDED_NUMERICAL_CHANGED_PI6_CAVITY_FRONTIER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        replacement_recipe=dict(frame=[7,8,52],norm=2,contacts=[2,2,'m'],m_values=REPLACEMENT_LABELS,
            white_profile='all right toS6',both_norm_branches=True),
        changed_or_deleted_old_walls=[12,47,53],deletion_cost=cost(),profiles=profiles,
        exact_inherited_redundant53_projection=redundant,
        old_control=dict(walls=64,finite_vertices=198,S6_open_ends=15),
        replacement_trials=trials,accepted_replacements=accepted,changed_seed_pools=records,
        scope='Only12-to52 labels3,4,5, retaining12-right7/8, both norm branches. Each viable branch is blocked only by53; deleting that zero-finite-vertex, timelike-S6-redundant wall leaves61 unchanged native contacts plus the new12. Subsequent endpoint profiles are recomputed in changed geometry, not imported from the rigid pi6 obstruction. Numerical only; no exact geometry or finite-volume claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=json.loads(json.dumps(compute()));p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert json.loads(p.read_text())==data
    else:p.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],'accepted',len(data['accepted_replacements']),
          'seed_counts',[len(r['seeds']) for r in data['changed_seed_pools']],flush=True)
