#!/usr/bin/env python3
"""A different four-duty cavity: delete47, retain52, coordinate new contacts.

Exact integer deletion ledger; bounded floating-point geometric scout only.
No inherited failed six-cap bank is extended or rerun.
"""
import argparse
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, product
from fractions import Fraction as F
import json
from pathlib import Path

import numpy as np
import fresh_nonarithmetic_replace47_b_pair_20261017 as C

ROOT=C.ROOT; B=C.B; M=B.METRIC
J=np.array([[2.,1.],[1.,2.]])/3
TAG='fresh_nonarithmetic_release47_cycle_patch_20261017'
CYCLE=[7,12,52,30]
ALPHABET=list(range(2,13))+[16,24,32]
INPUT=ROOT/'results/fresh_nonarithmetic_b_release52_patch_20261017.json'


def exact_left_domain():
    # The pi/6 pair12,52 has inverse diagonal2. Any nonzero Coxeter contact
    # has magnitude>=1, hence Schur cost>=2; finite extension forces right
    # contacts to both. The other corner then reduces to A3+one.
    # A3 inverse (ordering5,6,7) is 1/4 [[3,2,1],[2,4,2],[1,2,3]].
    # Nonzero contacts have magnitude>=1; inverse entries are nonnegative.
    rows=[]
    costs={2:(F(0),F(0)),3:(F(3,4),F(0)),4:(F(3,2),F(0)),5:(F(9,8),F(3,8))}
    for m in [2,3,4,5]:rows.append(dict(labels=[2,2,m,2,2],cost=costs[m]))
    for m in [3,4,5]:rows.append(dict(labels=[m,2,2,2,2],cost=costs[m]))
    rows.append(dict(labels=[2,3,2,2,2],cost=(F(1),F(0))))
    for r in rows:
        a,b=r['cost'];assert 2-a>0 and (2-a)**2>5*b*b
        r['cost']=[str(a),str(b)]
    assert F(3,4)*3>2  # Any isolated endpoint bond m>=6 already fails.
    assert F(3,4)+F(3,4)+2*F(1,4)==2  # Both endpoints contacted.
    assert F(1)+F(3,4)+2*F(1,2)>2  # Center and endpoint contacted.
    assert 2*F(1)==2  # A bond m>=4 to the A3 center fails.
    return dict(status='PASS_EXACT_COMPLETE_LEFT_FINITE_PROFILE_DOMAIN',
        label_order=[5,6,7,12,52],profiles=rows,
        inverse_A3=[['3/4','1/2','1/4'],['1/2','1','1/2'],['1/4','1/2','3/4']],
        field='cost pair a+b*sqrt5',branch_upper_bound=16,
        I2_6_inverse_diagonal='2',
        scope='Complete necessary finite-contact profiles for one actual-S6 cap terminating7/12 and12/52, within the previously proved seven A2 profiles. Inverse-diagonal2 forces right12/52; exact A3 Schur costs exclude all high-label tails. Native old-wall geometry is screened only numerically below.')


def exact_disk():
    old=list('abcd');new=list('LPRQ')
    triangles=[list(t) for t in ['abL','bcL','cdP','daQ','cLP','dPR','dRQ','aQL','LPR','LRQ']]
    counts=Counter(tuple(sorted(e)) for t in triangles for e in combinations(t,2))
    boundary={tuple(sorted(e)) for e in ['ab','bc','cd','da']}
    assert {e for e,c in counts.items() if c==1}==boundary
    assert all(c==(1 if e in boundary else 2) for e,c in counts.items())
    for v in old+new:
        links=defaultdict(set)
        for t in triangles:
            if v in t:
                a,b=[i for i in t if i!=v];links[a].add(b);links[b].add(a)
        pending=[next(iter(links))];seen=set()
        while pending:
            a=pending.pop()
            if a in seen:continue
            seen.add(a);pending.extend(links[a]-seen)
        assert len(seen)==len(links)
        degrees=sorted(map(len,links.values()))
        assert degrees==([1,1]+[2]*(len(links)-2) if v in old else [2]*len(links))
    a=sum((i in old)!=(j in old) for i,j in counts)
    b=sum(i in new and j in new for i,j in counts)
    xyz=Counter(sum(i in new for i in t) for t in triangles)
    assert (a,b,xyz[1],xyz[2],xyz[3])==(8,5,4,4,2)
    assert len(old)+len(new)-len(counts)+len(triangles)==1
    assert 2*(a+b)-3*len(triangles)==-4
    return dict(status='PASS_EXACT_ABSTRACT_FOUR_CAP_DISK',old_cycle=CYCLE,
        old_labels=dict(zip(old,CYCLE)),new_boundary_order=['L','P','R','Q'],
        triangles=triangles,counts=dict(a=a,b=b,x=4,y=4,z=2),delta_ends=-4,
        scope='Conditional simple no-ideal disk retaining old incidences, not geometric existence.')


def setup():
    z,labels,ns,h=C.setup()
    bank=json.loads(INPUT.read_text())
    profiles=bank['profiles']
    assert len(profiles)==7 and all(len(p['offsets'])==len(labels) for p in profiles)
    assert [i for i in labels if i!=52]==bank['retained_labels']
    assert max(np.max(abs(ns[int(i)]-v)) for i,v in bank['old_projected_normals'].items())<1e-12
    assert np.max(abs(h-bank['common_projected_center']))<1e-12
    return z,labels,ns,h,profiles


def deletion_cost():
    old=json.loads((ROOT/'results/fresh_nonarithmetic_b_release_costs_20261017.json').read_text())
    vertices=[tuple(v) for v in old['finite_vertices']]
    assert len(vertices)==198
    edges=defaultdict(list); faces=defaultdict(list)
    for i,v in enumerate(vertices):
        for e in combinations(v,8): edges[e].append(i)
        for k in range(8):
            for w in combinations(v,k): faces[w].append(i)
    lost={i for i,v in enumerate(vertices) if 47 in v}; touched=[]
    for w,vs in faces.items():
        if not lost.intersection(vs) or len(vs)<10-len(w):continue
        if all(all(len(edges[tuple(j for j in vertices[i] if j!=k)])==2
                   for k in vertices[i] if k not in w) for i in vs):
            touched.append(dict(support=list(w),dimension=9-len(w),finite_vertices=vs))
    return dict(removed_wall=47,lost_finite_supports=[list(vertices[i]) for i in sorted(lost)],
                old_compact_faces_touched=touched,
                counts_by_dimension=dict(Counter(t['dimension'] for t in touched)),
                scope='Exact incidence loss of prior certificates, not geometry of a promoted replacement.')


def native(y,pi,v,pj,profiles):
    return float(y@M@v+np.array(profiles[pi]['contacts'])@J@np.array(profiles[pj]['contacts']))


def solve_frame(ids, frame, fprofiles, needed, data):
    """ids>=0 are old walls; -1 denotes the selected new seed."""
    z,labels,ns,h,profiles=data;retained=[i for i in labels if i!=47]
    old=np.array([ns[i] for i in retained]);ix=[labels.index(i) for i in retained]
    gram=frame@M@frame.T;counts=Counter();out=[];seen=set()
    if abs(np.linalg.det(gram))<1e-10:return [],Counter(singular_frame=1)
    gi=np.linalg.inv(gram); k=np.linalg.svd(frame@M)[2][-1];kq=k@M@k
    if abs(kq)<1e-10:return [],Counter(null_kernel=1)
    mss=list(product(ALPHABET,repeat=3))
    desired=-np.array([[B.theta(m) for m in row] for row in mss])
    for pi,pr in enumerate(profiles):
        cp=np.array(pr['contacts']);off=np.array(pr['offsets']);q=pr['square']
        shifts=np.array([off[labels.index(i)] if i>=0 else cp@J@np.array(profiles[pj]['contacts'])
                         for i,pj in zip(ids,fprofiles)])
        cc=desired-shifts;pp=cc@gi@frame; good=np.ones(len(mss),dtype=bool)
        for ab in needed:
            gg=gram[np.ix_(ab,ab)]
            if np.linalg.eigvalsh(gg)[0]<=1e-9:good[:]=False;break
            good &= q-np.einsum('ij,jk,ik->i',cc[:,ab],np.linalg.inv(gg),cc[:,ab])>1e-8
        rem=(q-np.einsum('ij,jk,ik->i',pp,M,pp))/kq
        indices=np.where(good&(rem>1e-9))[0]
        counts['real_finite_endpoint_branches']+=2*len(indices)
        for branch in [-1,1]:
            yy=pp[indices]+branch*np.sqrt(rem[indices])[:,None]*k
            vv=yy@M@old.T+off[ix]
            for j,i in enumerate(ids):
                if i>=0:vv[:,retained.index(i)]=desired[indices,j]
            centers=yy@M@h+pr['center_offset']
            for pos in np.where((centers<-1e-8)&(vv.max(axis=1)<2e-7))[0]:
                rel=[B.relation(v) for v in vv[pos]]
                if any(r is None or r=='parallel_or_unresolved' for r in rel):
                    counts['native_angle_gap']+=1;continue
                y=yy[pos];key=(pi,)+tuple(np.round(y,7))
                if key in seen:continue
                seen.add(key)
                assert abs(y@M@y-q)<1e-5
                out.append(dict(index=len(out),profile=pi,frame=ids,labels=list(mss[indices[pos]]),
                    branch=branch,projection=y.tolist(),center_product=float(centers[pos]),
                    finite_old_contacts=[[i,r] for i,r in zip(retained,rel) if isinstance(r,int)]))
    return out,counts


def endpoint_mask(row,data):
    z,labels,ns,h,profiles=data;y=np.array(row['projection']);q=profiles[row['profile']]['square']
    out=[]
    for e,a in enumerate(CYCLE):
        b=CYCLE[(e+1)%4];f=np.array([ns[a],ns[b]]);g=f@M@f.T;c=f@M@y
        assert np.linalg.eigvalsh(g)[0]>1e-8
        if q-c@np.linalg.inv(g)@c>1e-8:out.append(e)
    return out


def compute():
    data=setup();z,labels,ns,h,profiles=data
    baseline=B.ledger(np.array([ns[i] for i in z['actual_sides']]),z['actual_sides'])
    actual=[i for i in z['actual_sides'] if i!=47]
    deletion=B.ledger(np.array([ns[i] for i in actual]),actual)
    assert (baseline['finite'],baseline['open_ends'])==(25,15)
    assert (deletion['finite'],deletion['open_ends'])==(23,15)
    costs=deletion_cost();seeds=[];records=[];batches=[];counts=Counter()
    leftdomain=exact_left_domain();disk=exact_disk()
    for rotation in range(4):
        cyc=CYCLE[rotation:]+CYCLE[:rotation];ids=cyc[:3]
        pool,cc=solve_frame(ids,np.array([ns[i] for i in ids]),[0]*3,[(0,1),(1,2)],data)
        counts.update(cc)
        print('MULTI_ENDPOINT_FRAME',ids,'seeds',len(pool),dict(cc),flush=True)
        for rec in pool:
            rec['rotation']=rotation;rec['index']=len(seeds);rec['endpoint_duties']=endpoint_mask(rec,data)
            seeds.append(rec)
    assert len(seeds)==3 and all(s['rotation']==0 and s['profile']==0 for s in seeds)
    assert [s['labels'] for s in seeds]==[[3,2,2],[4,2,2],[5,2,2]]
    assert np.max(abs(np.array(seeds[0]['projection'])-ns[47]))<1e-7
    for seed in seeds:
        rotation=seed['rotation'];cyc=CYCLE[rotation:]+CYCLE[:rotation]
        L=np.array(seed['projection']);lp=seed['profile']
        # Remaining boundary edges cyc2-cyc3 and cyc3-cyc0.
        pp,pc=solve_frame([cyc[2],cyc[3],-1],np.array([ns[cyc[2]],ns[cyc[3]],L]),
                         [0,0,lp],[(0,1),(0,2)],data)
        qq,qc=solve_frame([cyc[3],cyc[0],-1],np.array([ns[cyc[3]],ns[cyc[0]],L]),
                         [0,0,lp],[(0,1),(1,2)],data)
        paircounts=Counter();pairs=[];vmin=None;vmax=None;four_stage=[]
        for p in pp:
            pv=np.array(p['projection'])
            for q in qq:
                qv=np.array(q['projection']);val=native(pv,p['profile'],qv,q['profile'],profiles)
                vmin=val if vmin is None else min(vmin,val);vmax=val if vmax is None else max(vmax,val)
                rel=B.relation(val);paircounts[str(rel)]+=1
                if rel=='dotted':
                    # Replace the missing P--Q edge by R, with finite P--R,
                    # R--Q and L--R contacts. R has just one old anchor.
                    rr,rc=solve_frame([cyc[3],-1,-2],np.array([ns[cyc[3]],pv,qv]),
                        [0,p['profile'],q['profile']],[(0,1),(0,2)],data)
                    rcounts=Counter();joins=[];rmin=None;rmax=None;scored=[];best_partial=None
                    for r in rr:
                        rv=np.array(r['projection']);rl=native(rv,r['profile'],L,lp,profiles)
                        rmin=rl if rmin is None else min(rmin,rl);rmax=rl if rmax is None else max(rmax,rl)
                        lr=B.relation(rl);rcounts[str(lr)]+=1
                        if isinstance(lr,int) or lr=='dotted':
                            # A failed prescribed disk can still improve another
                            # region. Score every coherent four-cap collection.
                            ledger=B.ledger(np.array([ns[i] for i in actual]+[L,pv,rv,qv]),actual+[66,67,68,69])
                            score=dict(R_index=r['index'],profile=r['profile'],labels=r['labels'],
                                branch=r['branch'],LR_relation=lr,finite=ledger['finite'],
                                exterior=ledger['exterior'],open_ends=ledger['open_ends'])
                            scored.append(score)
                            if best_partial is None or (ledger['open_ends'],ledger['exterior']) < (best_partial['boundary']['open_ends'],best_partial['boundary']['exterior']):
                                best_partial=dict(R=r,LR_relation=lr,boundary=ledger)
                            if ledger['open_ends']<baseline['open_ends']:
                                batches.append(dict(seed=seed['index'],P=p['index'],Q=q['index'],R=r,
                                    LR_relation=lr,topology='complete_boundary_improvement',boundary=ledger))
                                print('BOUNDARY_IMPROVEMENT',seed['index'],p['index'],q['index'],r['index'],ledger['open_ends'],flush=True)
                        if not isinstance(lr,int):continue
                        if any(np.linalg.eigvalsh(f@M@f.T)[0]<=1e-8
                               for f in (np.array([L,pv,rv]),np.array([L,qv,rv]))):
                            rcounts['nonfinite_new_corner']+=1;continue
                        ledger=B.ledger(np.array([ns[i] for i in actual]+[L,pv,rv,qv]),actual+[66,67,68,69])
                        row=dict(seed=seed['index'],P=p['index'],Q=q['index'],R=r,
                            LR_label=lr,LR_product=rl,topology='four_cap_disk',boundary=ledger)
                        joins.append(row);batches.append(row)
                        print('FOUR_DISK_CANDIDATE',seed['index'],p['index'],q['index'],r['index'],ledger['open_ends'],flush=True)
                    four_stage.append(dict(P=p['index'],Q=q['index'],R_counts=dict(rc),
                        R_native_accepted=len(rr),LR_relation_counts=dict(rcounts),LR_product_range=[rmin,rmax],
                        coherent_full_boundary_scores=scored,best_coherent_partial=best_partial,joins=joins))
                    print('FOUR_FRAME_DONE',p['index'],q['index'],'coherent',len(scored),
                          'best',None if best_partial is None else best_partial['boundary']['open_ends'],flush=True)
                if not isinstance(rel,int):continue
                # Both the shared old-side corner and the all-new triangle are required.
                ff=np.array([ns[cyc[3]],pv,qv]);gg=ff@M@ff.T
                if np.linalg.eigvalsh(gg)[0]<=1e-8:paircounts['nonfinite_shared_old_corner']+=1;continue
                ff=np.array([L,pv,qv]);gg=ff@M@ff.T
                if np.linalg.eigvalsh(gg)[0]<=1e-8:paircounts['nonfinite_all_new_corner']+=1;continue
                ledger=B.ledger(np.array([ns[i] for i in actual]+[L,pv,qv]),actual+[66,67,68])
                row=dict(seed=seed['index'],P=p['index'],Q=q['index'],mutual_label=rel,
                         product=val,boundary=ledger)
                pairs.append(row);batches.append(row)
                print('DISK_CANDIDATE',seed['index'],p['index'],q['index'],ledger['open_ends'],flush=True)
        records.append(dict(seed=seed['index'],P=pp,Q=qq,P_counts=dict(pc),Q_counts=dict(qc),
                            pair_counts=dict(paircounts),product_range=[vmin,vmax],candidates=pairs,
                            four_cap_bridge_stage=four_stage))
        print('SEED_DONE',seed['index'],'rotation',rotation,'P',len(pp),'Q',len(qq),'closed',len(pairs),flush=True)
    batches.sort(key=lambda r:(r['boundary']['open_ends'],r['boundary']['exterior'],-r['boundary']['finite']))
    stages=[r for p in records for r in p['four_cap_bridge_stage']]
    scores=[s for r in stages for s in r['coherent_full_boundary_scores']]
    scorehist=Counter((s['finite'],s['exterior'],s['open_ends']) for s in scores)
    assert len(stages)==62 and sum(r['R_native_accepted'] for r in stages)==8473
    assert len(scores)==8113 and scorehist==Counter({(31,7,15):8113})
    assert not batches
    paths=[Path(__file__),Path(C.__file__),Path(B.__file__),INPUT,C.PARENT,
           ROOT/'results/fresh_nonarithmetic_replace47_close12_20261017.json',
           ROOT/'results/fresh_nonarithmetic_b_release_costs_20261017.json']
    return dict(status='BOUNDED_NUMERICAL_RELEASE47_FOUR_DUTY_DISK_SCOUT',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        deleted_wall=47,boundary_cycle=CYCLE,alphabet=ALPHABET,profiles=profiles,
        retained_labels=[i for i in labels if i!=47],baseline=baseline,after_deletion=deletion,
        deletion_cost=costs,exact_left_profile_domain=leftdomain,abstract_disk=disk,
        seed_counts=dict(counts),seeds=seeds,patches=records,batches=batches,
        numerical_summary=dict(coherent_PQ_frames=62,R_native_passes=8473,
            coherent_collections_scored=8113,complete_boundary_histogram=[dict(finite=f,exterior=e,open_ends=i,count=v) for (f,e,i),v in sorted(scorehist.items())],
            improving_collections=0),
        scope='Changed cavity with52 retained. First a three-cap triangulated disk is tested. For every coherent dotted P/Q pair, an additional R is then framed by their shared old side and P/Q; finite L/R plus two all-new finite corners would complete a four-cap disk without demanding P/Q meet. All seven A2 profiles and listed finite label values are tested. Every retained native contact and common old center are checked numerically. No all-label exclusion, exact promoted chamber or finite-volume conclusion follows from this scout.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=json.loads(json.dumps(compute()));p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert json.loads(p.read_text())==data
    else:p.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],'seeds',len(data['seeds']),'disks',len(data['batches']),flush=True)
