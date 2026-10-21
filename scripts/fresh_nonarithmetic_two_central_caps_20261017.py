#!/usr/bin/env python3
"""Two-center disk scout, changing topology rather than an angle cutoff.

Uses both interior-edge diagonals in the twenty frozen literal rings.
Exact abstract topology; bounded numerical geometry, no promotion.
"""
import argparse
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import numpy as np
import fresh_nonarithmetic_release47_central_patch_20261017 as C

F,A,B,M=C.F,C.A,C.B,C.M
ROOT=C.ROOT
TAG='fresh_nonarithmetic_two_central_caps_20261017'
INPUT=ROOT/'results'/f'{C.TAG}.json'


def disk():
    old=list('abcd');new=list('LPRQST')
    triangles=[list(t) for t in ['abL','bcL','cdP','daQ','cLP','dPR','dRQ','aQL',
                                'SLP','SPR','TLQ','TQR','SLT','SRT']]
    edges=Counter(tuple(sorted(e)) for t in triangles for e in combinations(t,2))
    boundary={tuple(sorted(e)) for e in ['ab','bc','cd','da']}
    assert {e for e,k in edges.items() if k==1}==boundary
    assert all(k==(1 if e in boundary else 2) for e,k in edges.items())
    for v in old+new:
        link=defaultdict(set)
        for tri in triangles:
            if v in tri:
                a,b=[i for i in tri if i!=v];link[a].add(b);link[b].add(a)
        seen=set();todo=[next(iter(link))]
        while todo:
            u=todo.pop()
            if u not in seen:seen.add(u);todo.extend(link[u]-seen)
        assert set(link)==seen
        assert sorted(map(len,link.values()))==([1,1]+[2]*(len(link)-2) if v in old else [2]*len(link))
    a=sum((i in old)!=(j in old) for i,j in edges)
    b=sum(i in new and j in new for i,j in edges)
    xyz=Counter(sum(i in new for i in tri) for tri in triangles)
    assert (a,b,xyz[1],xyz[2],xyz[3])==(8,11,4,4,6)
    assert len(old+new)-len(edges)+len(triangles)==1
    assert 2*(a+b)-3*len(triangles)==-4
    return dict(status='PASS_EXACT_ABSTRACT_TWO_CENTER_DISK',triangles=triangles,
        old_cycle=A.CYCLE,ring_order=['L','P','R','Q'],counts=dict(a=a,b=b,x=4,y=4,z=6),
        boundary_vertices=4,new_vertices=6,edges=23,triangles_count=14,delta_ends=-4,
        scope='Conditional retained genuine old-boundary/simple/no-ideal disk. No claim of geometric existence.')


def coherent_pool(frame_indices,ring,vec,newdata):
    frame=np.array([vec[i] for i in frame_indices])
    pool,counts=A.solve_frame([-1,-2,-3],frame,
        [ring[i]['profile'] for i in frame_indices],[(0,1),(1,2)],newdata)
    omitted=next(i for i in range(4) if i not in frame_indices)
    accepted=[];cross=Counter()
    for row in pool:
        val=A.native(np.array(row['projection']),row['profile'],vec[omitted],ring[omitted]['profile'],newdata[-1])
        rel=B.relation(val);cross[str(rel)]+=1
        if isinstance(rel,int) or rel=='dotted':
            accepted.append(dict(**row,omitted_ring=omitted,omitted_relation=rel,omitted_product=val))
    return accepted,dict(counts),dict(cross)


def run():
    saved=json.loads(INPUT.read_text());front=json.loads(C.FRONTIER.read_text())
    parents=json.loads(C.INPUT.read_text())
    data=A.setup();z,labels,ns,h,profiles=data
    actual=[i for i in z['actual_sides'] if i!=47]
    baseline=B.ledger(np.array([ns[i] for i in z['actual_sides']]),z['actual_sides'])
    assert (baseline['finite'],baseline['exterior'],baseline['open_ends'])==(25,7,15)
    stages=[];allcounts=Counter();hist=Counter();improved=[];complete=[];relaxed_hist=Counter();relaxed_best=None
    for saved_stage in saved['stages']:
        parent_index=saved_stage['parent_stage'];st=parents['stages'][parent_index]
        group=next(g for g in front['changed_seed_pools']
            if front['replacement_trials'][g['replacement']]['label_to52']==st['replacement_label'])
        replacement=front['replacement_trials'][group['replacement']]
        nmap={**ns,12:np.array(replacement['projection'])};newlabels=[i for i in labels if i!=53]
        newprofiles=[{**p,'offsets':[p['offsets'][labels.index(i)] for i in newlabels]} for p in profiles]
        newdata=(z,newlabels,nmap,h,newprofiles)
        L=group['seeds'][st['seed']];patch=group['patches'][st['seed']]
        P=patch['P'][st['P']];Q=patch['Q'][st['Q']];R=st['best_coherent']['R']
        ring=[L,P,R,Q];vec=[np.array(r['projection']) for r in ring]
        oldnormals=np.array([nmap[i] for i in actual]+vec);alllabels=actual+[66,67,68,69]
        base=B.ledger(oldnormals,alllabels)
        assert (base['finite'],base['exterior'],base['open_ends'])==(31,7,15)
        for diagonal in [0,1]:
            ix=[(diagonal+i)%4 for i in range(4)]
            si=[ix[0],ix[1],ix[2]];ti=[ix[0],ix[3],ix[2]]
            sp,sc,sx=coherent_pool(si,ring,vec,newdata)
            tp,tc,tx=coherent_pool(ti,ring,vec,newdata)
            print('TWO_CENTER_POOLS',parent_index,diagonal,len(sp),len(tp),flush=True)
            counts=Counter();scores=[];joins=[];product_range=None;closest=None;best=None;relaxed=[]
            for s in sp:
                sv=np.array(s['projection'])
                for t in tp:
                    tv=np.array(t['projection']);val=A.native(sv,s['profile'],tv,t['profile'],newprofiles)
                    rel=B.relation(val);counts[str(rel)]+=1
                    if product_range is None:product_range=[val,val]
                    else:product_range=[min(product_range[0],val),max(product_range[1],val)]
                    if -2<val<0:
                        continuous_label=np.pi/np.arccos(-val/2)
                        bracket={max(2,int(np.floor(continuous_label))),max(2,int(np.ceil(continuous_label)))}
                        near=min(sorted(bracket),key=lambda m:abs(val+B.theta(m)))
                        gap=abs(val+B.theta(near))
                        if closest is None or gap<closest['absolute_gap']:
                            closest=dict(S=s['index'],T=t['index'],product=val,label=near,absolute_gap=gap)
                        eigs=[float(np.linalg.eigvalsh(g@M@g.T)[0]) for g in
                            [np.array([vec[i],sv,tv]) for i in [ix[0],ix[2]]]]
                        if min(eigs)>1e-8:
                            rb=B.ledger(np.vstack([oldnormals,sv,tv]),alllabels+[70,71])
                            rr=dict(S=s['index'],T=t['index'],mutual_product=val,
                                nearest_label=near,angle_gap=gap,central_min_eigenvalues=eigs,
                                finite=rb['finite'],exterior=rb['exterior'],open_ends=rb['open_ends'],
                                near_ideal_rays=rb['rays']-rb['finite']-rb['exterior'],
                                degenerate_faces=rb['edges'].get('degenerate_face',0))
                            relaxed.append(rr);relaxed_hist[rb['finite'],rb['exterior'],rb['open_ends']]+=1
                            if relaxed_best is None or gap<relaxed_best['angle_gap']:
                                relaxed_best=dict(parent_stage=parent_index,diagonal=diagonal,
                                    replacement=replacement,ring=ring,S=s,T=t,angle_gap=gap,
                                    nearest_label=near,mutual_product=val,central_min_eigenvalues=eigs,boundary=rb)
                    if not isinstance(rel,int) and rel!='dotted':continue
                    if np.max(abs(sv-tv))<1e-7:
                        counts['duplicate_normal']+=1;continue
                    boundary=B.ledger(np.vstack([oldnormals,sv,tv]),alllabels+[70,71])
                    score=dict(S=s['index'],T=t['index'],mutual_relation=rel,product=val,
                        finite=boundary['finite'],exterior=boundary['exterior'],open_ends=boundary['open_ends'],
                        near_ideal_rays=boundary['rays']-boundary['finite']-boundary['exterior'],
                        degenerate_faces=boundary['edges'].get('degenerate_face',0))
                    scores.append(score);hist[(score['finite'],score['exterior'],score['open_ends'])]+=1
                    key=(score['open_ends'],score['exterior'],-score['finite'])
                    if best is None or key<best['key']:
                        best=dict(key=key,S=s,T=t,mutual_relation=rel,boundary=boundary)
                    if isinstance(rel,int):
                        eigs=[float(np.linalg.eigvalsh(g@M@g.T)[0]) for g in
                            [np.array([vec[i],sv,tv]) for i in [ix[0],ix[2]]]]
                        score['central_endpoint_min_eigenvalues']=eigs
                        if min(eigs)>1e-8:
                            hit=dict(parent_stage=parent_index,diagonal=diagonal,S=s,T=t,
                                mutual_relation=rel,boundary=boundary)
                            joins.append(hit);complete.append(hit)
                            print('TWO_CENTER_DISK',parent_index,diagonal,s['index'],t['index'],rel,boundary['open_ends'],flush=True)
                        else:counts['nonfinite_central_endpoint']+=1
                    if score['open_ends']<15:
                        improved.append(dict(parent_stage=parent_index,diagonal=diagonal,replacement=replacement,
                            ring=ring,S=s,T=t,mutual_relation=rel,boundary=boundary))
                        print('TWO_CENTER_IMPROVEMENT',parent_index,diagonal,s['index'],t['index'],rel,boundary['open_ends'],flush=True)
            allcounts.update(counts)
            stages.append(dict(parent_stage=parent_index,diagonal=diagonal,
                replacement_label=st['replacement_label'],ring=ring,S_frame=si,T_frame=ti,
                S_native_counts=sc,T_native_counts=tc,S_fourth_ring_counts=sx,T_fourth_ring_counts=tx,
                S=sp,T=tp,pair_counts=dict(counts),mutual_product_range=product_range,
                closest_finite_gap=closest,complete_boundary_scores=scores,best=best,complete_disks=joins))
            stages[-1]['continuous_angle_relaxations']=relaxed
            print('TWO_CENTER_DONE',parent_index,diagonal,dict(counts),'scores',len(scores),
                'best',None if best is None else best['boundary']['open_ends'],flush=True)
    assert len(stages)==40
    assert allcounts==Counter({'None':8781,'dotted':4215})
    assert hist==Counter({(35,7,15):4215}) and relaxed_hist==Counter({(37,5,11):316})
    assert not improved and not complete
    assert all(not r['near_ideal_rays'] and not r['degenerate_faces']
               for s in stages for r in s['complete_boundary_scores']+s['continuous_angle_relaxations'])
    paths=[Path(__file__),Path(C.__file__),Path(A.__file__),Path(F.__file__),INPUT,C.INPUT,C.FRONTIER,A.INPUT,A.C.PARENT]
    return dict(status='BOUNDED_NUMERICAL_TWO_CENTER_RING_DISK_SCOUT',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        abstract_disk=disk(),baseline=baseline,control=dict(walls=64,finite_vertices=198,S6_open_ends=15),
        changed_or_deleted_old_walls=[12,47,53],inherited_deletion_cost=front['deletion_cost'],
        alphabet=A.ALPHABET,profiles=profiles,selection='Both interior-edge diagonals in all20 frozen literal ring representatives.',
        stages=stages,pair_counts=dict(allcounts),
        boundary_histogram=[dict(finite=f,exterior=e,open_ends=n,count=c) for (f,e,n),c in sorted(hist.items())],
        continuous_angle_relaxation_histogram=[dict(finite=f,exterior=e,open_ends=n,count=c) for (f,e,n),c in sorted(relaxed_hist.items())],
        closest_finite_endpoint_relaxation=relaxed_best,
        complete_disks=complete,improvements=improved,
        scope='Exact abstract two-center topology, bounded numerical search with all7 A2profiles and unchanged finite alphabet. Both real norm branches, all retained native wall gates, all ring contacts, common center and every coherent pair full S6 ledger are checked. Separate continuous-angle relaxations keep all old Coxeter gates but permit a non-Coxeter S/T angle, explicitly not accepted. Finite/ideal/exterior tolerance results are discovery only. No all-label two-cap obstruction or new chamber/volume claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['pair_counts'],d['boundary_histogram'],'improvements',len(d['improvements']),flush=True)
