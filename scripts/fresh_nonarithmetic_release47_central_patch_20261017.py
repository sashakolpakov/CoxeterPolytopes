#!/usr/bin/env python3
"""Close the four-new-side ring with a central cap having no old anchor.

Exact abstract disk, bounded numerical geometry; predecessors remain frozen.
"""
import argparse
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import numpy as np
import fresh_nonarithmetic_release47_replace12_frontier_20261017 as F

A=F.A;ROOT=F.ROOT;B=F.B;M=F.M
TAG='fresh_nonarithmetic_release47_central_patch_20261017'
INPUT=ROOT/'results/fresh_nonarithmetic_release47_nonzero_a2_bridge_20261017.json'
FRONTIER=ROOT/'results'/f'{F.TAG}.json'


def disk():
    old=list('abcd');new=list('LPRQS')
    triangles=[list(t) for t in ['abL','bcL','cdP','daQ','cLP','dPR','dRQ','aQL','SLP','SPR','SRQ','SQL']]
    counts=Counter(tuple(sorted(e)) for t in triangles for e in combinations(t,2))
    boundary={tuple(sorted(e)) for e in ['ab','bc','cd','da']}
    assert {e for e,c in counts.items() if c==1}==boundary
    assert all(c==(1 if e in boundary else 2) for e,c in counts.items())
    for v in old+new:
        links=defaultdict(set)
        for tri in triangles:
            if v in tri:
                i,j=[r for r in tri if r!=v];links[i].add(j);links[j].add(i)
        seen=set();pending=[next(iter(links))]
        while pending:
            i=pending.pop()
            if i not in seen:seen.add(i);pending.extend(links[i]-seen)
        assert len(seen)==len(links)
        deg=sorted(map(len,links.values()))
        assert deg==([1,1]+[2]*(len(links)-2) if v in old else [2]*len(links))
    a=sum((i in old)!=(j in old) for i,j in counts)
    b=sum(i in new and j in new for i,j in counts)
    xyz=Counter(sum(i in new for i in tri) for tri in triangles)
    assert (a,b,xyz[1],xyz[2],xyz[3])==(8,8,4,4,4)
    assert len(old)+len(new)-len(counts)+len(triangles)==1
    assert 2*(a+b)-3*len(triangles)==-4
    return dict(status='PASS_EXACT_ABSTRACT_FIVE_CAP_CENTRAL_DISK',old_cycle=A.CYCLE,
        triangles=triangles,counts=dict(a=a,b=b,x=4,y=4,z=4),delta_ends=-4,
        central_cap='S: no old anchor, four new neighbors L,P,R,Q',
        scope='Conditional simple/no-ideal disk with unchanged genuine old boundary incidences; not a geometric existence assertion.')


def isolated_high_contacts(ring,vec,data):
    """All possible >=6 central bonds, using exact finite-corner isolation.

    A high bond with any other incident nonright bond has a nonpositive
    3x3 principal determinant. Thus its two ends must be isolated in
    each required finite corner. The saved ring has at most one high
    central bond: eligible opposite vertices never both occur here.
    Geometry remains floating-point; this is an exact label reduction.
    """
    z,labels,ns,h,profiles=data
    edge_labels=[ring[1]['labels'][2],ring[2]['labels'][1],ring[2]['labels'][2],ring[3]['labels'][2]]
    eligible=[i for i in range(4) if ring[i]['profile']==0 and edge_labels[i-1]==edge_labels[i]==2]
    assert not any(i in eligible and (i+2)%4 in eligible for i in range(4))
    retained=[i for i in labels if i!=47];old=np.array([ns[i] for i in retained])
    rows=[];hits=[];counts=Counter()
    for high in eligible:
        frame_ids=[i for i in range(4) if i!=high]
        for other_m in [2,3,4,5]:
            desired=[-B.theta(other_m) if i==(high+2)%4 else 0. for i in frame_ids]
            for branch,s in A.C.R.roots(np.array([vec[i] for i in frame_ids]),desired):
                val=float(s@M@vec[high]);relation=B.relation(val)
                row=dict(high_ring_vertex=high,opposite_label=other_m,branch=branch,
                    projection=s.tolist(),high_product=val,numerical_relation=relation,
                    center_product=float(s@M@h))
                if not (-2+2e-7 < val <= -np.sqrt(3)+2e-7):status='outside_finite_high_label_range'
                elif not all(np.linalg.eigvalsh(f@M@f.T)[0]>1e-8 for f in
                    [np.array([vec[i],vec[(i+1)%4],s]) for i in range(4)]):status='nonfinite_corner'
                elif s@M@h>=-1e-8:status='wrong_common_center'
                else:
                    values=old@M@s;relations=[B.relation(v) for v in values]
                    bad=[dict(wall=i,product=float(v),numerical_relation=t) for i,v,t in zip(retained,values,relations)
                         if t is None or t=='parallel_or_unresolved']
                    row['old_native_blockers']=bad
                    if bad:status='old_native_blocker'
                    elif not isinstance(relation,int) or relation<6:status='high_product_gap'
                    else:status='candidate';hits.append(row)
                row['status']=status;rows.append(row);counts[status]+=1
    return dict(ring_edge_labels=edge_labels,eligible_high_vertices=eligible,
        at_most_one_high_bond=True,rows=rows,counts=dict(counts),candidates=hits,
        scope='Complete high-label profile reduction for these literal ring labels; real roots and rejecting products are numerical.')


def compute():
    raw=json.loads(INPUT.read_text());front=json.loads(FRONTIER.read_text())
    data=A.setup();z,labels,ns,h,profiles=data
    actual=[i for i in z['actual_sides'] if i!=47]
    stages=[];hits=[];hist=Counter();allcounts=Counter()
    for index,st in enumerate(raw['stages']):
        if st['best_coherent'] is None:continue
        group=next(g for g in front['changed_seed_pools']
            if front['replacement_trials'][g['replacement']]['label_to52']==st['replacement_label'])
        replacement=front['replacement_trials'][group['replacement']]
        newns={**ns,12:np.array(replacement['projection'])};newlabels=[i for i in labels if i!=53]
        newprofiles=[{**p,'offsets':[p['offsets'][labels.index(i)] for i in newlabels]} for p in profiles]
        newdata=(z,newlabels,newns,h,newprofiles)
        L=group['seeds'][st['seed']];patch=group['patches'][st['seed']]
        P=patch['P'][st['P']];Q=patch['Q'][st['Q']];R=st['best_coherent']['R']
        ring=[L,P,R,Q];vec=[np.array(r['projection']) for r in ring]
        cap,cc=A.solve_frame([-1,-2,-3],np.array(vec[:3]),
            [r['profile'] for r in ring[:3]],[(0,1),(1,2)],newdata)
        high=isolated_high_contacts(ring,vec,newdata)
        counts=Counter();scored=[];best=None;closed=[]
        for s in cap:
            sv=np.array(s['projection']);value=A.native(sv,s['profile'],vec[3],Q['profile'],newprofiles)
            rel=B.relation(value);counts[str(rel)]+=1
            if not isinstance(rel,int) and rel!='dotted':continue
            boundary=B.ledger(np.array([newns[i] for i in actual]+vec+[sv]),actual+[66,67,68,69,70])
            row=dict(S=s['index'],Q_relation=rel,Q_product=value,
                finite=boundary['finite'],exterior=boundary['exterior'],open_ends=boundary['open_ends'])
            scored.append(row);hist[(row['finite'],row['exterior'],row['open_ends'])]+=1
            if best is None or (row['open_ends'],row['exterior'])<(best['boundary']['open_ends'],best['boundary']['exterior']):
                best=dict(S=s,Q_relation=rel,boundary=boundary)
            finite_corners=all(np.linalg.eigvalsh(frame@M@frame.T)[0]>1e-8 for frame in
                [np.array([vec[2],vec[3],sv]),np.array([vec[3],vec[0],sv])])
            if isinstance(rel,int) and finite_corners:closed.append(dict(S=s,Q_relation=rel,boundary=boundary))
            if row['open_ends']<15:
                hits.append(dict(stage=index,replacement=replacement,L=L,P=P,R=R,Q=Q,S=s,
                    Q_relation=rel,boundary=boundary,all_four_new_corners_finite=finite_corners))
                print('CENTRAL_IMPROVEMENT',index,s['index'],rel,row['open_ends'],flush=True)
        allcounts.update(counts)
        stages.append(dict(parent_stage=index,replacement_label=st['replacement_label'],
            seed=st['seed'],P=st['P'],Q=st['Q'],R=R,
            candidate_counts=dict(cc),S_native_accepted=len(cap),Q_counts=dict(counts),
            coherent_boundary_scores=scored,best=best,complete_disks=closed,
            isolated_high_label_screen=high))
        for row in high['candidates']:
            sv=np.array(row['projection'])
            boundary=B.ledger(np.array([newns[i] for i in actual]+vec+[sv]),actual+[66,67,68,69,70])
            row['boundary']=boundary
            if boundary['open_ends']<15:
                hits.append(dict(stage=index,high_label_candidate=row,replacement=replacement,
                    L=L,P=P,R=R,Q=Q,boundary=boundary))
        print('CENTRAL_DONE',index,'S',len(cap),'Q',dict(counts),'closed',len(closed),
            'best',None if best is None else best['boundary']['open_ends'],'high',high['counts'],flush=True)
    assert len(stages)==20
    paths=[Path(__file__),Path(F.__file__),Path(A.__file__),INPUT,FRONTIER,A.INPUT,A.C.PARENT]
    return dict(status='BOUNDED_NUMERICAL_ZERO_OLD_ANCHOR_CENTRAL_PATCH',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        abstract_disk=disk(),alphabet=A.ALPHABET,profiles=profiles,
        selection='For each of the20 parent stages with coherent R, take its saved deterministic best_coherent representative; no old-angle cutoff increase.',
        exact_high_label_reduction='If a>=sqrt3 is a finite high bond and b>=1 any other incident nonright bond, a 3x3 norm2 principal determinant is at most8-2(a^2+b^2)<=0. Therefore both high-bond endpoints are isolated in each finite corner. In each saved ring the eligible high vertices contain no opposite pair, so at most one central bond can be high. Its other three contacts are two right and one label2..5, which determines all norm branches without sampling the high label.',
        stages=stages,Q_counts=dict(allcounts),
        full_boundary_histogram=[dict(finite=f,exterior=e,open_ends=n,count=c) for (f,e,n),c in sorted(hist.items())],
        hits=hits,scope='Central cap is framed solely on three NEW ring walls, with no old finite anchor. All seven A2 profiles, listed finite labels, both norm branches, every retained native gate, center and fourth-ring contact are screened. Every coherent five-cap collection is boundary-scored. Numerical discovery only; no all-parameter exclusion or geometry promotion.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(compute()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['Q_counts'],out['full_boundary_histogram'],'hits',len(out['hits']),flush=True)
