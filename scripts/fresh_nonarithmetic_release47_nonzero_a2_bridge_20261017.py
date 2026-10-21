#!/usr/bin/env python3
"""One bounded fourth-cap bridge in the changed12, nonzero-A2 family.

Numerical discovery only. Every coherent four-cap collection is scored by
the complete intrinsic S6 boundary; no geometry is accepted from angles alone.
"""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
import fresh_nonarithmetic_release47_replace12_frontier_20261017 as F

A=F.A; ROOT=F.ROOT; B=F.B; M=F.M
TAG='fresh_nonarithmetic_release47_nonzero_a2_bridge_20261017'
INPUT=ROOT/'results'/f'{F.TAG}.json'


def compute():
    parent=json.loads(INPUT.read_text())
    assert parent['status']=='BOUNDED_NUMERICAL_CHANGED_PI6_CAVITY_FRONTIER'
    data=A.setup();z,labels,ns,h,profiles=data
    actual=[i for i in z['actual_sides'] if i!=47]
    stages=[];improved=[];total_counts=Counter();hist=Counter()
    for group in parent['changed_seed_pools']:
        replacement=parent['replacement_trials'][group['replacement']]
        m=replacement['label_to52'];newns={**ns,12:np.array(replacement['projection'])}
        newlabels=[i for i in labels if i!=53]
        newprofiles=[{**p,'offsets':[p['offsets'][labels.index(i)] for i in newlabels]} for p in profiles]
        newdata=(z,newlabels,newns,h,newprofiles)
        for seed,patch in zip(group['seeds'],group['patches']):
            assert seed['index']==patch['seed']
            if seed['profile']==0:continue
            rot=seed['rotation'];cyc=A.CYCLE[rot:]+A.CYCLE[:rot]
            L=np.array(seed['projection']);lp=seed['profile']
            for pair in patch['complete_boundary_scores']:
                assert pair['mutual_relation']=='dotted'
                p=patch['P'][pair['P']];q=patch['Q'][pair['Q']]
                pv=np.array(p['projection']);qv=np.array(q['projection'])
                rr,rc=A.solve_frame([cyc[3],-1,-2],np.array([newns[cyc[3]],pv,qv]),
                    [0,p['profile'],q['profile']],[(0,1),(0,2)],newdata)
                counts=Counter();scores=[];best=None;joins=[]
                for r in rr:
                    rv=np.array(r['projection'])
                    value=A.native(rv,r['profile'],L,lp,newprofiles);rel=B.relation(value)
                    counts[str(rel)]+=1
                    if not isinstance(rel,int) and rel!='dotted':continue
                    boundary=B.ledger(np.array([newns[i] for i in actual]+[L,pv,rv,qv]),actual+[66,67,68,69])
                    score=dict(R_index=r['index'],profile=r['profile'],labels=r['labels'],
                        branch=r['branch'],LR_product=value,LR_relation=rel,
                        finite=boundary['finite'],exterior=boundary['exterior'],open_ends=boundary['open_ends'])
                    scores.append(score)
                    hist[(score['finite'],score['exterior'],score['open_ends'])]+=1
                    if best is None or (score['open_ends'],score['exterior'])<(best['boundary']['open_ends'],best['boundary']['exterior']):
                        best=dict(R=r,LR_relation=rel,boundary=boundary)
                    if boundary['open_ends']<15:
                        hit=dict(replacement_label=m,seed=seed,p=p,q=q,R=r,LR_relation=rel,boundary=boundary)
                        improved.append(hit)
                        print('IMPROVED',m,seed['index'],p['index'],q['index'],r['index'],boundary['open_ends'],flush=True)
                    if isinstance(rel,int):
                        finite_new_corners=all(np.linalg.eigvalsh(frame@M@frame.T)[0]>1e-8
                            for frame in [np.array([L,pv,rv]),np.array([L,qv,rv])])
                        if finite_new_corners:joins.append(dict(R=r,LR_relation=rel,boundary=boundary))
                        else:counts['nonfinite_prescribed_new_corner']+=1
                total_counts.update(counts)
                stages.append(dict(replacement_label=m,seed=seed['index'],seed_profile=lp,
                    P=p['index'],Q=q['index'],frame=[cyc[3],-1,-2],
                    candidate_counts=dict(rc),R_native_accepted=len(rr),LR_counts=dict(counts),
                    coherent_full_boundary_scores=scores,best_coherent=best,joins=joins))
                print('FRAME_DONE',len(stages),m,seed['index'],p['index'],q['index'],
                    'R',len(rr),'LR',dict(counts),'coherent',len(scores),flush=True)
    assert len(stages)==27
    paths=[Path(__file__),Path(F.__file__),Path(A.__file__),INPUT,A.INPUT,A.C.PARENT]
    return dict(status='BOUNDED_NUMERICAL_CHANGED12_NONZERO_A2_FOUR_CAP_BRIDGE',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        selection='All27 coherent dotted P/Q frames whose seed has a nonzero A2 profile in the frozen changed12 source domain.',
        alphabet=A.ALPHABET,profiles=profiles,stages=stages,LR_counts=dict(total_counts),
        complete_boundary_histogram=[dict(finite=f,exterior=e,open_ends=n,count=c) for (f,e,n),c in sorted(hist.items())],
        improved=improved,
        scope='Fourth cap R has finite contacts to the shared old side,P,Q and finite projected corners old/P/R and old/Q/R. All seven A2 profiles and both real norm branches are tested over the stated finite alphabet. L/R is checked separately. Every coherent four-cap collection receives a full numerical S6 boundary. This is a bounded discovery test, not an all-label exclusion or an exactly certified replacement chamber.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(compute()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],len(out['stages']),out['LR_counts'],out['complete_boundary_histogram'],flush=True)
