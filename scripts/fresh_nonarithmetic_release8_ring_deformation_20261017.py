#!/usr/bin/env python3
"""Continuous ring deformation and its full deletion cost: numerical discovery.

A scaffold contact is released; the two central walls are made right by a
scalar root solve. Neither floating-point feasibility nor a projective cone
score is promoted to an exact chamber certificate.
"""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import fresh_nonarithmetic_two_central_caps_20261017 as C

A,B,M=C.A,C.B,C.M
ROOT=C.ROOT
TAG='fresh_nonarithmetic_release8_ring_deformation_20261017'
PARENT=ROOT/'results'/f'{C.TAG}.json'
COST=ROOT/'results/fresh_nonarithmetic_b_release_costs_20261017.json'


def run():
    parent=json.loads(PARENT.read_text());case=parent['closest_finite_endpoint_relaxation']
    assert case['parent_stage']==25 and case['diagonal']==0
    z,labels,ns,h,profiles=A.setup()
    assert np.allclose(profiles[4]['contacts'],[0,-1],rtol=0,atol=1e-15)
    ring=case['ring'];replacement=case['replacement']
    target_rows=[replacement]+ring+[case['S'],case['T']]
    names=[12,'L','P','R','Q','S','T'];pis=[0]+[r['profile'] for r in target_rows[1:]]

    def make(frame,profile,contacts,target):
        gram=frame@M@frame.T
        if abs(np.linalg.det(gram))<1e-10:raise ValueError('singular numerical frame')
        part=np.linalg.solve(gram,contacts)@frame
        kernel=np.linalg.svd(frame@M)[2][-1];den=float(kernel@M@kernel)
        if abs(den)<1e-10:raise ValueError('null numerical kernel')
        rem=(profiles[profile]['square']-part@M@part)/den
        if rem<=0:raise ValueError('no two real norm branches')
        choices=[part+s*np.sqrt(rem)*kernel for s in [-1,1]]
        return min(choices,key=lambda x:np.linalg.norm(x-target))

    def node(frame,fprofiles,row):
        cp=np.array(profiles[row['profile']]['contacts'])
        contacts=-np.array([B.theta(m) for m in row['labels']])-np.array([
            cp@A.J@np.array(profiles[p]['contacts']) for p in fprofiles])
        return make(np.array(frame),row['profile'],contacts,np.array(row['projection']))

    def build(t):
        n12=make(np.array([ns[i] for i in [7,8,52]]),0,
                 np.array([0.,t,-B.theta(replacement['label_to52'])]),
                 np.array(replacement['projection']))
        L=node([ns[7],n12,ns[52]],[4,0,0],ring[0])
        P=node([ns[52],ns[30],L],[0,0,ring[0]['profile']],ring[1])
        Q=node([ns[30],ns[7],L],[0,4,ring[0]['profile']],ring[3])
        R=node([ns[30],P,Q],[0,ring[1]['profile'],ring[3]['profile']],ring[2])
        vv=[L,P,R,Q]
        S=node([vv[i] for i in [0,1,2]],[ring[i]['profile'] for i in [0,1,2]],case['S'])
        T=node([vv[i] for i in [0,3,2]],[ring[i]['profile'] for i in [0,3,2]],case['T'])
        return A.native(S,pis[-2],T,pis[-1],profiles),[n12,L,P,R,Q,S,T]

    initial,iv=build(0.)
    error=max(float(np.linalg.norm(x-np.array(r['projection']))) for x,r in zip(iv,target_rows))
    assert error<1e-7 and abs(initial-case['mutual_product'])<1e-8
    bracket=[0.,0.0001];values=[build(t)[0] for t in bracket]
    assert values[0]<0<values[1]
    parameter=float(brentq(lambda t:build(t)[0],*bracket,xtol=1e-15,rtol=1e-14))
    mutual,v=build(parameter)
    assert abs(mutual)<1e-8
    baseline=B.ledger(np.array([ns[i] for i in z['actual_sides']]),z['actual_sides'])
    assert (baseline['finite'],baseline['exterior'],baseline['open_ends'])==(25,7,15)

    def old_gates(vector,profile,deleted):
        kept=[i for i in labels if i not in deleted]
        vals=np.array([ns[i] for i in kept])@M@vector+np.array([
            profiles[profile]['offsets'][labels.index(i)] for i in kept])
        rels=[B.relation(x) for x in vals]
        return [dict(wall=i,product=float(x),relation=r) for i,x,r in zip(kept,vals,rels)]

    pair_gates=[]
    for i in range(7):
        for j in range(i):
            val=A.native(v[i],pis[i],v[j],pis[j],profiles)
            pair_gates.append(dict(walls=[names[i],names[j]],product=val,relation=B.relation(val)))
    assert all(isinstance(r['relation'],int) or r['relation']=='dotted' for r in pair_gates)
    centers=[float(x@M@h+profiles[p]['center_offset']) for x,p in zip(v,pis)]
    assert max(centers)<-1e-8
    variants=[]
    finite_supports=json.loads(COST.read_text())['finite_vertices'];assert len(finite_supports)==198
    for deleted,include_scaffold in [([8,12,47,53],True),([12,47,53],False)]:
        ii=range(7) if include_scaffold else range(1,7)
        gates={str(names[i]):old_gates(v[i],pis[i],deleted) for i in ii}
        blockers=[dict(new=name,**r) for name,rows in gates.items() for r in rows
                  if not isinstance(r['relation'],int) and r['relation']!='dotted']
        assert not blockers
        actual=[i for i in z['actual_sides'] if i not in deleted]
        if include_scaffold:actual.append(12)
        mapping={**ns,12:v[0]}
        boundary=B.ledger(np.array([mapping[i] for i in actual]+v[1:]),actual+[66,67,68,69,70,71])
        lost=[s for s in finite_supports if set(s)&set(deleted)]
        variants.append(dict(deleted_old_walls=deleted,scaffold12_retained=include_scaffold,
            nominal_wall_count=64-len(deleted)+len(ii),native_old_products=gates,
            native_blockers=blockers,boundary=boundary,
            lost_old_finite_supports=lost,lost_old_finite_count=len(lost),
            loss_scope='Exact loss of old198 finite-support certificates; not a new ambient vertex census.'))
    assert [r['boundary']['open_ends'] for r in variants]==[18,15]
    assert [r['nominal_wall_count'] for r in variants]==[67,67]

    # Keep8 and discard the scaffold12. Attempt to fill the new left cavity
    # without forgetting native contacts to any of the six constructed walls.
    newlabels=[i for i in labels if i not in [12,47,53]]
    newprofiles=[{**p,'offsets':[p['offsets'][labels.index(i)] for i in newlabels]} for p in profiles]
    data=(z,newlabels,ns,h,newprofiles)
    def cap_pool(ids,frame,fprofiles,needed):
        candidates,counts=A.solve_frame(ids,np.array(frame),fprofiles,needed,data)
        accepted=[];rejected=[]
        for row in candidates:
            vec=np.array(row['projection'])
            vals=[A.native(vec,row['profile'],x,p,profiles) for x,p in zip(v[1:],pis[1:])]
            rels=[B.relation(x) for x in vals]
            bad=[dict(wall=name,product=float(x),relation=r) for name,x,r in zip(names[1:],vals,rels)
                 if not isinstance(r,int) and r!='dotted']
            if bad:rejected.append(dict(candidate=row,blockers=bad))
            else:accepted.append(dict(**row,new_products=vals,new_relations=rels))
        return dict(frame=ids,needed=needed,raw_count=len(candidates),producer_counts=dict(counts),
                    accepted=accepted,rejected=rejected)

    restore=cap_pool([7,8,52],[ns[7],ns[8],ns[52]],[4,0,0],[(0,1),(1,2)])
    assert restore['raw_count']==13 and not restore['accepted']
    left=cap_pool([7,8,-1],[ns[7],ns[8],v[1]],[4,0,pis[1]],[(0,1),(0,2)])
    right=cap_pool([8,52,-1],[ns[8],ns[52],v[1]],[0,0,pis[1]],[(0,1),(1,2)])
    pair_rows=[];pair_counts=Counter()
    for u in left['accepted']:
        for w in right['accepted']:
            value=A.native(np.array(u['projection']),u['profile'],np.array(w['projection']),w['profile'],profiles)
            rel=B.relation(value);pair_counts[str(rel)]+=1
            pair_rows.append(dict(left=u['index'],right=w['index'],product=value,relation=rel))
    assert len(left['accepted'])==2 and len(right['accepted'])==38
    assert pair_counts==Counter({'None':53,'dotted':23})
    three_seed=cap_pool([7,-1,52],[ns[7],v[1],ns[52]],[4,pis[1],0],[(0,1),(1,2)])
    assert three_seed['raw_count']==0
    paths=[Path(__file__),Path(C.__file__),Path(A.__file__),Path(B.__file__),PARENT,A.INPUT,A.C.PARENT,COST]
    return dict(status='BOUNDED_NUMERICAL_CONTINUOUS_RING_CLOSURE_WITH_FULL_DELETION_COSTS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        parent_stage=25,parent_diagonal=0,baseline_reconstruction_error=error,
        free_parameter='Native product of scaffold12 with old8; products with7 and52 remain fixed.',
        bracket=bracket,bracket_mutual_products=values,parameter=parameter,
        solved_mutual_product=mutual,target_mutual_label=2,
        roots=[dict(name=name,profile=p,projection=x.tolist(),center_product=ct)
               for name,p,x,ct in zip(names,pis,v,centers)],
        pair_products=pair_gates,baseline_boundary=baseline,variants=variants,
        restoration_alphabet=A.ALPHABET,restoration_one_cap=restore,
        restoration_two_caps=dict(left=left,right=right,pair_counts=dict(pair_counts),pairs=pair_rows),
        restoration_three_cap_seed=three_seed,
        optimizer_invocations=0,scalar_root_solves=1,new_exact_chamber=False,new_finite_volume_example=False,
        scope='One explicit continuous recipe near parentstage25; numerical branch selection/root solving/native gates/complete intrinsic S6 ledgers only. Releasing12-right8 fixes the mutual right angle. The two globally native-coherent wall selections have18or15 open ends after full deletion costs, not the invalid11-end selection retaining both8andscaffold12. Restoration searches use all7 inherited A2profiles but only the declared finite label alphabet. No all-parameter exclusion, interval-certified root, exact chamber, ambient vertex census or finite-volume claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=json.loads(json.dumps(run()));out=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],'parameter',result['parameter'],'ends',
          [r['boundary']['open_ends'] for r in result['variants']],
          'restore_pairs',result['restoration_two_caps']['pair_counts'],flush=True)
