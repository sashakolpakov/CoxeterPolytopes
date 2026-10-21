#!/usr/bin/env python3
"""Spectator-gated D17 node3 transfer: exact local and retention interfaces.

The accompanying report supplies the universal guard-persistence argument.
This source reads only small certificates; all hashes are streamed.
"""
import argparse
import json
import signal
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path

import fresh_order17_d5_alternate_half_zero_transfer_20261018 as O
import fresh_finite_d17_a4_d13_half_transfer_20261018 as L
import fresh_d17_node3_large_d_guard_upper_20261018 as U

ROOT=O.ROOT; T=O.T; A=O.A; Q=O.Q; D=O.D; P=O.P; B=P.B
TAG='fresh_d17_node3_large_d_guard_transfer_20261018'
CHANGE={(3,13),(3,14)}
CASES=[1293,2019,2049,2508]


def load(stem):
    path=ROOT/f'results/{stem}.json'
    assert path.stat().st_size<=20*1024*1024,(stem,'stream instead')
    with path.open() as stream:
        out=json.load(stream)
    U.memory_guard()
    assert out['status'].startswith('PASS')
    return out


def has_guard(fs):
    return any(kind=='D' and rank>=18 for kind,rank in fs)


def context_delta(fs):
    if not has_guard(fs):
        return 0
    co=Counter(fs)
    return -64*co['D',17]+128*co['A',4]*co['A',13]+(
        128*co['A',2]*co['A',14] if co['A',1] else 0)


def value(fs):
    return O.value(fs)+context_delta(fs)


def change(g,white,marks):
    """Delta in the formal state with a preserved external whole-D guard."""
    ids=sorted(set(white)|set(marks))
    if Q.inertia(A.sub(g,ids))!=[len(ids),0,0]:
        return dict(delta=F(),roles=[],components=[])
    blocks=D.finite_blocks(g,ids)
    aa,bb=[next(b for b in blocks if m in b['nodes'])for m in marks]
    roles=[]; delta=F()
    if aa is bb and aa['type']=='D17':
        images=T.mark_pairs(g,aa['nodes'],marks,17,'D')
        if set(images)&CHANGE:
            assert set(images)<=CHANGE
            delta=-F(1,2)
            roles.append(dict(role='guarded_D17_node3_zero',images=images))
    elif aa is not bb:
        for sm,la,m,lm in((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
            if(sm['type'],la['type'])==('A4','A13'):
                si=T.mark_pairs(g,sm['nodes'],[m],4,'A')
                li=T.mark_pairs(g,la['nodes'],[lm],13,'A')
                if si==[(0,),(3,)] and li==[(3,),(9,)]:
                    delta=F(1,2)
                    roles.append(dict(role='guarded_A4end_A13fourth_half',small_images=si,large_images=li))
            if(sm['type'],la['type'])==('A2','A14'):
                si=T.mark_pairs(g,sm['nodes'],[m],2,'A')
                li=T.mark_pairs(g,la['nodes'],[lm],14,'A')
                witnesses=[b for b in blocks if b['type']=='A1' and not set(marks)&set(b['nodes'])]
                if si==[(0,),(1,)] and li==[(3,),(10,)] and witnesses:
                    delta=F(1,2)
                    roles.append(dict(role='guarded_A2_A14fourth_A1_half',small_images=si,large_images=li,whole_A1=witnesses))
    assert len(roles)<=1
    return dict(delta=delta,roles=roles,components=blocks)


def local(prior,old):
    support=prior['local_proof']['complete_selected_support']
    assert support['complete_eligible_patterns']==140 and support['complete_maps']==168
    hits=[r for r in support['rows']if set(map(tuple,r['marked_images']))&CHANGE]
    assert [r['index']for r in hits]==CASES
    assert all(len(r['base_nodes'])==17 and set(map(tuple,r['marked_images']))<=CHANGE for r in hits)
    rows=[]
    for index in CASES:
        p=next(p for p in old['selected_proofs']if p['index']==index)
        g=A.matrix(p['gram']);z=list(map(F,p['radical']));white=p['white']
        black=sorted(set(range(len(g)))-set(white))
        assert len(g)==19 and Q.inertia(g)==[17,1,1]
        assert Q.inertia(A.sub(g,white))==[15,0,0]
        assert all(Q.inertia(A.sub(g,white+[b]))==[16,0,0]for b in black)
        assert all(sum(g[i][j]*z[j]for j in range(19))==0 for i in range(19))
        a=black[0]
        adjacent={m for c in p['corners']if a in c['marks']for m in c['marks']if m!=a}
        u=[a,next(m for m in black if m!=a and m not in adjacent)]
        v=[m for m in black if m not in u]
        selected=next(c for c in p['corners']if c['selected'])
        inspected=Q.inspect_corner(g,[1]*19,u+v,*selected['marks'])
        assert inspected['possible_contacts']==p['allowed_contacts']
        allowed=p['allowed_contacts'];assert all(z[i]for i in allowed)
        du=Q.graph_dist(g,u);dv=Q.graph_dist(g,v);short=[]
        for i,j in combinations(allowed,2):
            if z[i]*z[j]<0:
                length=min(du[i]+dv[j],du[j]+dv[i])+2
                assert length<=8
                short.append(dict(nodes=[i,j],length=length))
        corners=[]
        for c in p['corners']:
            ch=change(g,white,c['marks'])
            assert sorted(b['nodes']for b in ch['components'])==sorted(b['nodes']for b in c['components'])
            pay=F(c['payment'])+ch['delta'];assert pay>=0
            corners.append(dict(marks=c['marks'],previous=c['payment'],**ch,payment=pay))
        assert sorted(c['payment']for c in corners)==[0,0,F(1,2),F(1,2)]
        assert sorted(c['delta']for c in corners)==[-F(1,2),0,0,F(1,2)]
        rows.append(dict(index=index,gram=p['gram'],white=white,radical=z,
            allowed_contacts=allowed,mixed_sign_shortcuts=short,corners=corners,
            guard_false_total='1',guard_true_total='1',whole_D_guard_preserved_by_core_orthogonality=True))
    return dict(complete_patterns=140,complete_maps=168,affected_contexts=CASES,
        affected_maps=sum(len(r['maps'])for r in hits),proofs=rows)


def literals(current,old):
    # The only overlap datum T.literal_inputs needs is this physical Gram.
    gram=old['selected_proofs'][0]['gram']
    docs={s:load(s)for s in O.LITERAL if s!=A.OVER}
    docs[A.OVER]={'isolation':{'cores':{'(25, 25, ((9, 9),))':{'gram':gram}}}}
    rows=[];hits=Counter();negative=[]
    for prev,(stem,case,g,white,cs)in zip(current['current_literal_interface']['rows'],T.literal_inputs(docs),strict=True):
        assert(prev['source'],prev['case'])==(stem,O.enc(case))
        corners=[]
        for q,(marks,*_)in zip(prev['corners'],cs,strict=True):
            assert q['marks']==marks
            ch=change(g,white,marks);pay=F(q['lower'])+ch['delta'];assert pay>=0
            hits.update(r['role']for r in ch['roles'])
            if ch['delta']<0:
                assert stem==T.A.D13
                negative.append(dict(source=stem,case=case,marks=marks,core_order=len(g)))
            corners.append(dict(marks=marks,previous=q['lower'],**ch,lower=pay))
        total=sum(c['lower']for c in corners);assert total>=1
        rows.append(dict(source=stem,case=case,corners=corners,guard_true_total=total,guard_false_total=prev['total']))
    assert len(rows)==216 and len(negative)==8
    assert hits['guarded_A2_A14fourth_A1_half']==8 and hits['guarded_A4end_A13fourth_half']==0,dict(hits)
    return dict(configurations=216,rows=rows,changed_roles=dict(hits),
        negative_outputs=negative,affine_D13_equality_guard_preservation_required=True)


def shared_images(g,blocks,marks):
    block=next((b for b in blocks if b['type']=='D17'and set(marks)<=set(b['nodes'])),None)
    if block is None:
        return []
    images=T.mark_pairs(A.matrix(g),block['nodes'],marks,17,'D')
    assert not set(images)&CHANGE or set(images)<=CHANGE
    return images if set(images)&CHANGE else []


def interfaces():
    packets=[];hits=[];count=0
    for stem,keys in B.FINITE.items():
        doc=load(stem);n=0
        for section in keys:
            for path,p,gg,white,blocks in B.corners(doc[section],section):
                images=shared_images(gg,blocks,p['marks'])
                if images:
                    assert stem==L.TAG
                    index=doc['selected_proofs'][int(path.split('/')[1])]['index']
                    assert index in CASES
                    hits.append(dict(source=stem,index=index,path=path,images=images))
                n+=1;count+=1
        packets.append(dict(source=stem,corners=n))
    assert count==564 and len(hits)==4
    # These current proof packets were introduced after the 22-packet inventory.
    newer=[]
    collections=[(P.R.TAG,('local_proof','whole_proofs')),
        (P.H.TAG,('local_proof','whole_proofs')),
        (P.S.TAG,('local_proof','whole_proofs')),
        (P.N.TAG,('local_proof','proofs')),
        (P.R.R.TAG,('local_seals',)),
        (P.H.Z.U.TAG,('local_proof',)),
        (P.TAG,('local_proof','old_ideal17_seals')),
        (P.TAG,('local_proof','A16_seals'))]
    for stem,keys in collections:
        doc=load(stem);xs=doc
        for key in keys:xs=xs[key]
        if isinstance(xs,dict):xs=[xs]
        n=0
        for p in xs:
            physical=p.get('old_complete_seal',p)
            for c in p.get('corners',[]):
                g=physical.get('gram',physical.get('old_seal',{}).get('gram'))
                blocks=c.get('components',[])
                if blocks:
                    assert g is not None
                    assert not shared_images(g,blocks,c['marks'])
                n+=1
        newer.append(dict(source=stem,path=keys,corners=n))
    chosen=[]
    for stem in(B.PART,B.REST):
        doc=load(stem)
        rows=[(path,c)for path,c in P.W.walk(doc)if F(c['payment'])>0]
        assert len(rows)==(280 if stem==B.PART else 52)
        assert all(not any(b['type']=='D17'for b in c['roles'])for _,c in rows)
        chosen.append(dict(source=stem,positive_outputs=len(rows),D17_outputs=0))
    return dict(later_packets=packets,later_corners=count,replaced_own_selected_outputs=hits,
        newer_seals=newer,chosen_D12_outputs=chosen)


def run():
    current=load(O.TAG);prior=load(P.TAG);old=load(L.TAG);upper=load(U.TAG)
    local_proof=local(prior,old)
    literal=literals(current,old)
    retained=interfaces()
    assert F(upper['increasing_family_upper'])==365<372
    witness=[('A',1),('A',16),('A',16),('D',24)]
    changed=[('D',17),('D',17),('D',23)]
    assert O.value(witness)==value(witness)==64*372
    assert O.value(changed)==64*372 and value(changed)==64*370
    paths={Path(__file__),Path(U.__file__),Path(O.__file__),Path(L.__file__),
        ROOT/f'results/{U.TAG}.json',ROOT/f'results/{O.TAG}.json',ROOT/f'results/{P.TAG}.json',ROOT/f'results/{L.TAG}.json'}
    paths.update(ROOT/f'results/{s}.json'for s in set(O.LITERAL)|set(B.FINITE)|{r['source']for r in retained['newer_seals']}|{B.PART,B.REST})
    return O.enc(dict(status='PASS_D17_NODE3_LARGE_D_GUARD_TRANSFER',
        sources={str(p.relative_to(ROOT)):U.digest(p)for p in sorted(paths)},
        changed_D17_pairs=sorted(CHANGE),guard='presence of a separate whole finite D_s with s>=18',
        local_proof=local_proof,current_literal_interface=literal,current_retention=retained,
        finite=dict(complete_increasing_family_upper='365',inherited_complete_upper='372',
            new_exact_global_formal_maximum='372',unchanged_attainer=witness,
            improved_old_attainer=changed,improved_cost='370',new_global_DP=False),
        ideal_weights_unchanged=True,triangle_minimum='2',quadrilateral_minimum='1',
        completed_dimension_bound=57,new_dimension_bound=None,
        scope='Four full17 seals, preserved inert large-D guard, eight affine-D13 equality repairs, current output retention and complete conditional arithmetic upper. No dimension56 theorem.'))


if __name__=='__main__':
    signal.alarm(59)
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists()
    out=run();U.memory_guard()
    if args.verify:
        assert load(TAG)==out
    else:
        path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['finite']['new_exact_global_formal_maximum'],len(out['current_literal_interface']['negative_outputs']))
