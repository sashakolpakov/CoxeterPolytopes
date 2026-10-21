#!/usr/bin/env python3
"""Zero two D16 halves using distinct whole A16 or A1/A15 alternatives.

Current216 literal and typed domination preflight precedes eight exact
whole-core reconstructions. Complete saved support is validated, not searched.
No finite forest optimization.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_d5_d12_alternate_coupon_rebalance_20261018 as C
import fresh_rank16_existing_payment_exchange_20261018 as E
import fresh_center_reversion_rank16_balance_ambient60_frontier_20261018 as B

ROOT=C.ROOT;T=E.T;A=E.A;Q=E.Q;D=E.D;P=E.P
TAG='fresh_finite_d16_alternate_half_zero_transfer_20261018'
CHANGE={(2,14),(2,15)};AHALF={(0,13),(2,15)}
EXPECTED=[1035,1049,1905,1964,2298,2357,2421,2429]
NEW_COUPON=[1049,2298,2357,2429]
R15='fresh_finite_rank15_existing_half_transfer_20261018'
L15='fresh_finite_d15_three_singleton_half_transfer_20261018'
F14='fresh_finite_rank14_existing_payment_transfer_20261018'
STEMS=tuple(dict.fromkeys((*C.STEMS,C.TAG,*E.STEMS,E.TAG,B.TAG,R15,L15,F14)))
enc=C.enc


def actual_base(t):return C.actual_base(t)-(16 if t==('D',16) else 0)
def upper_base(t):return C.upper_base(t)-(8 if t==('D',16) else 0)
def actual_cross(a,b):return C.actual_cross(a,b)+(16 if {a,b}=={('A',1),('A',15)} else 0)
def upper_cross(a,b):return C.upper_cross(a,b)


def coupon(g,ids,marks):
    blocks=D.finite_blocks(g,ids);aa,bb=[next(b for b in blocks if m in b['nodes']) for m in marks]
    if aa is bb:return []
    out=[]
    for sm,la,m,n in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
        if sm['type']=='A1' and la['type']=='A15':
            images=T.mark_pairs(g,la['nodes'],[n],15,'A')
            if set(images)=={(2,),(12,)}:
                out.append(dict(small=sm,small_mark=m,large=la,large_mark=n,large_images=images,weight='1/2'))
    assert len(out)<=1
    return out


def literal(docs):
    old={C.R.key(r['source'],r['case']):r for r in docs[C.TAG]['combined_literal_interface']['rows']}
    assert len(old)==len(docs[C.TAG]['combined_literal_interface']['rows'])==216
    rows=[];counts=Counter()
    for stem,index,g,white,cs in T.literal_inputs(docs):
        prior=old[C.R.key(stem,index)];corners=[];assert len(prior['corners'])==len(cs)
        for j,(marks,*_) in enumerate(cs):
            assert marks==prior['corners'][j]['marks'];ids=sorted(set(white)|set(marks))
            debit=F();roles=[];images=[]
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                blocks=D.finite_blocks(g,ids);shared=next((b for b in blocks if set(marks)<=set(b['nodes'])),None)
                if shared and shared['type']=='D16':
                    images=T.mark_pairs(g,shared['nodes'],marks,16,'D')
                    if set(images)&CHANGE:assert set(images)<=CHANGE;debit=F(1,2)
                roles=coupon(g,ids,marks)
            credit=F(len(roles),2);before=F(prior['corners'][j]['retained_payment_lower_bound']);after=before-debit+credit
            assert after>=0;counts['debits']+=bool(debit);counts['credits']+=bool(credit)
            corners.append(dict(marks=marks,canonical_images=images,debit=debit,coupon_increment=credit,
                coupons=roles,previous_lower=before,retained_payment_lower_bound=after))
        total=sum(r['retained_payment_lower_bound'] for r in corners);assert total>=1
        rows.append(dict(source=stem,case=index,corners=corners,total_lower_bound=total))
    assert len(rows)==216
    return dict(configurations=216,counts=dict(counts),rows=rows,
        changed=[r for r in rows if any(c['debit'] or c['coupon_increment'] for c in r['corners'])],
        minimum_total=min(r['total_lower_bound'] for r in rows))


def weights_and_interfaces(docs):
    prior=docs[B.TAG]['local_proofs']['finite_weighted_pairs']
    weights={typ:{tuple(r['pair']):F(r['weight']) for r in prior[typ]} for typ in ('A16','D16')}
    assert sum(weights['A16'].values())==F(C.actual_base(('A',16)),16)==F(441,4)
    assert sum(weights['D16'].values())==F(C.actual_base(('D',16)),16)==F(443,4)
    assert all(weights['D16'][p]==F(1,2) for p in CHANGE)
    assert all(weights['A16'][p]==F(1,2) for p in AHALF)
    for p in CHANGE:weights['D16'][p]=F()
    assert sum(weights['D16'].values())==F(actual_base(('D',16)),16)==F(439,4)
    keep={(1,12),(1,13),(1,14),(1,15),(2,13)}
    assert not keep&CHANGE and all(weights['D16'][p]==F(1,2) for p in keep)
    assert set(map(tuple,docs[R15]['preserved_D16_alternative_pairs']))=={(1,14),(1,15)}
    assert set(map(tuple,docs[L15]['existing_D16_half_pairs']))=={(1,12),(1,13)}
    assert docs[F14]['existing_D16_half']==[2,13]
    return dict(A16_weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(weights['A16'].items())],
        D16_weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(weights['D16'].items())],
        preserved_later_D16_half_pairs=sorted(keep),existing_whole_A16_half_pairs=sorted(AHALF),
        rank15_older_pairs=docs[R15]['preserved_D16_alternative_pairs'],
        rank15_three_singleton_pairs=docs[L15]['existing_D16_half_pairs'],
        finite14_pair=docs[F14]['existing_D16_half'])


def typed_costs():
    rep=C.P.Z.B.B.representative;types=[t for t in C.P.Z.B.B.types() if t[1]<=59]
    for a in types:
        assert upper_base(rep(a))>=actual_base(a)>=0
        for b in types:assert upper_cross(rep(a),rep(b))>=actual_cross(a,b)==actual_cross(b,a)>=0
    assert C.actual_cross(('A',1),('A',15))==8 and actual_cross(('A',1),('A',15))==24
    assert upper_cross(('A',1),('D',15))==24
    assert upper_base(('D',16))==actual_base(('A',16))==1764 and actual_base(('D',16))==1756
    return dict(unit='1/16',types=len(types),simultaneous_type_pairs=len(types)**2,
        old_actual_D16='443/4',new_actual_D16='439/4',actual_A16_unchanged='441/4',
        old_rank16_representative_upper='443/4',new_rank16_representative_upper='441/4',
        A1_A15_old_actual='1/2',A1_A15_new_actual='3/2',
        A1_D15_representative_allowance_unchanged='3/2',new_residual_key_or_coupon_upper_cost=False,
        all_other_actual_and_representative_costs_unchanged=True,finite_forest_recomputed=False)


def seals(docs,current):
    half=docs[T.TAG];over=docs[A.OVER];sup=E.support('D',half)
    target=next(r for r in half['records'] if r['rank']==16)['target']
    assert E.validate_maps(sup,target,over)==578 and sup['complete_eligible_patterns']==394
    _,dist,_=P.distance_data(target);assert all(dist[i][j]>=10 for i,j in CHANGE)
    assert all({tuple(sorted(m[i] for i in p)) for p in CHANGE}==CHANGE for m in P.embeddings(target,target))
    old={r['index']:r for r in docs[E.TAG]['selected_proofs'] if r['kind']=='D'}
    aw={tuple(r['pair']):F(r['weight']) for r in current['A16_weighted_pairs']}
    proofs=[];hits=[]
    for row in sup['rows']:
        maps=[m for m in row['maps'] if tuple(m['marked_images']) in CHANGE]
        if not maps:continue
        c=over['isolation']['cases'][row['index']];core=over['isolation']['cores'][c['core_key']]
        assert len(c['marked_component'])==16
        proof=E.core_proof(c,core,'D',E.DHALF);assert E.enc(proof)==old[c['index']]
        g=A.matrix(core['gram']);white=proof['white'];corners=[]
        for previous in proof['corners']:
            marks=previous['marks'];ids=sorted(set(white)|set(marks));blocks=D.finite_blocks(g,ids)
            shared=next((b for b in blocks if set(marks)<=set(b['nodes'])),None)
            value=F();reason=None;images=[];new=coupon(g,ids,marks)
            if previous['reason']=='retained_selected_half':
                assert set(previous['canonical_images'])<=CHANGE;reason='selected_D16_zero'
            elif previous['reason']=='existing_whole_small_D13_fork_half':
                small,large=previous['roles'];assert large['type']=='D13'
                # Check intrinsic leaf role rather than preserving a now
                # reduced center coupon by a broad component allowance.
                m=next(m for m in marks if m in large['nodes'])
                assert set(T.mark_pairs(g,large['nodes'],[m],13,'D'))=={(11,),(12,)}
                assert small['type'] in ('A1','A4')
                value=F(1,2);reason='existing_whole_A1_or_A4_D13_shortleaf_half'
            elif shared and shared['type']=='A16':
                images=T.mark_pairs(g,shared['nodes'],marks,16,'A');assert set(images)==AHALF
                assert all(aw[p]==F(1,2) for p in images)
                value=F(1,2);reason='existing_whole_A16_half'
            elif new:
                value=F(1,2);reason='new_whole_A1_A15_third_half'
            corners.append(dict(marks=marks,components=blocks,payment=value,reason=reason,
                canonical_images=images,new_coupons=new))
        assert sum(c['payment'] for c in corners)==1 and sum(bool(c['payment']) for c in corners)==2
        assert sum(c['reason']=='new_whole_A1_A15_third_half' for c in corners)==int(c['index'] in NEW_COUPON)
        assert sum(c['reason']=='existing_whole_A16_half' for c in corners)==int(c['index'] not in NEW_COUPON)
        proofs.append(dict(index=c['index'],old_seal=proof,corners=corners,total='1'))
        hits.append(dict(index=c['index'],maps=maps))
    assert [r['index'] for r in proofs]==EXPECTED
    return dict(proofs=proofs,complete_eligible_patterns=394,validated_saved_maps=578,affected_contexts=hits,
        affected_maps=sum(len(r['maps']) for r in hits),all_selected_components_whole_D16=True,
        proper_displayed_subcomponent_cases=0,geometry_census_repeated=False)


def run(preflight=False):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem,d in docs.items():
        assert d['status'].startswith('PASS'),stem
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    current=weights_and_interfaces(docs);lit=literal(docs);cost=typed_costs()
    if preflight:return enc(dict(status='PASS_D16_ALTERNATE_HALF_ZERO_PREFLIGHT',current_interfaces=current,
        literal_counts=lit['counts'],literal_changed=lit['changed'],typed_costs=cost))
    local=seals(docs,current)
    paths={Path(__file__),Path(C.__file__),Path(E.__file__),Path(E.X.__file__),Path(B.__file__),
        Path(T.__file__),Path(T.T.__file__),Path(A.__file__),Path(Q.__file__),Path(D.__file__),Path(P.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_FINITE_D16_ALTERNATE_HALF_ZERO_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        target=next(r for r in docs[T.TAG]['records'] if r['rank']==16)['target'],changed_pairs=sorted(CHANGE),
        old_weight='1/2',new_weight='0',old_cost='443/4',new_cost='439/4',**current,
        selected_proofs=local,combined_literal_interface=lit,typed_costs=cost,
        new_coupon=dict(small='whole_A1',large='whole_A15',large_role='third node from either endpoint',
            large_marks=[2,12],weight='1/2',eligible_marked_pairs=2,component_pair_increment='1',
            old_actual_component_cost='1/2',new_actual_component_cost='3/2',
            unchanged_representative_pair=['A1','D15'],unchanged_representative_allowance='3/2'),
        all_ideal_weights_unchanged=True,all_other_finite_weights_unchanged=True,
        triangle_minimum=2,quadrilateral_minimum=1,finite_forest_recomputed=False,new_dimension_bound=None,
        scope='All216 current literal interfaces and simultaneous180-type domination; complete saved394/578 D16 map validation and eight freshly reconstructed whole-core seals. Arbitrary-white and current-route proof is written; historical census/complete map enumeration inherited. No finite forest optimization or dimension bound.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--preflight',action='store_true');a=ap.parse_args()
    assert not(a.verify and a.preflight)
    result=run(a.preflight)
    if a.preflight:print(json.dumps(result,indent=2),flush=True)
    else:
        path=ROOT/f'results/{TAG}.json'
        if a.verify:assert result==json.loads(path.read_text())
        else:assert not path.exists();path.write_text(json.dumps(result,indent=2)+'\n')
        print(result['status'],'D16',result['new_cost'],'contexts',len(result['selected_proofs']['proofs']),
            'literal',result['combined_literal_interface']['counts'],flush=True)
