#!/usr/bin/env python3
"""Compose D12 halves with the complete current actual-face proof.

The report supplies the universal output inventory. This binder checks its
finite terminal interfaces, full routing partition and frozen local inputs;
it does not replace arbitrary-support proofs by a finite list of faces.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_order10_low8_ambient58_frontier_20261018 as B
import fresh_finite_d12_long_leaf_support_scout_20261018 as S

ROOT=B.ROOT
TAG='fresh_d12_half_face_composition_20261018'
PART='fresh_d12_radical_saturated_half_replacements_20261018'
REST='fresh_d12_remaining_selected_half_replacements_20261018'
LOW='fresh_d12_low_protected_adaptive_interfaces_20261018'
IDEAL='fresh_d12_universal_ideal_output_repairs_20261018'
FINITE={
 'fresh_finite_a16_d16_d17_sparse_coupon_transfer_20261018':('radical_cases','proper_completion_cases'),
 'fresh_finite_d17_five_unit_reduction_20261018':('radical_transfers',),
 'fresh_rank16_existing_payment_exchange_20261018':('selected_proofs',),
 'fresh_finite_rank15_existing_half_transfer_20261018':('selected_proofs',),
 'fresh_finite_b15_d15_existing_payment_transfer_20261018':('selected_proofs',),
 'fresh_finite_d15_three_singleton_half_transfer_20261018':('selected_proofs',),
 'fresh_finite_d16_fork_exchange_20261018':('selected_proofs','proper_completion_proofs'),
 'fresh_finite_d16_second_pair_zero_transfer_20261018':('radical_proofs','proper_completion_proofs'),
 'fresh_d16_d17_eighth_balance_20261018':('reweighted_proofs',),
 'fresh_finite_d16_alternate_half_zero_transfer_20261018':('selected_proofs',),
 'fresh_finite_d17_leaf_a16_half_transfer_20261018':('selected_proofs',),
 'fresh_finite_d17_a4_d13_half_transfer_20261018':('selected_proofs',),
 'fresh_finite_d17_d5_d12_unit_transfer_20261018':('selected_proofs',),
 'fresh_a16_a17_full_payment_exchange_20261018':('selected_proofs',),
 'fresh_a17_five_pair_alternate_transfer_20261018':('selected_proofs',),
 'fresh_finite_a16_endpoint_sparse_half_transfer_20261018':('selected_proofs',),
 'fresh_d17_existing_three_corner_transfer_20261018':('selected_proofs',),
 'fresh_finite_rank14_existing_payment_transfer_20261018':('selected_proofs',),
 'fresh_finite_a14_near_endpoint_zero_transfer_20261018':('selected_proofs',),
 'fresh_finite_a15_three_coupon_zero_transfer_20261018':('selected_proofs',),
 'fresh_finite_d15_a4_a11_half_zero_transfer_20261018':('selected_proofs',),
 'fresh_finite_d16_d17_two_unit_exchange_20261018':('radical_exchange_proofs',),
}
INVENTORY=(
 'fresh_finite8_ideal5_combined_payments_20261017',
 'fresh_finite_a17_endpoint_payment_support_20261017',
 'fresh_finite_d17_long_pair_support_20261017',
 'fresh_finite_d16_d17_long_pair_half_transfer_20261017',
 'fresh_large_finite_d_leaf_half_transfer_20261017',
 'fresh_small_d16_coupon_and_ideal17_rebalance_20261018',
 'fresh_d5_d12_alternate_coupon_rebalance_20261018',
 'fresh_finite_a15_endpoint_coupon_zero_20261018',
 'fresh_rank15_sparse_payment_exchange_20261018',
 'fresh_rank13_contextual_coupon_rebalance_20261018',
 'fresh_depth6_finite6_zero_transfer_20261018',
 'fresh_depth6_finite7_adaptive_transfer_20261018',
 'fresh_depth6_finite8_zero_transfer_20261018',
)
STEMS=tuple(dict.fromkeys((B.TAG,S.TAG,S.A.OVER,PART,REST,LOW,IDEAL,*FINITE,*INVENTORY)))


def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def actual_base(t):return B.actual_base(t)-(32 if t==('D',12) else 0)
def actual_cross(a,b):
    old=B.actual_cross(a,b)
    pair=sorted((a,b),key=lambda t:(t[1],t[0]))
    if pair==[('A',1),('A',11)]:return old+32
    if pair==[('A',2),('A',11)]:return old+64
    return old


def corners(x,path='',g=None,white=None):
    if isinstance(x,dict):
        inherited=x.get('old_seal',{})
        g=x.get('gram',inherited.get('gram',g))
        white=x.get('white',inherited.get('white',white))
        blocks=x.get('components',x.get('whole_components'))
        if 'marks' in x and blocks is not None:
            assert len(x['marks'])==2 and isinstance(blocks,list)
            assert g is not None and white is not None,(path,list(x))
            yield path,x,g,white,blocks
        for k,v in x.items():yield from corners(v,path+'/'+k,g,white)
    elif isinstance(x,list):
        for i,v in enumerate(x):yield from corners(v,path+'/'+str(i),g,white)


def finite_outputs(docs):
    out=[];target=S.T.target(12,'D')
    for stem,keys in FINITE.items():
        rows=[]
        for key in keys:
            assert key in docs[stem]
            for path,c,gg,white,blocks in corners(docs[stem][key],key):
                g=S.A.matrix(gg);marks=c['marks'];ids=sorted(set(white)|set(marks))
                parts=S.A.components(g,ids)
                assert sorted(parts)==sorted(b['nodes']for b in blocks),(stem,path)
                shared=next((p for p in parts if set(marks)<=set(p)),None)
                images=[]
                if shared is not None and len(shared)==12:
                    images=[m['marked_images']for m in S.T.T.mapped(g,shared,marks,target)]
                assert not any(tuple(p)in S.PAIRS for p in images),(stem,path)
                rows.append(dict(path=path,marks=marks,whole_components=blocks,
                    shared_component=shared,whole_D12_images=images,new_D12_debit='0'))
        assert rows,stem
        out.append(dict(source=stem,proof_sections=list(keys),checked_corners=len(rows),corners=rows))
    return out


def run():
    docs={stem:json.loads((ROOT/f'results/{stem}.json').read_text())for stem in STEMS}
    for stem,d in docs.items():
        assert d['status'].startswith('PASS'),stem
        for p,h in d['sources'].items():assert digest(p)==h,(stem,p)
    assert docs[B.TAG]['completed_unrestricted_bound']==58
    over=docs[S.A.OVER]['isolation']['cases'];partition=Counter();eligible=[]
    for c in over:
        bound=c.get('component_order_bound')
        if bound is None:
            assert c['method']in('finite_leaf9_or_ideal_feature5','nine_context_adaptive_finite_payment')
            partition[c['method']]+=1
        elif bound<12:partition['bound_below12']+=1
        elif len(c['marked_component'])>12:partition['display_above12']+=1
        else:partition['eligible']+=1;eligible.append(c['index'])
    assert partition==dict(bound_below12=259,display_above12=1514,eligible=758,
        finite_leaf9_or_ideal_feature5=74,nine_context_adaptive_finite_payment=9)
    scout=docs[S.TAG]
    assert eligible==[r['index']for r in scout['selected_patterns']]
    assert len(eligible)==scout['complete_eligible_patterns']==758 and scout['complete_maps']==938
    assert len(scout['affected_contexts'])==137 and scout['affected_maps']==288
    first={r['index']for r in docs[PART]['proofs']}
    last={r['index']for r in docs[REST]['new_proofs']}
    affected={r['index']for r in scout['affected_contexts']}
    assert len(first)==102 and len(last)==35 and not first&last and first|last==affected
    assert docs[REST]['selected_contexts_remaining']==0
    assert docs[REST]['selected_patterns']==[
        {k:r[k] for k in ('index','maps','affected_maps')}
        for r in scout['selected_patterns']]
    assert len(scout['literal_interfaces'])==216
    assert all(not c['candidate_debit']for r in scout['literal_interfaces']for c in r['corners'])
    assert all(not r['affected'] and not r['violations']for r in scout['literal_trials'])
    low=docs[LOW]['summary']
    assert(low['low_contexts'],low['unchanged_low'],low['repaired_low'],low['protected_paths'],low['adaptive_paths'])==(66,62,4,74,9)
    ideal=docs[IDEAL]
    assert len(ideal['affine_D11']['rows'])==32 and len(ideal['affine_D13']['rows'])==12
    assert ideal['ideal_weights_changed'] is False and ideal['additional_coupon_beyond_proposal']=='0'
    outputs=finite_outputs(docs)
    assert len(outputs)==22 and sum(r['checked_corners']for r in outputs)==564
    types=[t for t in B.C.C.C.P.Z.B.B.types()if t[1]<=58]
    rep=B.C.C.C.P.Z.B.B.representative
    for a in types:
        assert actual_base(a)>=0
        for b in types:
            assert 0<=actual_cross(a,b)==actual_cross(b,a)<=B.upper_cross(rep(a),rep(b))
    assert actual_base(('D',12))==32*65
    assert actual_base(('A',12))==actual_base(('B',12))==32*66
    assert actual_cross(('A',1),('A',11))==32*2
    assert actual_cross(('A',2),('A',11))==32*3
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{stem}.{ext}'for stem in STEMS
        for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return dict(status='PASS_D12_HALF_CURRENT_FACE_COMPOSITION',
        sources={p:digest(p)for p in sorted(paths)},baseline=B.TAG,
        finite_pair_rule=dict(type='whole D12',pairs=sorted(map(list,S.PAIRS)),old='1',new='1/2',old_cost='66',new_cost='65'),
        added_coupons=[dict(small='whole A1',large='whole A11 endpoint',weight='1/2',pair_increment='1'),
            dict(small='whole A2 any mark',large='whole A11 endpoint',weight='1/2',pair_increment='2')],
        complete_high_partition=dict(partition),eligible_indices=eligible,
        selected_contexts=137,selected_maps=288,low_interface=low,
        universal_ideal_repair_counts=dict(D11=32,D13=12),
        later_finite_output_interfaces=outputs,finite_output_corners=sum(r['checked_corners']for r in outputs),
        inherited_literal_interfaces=216,literal_D12_debits=0,
        simultaneous_type_pairs=len(types)**2,old_representative_cross_allowances_suffice=True,
        ideal_weights_changed=False,contextual_coupons_changed=False,
        triangle_minimum='2',quadrilateral_minimum='1',finite_forest_count_performed=False,
        completed_dimension_bound=58,new_dimension_bound=None,
        scope='Composition of complete local arbitrary-support theorems and explicit universal output inventory. Fresh full high routing and22 later finite terminal-output interfaces; frozen137/66/74/9 and universal ideal repairs. No new forest maximum, ambient bound or example.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists()
    out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['finite_output_corners'],out['simultaneous_type_pairs'],flush=True)
