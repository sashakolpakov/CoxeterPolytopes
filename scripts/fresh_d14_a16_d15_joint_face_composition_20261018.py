#!/usr/bin/env python3
"""Compose the D14/A16 exchange and the small/D15 center retention.

The universal proof-output replay is in the report. This binder checks the
complete routing, frozen selected/interface domains, exact weight tables
and simultaneous actual costs; it does not perform a forest count.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import fresh_d12_half_face_composition_20261018 as B
import fresh_d12_type_split_ambient58_forest_frontier_20261018 as C
import fresh_d14_long_leaf_a16_half_selected_transfer_20261018 as S
import fresh_small_d15_center_three_eighth_retention_20261018 as U

ROOT=B.ROOT
TAG='fresh_d14_a16_d15_joint_face_composition_20261018'
EXTRA=(
 'fresh_order11_bcd_half_pair_transfer_20261018',
 'fresh_order14_d_near_fork_half_transfer_20261018',
 'fresh_order15_bd_basic_fractional_transfer_20261018',
 'fresh_order10_d_opposite_fork_fractional_transfer_20261018',
 'fresh_order10_b_terminal_half_transfer_20261018',
 'fresh_order10_c_endpoint_zero_transfer_20261018',
 'fresh_order10_cycle_half_transfer_20261018',
 'fresh_finite_a17_endpoint_payment_support_20261017',
 'fresh_small_d15_fork_coupon_provenance_20261018',
 'fresh_affine_d11_leaf_pair_transfer_20261017',
 'fresh_cycle16_half_coupon_transfer_20261017',
 'fresh_b14_b15_singleton_completion_transfer_20261017',
)
STEMS=tuple(dict.fromkeys((B.TAG,C.TAG,S.TAG,U.TAG,*B.STEMS,*EXTRA,
 *(s for _,ss,_ in U.I.ROUTES for s in ss),
 'fresh_affine_d12_leaf_pair_transfer_20261017',
 'fresh_order15_affine_b_pair_reduction_20261017',S.H.TAG,S.AW)))


def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def actual_base(t):return C.actual_base(t)-(32 if t==('D',14)else 0)+(16 if t==('A',16)else 0)
def actual_cross(a,b):return U.actual_cross(a,b)
def context(fs,upper=False):return C.context(fs,upper)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    selected=docs[S.TAG];center=docs[U.TAG];baseline=docs[B.TAG]
    assert selected['baseline']==center['baseline']==B.TAG
    assert center['triangle_minimum']=='2'and center['quadrilateral_minimum']=='1'
    assert center['same_component_weights_changed'] is False and center['ideal_weights_changed'] is False
    assert selected['no_new_cross_coupon'] is True and selected['ideal_weights_changed'] is False
    assert baseline['triangle_minimum']=='2'and baseline['quadrilateral_minimum']=='1'
    over=docs[S.A.OVER]['isolation']['cases'];partition=Counter();eligible=[]
    for c in over:
        bound=c.get('component_order_bound')
        if bound is None:
            assert c['method']in('finite_leaf9_or_ideal_feature5','nine_context_adaptive_finite_payment')
            partition[c['method']]+=1
        elif bound<14:partition['bound_below14']+=1
        elif len(c['marked_component'])>14:partition['display_above14']+=1
        else:partition['eligible']+=1;eligible.append(c['index'])
    assert partition==dict(bound_below14=1378,display_above14=553,eligible=600,
        finite_leaf9_or_ideal_feature5=74,nine_context_adaptive_finite_payment=9)
    assert len(over)==2614 and eligible==[r['index']for r in selected['selected_patterns']]
    assert sum(len(r['maps'])for r in selected['selected_patterns'])==1104
    assert sum(len(r['affected_maps'])for r in selected['selected_patterns'])==36
    affected=[r['index']for r in selected['selected_patterns']if r['affected_maps']]
    assert affected==[r['index']for r in selected['selected_proofs']]==S.EXPECTED
    assert selected['whole_variants']==24 and selected['proper_projection_exclusions']==4
    for p in selected['selected_proofs']:
        for v in p.get('variants',[]):
            assert F(v['total'])==sum(F(c['payment'])for c in v['corners'])>=1
            # No changed small/D15 center coupon is used by a new alternative.
            assert all(not(c['roles']and any(r['type']=='D15'for r in c['roles']))for c in v['corners'])
    literal=selected['literal_interfaces'];assert len(literal['records'])==216 and literal['changed_corners']==16
    assert all(F(r['total'])>=1 for r in literal['records'])
    later=selected['later_finite_interfaces'];assert len(later)==22 and sum(len(r['corners'])for r in later)==564
    assert all(c['new_debit']=='0'for r in later for c in r['corners'])
    assert len(selected['affine_D12_two_A1_interfaces'])==16
    assert all(r['total']=='1'and len(r['corners'])==2 for r in selected['affine_D12_two_A1_interfaces'])
    dw={tuple(r['pair']):F(r['weight'])for r in selected['D14_weighted_pairs']}
    aw={tuple(r['pair']):F(r['weight'])for r in selected['A16_weighted_pairs']}
    olda={tuple(r['pair']):F(r['weight'])for r in docs[S.AW]['A16_weighted_pairs']}
    oldd={tuple(r['pair']):F(r['weight'])for r in next(r for r in docs[S.H.TAG]['rules']if r['kind']=='D')['weighted_pairs']}
    assert set(dw)==set(oldd)and set(aw)==set(olda)
    assert {p for p in dw if dw[p]!=oldd[p]}==S.CHANGES
    assert {p for p in aw if aw[p]!=olda[p]}=={S.RESTORE}
    assert all(dw[p]==F(1,2)and oldd[p]==1 for p in S.CHANGES)
    assert aw[S.RESTORE]==F(1,2)and olda[S.RESTORE]==0
    assert min(dw.values())==F(1,2)and dw[0,10]==1 and dw[0,11]==F(1,2)
    assert dw[1,12]==dw[1,13]==1 and aw[0,13]==aw[2,15]==F(1,2)
    assert sum(dw.values())==F(actual_base(('D',14)),32)==F(179,2)
    assert sum(aw.values())==F(actual_base(('A',16)),32)==F(443,4)
    unchanged=0
    for a in C.TYPES:
        assert actual_base(a)>=0
        if a not in(('D',14),('A',16)):
            assert actual_base(a)==C.actual_base(a);unchanged+=1
        for b in C.TYPES:
            assert 0<=actual_cross(a,b)==actual_cross(b,a)==C.actual_cross(a,b)-U.debit(a,b)
    assert unchanged==175
    low=docs[B.LOW]['summary']
    assert(low['low_contexts'],low['protected_paths'],low['adaptive_paths'])==(66,74,9)
    assert baseline['selected_contexts']==137 and baseline['selected_maps']==288
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in
        [('scripts','py'),('results','json'),('results','md')]}
    return S.enc(dict(status='PASS_D14_A16_SMALL_D15_JOINT_FACE_COMPOSITION',
        sources={p:digest(p)for p in sorted(paths)},baseline=B.TAG,
        selected_exchange=S.TAG,center_retention=U.TAG,complete_high_partition=dict(partition),
        eligible_indices=eligible,selected_contexts=17,selected_maps=36,whole_variants=24,
        proper_completion_exclusions=4,historical_literal_records=216,historical_literal_debits=16,
        later_terminal_corners=564,later_terminal_debits=0,affine_D12_two_A1_cases=16,
        low_contexts=66,protected_contexts=74,adaptive_contexts=9,
        D12_selected_repairs_retained=137,D12_alternative_D14_minimum='1/2',
        D14_cost='179/2',A16_cost='443/4',small_D15_center_weight='3/8',
        D14_new_pairs=sorted(S.CHANGES),A16_restored_pair=S.RESTORE,
        small_D15_total_costs={str(r):F(11*r,8)for r in(3,4,5)},
        unchanged_actual_bases=unchanged,simultaneous_cross_identities=len(C.TYPES)**2,
        ideal_weights_changed=False,contextual_rules_changed=False,new_cross_coupons=[],
        triangle_minimum='2',quadrilateral_minimum='1',finite_forest_count_performed=False,
        completed_unrestricted_bound=58,new_dimension_bound=None,
        scope='Whole-face replay of complete current universal and selected proof with both local rules. Complete selected routing and interfaces supplement arbitrary-support proofs; fresh finite forest optimization and any dimension consequence are separate.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists()
    out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['selected_contexts'],out['simultaneous_cross_identities'],flush=True)
