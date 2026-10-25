#!/usr/bin/env python3
"""CURRENT all24 D14 receipts with their source-proved whole witnesses.

The frozen old variant models omit witness declarations at later quarter
corners. The accepted all24 full-core seal supplies those declarations, at
its complete radical/inverse/saturation premises and modulo independent S8.
No old A16 half is restored and no selected or global coverage is claimed.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import types

ROOT = Path(__file__).resolve().parents[1]
LOADER = 'scripts/index_finite57_funded_selected_whole_providers_20261025.py'
LOADER_SHA = '0a245deec21b283f0d47af0323cf9557fd6ba4465f6e25a2b7c47897f2a46e86'
DIRECT = 'scripts/check_finite57_terminal_collectors_direct_receivers_20261025.py'
DIRECT_SHA = 'a2a9d61a7637e45d5a4d7af1c3283e120ca7685d07e9f686ae7b6ee585520eca'
VARIANTS = 'fresh_d14_long_leaf_a16_half_selected_transfer_20261018'
GUARD = 'fresh_d14_small_component_context_retention_20261018'
QUARTERS = 'fresh_rank15_a1_context_composed_ambient58_frontier_20261018'
MAX_OUTPUT = 500 * 1024
EXTRA_PINS = {
    'scripts/'+QUARTERS+'.py':'fad44f22724b8f65131a9f3785f4f1281a0b15e66d20710ede6dd4f016a9d850',
    'results/'+QUARTERS+'.json':'4d04bee005a81bde4631dbd9ce8cdf0cb24a67be17f1cb6582fca47779509256',
    'results/'+QUARTERS+'.md':'f0d7201ff57e4c23969b20db72acfce6707457124f45babedbddbd2025092bd9',
    'results/'+GUARD+'.md':'10fb43eb7646aeb36b1f14c776d3defc76a6f63923d8d96c2f29bbbe01170f02',
    'results/'+VARIANTS+'.md':'ddd6fb81855bc30c6d95da353bc4629cb75f1c2126a85e674d37d573f4af51a2',
    'results/audit_finite57_earlier_exchanges_and_selected_d12_continuation_20261024.md':
        '3fefc45fb7d7e73de0dd7f156c36274ef7b02477dcad13dfab650ea371f9068a',
    'results/review_finite57_bounded_receipts_and_selected_d12_root_20261024.md':
        '201cb6385f376305eb709a12b17191f3b136bda843a1df0860ab93df3ac7e3db',
}
FULL_SEAL = dict(
    kind='accepted_source_specific_full_core_seal',runtime_verified=False,
    source='results/'+GUARD+'.md',section='Complete selected outputs',
    original_source='results/'+VARIANTS+'.md',
    premise='Only the exact24 saved whole variants of whole_D14_transfer: outside an independently paid same-face short8 escape, the radical or saturated-zero-residual source proof seals the ENTIRE old core and any admitted saturated root. Every displayed finite component at EVERY corner, including each positive unmarked witness, is therefore whole; no indirect tail has a first contact. Complete native/contact/inverse/radical/saturation and physical-variant predicates remain accepted inputs.',
    conclusion='Every unmarked component explicitly classified in THIS exact variant/corner is an actual whole factor, not merely a displayed subdiagram.',
    independent_S8_runtime_verified_by_this_batch=False)


def need(ok,message):
    if not ok:
        raise ValueError(message)


def raw(name):
    path = (ROOT/name).resolve()
    need(path.is_relative_to(ROOT) and path.is_file(), 'existing repository source')
    with path.open('rb') as stream:
        data = stream.read(2*1024**2+1)
    need(0 < len(data) <= 2*1024**2, 'bounded complete source')
    return data


def unique_domain(rows,key,description):
    need(type(rows) is list and rows, 'complete nonempty '+description)
    out = {}
    for row in rows:
        identity = key(row)
        need(identity not in out, 'unique exact '+description)
        out[identity] = row
    return out


def prepare():
    encoded = raw(LOADER)
    need(sha256(encoded).hexdigest() == LOADER_SHA, 'unchanged helper-only whole loader')
    helper = types.ModuleType('finite57_d14_witness_loader')
    helper.__file__ = str(ROOT/LOADER)
    exec(compile(encoded,helper.__file__,'exec'),helper.__dict__)
    routes,c = helper.prepare()
    need(helper.COLLECTOR == DIRECT and helper.COLLECTOR_SHA == DIRECT_SHA,
         'exact DIRECT current receiver collector')
    pins = {LOADER:LOADER_SHA,DIRECT:DIRECT_SHA,helper.ROUTES:helper.ROUTES_SHA,
            c.MANIFEST:c.MANIFEST_SHA,**EXTRA_PINS}
    for name,pin in pins.items():
        need(name not in c.PINS or c.PINS[name] == pin, 'accepted source pin collision')
        c.PINS[name] = pin
        need(c.digest(name) == pin, 'source/proof byte identity')
        c.READ.add(name)
    return helper,routes,c


def run(c):
    proofs = c.field(VARIANTS,('selected_proofs',))
    guards = c.field(GUARD,('selected_whole_variants',))
    quarters = c.field(QUARTERS,('local_proof','whole_variants'))
    final = c.field(c.B13,('current_interfaces','D14_variants'))
    identity = lambda p:(p['index'],c.canonical(p['extra_contact']))
    by_guard = unique_domain(guards,identity,'D14 whole-small guards')
    by_quarter = unique_domain(quarters,identity,'latest guarded quarter/current increments')
    by_final = unique_domain(final,identity,'final B13 current corner debits')
    physical = {}
    for pi,proof in enumerate(proofs):
        for vi,variant in enumerate(proof.get('variants',[])):
            need(proof['method'] == 'whole_D14_transfer' and proof['sealing'],
                 'exact source whole-D14 arbitrary-extension proof')
            key = proof['index'],c.canonical(variant['extra_contact'])
            need(key not in physical, 'one exact physical source variant')
            physical[key] = pi,vi,proof,variant
    need(len(physical) == len(by_guard) == len(by_quarter) == len(by_final) == 24
         and set(physical) == set(by_guard) == set(by_quarter) == set(by_final),
         'complete equality of all24 physical/guard/quarter/final domains')
    receipts = [];joins = [];domain = []
    for key,(pi,vi,proof,v) in physical.items():
        guard,latest,debits = by_guard[key],by_quarter[key],by_final[key]
        corners = v['corners'];g = v['gram'];white = v['white']
        need(len(corners) == len(latest['corners']) == len(debits['corners']) == 4,
             'all FOUR original ordered current corners')
        need(len({tuple(sorted(p['marks'])) for p in corners}) == 4,
             'four distinct actual source corners')
        # This mask is justified ONLY by the exact source full-core proof,
        # including saturated variants. It is never inferred from a model.
        seal = dict(FULL_SEAL,source_index=proof['index'],source_proof_ordinal=pi,
            source_variant_ordinal=vi,extra_contact=v['extra_contact'],
            exact_saved_sealing=proof['sealing'])
        witnesses = c.all_witnesses(g,white,corners,seal['premise'])
        got = c.receipt(VARIANTS,[key[0],key[1]],g,white,corners,
            dict(kind='ALL24_CURRENT_WITH_SOURCE_PROVED_FULL_CORE_WITNESSES',
                 source_selector=['selected_proofs',pi,'variants',vi],full_core_seal=seal),
            sealed_witnesses=witnesses)
        need(got['direct_current_receiver_condition'] is True and got['legacy_pointwise_bounds_used'] is False,
             'DIRECT CURRENT total, never old pointwise lower credit')
        per_corner = [];required_total = F()
        for old,paid,delta,actual in zip(corners,latest['corners'],debits['corners'],got['corners'],strict=True):
            marks = old['marks']
            need(marks == paid['marks'] == delta['marks'] == actual['marks'],
                 'same ordered physical corner in every current source ledger')
            observed = c.same_images(g,white,marks,'B',13)
            need(observed == delta['images'], 'complete final finite-B13 images at this SAME corner')
            debit = F(delta['debit']);joint = F(paid['joint_lower'])
            need(debit >= 0 and joint >= 0 and joint == F(paid['guarded_quarter_lower'])+F(paid['rank16_transfer_increment']),
                 'exact latest quarter/current-increment ledger identity')
            required = max(F(),joint-debit)
            need(F(actual['weight']) >= required, ('CURRENT corner meets its latest source-bound lower',
                dict(variant_key=list(key),marks=marks,actual_weight=actual['weight'],
                     required=str(required),joint=str(joint),debit=str(debit),
                     paid_source_row=paid,actual_forest=actual.get('forest'),actual_corner=actual)))
            required_total += required
            mask = witnesses[tuple(sorted(marks))]
            if paid['role'] in ('A15','D15') and F(paid['guarded_quarter_lower']) > 0:
                need(paid['has_unmarked_whole_A1'] is True and mask,
                     'retained rank15 quarter has actual source-sealed unmarked witness')
            per_corner.append(dict(marks=marks,rank15_joint_lower=str(joint),B13_debit=str(debit),
                required_current_lower=str(required),actual_current_lower=actual['weight'],
                whole_unmarked_witnesses=mask,rank15_role=paid['role'],
                has_unmarked_whole_A1=paid['has_unmarked_whole_A1']))
        need(required_total >= 1 and F(latest['total']) == sum((F(p['joint_lower']) for p in latest['corners']),F())
             and F(guard['current_total']) == F(latest['total']),
             'complete latest current source receipt remains at least one after final debit')
        need(len(guard['guarded_outputs']) == 1, 'one exact selected guarded D14 output')
        selected = guard['guarded_outputs'][0];marks = selected['marks']
        need(c.same_images(g,white,marks,'D',14) == selected['images']
             and set(map(tuple,selected['images'])) <= {(0,12),(0,13)},
             'same complete long-leaf D14 native images')
        mask = witnesses[tuple(sorted(marks))]
        for block in selected['whole_witnesses']:
            need(block['type'] in ('A1','A2') and sorted(block['nodes']) in [sorted(p) for p in mask],
                 'selected D14 guard uses THIS actual wholeA1/A2 factor')
        selected_current = next(p for p in got['corners'] if p['marks'] == marks)
        need(F(selected_current['weight']) >= F(selected['retained_half']) == F(1,2),
             'guarded selected D14 current half, no unconditional unit')
        domain.append([key[0],v['extra_contact']])
        receipts.append(got)
        joins.append(dict(index=key[0],extra_contact=v['extra_contact'],source_proof_ordinal=pi,
            source_variant_ordinal=vi,full_core_seal_premise=seal,per_corner=per_corner,
            current_total=got['total'],required_current_total=str(required_total),
            selected_guard=selected,legacy_restored_A16_credit_used=False,
            independent_S8_runtime_verified_by_this_batch=False))
    need(len(receipts) == len(joins) == len(domain) == 24, 'complete all24 result, never truncate')
    return dict(domain=domain,receipts=receipts,source_lower_joins=joins,
        complete_all24_at_declared_premises=True,old_restored_A16_half_used=False,
        full_core_seal_and_numerical_predicates_regenerated=False,
        selected_full_map_coverage_verified=False,independent_S8_runtime_verified=False,
        current_terminal_bridge_complete=False,global_finite57_verified=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',choices=('d14-current24',),required=True)
    parser.add_argument('--parent-lock',action='store_true',required=True)
    args = parser.parse_args()
    own = str(Path(__file__).resolve().relative_to(ROOT));before = raw(own)
    helper,routes,c = prepare()
    result = run(c)
    for name in c.READ:need(c.digest(name) == c.PINS[name], 'consumed D14 source changed')
    need(raw(own) == before and sha256(raw(LOADER)).hexdigest() == LOADER_SHA
         and sha256(raw(DIRECT)).hexdigest() == DIRECT_SHA
         and sha256(raw(helper.ROUTES)).hexdigest() == helper.ROUTES_SHA
         and sha256(c.raw(c.MANIFEST)).hexdigest() == c.MANIFEST_SHA,
         'worker/frozen direct preparation chain changed')
    payload = dict(status='PASS_CURRENT_D14_ALL24_WITNESS_RECEIVERS_AT_EXPLICIT_PREMISES',
        batch=args.batch,source_sha256=sha256(before).hexdigest(),result=result,
        consumed_inputs={name:c.PINS[name] for name in sorted(c.READ)},fuel_used=c.FUEL[0],
        scope='Only all24 exact source physical D14 variants with source-proved whole witnesses at every current corner, latest guarded quarters, D14 guard and final B13 debit joins. Full-core/variant/numerical geometry remains an explicit accepted premise modulo independent S8; no selected full-map or global coverage.',
        independent_S8_verified=False,selected_provider_union_complete=False,
        full_transitive_regeneration=False,geometric_formalization=False,global_finite57_verified=False)
    output = json.dumps(payload,sort_keys=True)
    need(len(output.encode())+1 <= MAX_OUTPUT, 'complete current24 output cap; refuse, never omit corners')
    print(output)


if __name__ == '__main__':
    main()
