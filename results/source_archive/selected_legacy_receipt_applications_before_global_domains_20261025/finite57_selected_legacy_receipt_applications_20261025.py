#!/usr/bin/env python3
"""Exact selected-key applications of source receipt and exclusion arguments.

This pure helper opens no files and imports no repository worker. Its caller
must inject byte-verifying SOURCE_PINNER and complete parsed SOURCE_READER,
load the frozen current literal W, and supply ONLY F.actual-validated units.
No standalone execution, native-map enumeration or new universal role is
provided. All geometry and independent same-face S8 premises stay explicit.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

from fractions import Fraction as F
from hashlib import sha256
import json
import re

SOURCE_READER = None
SOURCE_PINNER = None
CACHE = {}
SHA = re.compile(r'[0-9a-f]{64}\Z')
FIELDS = {'high_index','core_key','base_nodes','selected_marks','target_kind',
          'target_order','affine','induced_map','marked_image'}
TEX = 'manuscripts/finite_volume57.tex'
TEX_SHA = 'c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0'
LITERAL = 'scripts/finite57_literal_weights_audit_20261024.py'
LITERAL_SHA = 'f6eba7af77af50b821e9e99c589388c1352f9ba5a99052aaeee64393f5ada1e0'
D12 = 'results/fresh_d12_radical_saturated_half_replacements_20261018.json'
D12_DOMAIN = 'results/fresh_d12_remaining_selected_half_replacements_20261018.json'
D14 = 'results/fresh_d14_long_leaf_a16_half_selected_transfer_20261018.json'
GUARD = 'results/fresh_d14_small_component_context_retention_20261018.json'
D14_REPORT = 'results/fresh_d14_long_leaf_a16_half_selected_transfer_20261018.md'
GUARD_REPORT = 'results/fresh_d14_small_component_context_retention_20261018.md'
QUARTERS = 'results/fresh_rank15_a1_context_composed_ambient58_frontier_20261018.json'
B13 = 'results/fresh_b13_endpoint_retention_20261018.json'
A14 = 'results/fresh_a14_endpoint_a1_context_half_transfer_20261018.json'
A14_REPORT = 'results/fresh_a14_endpoint_a1_context_half_transfer_20261018.md'
BC12 = 'results/fresh_order12_bc_unit_coupon_transfer_20261017.json'
FINITE_PROOF = 'results/audit_finite57_earlier_exchanges_and_selected_d12_continuation_20261024.md'
FINITE_ROOT = 'results/review_finite57_bounded_receipts_and_selected_d12_root_20261024.md'
BC_PROOF = 'results/fresh_order12_bc_unit_coupon_transfer_20261017.md'
G4_ROOT = 'results/finite57_g4_selected_terminal_source_clearance_root_20261025.json'
HIGH = 'results/fresh_distance8_high_overlap_20261017.json'
B17_NOTE = 'results/finite57_selected512_affine_b17_forced_entry_source_proof_20261025.md'
D12_WORKER = 'scripts/check_finite57_current_d12_terminal_bindings_guard_safe_20261025.py'
D12_SHA = 'cb24da35f459dc3b54240ddbf84048c0cfcbeec21f58f83152c8edf8ea9a5b79'
G4_WORKER = 'scripts/check_finite57_g4_lower_cycle_terminal_bindings_20261025.py'
G4_SHA = 'f17014594587a68ef766ae245db372550fe89adc02bedf9e44c477a8d79eb242'
D14_WORKER = 'scripts/check_finite57_d14_current24_witness_receivers_20261025.py'
D14_SHA = '8f606376b015f7b1571635a30d5d2b23e0fb34bc3dfc5efe88484e1008ad6b8b'
D14_STATUS = 'PASS_CURRENT_D14_ALL24_WITNESS_RECEIVERS_AT_EXPLICIT_PREMISES'
PINS = {
    TEX:TEX_SHA,LITERAL:LITERAL_SHA,
    D12:'2d39a231e9d65d295017695913c32951df792edcf1cfc876ab0c272f199b379c',
    D12_DOMAIN:'931c5abc3ff4807a8bc31afb6a8c54ad8f8aeabf4c4b7cfe24288c715fc755a3',
    D14:'fbf2be51d353ad9eef59a29f8017be59f92057cc05bbece5cce8900c33a14ac2',
    GUARD:'fb484aa601a397838fa63ac61c48c2a8d5f7fdc32c5d1de29505919126fee0e5',
    D14_REPORT:'ddd6fb81855bc30c6d95da353bc4629cb75f1c2126a85e674d37d573f4af51a2',
    GUARD_REPORT:'10fb43eb7646aeb36b1f14c776d3defc76a6f63923d8d96c2f29bbbe01170f02',
    QUARTERS:'4d04bee005a81bde4631dbd9ce8cdf0cb24a67be17f1cb6582fca47779509256',
    B13:'2bb463c567b393d7c95f211f8199c013bc0c34716a4dc98c98251fad7db2b2da',
    A14:'14af340c03fe54cdd6a39c4110e45b701cb2839f941539c3fac116097777e6da',
    A14_REPORT:'bdfacca37a915bd77a234ed2dd065971d0d1a4ebf1048f0fd520e37618af49dc',
    BC12:'344345cba676d15466d8809d4cc8f7d828442b396d2dd38709d5014aca9f9a84',
    FINITE_PROOF:'3fefc45fb7d7e73de0dd7f156c36274ef7b02477dcad13dfab650ea371f9068a',
    FINITE_ROOT:'201cb6385f376305eb709a12b17191f3b136bda843a1df0860ab93df3ac7e3db',
    BC_PROOF:'079467a45fdd8d2715899de2adce183c7fd9cccf67aaed269c2586024d3b691c',
    G4_ROOT:'29e14b3c73d11bba48a97c82e5d185433a74270486c4c605413979a80bc23a06',
    HIGH:'f8a2ed3bdef6064febf7336bbe4ad75d71c7d35b5747d9f57b0ebf26b03d0cc9',
    B17_NOTE:'99b51ec8f1bd696b47d617f840129072ab7f9740e66408a65f79a84af8141169',
}
D12_CASES = {
    448:dict(proof=27,domain=170,selected=2,paying=3,core='(21, 25, ((8, 9),))'),
    1600:dict(proof=68,domain=520,selected=1,paying=3,core='(25, 24, ((9, 9),))'),
}
D14_CASES = {
    1920:dict(proof=12,variant=0,guard=15,selected=1,core='(25, 25, ((9, 9),))'),
    2240:dict(proof=14,variant=0,guard=19,selected=0,core='(25, 22, ((9, 8),))'),
}
B17_MAPS = [
    [9,8,7,6,10,11,12,13,0,15,16,1,2,3,4,5],
    [9,8,7,6,10,11,12,13,0,16,15,1,2,3,4,5],
]


def need(ok,message):
    if not ok:
        raise ValueError(message)


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=True)


def digest(value):
    need(type(value) is str and SHA.fullmatch(value), 'exact source/evidence SHA256')
    return value


def rational(value):
    need(type(value) in (str,int) and len(str(value)) <= 128, 'exact bounded receipt rational')
    return F(value)


def pin(name,expected=None):
    need(callable(SOURCE_PINNER), 'mandatory injected source/proof/evidence byte pinner')
    expected = PINS[name] if expected is None else digest(expected)
    SOURCE_PINNER(name,expected)
    return dict(path=name,sha256=expected)


def source(name):
    pin(name)
    need(callable(SOURCE_READER), 'mandatory injected complete parsed-source reader')
    if name not in CACHE:
        CACHE[name] = SOURCE_READER(name,PINS[name])
    need(type(CACHE[name]) is dict, 'complete original source object')
    return CACHE[name]


def one(rows,predicate,message):
    need(type(rows) is list, 'complete named source/receipt list')
    found = [row for row in rows if predicate(row)]
    need(len(found) == 1,message)
    return found[0]


def full_key(key):
    need(type(key) is dict and set(key) == FIELDS, 'complete exact nine-field selected key')
    need(type(key['high_index']) is int and 0 <= key['high_index'] < 2614
         and type(key['core_key']) is str and 0 < len(key['core_key']) <= 256
         and type(key['affine']) is bool
         and key['target_kind'] in (('A','B','C','D') if key['affine'] else ('A','B','D'))
         and type(key['target_order']) is int and 10 <= key['target_order'] <= 17,
         'exact selected target/high schema')
    nodes,marks,mp,pair = (key[n] for n in ('base_nodes','selected_marks','induced_map','marked_image'))
    for values in (nodes,marks,mp,pair):
        need(type(values) is list and all(type(v) is int and v >= 0 for v in values)
             and len(values) == len(set(values)), 'ordered distinct integer selected-key lists')
    need(10 <= len(nodes) <= key['target_order'] and max(nodes) < 64
         and len(mp) == len(nodes) and set(mp) <= set(range(key['target_order']))
         and len(marks) == len(pair) == 2 and set(marks) <= set(nodes)
         and pair == sorted(mp[nodes.index(m)] for m in marks),
         'full native-map/ordered physical mark association')


def member(key,rows):
    need(type(rows) is list and rows and all(type(p) is dict for p in rows),
         'complete original affected-map list')
    matching = [i for i,p in enumerate(rows) if p.get('images') == key['induced_map']
                and p.get('marked_images') == key['marked_image']]
    need(len(matching) == 1, 'THIS full native map in the complete original affected-map list')
    return matching[0]


def unit(units,worker,batch,source_sha,status,pin_field,required):
    digest(source_sha);pin(worker,source_sha)
    need(type(units) is dict and (worker,batch) in units, 'required exact actual receipt is missing')
    item = units[worker,batch]
    need(type(item) in (tuple,list) and len(item) == 2, 'validated actual unit pair')
    entry,doc = item
    need(type(entry) is dict and type(doc) is dict
         and entry['worker'] == worker and entry['batch'] == batch
         and entry['inner_status'] == doc['status'] == status
         and doc['batch'] == batch and doc['source_sha256'] == source_sha,
         'exact actual worker/source/inner-status/batch association')
    evidence = dict(worker=worker,batch=batch,source_sha256=source_sha)
    for field,hfield,sizefield in (('path','sha256','bytes'),('stdout_path','stdout_sha256','stdout_bytes')):
        name = entry[field]
        need(type(name) is str and name.startswith('results/')
             and type(entry[sizefield]) is int and entry[sizefield] > 0,
             'complete actual ledger/stdout identity')
        pin(name,entry[hfield])
        evidence[field] = name;evidence[hfield] = entry[hfield];evidence[sizefield] = entry[sizefield]
    pins = doc[pin_field]
    need(type(pins) is dict and all(pins.get(name) == PINS[name] for name in required),
         'actual receipt consumed THIS exact physical/literal source')
    need(doc.get('global_finite57_verified',False) is False
         and doc.get('current_terminal_bridge_complete',False) is False,
         'local receipt scope')
    return doc,evidence


def context(gram,white):
    need(type(gram) is list and 0 < len(gram) <= 64
         and all(type(row) is list and len(row) == len(gram) for row in gram)
         and type(white) is list and len(white) == len(set(white))
         and all(type(v) is int and 0 <= v < len(gram) for v in white),
         'complete original physical Gram/white source schema')
    return dict(gram=gram,white=white,
        context_sha256=sha256(canonical([gram,white]).encode()).hexdigest())


def d12(key,units):
    case = D12_CASES[key['high_index']]
    original = source(D12);domain = source(D12_DOMAIN)
    proof = original['proofs'][case['proof']]
    saved = domain['selected_patterns'][case['domain']]
    need(proof['index'] == saved['index'] == key['high_index']
         and proof['core_key'] == key['core_key'] == case['core']
         and proof['method'] == 'radical_annihilation',
         'exact D12 source proof/high/core/radical method')
    selected = proof['corners'][case['selected']]
    need(selected['marks'] == key['selected_marks'] and selected['reason'] == 'proposed_D12_half'
         and rational(selected['payment']) == F(1,2)
         and len(selected['roles']) == 1 and selected['roles'][0]['type'] == 'D12'
         and selected['roles'][0]['nodes'] == key['base_nodes'],
         'same whole D12 selected physical corner')
    mi = member(key,saved['affected_maps'])
    need(any(p['images'] == key['induced_map'] and p['marked_images'] == key['marked_image']
             for p in saved['maps']), 'same complete saved D12 map domain')
    need(proof['radical'] and proof['remaining_contact_positions'] and proof['mixed_sign_shortcuts']
         and all(type(p['length']) is int and 0 < p['length'] <= 8 for p in proof['mixed_sign_shortcuts']),
         'exact retained D12 radical/contact/shortpath predicates')
    doc,evidence = unit(units,D12_WORKER,'d12-first',D12_SHA,
        'PASS_GUARD_SAFE_CURRENT_D12_TERMINAL_LITERAL_BINDING','inputs',(D12,LITERAL,TEX))
    paid = []
    for ci in (case['selected'],case['paying']):
        path = '/proofs/'+str(case['proof'])+'/corners/'+str(ci)
        actual = one(doc['rows'],lambda p:p['path'] == path,'one exact actual D12 terminal path')
        original_corner = proof['corners'][ci]
        need(actual['marks'] == original_corner['marks']
             and actual['gram_path'] == '/proofs/'+str(case['proof'])+'/gram'
             and actual['white_path'] == '/proofs/'+str(case['proof'])+'/white'
             and actual['possible_negative_D13_guard'] is False
             and rational(actual['chosen']) == rational(original_corner['payment'])
             and rational(actual['current_lower']) >= rational(actual['chosen']),
             'same physical D12 source corner/current terminal receiver')
        paid.append(dict(source_selector=['proofs',case['proof'],'corners',ci],
            actual_path=path,marks=actual['marks'],current_lower=actual['current_lower']))
    need(paid[0]['marks'] != paid[1]['marks']
         and sum((rational(p['current_lower']) for p in paid),F()) >= 1,
         'distinct current D12-source corners pay THIS same face')
    return dict(source=D12,proof_selector=['proofs',case['proof']],
        map_source=D12_DOMAIN,map_selector=['selected_patterns',case['domain'],'affected_maps',mi],
        physical=context(proof['gram'],proof['white']),source_sealing=dict(
            method=proof['method'],radical=proof['radical'],gauge=proof['gauge'],
            possible_contacts=proof['remaining_contact_positions'],
            mixed_sign_shortcuts=proof['mixed_sign_shortcuts']),
        paying_corners=paid,actual_receipts=[evidence],proof_references=[pin(FINITE_PROOF),pin(FINITE_ROOT)])


def d14(key,units):
    case = D14_CASES[key['high_index']]
    original = source(D14);guards = source(GUARD)
    proof = original['selected_proofs'][case['proof']];v = proof['variants'][case['variant']]
    guard = guards['selected_whole_variants'][case['guard']]
    need(proof['index'] == guard['index'] == key['high_index']
         and proof['core_key'] == key['core_key'] == case['core']
         and proof['method'] == 'whole_D14_transfer'
         and proof['sealing']['method'] == 'radical_annihilation'
         and v['extra_contact'] is None and guard['extra_contact'] is None,
         'exact D14 original radical physical variant')
    mi = member(key,proof['affected_maps'])
    chosen = one(guard['guarded_outputs'],lambda p:p['marks'] == key['selected_marks'],
        'exact selected D14 whole-small guarded output')
    selected = v['corners'][case['selected']]
    factor = one(selected['components'],
        lambda p:p['type'] == 'D14' and set(key['selected_marks']) <= set(p['nodes']),
        'THIS source selected D14 whole factor')
    need(selected['marks'] == key['selected_marks'] and factor['nodes'] == key['base_nodes']
         and key['marked_image'] in chosen['images'] and rational(chosen['retained_half']) == F(1,2),
         'same source D14 base/ordered marks/native image/guarded half')
    doc,evidence = unit(units,D14_WORKER,'d14-current24',D14_SHA,D14_STATUS,
        'consumed_inputs',(D14,GUARD,QUARTERS,A14,B13,LITERAL,TEX,D14_REPORT,GUARD_REPORT,A14_REPORT))
    result = doc['result']
    expected = [[p['index'],b['extra_contact']] for p in original['selected_proofs'] for b in p.get('variants',[])]
    need(len(expected) == 24 and result['domain'] == expected
         and result['complete_all24_at_declared_premises'] is True,
         'complete EXACT source physical24 domain, not count alone')
    expected_receivers = [[high,canonical(extra)] for high,extra in expected]
    need(type(result['receipts']) is list and len(result['receipts']) == 24
         and [p['index'] for p in result['receipts']] == expected_receivers
         and all(p['source'] == D14[8:-5] for p in result['receipts'])
         and type(result['source_lower_joins']) is list
         and [[p['index'],p['extra_contact']] for p in result['source_lower_joins']] == expected,
         'complete ordered actual24 receiver/lower-join identities')
    actual = one(result['receipts'],lambda p:p['source'] == D14[8:-5]
        and p['index'] == [key['high_index'],'null'],'one actual current D14 physical receiver')
    join = one(result['source_lower_joins'],
        lambda p:p['index'] == key['high_index'] and p['extra_contact'] is None,
        'one exact D14 latest-source lower join')
    physical = context(v['gram'],v['white'])
    seal = join['full_core_seal_premise']
    need(join['source_proof_ordinal'] == case['proof']
         and join['source_variant_ordinal'] == case['variant']
         and type(seal) is dict and seal['kind'] == 'accepted_source_specific_full_core_seal'
         and seal['runtime_verified'] is False
         and seal['source'] == GUARD_REPORT and seal['original_source'] == D14_REPORT
         and seal['source_index'] == key['high_index']
         and seal['source_proof_ordinal'] == case['proof']
         and seal['source_variant_ordinal'] == case['variant']
         and seal['extra_contact'] is None and seal['exact_saved_sealing'] == proof['sealing']
         and seal['independent_S8_runtime_verified_by_this_batch'] is False,
         'same exact source-specific conditional FULL-core seal, not a generic witness flag')
    need(actual['context_sha256'] == physical['context_sha256']
         and actual['direct_current_receiver_condition'] is True
         and actual['legacy_pointwise_bounds_used'] is False
         and actual['provenance']['source_selector'] == ['selected_proofs',case['proof'],'variants',case['variant']]
         and actual['provenance']['full_core_seal'] == seal
         and len(actual['corners']) == len(join['per_corner']) == len(v['corners']) == 4,
         'actual D14 source/context/four-corner DIRECT binding')
    paid = []
    for old,got,lower in zip(v['corners'],actual['corners'],join['per_corner'],strict=True):
        joint = rational(lower['rank15_joint_lower'])
        a14_debit = rational(lower['A14_endpoint_debit'])
        a14_retained = rational(lower['A14_retained_lower'])
        b13_debit = rational(lower['B13_debit'])
        need(old['marks'] == got['marks'] == lower['marks']
             and joint >= 0 and a14_debit >= 0 and b13_debit >= 0
             and a14_retained == joint-a14_debit and a14_retained >= 0
             and rational(got['weight']) == rational(lower['actual_current_lower'])
             and rational(got['weight']) >= rational(lower['required_current_lower'])
             and rational(lower['required_current_lower']) ==
                 max(F(),a14_retained-b13_debit)
             and type(lower['whole_unmarked_witnesses']) is list,
             'same ordered D14 corner/current lower/debit/witness association')
        paid.append(dict(marks=got['marks'],current_lower=got['weight'],
            whole_unmarked_witnesses=lower['whole_unmarked_witnesses']))
    chosen_actual = one(paid,lambda p:p['marks'] == key['selected_marks'],'actual D14 selected corner')
    need(chosen['whole_witnesses'] and all(p['type'] in ('A1','A2')
         and p['nodes'] in chosen_actual['whole_unmarked_witnesses'] for p in chosen['whole_witnesses'])
         and rational(chosen_actual['current_lower']) >= F(1,2),
         'whole source A1/A2 witnesses at THIS selected guarded D14 corner')
    total = sum((rational(p['current_lower']) for p in paid),F())
    need(len({tuple(sorted(p['marks'])) for p in paid}) == 4
         and total == rational(actual['total']) == rational(join['current_total']) >= 1,
         'complete current D14 DISTINCT-corner sum pays THIS same face')
    for name in (QUARTERS,A14,B13):pin(name)
    return dict(source=D14,proof_selector=['selected_proofs',case['proof'],'variants',case['variant']],
        map_source=D14,map_selector=['selected_proofs',case['proof'],'affected_maps',mi],
        physical=physical,source_sealing=seal,selected_whole_guard=chosen,paying_corners=paid,
        actual_receipts=[evidence],proof_references=[
            pin(FINITE_PROOF),pin(FINITE_ROOT),pin(D14_REPORT),pin(GUARD_REPORT),pin(A14_REPORT)])


def bc12(key,units):
    original = source(BC12);row = original['B']['cases'][8];core = row['inherited_geometry']
    need(row['index'] == core['index'] == key['high_index'] == 1728
         and core['core_key'] == key['core_key'] == '(25, 9, ((9, 5),))'
         and core['base_nodes'] == key['base_nodes']
         and core['preferred_corner'] == key['selected_marks']
         and core['method'] == 'radical' and core['outcome'] == 'sealed_radical_core'
         and core['mixed_sign_shortcut'] == 3,
         'exact original BC12 high1728 full-core radical source')
    need(type(row['embeddings']) is list and row['embeddings'].count(key['induced_map']) == 1,
         'THIS complete BC12 native map in the original exact source domain')
    need(len(row['payment_variants']) == 1 and row['payment_variants'][0]['variant'] == 'sealed_payment',
         'sole exact source BC12 radical physical variant')
    variant = row['payment_variants'][0];prefix = 'BC12/selected/1728/sealed_payment'
    doc,evidence = unit(units,G4_WORKER,'order12-selected',G4_SHA,
        'PASS_G4_CURRENT_TERMINALS_DECLARED_SCOPE','read_inputs',(BC12,LITERAL,TEX))
    got = one(doc['physical_receipts'],lambda p:p['key'] == prefix,'exact actual G4 BC12 physical receipt')
    rows = [p for p in doc['rows'] if p['key'].startswith(prefix+'/corner/')]
    need(len(rows) == len(variant['corners']) == got['distinct_corners'] == 4
         and {p['key'] for p in rows} == {prefix+'/corner/'+str(i) for i in range(4)},
         'complete exact current G4 BC12 ordered physical corner domain')
    paid = []
    for ci,old in enumerate(variant['corners']):
        current = one(rows,lambda p:p['key'] == prefix+'/corner/'+str(ci),'one exact G4 BC12 corner')['current']
        need(current['marks'] == old['marks'] and current['guard_possible'] is False
             and rational(current['weight']) >= rational(old['new_weight']),
             'same actual BC12 source/current marked corner and whole marked-factor premise')
        paid.append(dict(source_selector=['B','cases',8,'payment_variants',0,'corners',ci],
            actual_path=prefix+'/corner/'+str(ci),marks=current['marks'],kind=current['kind'],
            current_lower=current['weight']))
    selected = one(rows,lambda p:p['current']['marks'] == key['selected_marks'],
        'same source BC12 selected affine corner')['current']
    need(selected['kind'] == 'ideal' and rational(selected['weight']) == 0
         and key['marked_image'] in selected['marked_images'],
         'BC12 selected intrinsic zero stays zero; only the same-face receipt pays')
    total = sum((rational(p['current_lower']) for p in paid),F())
    need(len({tuple(sorted(p['marks'])) for p in paid}) == 4
         and total == rational(got['current_total']) >= rational(got['required_total']) == 1,
         'complete actual G4 CURRENT distinct-corner same-face sum')
    return dict(source=BC12,proof_selector=['B','cases',8,'inherited_geometry'],
        map_source=BC12,map_selector=['B','cases',8,'embeddings',row['embeddings'].index(key['induced_map'])],
        physical=context(core['gram'],core['white']),source_sealing=dict(
            method=core['method'],radical=core['radical'],gauge_squared=core['gauge_squared'],
            allowed_contacts=core['allowed_contacts'],positive_contact_domain=core['positive_contact_domain'],
            negative_contact_domain=core['negative_contact_domain'],mixed_sign_shortcut=core['mixed_sign_shortcut']),
        paying_corners=paid,actual_receipts=[evidence],proof_references=[pin(BC_PROOF),pin(G4_ROOT)])


def b17_exclusion(key,W):
    original = source(HIGH)['isolation']
    case = original['cases'][512];core = original['cores'][key['core_key']]
    need(case['index'] == key['high_index'] == 512
         and case['core_key'] == key['core_key'] == '(22, 25, ((4, 9),))'
         and case['method'] == 'lorentz_residual'
         and case['component_order_bound'] == 17
         and case['marked_component'] == key['base_nodes'] ==
             [0,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17]
         and case['preferred_corner'] == key['selected_marks'] == [0,11],
         'exact affine B17 source high/core/ordered base/marks/bound')
    need(key['induced_map'] in B17_MAPS
         and set(range(17))-set(key['induced_map']) == {14}
         and key['induced_map'][key['base_nodes'].index(8)] == 13,
         'THIS exact saved B17 full map omits only14 and sends original8 to13')
    target = W.classical_graph('B',17,affine=True)
    need(target[13].get(14) == 4 and target[14] == {13:4},
         'frozen affine B17 missing terminal forces the native label4 entry')
    g,J,gauge = core['gram'],core['inverse'],core['gauge']
    original_black = set(case['preferred_corner']) | set(case['omitted_black'])
    physical = context(g,[i for i in range(18) if i not in original_black])
    need(len(g) == 18 and core['signature'] == [17,1,0] and gauge == [1]*18
         and all(rational(g[i][i]) == 2 for i in range(18))
         and type(J) is list and len(J) == 18
         and all(type(row) is list and len(row) == 18 for row in J),
         'exact original B17 norm-two/gauge/nondegenerate-Lorentz/inverse schema')
    end = [1,7,8];block = [[rational(J[i][j]) for j in end] for i in end]
    need(case['possible_contacts'] == [1,7,8,12]
         and case['end_sets'] == [[1,7,8],[12]]
         and case['omitted_black'] == [1,12]
         and set(key['base_nodes']) | set(case['omitted_black']) == set(range(18))
         and not set(key['base_nodes']) & set(case['omitted_black'])
         and block == [[F(8),F(6),F(3)],[F(6),F(6),F(3)],[F(3),F(3),F(2)]]
         and all(value > 0 for row in block for value in row)
         and rational(J[8][8]) == 2,
         'exact saved B17 complete contacts/ends/new-root exclusion and positive inverse block')
    lower = F(2)*gauge[8]*rational(J[8][8])
    need(lower == 4 and lower > 2, 'forced doubled-entry projection exceeds new-root norm two')
    return dict(application_kind='EXCLUDED_COMPLETION_UNDER_NO_S8',
        source=HIGH,proof_selector=['isolation','cases',512],map_source=B17_NOTE,
        map_selector=dict(format='markdown',locator='first keyed map' if B17_MAPS.index(key['induced_map']) == 0
            else 'second keyed map',literal_induced_map=key['induced_map']),
        physical=physical,source_sealing=dict(
            core_selector=['isolation','cores',key['core_key']],
            inverse_selector=['isolation','cores',key['core_key'],'inverse'],
            signature=core['signature'],gauge=gauge,possible_contacts=case['possible_contacts'],
            end_sets=case['end_sets'],omitted_black=case['omitted_black'],
            positive_end_nodes=end,positive_end_inverse_block=[[str(v) for v in row] for row in block],
            forced_entry=dict(missing_target_node=14,old_target_neighbor=13,original_node=8,
                native_label=4,physical_contact_square='2',inverse_diagonal=J[8][8],
                projection_lower=str(lower),new_root_norm_square='2',
                genuinely_new_supporting_white_not_an_original_root=True),
            complete_inverse_and_contact_end_predicates_runtime_regenerated=False),
        paying_corners=[],actual_receipts=[],
        source_branch_conclusion=dict(kind='THIS_EXACT_SELECTED_COMPLETION_IMPOSSIBLE_UNDER_NO_S8',
            reason='The mandatory genuinely new terminal has projection square>=4, exceeding its norm square2. Same-end additions increase the projection; other-end contact supplies the separate same-face S8 branch.'),
        proof_references=[pin(B17_NOTE),pin(FINITE_PROOF),pin(FINITE_ROOT)])


def apply_selected(key,W,units):
    full_key(key)
    family = None
    if not key['affine'] and key['target_kind'] == 'D':
        if key['target_order'] == 12 and key['high_index'] in D12_CASES:family = 'selected_D12_radical'
        elif key['target_order'] == 14 and key['high_index'] in D14_CASES:family = 'selected_D14_radical'
    elif key['affine'] and key['target_kind'] == 'B' and key['target_order'] == 12 and key['high_index'] == 1728:
        family = 'selected_BC12_radical'
    elif key['affine'] and key['target_kind'] == 'B' and key['target_order'] == 17 and key['high_index'] == 512:
        family = 'selected_B17_forced_terminal'
    if family is None:
        return None
    need(callable(SOURCE_READER) and callable(SOURCE_PINNER), 'both source readers must be injected')
    pin(TEX);pin(LITERAL)
    need(getattr(W,'TEX_SHA256',None) == TEX_SHA, 'frozen current literal W/canonical manuscript')
    if family == 'selected_D12_radical':
        need(key['marked_image'] in ([0,10],[0,11])
             and W.finite_same_weight(('D',12),*key['marked_image']) == F(1,2),
             'exact source D12 changed long-leaf role/current intrinsic half')
        bound = d12(key,units);intrinsic = '1/2'
    elif family == 'selected_D14_radical':
        need(key['marked_image'] in ([0,12],[0,13])
             and W.finite_same_weight(('D',14),*key['marked_image']) == 0,
             'exact source D14 long-leaf role/current unguarded intrinsic zero')
        bound = d14(key,units);intrinsic = '0'
    elif family == 'selected_BC12_radical':
        need(key['marked_image'] in ([9,10],[9,11])
             and W.affine_weight('B',12,*key['marked_image']) == 0,
             'exact source BC12 removed affine role/current intrinsic zero')
        bound = bc12(key,units);intrinsic = '0'
    else:
        need(key['marked_image'] in ([9,15],[9,16])
             and W.affine_weight('B',17,*key['marked_image']) == 0,
             'exact source B17 removed affine role/current intrinsic zero')
        bound = b17_exclusion(key,W);intrinsic = '0'
    frozen_key = json.loads(canonical(key))
    excluded = bound.get('application_kind') == 'EXCLUDED_COMPLETION_UNDER_NO_S8'
    return dict(status='BOUND_CONDITIONAL_WHOLE_FACE_APPLICATION',
        application_kind=bound.get('application_kind','EXACT_CURRENT_RECEIVER_SOURCE_BINDING'),
        theorem_id='finite57.source_application.'+family+
            ('.no_S8_forced_entry_exclusion.v1' if excluded else '.no_S8_current_receivers.v1'),
        application_key=dict(selected_full_key=frozen_key,source=bound['source'],
            source_proof_selector=bound['proof_selector'],source_map_selector=bound['map_selector']),
        rule_identity=dict(canonical_manuscript=pin(TEX),literal_source=pin(LITERAL)),
        physical_context=bound['physical'],source_sealing=bound['source_sealing'],
        source_references=[pin(bound['source']),pin(bound['map_source'])],
        proof_references=bound['proof_references'],actual_terminal_receipts=bound['actual_receipts'],
        source_branch_conclusion=bound.get('source_branch_conclusion',
            dict(kind='CURRENT_DISTINCT_CORNER_TOTAL_AT_LEAST_ONE_UNDER_NO_S8')),
        selected_whole_guard=bound.get('selected_whole_guard'),
        current_paying_corners=bound['paying_corners'],intrinsic_selected_current_lower=intrinsic,
        assumptions=[
            dict(kind='actual_same_face_source_occurrence_and_whole_selected_target',
                selected_full_key=frozen_key,source=bound['source'],selector=bound['proof_selector'],
                statement='The exact original high/core/ordered base/marks/full native map represents this SAME actual face and its stated whole finite/affine target.'),
            dict(kind='accepted_complete_source_numerical_and_geometric_predicates',runtime_regenerated=False,
                proof_references=bound['proof_references'],
                statement='Retain the exact source complete contact/end/native-map domain, physical gauge/kernel/inverse/signature, corner partitions and accepted source-specific full-core seal. Displayed components alone are not whole.'),
            dict(kind='NO_ADJACENT_WHITE_INTERIOR_PATH_AT_MOST_8_IN_THIS_SAME_FACE',
                runtime_verified=False,
                statement='This source branch pays under NO-S8. Independent conditional S8 is a separate premise and no short path is observed by this helper.')],
        provider_composition_premise=dict(
            status='EXACT_SOURCE_NO_S8_IMPLICATION_BOUND_TO_ACTUAL_CURRENT_RECEIVERS',
            same_full_key=frozen_key,
            implication=('For THIS exact B17 source/key and complete physical/native/contact/end/inverse predicates: NO same-face adjacent white-interior path<=8 excludes this actual completion.' if excluded else
                'For THIS source occurrence and complete physical/native/contact predicates: NO same-face adjacent white-interior path<=8 implies CURRENT total of DISTINCT corners>=1.'),
            independent_S8_composition='If this same face were unpaid, independent conditional S8 forbids every such short path. The accepted exhaustive source implication then contradicts unpaidness.',
            requiresS8_flag_alone_is_not_a_disjunction=True),
        conclusion=dict(quantity='total_current_charge_of_same_actual_face',lower='1',
            route='unpaid-face NO-S8 and impossible exact completion' if excluded else 'actual current distinct-corner receivers under NO-S8'),
        intrinsic_corner_unit=False,theorem_execution=False,
        selected_provider_union_complete=False,current_terminal_bridge_complete=False,
        full_transitive_regeneration=False,geometric_formalization=False,global_finite57_verified=False)


if __name__ == '__main__':
    raise SystemExit('Application helper only: root must load pinned definitions into the final consumer')
