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
SOURCE_ACCEPTANCE = None
APPLICATION_SOURCE_SHA = None
APPLICATION_SOURCE = 'scripts/finite57_selected_legacy_receipt_applications_20261025.py'
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
GENERIC_NOTE = 'results/finite57_generic_forced_entry_projection_source_assessment_20261025.md'
EXTENSION_PROOF = 'results/audit_finite57_extension_schemas_and_weight_binding_fresh_20261024.md'
EXTENSION_ROOT = 'results/review_finite57_literal_rule_and_extension_logic_root_20261024.md'
A13_SOURCE = 'results/fresh_a13_endpoint_a1_context_retention_20261018.json'
A15_SOURCE = 'results/fresh_rank15_sparse_payment_exchange_20261018.json'
D_LONG_SOURCE = 'results/fresh_finite_d16_d17_long_pair_half_transfer_20261017.json'
LATE_PROOF = 'results/audit_finite57_current_bounded_high_and_final_exchanges_20261024.md'
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
    GENERIC_NOTE:'c9a553ba04b42c3729ec7d79b8e8e60264d35ebfe83545223eacf2242baab44f',
    EXTENSION_PROOF:'1e6a165b957069b4d618a17361900aa0a92ee57f52bc7f7a22362ef05e57a847',
    EXTENSION_ROOT:'6ea02657eea8cecc37f37b8e024ffbd0be700b98e3bcb5926824e1b80cca4cc8',
    A13_SOURCE:'0512c3ef32c183779895876c5cb9a2164118f3abb33dfc441b14c78219e1a4c8',
    A15_SOURCE:'50d7949f65aa549a90048b5cd1badf57572a7065746ff3f93770d22aec4d2770',
    D_LONG_SOURCE:'c883a7afadcaeed4340ebcdaa160c9058da223f7a39fd27597e2b7e51c872428',
    LATE_PROOF:'602da0b68727129020708c11e8c19cbd6721fb2745d9ba98c45cba90aab3dec9',
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


def accepted_premises(ids):
    need(type(SOURCE_ACCEPTANCE) is dict and set(SOURCE_ACCEPTANCE) == {'path','sha256'},
         'mandatory separately curated root source-application acceptance identity')
    path = SOURCE_ACCEPTANCE['path'];sha = digest(SOURCE_ACCEPTANCE['sha256'])
    need(type(path) is str and path.startswith('results/'), 'root acceptance repository path')
    own_sha = digest(APPLICATION_SOURCE_SHA)
    pin(APPLICATION_SOURCE,own_sha);reference = pin(path,sha)
    cache_key = (path,sha)
    if cache_key not in CACHE:
        CACHE[cache_key] = SOURCE_READER(path,sha)
    doc = CACHE[cache_key]
    need(type(doc) is dict and doc['status'] == 'CURRENT_SOURCE_APPLICATION_PREMISES_ROOT_ACCEPTED'
         and type(doc['accepted_application_sources']) is dict
         and doc['accepted_application_sources'].get(APPLICATION_SOURCE) == own_sha,
         'root acceptance binds THIS exact helper source, not an ordinary-premise broadening')
    accepted = doc['accepted_supplemental_premise_ids']
    need(type(accepted) is list and all(type(s) is str for s in accepted)
         and len(accepted) == len(set(accepted)) and set(ids) <= set(accepted),
         'every relevant namespaced supplemental source/geometry premise explicitly root-accepted')
    return reference


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



CORE_CACHE = {}
GENERIC_PRODUCTS = 0
MAX_GENERIC_PRODUCTS = 2_000_000


def at(value,selector):
    for part in selector:
        value = value[part]
    return value


def authoritative(key):
    doc = source(HIGH)['isolation'];p = doc['cases'][key['high_index']]
    need(p['index'] == key['high_index'] and p['core_key'] == key['core_key']
         and p['marked_component'] == key['base_nodes']
         and p['preferred_corner'] == key['selected_marks']
         and type(p.get('component_order_bound')) is int
         and len(key['base_nodes']) <= key['target_order'] <= p['component_order_bound'],
         'same authoritative HIGH/core/ordered base/marks/whole-order bound')
    return p,doc['cores'][key['core_key']]


def opposite_pairs(p,n):
    left,right = p['left'][1],p['right'][1];mapping = p['right_copy_map']
    need(type(left) is list and type(right) is list and len(left) == len(right) == 2,
         'complete original opposite-pair selectors')
    if type(mapping) is list:
        values = mapping
    else:
        need(type(mapping) is dict and mapping, 'original copy-map list or canonical dictionary')
        keys = []
        for k in mapping:
            need(type(k) is int or (type(k) is str and k.isdecimal() and str(int(k)) == k),
                 'canonical original copy-map node')
            keys.append(int(k))
        need(len(keys) == len(set(keys)) and set(keys) == set(range(len(keys))),
             'complete integer original copy-map domain')
        values = [mapping[i] if i in mapping else mapping[str(i)] for i in range(len(keys))]
    need(0 < len(values) <= n and all(type(v) is int and 0 <= v < n for v in values)
         and len(values) == len(set(values))
         and all(type(v) is int and 0 <= v < len(values) for v in right),
         'original right-node to actual-core injection')
    u,v = left,[values[i] for i in right]
    need(all(type(i) is int and 0 <= i < n for i in u)
         and len(set(u+v)) == 4 and p['preferred_corner'][0] in u and p['preferred_corner'][1] in v,
         'four distinct original black roots and same ordered actual corner')
    return u,v


def matrix_values(gram):
    need(type(gram) is list and 0 < len(gram) <= 64
         and all(type(r) is list and len(r) == len(gram) for r in gram),
         'bounded complete rational source matrix')
    return [[rational(v) for v in row] for row in gram]


def source_physical(key,gram,white,unit_gauge=False):
    p,core = authoritative(key);g = matrix_values(gram);old = matrix_values(core['gram']);n = len(old)
    need(len(g) == n and g == old and len(core['gauge']) == n,
         'exact original physical Gram/source identity')
    scales = [rational(v) for v in core['gauge']]
    need(all(s > 0 for s in scales) and all(g[i][i] == 2*scales[i] for i in range(n))
         and (not unit_gauge or all(s == 1 for s in scales)),
         'actual positive-diagonal physical gauge, not a native label relaxation')
    u,v = opposite_pairs(p,n)
    need(white == sorted(set(range(n))-set(u+v)), 'same complete original white support')
    return p,core,context(gram,white)


def graph_distances(g,starts):
    d = {i:0 for i in starts};todo = list(starts)
    for i in todo:
        for j in range(len(g)):
            if i != j and g[i][j] and j not in d:
                d[j] = d[i]+1;todo.append(j)
    return d


def generic_forced_entry(key,W):
    """One GIVEN full map; no native-map enumeration or new-root search."""
    global GENERIC_PRODUCTS
    p,core = authoritative(key)
    if len(key['base_nodes']) >= key['target_order'] or p['method'] not in (
            'lorentz_residual','relevant_component','marked_component_isolated_by_end_inverses'):
        return None
    n = len(core['gram'])
    if core.get('inverse') is None or core['signature'] != [n-1,1,0]:
        return None
    if key['core_key'] not in CORE_CACHE:
        g,J = matrix_values(core['gram']),matrix_values(core['inverse'])
        scales = [rational(v) for v in core['gauge']]
        if len(J) != n or len(scales) != n or not all(s > 0 for s in scales):
            return None
        if not all(g[i][i] == 2*scales[i] and
                   all(g[i][j] == g[j][i] <= 0 and J[i][j] == J[j][i]
                       for j in range(i)) for i in range(n)):
            return None
        for i in range(n):
            for j in range(n):
                total = F()
                for k in range(n):
                    if g[i][k]:
                        GENERIC_PRODUCTS += 1
                        need(GENERIC_PRODUCTS <= MAX_GENERIC_PRODUCTS,
                             'complete cached GJ identity budget; refuse, never truncate')
                        total += g[i][k]*J[k][j]
                if total != int(i == j):
                    return None
        CORE_CACHE[key['core_key']] = g,J,scales
    g,J,scales = CORE_CACHE[key['core_key']]
    u,v = opposite_pairs(p,n);du,dv = graph_distances(g,u),graph_distances(g,v)
    if len(du) != n or len(dv) != n or min(du[i] for i in v) <= 8:
        return None
    allowed = p['possible_contacts'];ends = p.get('end_sets')
    if (type(allowed) is not list or len(allowed) != len(set(allowed))
            or any(type(i) is not int or not 0 <= i < n for i in allowed)
            or type(ends) is not list or len(ends) != 2
            or any(type(e) is not list or len(e) != len(set(e)) for e in ends)):
        return None
    left,right = map(set,ends)
    if left & right or left | right != set(allowed):
        return None
    splits = [k for k in range(7) if all(du[i] <= k for i in left)
              and all(dv[i] <= 6-k for i in right)]
    if not splits or any(J[i][j] < 0 for end in ends for i in end for j in end):
        return None
    base,mp = key['base_nodes'],key['induced_map'];marks = key['selected_marks']
    actual_corner = (set(range(n))-set(u+v)) | set(marks)
    if (not set(base) <= actual_corner or p['omitted_black'] != sorted(set(u+v)-set(marks))
            or any(g[i][j] for i in actual_corner-set(base) for j in base)):
        return None
    target = W.classical_graph(key['target_kind'],key['target_order'],affine=key['affine'])
    for ai,i in enumerate(base):
        for aj,j in enumerate(base[:ai]):
            label = target[mp[ai]].get(mp[aj],0)
            if ((label == 0 and g[i][j] != 0) or
                    (label in (3,4) and (g[i][j] >= 0 or
                     g[i][j]**2 != (1 if label == 3 else 2)*scales[i]*scales[j])) or
                    label not in (0,3,4)):
                return None
    image = dict(zip(mp,base))
    for new in sorted(set(range(key['target_order']))-set(mp)):
        forced = [(image[old],label) for old,label in target[new].items() if old in image]
        if not forced or any(label not in (3,4) for _,label in forced):
            continue
        nodes = {i for i,_ in forced};reason = None;lower = None;cross = None
        if not nodes <= set(allowed):
            reason = 'required native neighbor outside complete necessary-contact domain'
        elif nodes & left and nodes & right:
            reason = 'required neighbors meet both original ends, giving same-face S8'
        else:
            end = left if nodes <= left else right
            if not nodes <= end:
                continue
            lower = sum(((1 if label == 3 else 2)*scales[i]*J[i][i] for i,label in forced),F())
            if lower > 2:
                reason = 'mandatory boundary-root projection lower exceeds norm two'
            elif lower == 2:
                cross = next(([i,j,str(J[i][j])] for i in nodes for j in nodes if i < j and J[i][j] > 0),None)
                if cross:
                    reason = 'diagonal projection reaches norm two and a forced positive cross term is strict'
        if reason is None:
            continue
        _,_,physical = source_physical(key,core['gram'],sorted(set(range(n))-set(u+v)))
        return dict(application_kind='EXCLUDED_COMPLETION_UNDER_NO_S8',source=HIGH,
            proof_selector=['isolation','cases',key['high_index']],map_source=GENERIC_NOTE,
            map_selector=dict(format='given_full_native_map_verified_against_frozen_W',induced_map=mp),
            physical=physical,source_sealing=dict(
                kind='strict_generic_proper_boundary_root',method=p['method'],
                core_selector=['isolation','cores',key['core_key']],source_gauge=core['gauge'],
                source_signature=core['signature'],inverse_identity_checked_in_this_call=True,
                original_opposite_pairs=[u,v],complete_possible_contacts=allowed,end_sets=ends,
                distance_split=splits[0],full_end_inverse_blocks=[
                    [[str(J[i][j]) for j in e] for i in e] for e in ends],
                missing_target_node=new,forced_original_neighbors=[
                    dict(node=i,native_label=label,contact_square_lower=str((1 if label == 3 else 2)*scales[i]))
                    for i,label in forced],projection_diagonal_lower=None if lower is None else str(lower),
                forced_strict_cross=cross,genuinely_new_supporting_white_not_an_original_root=True),
            paying_corners=[],actual_receipts=[],source_branch_conclusion=dict(
                kind='THIS_EXACT_SELECTED_COMPLETION_IMPOSSIBLE_UNDER_NO_S8',reason=reason),
            proof_references=[pin(GENERIC_NOTE),pin(EXTENSION_PROOF),pin(EXTENSION_ROOT),
                              pin(FINITE_PROOF),pin(FINITE_ROOT)])
    return None


def bound_exclusion(key,W,proof,record,src,selector,map_selector):
    need(rational(record['projection_lower_bound']) > 2,
         'exact source forced-contact projection exceeds norm two')
    certificate = generic_forced_entry(key,W)
    need(certificate is not None, 'named proper source map needs its complete physical/end/forced-entry certificate')
    certificate['source'] = src;certificate['proof_selector'] = selector
    certificate['map_source'] = src;certificate['map_selector'] = map_selector
    certificate['source_sealing']['named_source_projection_record'] = record
    return certificate


def d12_current_row(doc,source_doc,path,gram_selector,white_selector):
    old = at(source_doc,path);pointer = '/'+('/'.join(map(str,path)))
    actual = one(doc['rows'],lambda p:p['path'] == pointer,'one exact actual D12 terminal path')
    need(actual['marks'] == old['marks'] and
         actual['gram_path'] == '/'+('/'.join(map(str,gram_selector))) and
         actual['white_path'] == '/'+('/'.join(map(str,white_selector))) and
         actual['possible_negative_D13_guard'] is False and
         rational(actual['chosen']) == rational(old['payment']) and
         rational(actual['current_lower']) >= rational(actual['chosen']),
         'same source D12 physical corner/current guarded-safe terminal binding')
    return dict(source_selector=path,actual_path=pointer,marks=actual['marks'],
                current_lower=actual['current_lower'])


def d12(key,W,units):
    domain = source(D12_DOMAIN);first = source(D12)
    saved = one(domain['selected_patterns'],lambda p:p['index'] == key['high_index'],
                'complete source D12 selected-domain row')
    si = domain['selected_patterns'].index(saved);mi = member(key,saved['affected_maps'])
    need(any(p['images'] == key['induced_map'] and p['marked_images'] == key['marked_image']
             for p in saved['maps']), 'same complete saved D12 native map')
    candidates = [(D12,['proofs',i],p,'d12-first') for i,p in enumerate(first['proofs']) if p['index'] == key['high_index']]
    candidates += [(D12_DOMAIN,['new_proofs',i],p,'d12-rest')
                   for i,p in enumerate(domain['new_proofs']) if p['index'] == key['high_index']]
    need(len(candidates) == 1, 'exact disjoint first/rest D12 source-proof partition')
    src,selector,p,batch = candidates[0];doc_source = first if src == D12 else domain
    need(p['core_key'] == key['core_key'], 'D12 same source core')
    source_physical(key,p['gram'],p['white'])
    if p['method'] == 'proper_completion_projection_exclusion':
        ri = member(key,p['maps']);record = p['maps'][ri]
        need(record['outcome'] == 'proper_completion_exceeds_norm2' and record['required_label'] == 3,
             'exact source proper D12 entry exclusion')
        return bound_exclusion(key,W,p,record,src,selector,selector+['maps',ri])
    need(len(key['base_nodes']) == 12, 'whole source D12 target')
    doc,evidence = unit(units,D12_WORKER,batch,D12_SHA,
        'PASS_GUARD_SAFE_CURRENT_D12_TERMINAL_LITERAL_BINDING','inputs',(src,LITERAL,TEX))
    branches = []
    selected_half = dict(marks=key['selected_marks'],current_lower='1/2',
        source_selector=['isolation','cases',key['high_index'],'preferred_corner'],
        basis='exact current W on THIS whole D12 key, not a fabricated terminal receipt')
    def branch(condition,corner_paths,gram_path,white_path,include_selected):
        paid = [d12_current_row(doc,doc_source,q,gram_path,white_path) for q in corner_paths]
        if include_selected:
            need(all(set(q['marks']) != set(key['selected_marks']) for q in paid),
                 'distinct protected alternative, never count selected corner twice')
            paid = [selected_half]+paid
        need(paid and len({tuple(sorted(q['marks'])) for q in paid}) == len(paid)
             and sum((rational(q['current_lower']) for q in paid),F()) >= 1,
             'EVERY complete D12 source branch has a distinct current same-face receipt')
        branches.append(dict(condition=condition,gram_selector=gram_path,white_selector=white_path,
            physical=context(at(doc_source,gram_path),at(doc_source,white_path)),
            paying_corners=paid,current_total_lower=str(sum((rational(q['current_lower']) for q in paid),F()))))
    method = p['method']
    if method == 'radical_annihilation':
        need(p['radical'] and type(p['remaining_contact_positions']) is list
             and type(p['mixed_sign_shortcuts']) is list
             and all(0 < q['length'] <= 8 for q in p['mixed_sign_shortcuts']),
             'complete source radical/contact/shortpath predicates')
        branch(dict(kind='all further whites core-orthogonal under NO-S8'),
               [selector+['corners',i] for i,c in enumerate(p['corners']) if rational(c['payment']) > 0],
               selector+['gram'],selector+['white'],False)
    elif method == 'saturated_lorentz_extensions':
        need([b['extra_contact'] for b in p['branches']] == [None]+p['saturated_positions']
             and p['orthogonal_whites_cannot_attach_saturated_root'] is True,
             'exact exhaustive D12 no-extra/unique-saturated-root branch domain')
        for bi,b in enumerate(p['branches']):
            bs = selector+['branches',bi]
            branch(dict(kind='no extra root or unique saturated root',extra_contact=b['extra_contact']),
                [bs+['corners',i] for i,c in enumerate(b['corners']) if rational(c['payment']) > 0],
                bs+['gram'],bs+['white'],False)
    elif method == 'protected_distinct_corner':
        need(p['selected_marks'] == key['selected_marks'] and p['selected_factor'] == key['base_nodes']
             and p['all_possible_contact_positions_disjoint_from_paying_factors'] is True
             and all(not set(r['nodes']) & set(p['possible_single_contacts']) for r in p['paying_corner']['roles']),
             'complete contact set disjoint from THIS paying factor')
        branch(dict(kind='arbitrary unprotected extensions; paying factors have no first contact'),
               [selector+['paying_corner']],selector+['gram'],selector+['white'],True)
    elif method == 'one_saturated_contact_at_paying_factor':
        need(p['selected_marks'] == key['selected_marks'] and p['selected_factor'] == key['base_nodes']
             and p['before']['marks'] == p['after']['marks']
             and rational(p['new_root_projection_square']) == 2
             and rational(p['new_root_residual_square']) == 0
             and p['core_orthogonal_whites_are_orthogonal_to_new_root'] is True,
             'complete source unique saturated paying-factor alternative')
        for name,gname,wname in (('before','gram','white'),('after','extended_gram','extended_white')):
            branch(dict(kind='no contact at q' if name == 'before' else 'unique saturated root at q',q=p['contact']),
                   [selector+[name]],selector+[gname],selector+[wname],True)
    elif method == 'both_ends_shared_unit_residual':
        need([b['active_end'] for b in p['one_sided_branches']] == [0,1]
             and p['selected_marks'] == key['selected_marks'] and p['selected_factor'] == key['base_nodes']
             and p['residual_Gram'] == [['1','1'],['1','1']]
             and p['all_other_whites_orthogonal_to_extended_core_in_no_short8_branch'] is True
             and 0 < p['common_tail_shortcut_length'] <= 8,
             'exact source complete one-sided/shared-residual full-star alternatives')
        for bi,b in enumerate(p['one_sided_branches']):
            branch(dict(kind='only this end active (including neither)',active_end=b['active_end']),
                [selector+['one_sided_branches',bi,'paying_corner']],selector+['gram'],selector+['white'],True)
        branch(dict(kind='both ends active; unique shared unit residual; all other whites orthogonal'),
            [selector+['extended_corners',i] for i,c in enumerate(p['extended_corners']) if rational(c['payment']) > 0],
            selector+['extended_gram'],selector+['extended_white'],False)
    else:
        raise ValueError('unsupported exact D12 source-proof method')
    return dict(source=src,proof_selector=selector,map_source=D12_DOMAIN,
        map_selector=['selected_patterns',si,'affected_maps',mi],
        physical=context(p['gram'],p['white']),
        source_sealing=dict(kind='accepted_complete_source_D12_method',method=method,proof_selector=selector,
            all_source_specific_contact_matrix_branch_predicates_retained=True),
        paying_corners=branches[0]['paying_corners'] if len(branches) == 1 else [],
        receipt_branches=branches,actual_receipts=[evidence],
        consumer_literal_evaluations=[dict(function='finite_same_weight',type=['D',12],
            marked_image=key['marked_image'],result='1/2',used_only_with_distinct_source_alternative=True)],
        proof_references=[pin(FINITE_PROOF),pin(FINITE_ROOT)])

def d14(key,W,units):
    original = source(D14);guards = source(GUARD)
    matches = [(i,p) for i,p in enumerate(original['selected_proofs']) if p['index'] == key['high_index']]
    if not matches:
        return None
    need(len(matches) == 1,'unique exact complete D14 selected source proof')
    pi,proof = matches[0]
    need(proof['core_key'] == key['core_key'],'same D14 source core')
    mi = member(key,proof['affected_maps'] if 'affected_maps' in proof else proof['maps'])
    source_physical(key,proof['gram'],proof['white'])
    if proof['method'] == 'proper_completion_exclusion':
        record = proof['maps'][member(key,proof['maps'])]
        need(record['outcome'] == 'arbitrary_contact_projection_exceeds_norm2',
             'exact original proper D14 boundary-root exclusion')
        return bound_exclusion(key,W,proof,record,D14,['selected_proofs',pi],['selected_proofs',pi,'maps',mi])
    need(proof['method'] == 'whole_D14_transfer' and len(key['base_nodes']) == 14
         and proof['sealing']['method'] in ('radical_annihilation','complete_saturated_lorentz_extensions'),
         'complete accepted source-specific D14 full-core method')
    doc,evidence = unit(units,D14_WORKER,'d14-current24',D14_SHA,D14_STATUS,
        'consumed_inputs',(D14,GUARD,QUARTERS,A14,B13,LITERAL,TEX,D14_REPORT,GUARD_REPORT,A14_REPORT))
    result = doc['result']
    expected = [[p['index'],b['extra_contact']] for p in original['selected_proofs'] for b in p.get('variants',[])]
    need(len(expected) == 24 and result['domain'] == expected
         and result['complete_all24_at_declared_premises'] is True,
         'complete EXACT source physical24 domain, not count alone')
    need([p['index'] for p in result['receipts']] == [[i,canonical(x)] for i,x in expected]
         and all(p['source'] == D14[8:-5] for p in result['receipts'])
         and [[p['index'],p['extra_contact']] for p in result['source_lower_joins']] == expected,
         'complete ordered actual24 receiver and latest-source lower-join identities')
    branches = []
    selected_guards = []
    for vi,v in enumerate(proof['variants']):
        extra = v['extra_contact'];selector = ['selected_proofs',pi,'variants',vi]
        guard = one(guards['selected_whole_variants'],
            lambda p:p['index'] == key['high_index'] and p['extra_contact'] == extra,
            'same source D14 physical variant guard')
        chosen = one(guard['guarded_outputs'],lambda p:p['marks'] == key['selected_marks'],
                     'same selected D14 whole-small guarded output')
        selected = one(v['corners'],lambda p:p['marks'] == key['selected_marks'],
                       'same selected D14 physical source corner')
        factor = one(selected['components'],
            lambda p:p['type'] == 'D14' and set(key['selected_marks']) <= set(p['nodes']),
            'THIS source selected D14 whole factor')
        need(factor['nodes'] == key['base_nodes'] and key['marked_image'] in chosen['images']
             and rational(chosen['retained_half']) == F(1,2),
             'same D14 base/ordered marks/full native image and guarded half')
        actual = one(result['receipts'],
            lambda p:p['index'] == [key['high_index'],canonical(extra)],'one exact current D14 variant')
        join = one(result['source_lower_joins'],
            lambda p:p['index'] == key['high_index'] and p['extra_contact'] == extra,
            'one exact D14 latest A14/B13/rank15 lower join')
        physical = context(v['gram'],v['white']);seal = join['full_core_seal_premise']
        expected_seal = dict(kind='accepted_source_specific_full_core_seal',runtime_verified=False,
                source=GUARD_REPORT,original_source=D14_REPORT,source_index=key['high_index'],
                source_proof_ordinal=pi,source_variant_ordinal=vi,extra_contact=extra,
                exact_saved_sealing=proof['sealing'],independent_S8_runtime_verified_by_this_batch=False)
        need(join['source_proof_ordinal'] == pi and join['source_variant_ordinal'] == vi
             and type(seal) is dict and all(seal.get(k) == value for k,value in expected_seal.items())
             and seal['section'] == 'Complete selected outputs'
             and type(seal['premise']) is str and bool(seal['premise'])
             and type(seal['conclusion']) is str and bool(seal['conclusion']),
             'exact conditional FULL-core seal and variant selector')
        need(actual['context_sha256'] == physical['context_sha256']
             and actual['direct_current_receiver_condition'] is True
             and actual['legacy_pointwise_bounds_used'] is False
             and actual['provenance']['source_selector'] == selector
             and actual['provenance']['full_core_seal'] == seal
             and len(actual['corners']) == len(join['per_corner']) == len(v['corners']) == 4,
             'actual D14 source/context/four-corner DIRECT binding')
        paid = []
        for old,got,lower in zip(v['corners'],actual['corners'],join['per_corner'],strict=True):
            joint = rational(lower['rank15_joint_lower'])
            debit = rational(lower['A14_endpoint_debit']);retained = rational(lower['A14_retained_lower'])
            b13 = rational(lower['B13_debit']);required = rational(lower['required_current_lower'])
            need(old['marks'] == got['marks'] == lower['marks']
                 and joint >= 0 and debit >= 0 and b13 >= 0 and retained == joint-debit >= 0
                 and rational(got['weight']) == rational(lower['actual_current_lower']) >= required
                 and required == max(F(),retained-b13)
                 and type(lower['whole_unmarked_witnesses']) is list,
                 'same physical D14 current corner with all later A14/B13 debits')
            paid.append(dict(marks=got['marks'],current_lower=got['weight'],
                rank15_joint_lower=lower['rank15_joint_lower'],A14_endpoint_debit=lower['A14_endpoint_debit'],
                A14_retained_lower=lower['A14_retained_lower'],B13_debit=lower['B13_debit'],
                whole_unmarked_witnesses=lower['whole_unmarked_witnesses']))
        selected_actual = one(paid,lambda p:p['marks'] == key['selected_marks'],'actual selected D14 corner')
        need(chosen['whole_witnesses'] and all(p['type'] in ('A1','A2')
             and p['nodes'] in selected_actual['whole_unmarked_witnesses'] for p in chosen['whole_witnesses'])
             and rational(selected_actual['current_lower']) >= F(1,2),
             'same guarded selected whole A1/A2 witness, not a displayed component')
        total = sum((rational(p['current_lower']) for p in paid),F())
        need(len({tuple(sorted(p['marks'])) for p in paid}) == 4
             and total == rational(actual['total']) == rational(join['current_total']) >= 1,
             'EVERY exclusive D14 physical variant has a current distinct-corner unit')
        selected_guards.append(chosen)
        branches.append(dict(condition=dict(extra_contact=extra,source_variant_ordinal=vi),
            physical=physical,source_selector=selector,source_full_core_seal=seal,
            paying_corners=paid,current_total_lower=str(total)))
    need(branches,'nonempty complete source D14 variant domain')
    return dict(source=D14,proof_selector=['selected_proofs',pi],map_source=D14,
        map_selector=['selected_proofs',pi,'affected_maps',mi],physical=context(proof['gram'],proof['white']),
        source_sealing=proof['sealing'],selected_whole_guard=selected_guards,
        paying_corners=branches[0]['paying_corners'] if len(branches) == 1 else [],
        receipt_branches=branches,actual_receipts=[evidence],
        proof_references=[pin(FINITE_PROOF),pin(FINITE_ROOT),pin(D14_REPORT),pin(GUARD_REPORT),pin(A14_REPORT)])


def bc12(key,W,units):
    original = source(BC12)
    candidates = [(i,p) for i,p in enumerate(original['B']['cases']) if p['index'] == key['high_index']]
    if not candidates:
        return None
    need(len(candidates) == 1,'unique complete order12 B/C source occurrence')
    ri,row = candidates[0];core = row['inherited_geometry']
    need(core['index'] == key['high_index'] and core['core_key'] == key['core_key']
         and core['base_nodes'] == key['base_nodes'] and core['preferred_corner'] == key['selected_marks']
         and core['method'] in ('radical','lorentz_residual'),
         'exact original order12 B/C physical source predicates')
    source_physical(key,core['gram'],core['white'])
    if key['target_kind'] == 'B':
        need(row['embeddings'].count(key['induced_map']) == 1,
             'THIS complete B12 native map in the original complete source domain')
        map_selector = ['B','cases',ri,'embeddings',row['embeddings'].index(key['induced_map'])]
    else:
        ci,crow = one([[i,p] for i,p in enumerate(original['C']['cases'])],
            lambda t:t[1]['index'] == key['high_index'],'exact complete C12 selected source occurrence')
        mi,mp = one([[i,p] for i,p in enumerate(crow['maps'])],
            lambda t:t[1]['map'] == key['induced_map'],'exact C12 full native map')
        need(crow['end_projection_certificate'] == core
             and rational(crow['retained_finite_unit']) == 1 and not row['payment_variants']
             and len(mp['forced_entries']) == 1
             and all(p['native_label'] == 4 and rational(p['forced_projection_lower_bound']) > 2
                     for p in mp['forced_entries']),
             'C12 missing terminal uses THIS same sealed finite B11/projection proof')
        map_selector = ['C','cases',ci,'maps',mi]
    doc,evidence = unit(units,G4_WORKER,'order12-selected',G4_SHA,
        'PASS_G4_CURRENT_TERMINALS_DECLARED_SCOPE','read_inputs',(BC12,LITERAL,TEX))
    branches = []
    if not row['payment_variants']:
        need(core['outcome'] == 'whole_finite_B11_before_any_affine_completion',
             'exact pre-affine whole finite B11 source alternative')
        name = 'BC12/selected/'+str(key['high_index'])
        actual = one(doc['rows'],lambda p:p['key'] == name,'same actual whole finite B11 corner')['current']
        need(actual['marks'] == key['selected_marks'] and actual['kind'] == 'finite'
             and actual['guard_possible'] is False and rational(actual['weight']) >= 1,
             'actual same-face pre-affine B11 unit, not a fictitious affine intrinsic unit')
        branches.append(dict(condition=dict(kind='the selected component is already whole finite B11'),
            physical=context(core['gram'],core['white']),paying_corners=[dict(
                actual_path=name,marks=actual['marks'],current_lower=actual['weight'])],
            current_total_lower=actual['weight']))
    else:
        expected = [name for name in ('sealed_payment','no_extra_root_payment','one_extra_root_payment') if name in core]
        need([p['variant'] for p in row['payment_variants']] == expected
             and core['outcome'] in ('sealed_radical_core','zero_or_one_sealed_unit_contact_root')
             and key['target_kind'] == 'B',
             'complete source radical or zero/one saturated extra-root branch domain')
        for vi,variant in enumerate(row['payment_variants']):
            name = variant['variant'];prefix = 'BC12/selected/'+str(key['high_index'])+'/'+name
            actual = one(doc['physical_receipts'],lambda p:p['key'] == prefix,
                         'same actual G4 B12 physical receiver')
            rows = [p for p in doc['rows'] if p['key'].startswith(prefix+'/corner/')]
            need(len(rows) == len(variant['corners']) == actual['distinct_corners'] == 4
                 and {p['key'] for p in rows} == {prefix+'/corner/'+str(i) for i in range(4)},
                 'complete current G4 B12 exact ordered corner domain')
            paid = []
            for ci,old in enumerate(variant['corners']):
                current = one(rows,lambda p:p['key'] == prefix+'/corner/'+str(ci),'one G4 B12 source corner')['current']
                need(current['marks'] == old['marks'] and current['guard_possible'] is False
                     and rational(current['weight']) >= rational(old['new_weight']),
                     'same physical B12 current corner and complete source whole-factor premise')
                paid.append(dict(source_selector=['B','cases',ri,'payment_variants',vi,'corners',ci],
                    actual_path=prefix+'/corner/'+str(ci),marks=current['marks'],current_lower=current['weight']))
            selected = one(rows,lambda p:p['current']['marks'] == key['selected_marks'],
                           'same B12 selected physical affine corner')['current']
            need(selected['kind'] == 'ideal' and rational(selected['weight']) == 0
                 and key['marked_image'] in selected['marked_images'],
                 'THIS selected affine zero remains zero')
            total = sum((rational(p['current_lower']) for p in paid),F())
            need(len({tuple(sorted(p['marks'])) for p in paid}) == 4
                 and total == rational(actual['current_total']) >= rational(actual['required_total']) == 1,
                 'EVERY source B12 branch has a current distinct-corner unit')
            extended = name == 'one_extra_root_payment'
            branches.append(dict(condition=dict(variant=name),
                physical=context(core['one_extra_root_gram'] if extended else core['gram'],
                                 core['white']+[len(core['gram'])] if extended else core['white']),
                paying_corners=paid,current_total_lower=str(total)))
    seal_selector = ['B','cases',ri,'inherited_geometry']
    if core['method'] == 'radical':
        seal = dict(kind='accepted_complete_order12_radical_seal',method=core['method'],
            outcome=core['outcome'],source_selector=seal_selector,gauge_squared=core['gauge_squared'],
            complete_allowed_contacts=core['allowed_contacts'],radical=core['radical'],
            positive_contact_domain=core['positive_contact_domain'],
            negative_contact_domain=core['negative_contact_domain'],mixed_sign_shortcut=core['mixed_sign_shortcut'])
    else:
        seal = dict(kind='accepted_complete_order12_saturated_lorentz_alternatives',method=core['method'],
            outcome=core['outcome'],source_selector=seal_selector,gauge_squared=core['gauge_squared'],
            complete_allowed_contacts=core['allowed_contacts'],complete_end_sets=core['end_sets'],
            both_ends_shortcut=core['both_ends_shortcut'],
            complete_normalized_singleton_projections=core['normalized_singleton_projections'],
            complete_same_end_pair_exclusions=core['same_end_pair_exclusions'],
            unique_possible_contact=core['unique_possible_contact'],contact_native_label=core['contact_native_label'],
            contact_normalized_projection=core['contact_normalized_projection'],contact_residual=core['contact_residual'])
    return dict(source=BC12,proof_selector=seal_selector,
        map_source=BC12,map_selector=map_selector,physical=context(core['gram'],core['white']),
        source_sealing=seal,paying_corners=branches[0]['paying_corners'] if len(branches) == 1 else [],
        receipt_branches=branches,actual_receipts=[evidence],
        proof_references=[pin(BC_PROOF),pin(G4_ROOT),pin(EXTENSION_PROOF),pin(EXTENSION_ROOT)])


def component(parts,kind,nodes,marks=()):
    return one(parts,lambda p:p['type'] == kind and p['nodes'] == nodes
               and set(marks) <= set(p['nodes']),'exact whole source component/nodes/marks')


def a13(key,W,units):
    doc = source(A13_SOURCE);local = doc['local_proof']
    candidates = [(i,p) for i,p in enumerate(local['selected_patterns']) if p['index'] == key['high_index']]
    if not candidates:
        return None
    si,saved = one(candidates,lambda t:True,'unique complete selected A13 source row')
    if not saved['affected_maps']:
        return None
    mi = member(key,saved['affected_maps'])
    need(local['all_affected_components_already_whole_A13'] is True
         and len(key['base_nodes']) == 13 and key['marked_image'] == [0,12],
         'whole affected A13 endpoint source domain')
    proofs = [(i,p) for i,p in enumerate(local['radical_proofs']) if p['index'] == key['high_index']]
    saturated = local['saturated_two_witness_proof']
    if proofs:
        pi,p = one(proofs,lambda t:True,'unique radical A13 source proof')
        selector = ['local_proof','radical_proofs',pi]
        need(p['all_further_whites_core_orthogonal_unless_short5'] is True,
             'source radical first-contact seal retained')
        variants = [(None,p,p['selected_components'],selector)]
    elif saturated['index'] == key['high_index']:
        p = saturated;selector = ['local_proof','saturated_two_witness_proof']
        need(p['maximum_further_core_contact_whites'] == 1
             and p['core_orthogonal_whites_orthogonal_to_forced_roots'] is True
             and p['at_least_one_original_A1_remains_whole'] is True
             and [v['extra_contact'] for v in p['variants']] ==
                 [None]+[r['contact'] for r in p['forced_roots']],
             'complete zero/one saturated root alternatives and persistent whole A1')
        variants = [(v['extra_contact'],v,v['components'],selector+['variants',i])
                    for i,v in enumerate(p['variants'])]
    else:
        raise ValueError('affected A13 map has no exact source sealing proof')
    need(p['core_key'] == key['core_key'] and p['selected_marks'] == key['selected_marks'],
         'same A13 source core/ordered marks')
    source_physical(key,p['gram'],p['white'],True)
    value = W.finite_corner_weight([('A',13),('A',1)],(0,0),(0,12))
    need(value >= 1 and W.finite_same_weight(('A',13),0,12) == 0,
         'CURRENT A13 whole-A1 context unit, no intrinsic endpoint payment')
    branches = []
    for extra,v,parts,vs in variants:
        component(parts,'A13',key['base_nodes'],key['selected_marks'])
        witnesses = v['whole_A1_witnesses']
        need(witnesses and all(w['type'] == 'A1' and not set(w['nodes']) & set(key['selected_marks'])
             and any(c['type'] == 'A1' and c['nodes'] == w['nodes'] for c in parts) for w in witnesses),
             'same whole source unmarked A1 witness, not an arbitrary displayed leaf')
        if extra is not None:
            n = len(p['gram']);gg = matrix_values(v['gram']);old = matrix_values(p['gram'])
            need(len(gg) == n+1 and [r[:n] for r in gg[:n]] == old
                 and v['white'] == p['white']+[n] and gg[n][n] == 2
                 and all(gg[n][j] == gg[j][n] == (-1 if j == extra else 0) for j in range(n)),
                 'exact unique norm2 saturated root physical extension')
        else:
            need(v['gram'] == p['gram'] and v['white'] == p['white'],'same no-extra-root A13 physical model')
        branches.append(dict(condition=dict(extra_contact=extra),source_selector=vs,
            physical=context(v['gram'],v['white']),whole_unmarked_A1_witnesses=witnesses,
            paying_corners=[dict(marks=key['selected_marks'],current_lower=str(value))],
            current_total_lower=str(value)))
    return dict(source=A13_SOURCE,proof_selector=selector,map_source=A13_SOURCE,
        map_selector=['local_proof','selected_patterns',si,'affected_maps',mi],
        physical=context(p['gram'],p['white']),source_sealing=dict(
            kind='accepted_source_complete_A13_radical_or_saturated_whole_A1_seal',
            proof_selector=selector),
        paying_corners=branches[0]['paying_corners'] if len(branches) == 1 else [],
        receipt_branches=branches,actual_receipts=[],
        consumer_literal_evaluations=[dict(function='finite_corner_weight',
            forest=[['A',13],['A',1]],first=[0,0],second=[0,12],result=str(value),
            scope='CURRENT scalar executed by THIS actual selected consumer; not a separate historic receipt')],
        proof_references=[pin(FINITE_PROOF),pin(FINITE_ROOT),pin(LATE_PROOF)])


def a15(key,W,units):
    doc = source(A15_SOURCE)
    ri,support = one([[i,r] for i,r in enumerate(doc['complete_selected_support'])],
        lambda t:t[1]['kind'] == 'A','complete A15 support family')
    candidates = [(i,p) for i,p in enumerate(support['affected']) if p['index'] == key['high_index']]
    if not candidates:
        return None
    si,saved = one(candidates,lambda t:True,'unique A15 affected source occurrence')
    need(saved['maps'].count(key['induced_map']) == 1,'THIS complete source A15 affected native map')
    pi,p = one([[i,p] for i,p in enumerate(doc['selected_proofs'])],
        lambda t:t[1]['index'] == key['high_index'],'same exact A15 source geometry proof')
    need(p['core_key'] == key['core_key'] and p['sealing']['method'] == 'actual_radical'
         and len(key['base_nodes']) == 15 and key['marked_image'] == [0,14]
         and key['marked_image'] in p['selected_images'],
         'whole selected A15 endpoints in the accepted complete radical source domain')
    source_physical(key,p['gram'],p['white'],True)
    selected = one(p['corners'],lambda c:c['marks'] == key['selected_marks'],'same A15 selected physical corner')
    component(selected['components'],'A15',key['base_nodes'],key['selected_marks'])
    need(selected['canonical_images'] == p['selected_images']
         and W.finite_same_weight(('A',15),*key['marked_image']) == 0,'CURRENT selected A15 zero remains zero')
    ci,payer = one([[i,c] for i,c in enumerate(p['corners'])],
        lambda t:t[1].get('reason') == 'restored_whole_A17_unit','distinct whole A17 source receiver')
    large = one(payer['components'],lambda c:c['type'] == 'A17' and set(payer['marks']) <= set(c['nodes']),
                'same whole A17 marked source component')
    pairs = payer['canonical_images']
    need(pairs and len(large['nodes']) == 17 and set(payer['marks']) != set(key['selected_marks'])
         and all(type(pair) is list and len(pair) == 2 for pair in pairs),
         'complete marked A17 canonical image orbit at a distinct actual corner')
    values = [W.finite_same_weight(('A',17),*pair) for pair in pairs]
    lower = min(values);need(lower >= 1,'EVERY exact source A17 canonical marked image is CURRENT unit')
    paid = [dict(marks=payer['marks'],current_lower=str(lower),
                 source_selector=['selected_proofs',pi,'corners',ci])]
    return dict(source=A15_SOURCE,proof_selector=['selected_proofs',pi],map_source=A15_SOURCE,
        map_selector=['complete_selected_support',ri,'affected',si,'maps',saved['maps'].index(key['induced_map'])],
        physical=context(p['gram'],p['white']),source_sealing=p['sealing'],paying_corners=paid,
        actual_receipts=[],consumer_literal_evaluations=[
            dict(function='finite_same_weight',type=['A',17],marked_image=pair,result=str(value),
                 scope='CURRENT scalar inside this actual selected consumer, not a prior receiver job')
            for pair,value in zip(pairs,values,strict=True)],
        proof_references=[pin(FINITE_PROOF),pin(FINITE_ROOT),pin(LATE_PROOF)])


def d_long(key,W,units):
    doc = source(D_LONG_SOURCE);rank = key['target_order']
    ri,record = one([[i,r] for i,r in enumerate(doc['records'])],
        lambda t:t[1]['rank'] == rank and t[1]['selected_support']['kind'] == 'D',
        'complete finite D16/D17 source family')
    support = record['selected_support']
    candidates = [(i,p) for i,p in enumerate(support['affected']) if p['index'] == key['high_index']]
    if not candidates:
        return None
    si,saved = one(candidates,lambda t:True,'unique complete affected finite D long source occurrence')
    mi = member(key,saved['maps'])
    need(one(support['rows'],lambda r:r['index'] == key['high_index'],'same full selected-support row') == saved
         and saved['core_key'] == key['core_key'] and saved['base_nodes'] == key['base_nodes']
         and saved['marks'] == key['selected_marks'] and saved['method'] == 'radical',
         'same complete finite D source core/base/ordered marks/full-map support')
    pi,p = one([[i,p] for i,p in enumerate(record['proofs'])],
        lambda t:t[1]['index'] == key['high_index'],'same exact finite D radical geometry proof')
    need(p['core_key'] == key['core_key'] and p['whole_core_sealed_unless_short8'] is True
         and p['all_mixed_sign_shortcuts']
         and all(0 < q['length'] <= 8 for q in p['all_mixed_sign_shortcuts'])
         and len(key['base_nodes']) == rank,'complete source radical all-contact whole-core seal')
    source_physical(key,p['core_gram'],p['white'],True)
    selected = one(p['corners'],lambda c:c['marks'] == key['selected_marks'],'same finite D selected physical corner')
    component(selected['components'],'D'+str(rank),key['base_nodes'],key['selected_marks'])
    own = W.finite_same_weight(('D',rank),*key['marked_image'])
    if own != F(1,2):
        return None  # In particular the later D17(2,14) zero needs its own successor theorem.
    ci,payer = one([[i,c] for i,c in enumerate(p['corners'])],
        lambda t:bool(t[1].get('coupon')),'one exact distinct finite A1/A2,D14 source coupon')
    need(len(payer['coupon']) == 1 and set(payer['marks']) != set(key['selected_marks']),
         'distinct actual physical coupon corner, not a second use of selected half')
    coupon = payer['coupon'][0];small = coupon['small_type']
    need(small in ('A1','A2') and coupon['large_type'] == 'D14' and coupon['whole_components'] is True
         and {coupon['small_mark'],coupon['large_mark']} == set(payer['marks']),
         'whole exactly typed source small/D14 coupon and actual ordered marks')
    component(payer['components'],small,coupon['small_nodes'],[coupon['small_mark']])
    component(payer['components'],'D14',coupon['large_nodes'],[coupon['large_mark']])
    images = coupon['large_canonical_mark_images'];smallrank = int(small[1:])
    need(images and all(type(q) is list and len(q) == 1 and type(q[0]) is int for q in images),
         'complete saved one-mark D14 image orbit')
    evaluations = [dict(function='finite_cross_weight',small_type=['A',smallrank],small_mark=sm,
        large_type=['D',14],large_mark=pair[0],
        result=str(W.finite_cross_weight(('A',smallrank),sm,('D',14),pair[0])))
        for sm in range(smallrank) for pair in images]
    alternative = min(rational(q['result']) for q in evaluations)
    need(own+alternative >= 1,'CURRENT selected half plus exact source DISTINCT coupon pays the same face')
    paid = [dict(marks=key['selected_marks'],current_lower=str(own)),
            dict(marks=payer['marks'],current_lower=str(alternative))]
    return dict(source=D_LONG_SOURCE,proof_selector=['records',ri,'proofs',pi],map_source=D_LONG_SOURCE,
        map_selector=['records',ri,'selected_support','affected',si,'maps',mi],
        physical=context(p['core_gram'],p['white']),source_sealing=dict(
            kind='accepted_complete_radical_source_core_seal',radical=p['radical'],
            allowed_contacts=p['allowed_contacts'],all_mixed_sign_shortcuts=p['all_mixed_sign_shortcuts']),
        paying_corners=paid,actual_receipts=[],consumer_literal_evaluations=[
            dict(function='finite_same_weight',type=['D',rank],marked_image=key['marked_image'],result=str(own))
            ]+evaluations,
        proof_references=[pin(FINITE_PROOF),pin(FINITE_ROOT),pin(LATE_PROOF)])


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
    need(callable(SOURCE_READER) and callable(SOURCE_PINNER), 'both source readers must be injected')
    pin(TEX);pin(LITERAL)
    need(getattr(W,'TEX_SHA256',None) == TEX_SHA, 'frozen current literal W/canonical manuscript')
    authoritative(key)
    intrinsic = (W.affine_weight(key['target_kind'],key['target_order'],*key['marked_image'])
                 if key['affine'] else
                 W.finite_same_weight((key['target_kind'],key['target_order']),*key['marked_image']))
    bound = None;family = None
    if not key['affine']:
        kind,rank,pair = key['target_kind'],key['target_order'],key['marked_image']
        if kind == 'D' and rank == 12 and pair in ([0,10],[0,11]):
            need(intrinsic == F(1,2),'CURRENT D12 long-leaf half')
            bound = d12(key,W,units);family = 'selected_D12_full_source_domain'
        elif kind == 'D' and rank == 14 and pair in ([0,12],[0,13]):
            need(intrinsic == 0,'CURRENT unguarded D14 long-leaf zero')
            bound = d14(key,W,units);family = 'selected_D14_full_source_domain'
        elif kind == 'A' and rank == 13 and pair == [0,12]:
            bound = a13(key,W,units);family = 'selected_A13_whole_A1_context'
        elif kind == 'A' and rank == 15 and pair == [0,14]:
            bound = a15(key,W,units);family = 'selected_A15_whole_A17_alternative'
        elif kind == 'D' and rank in (16,17) and pair[0] == rank-15 and intrinsic == F(1,2):
            bound = d_long(key,W,units);family = 'selected_D'+str(rank)+'_whole_small_D14_coupon'
    elif key['target_order'] == 12 and (
            key['target_kind'] == 'B' and key['marked_image'] in ([9,10],[9,11]) or
            key['target_kind'] == 'C' and key['marked_image'] in ([0,10],[1,11])):
        need(intrinsic == 0,'CURRENT removed affine B/C12 pair remains zero')
        bound = bc12(key,W,units);family = 'selected_BC12_full_source_domain'
    if bound is None:
        bound = generic_forced_entry(key,W);family = 'generic_strict_proper_forced_entry'
    if bound is None and key['affine'] and key['target_kind'] == 'B' and key['target_order'] == 17 and key['high_index'] == 512:
        need(key['marked_image'] in ([9,15],[9,16]) and intrinsic == 0,
             'exact independently reviewed fixed B17 fallback domain')
        bound = b17_exclusion(key,W);family = 'selected_B17_fixed_forced_terminal'
    if bound is None:
        return None
    frozen_key = json.loads(canonical(key))
    excluded = bound.get('application_kind') == 'EXCLUDED_COMPLETION_UNDER_NO_S8'
    supplemental = [
        dict(id='finite57.legacy.exact_source_occurrence',
            statement='This exact ordered high/core/base/marks/full native map describes the SAME actual source face and stated whole target.',
            facts=dict(source=pin(bound['source']),proof_selector=bound['proof_selector'],
                same_full_key_field='/application_key/selected_full_key',
                physical_gram_white_field='/physical_context')),
        dict(id='finite57.legacy.complete_native_and_contact_domains',
            statement='The frozen source complete native-map orbit and necessary contact domain, original opposite-wall data, gauge/kernel/inverse/signature and acute original-wall hypotheses apply to THIS occurrence; no supporting-white contact outside that complete domain is omitted.',
            facts=dict(map_source=pin(bound['map_source']),map_selector=bound['map_selector'],
                given_full_native_map=key['induced_map'],marked_image=key['marked_image'],
                target=dict(affine=key['affine'],kind=key['target_kind'],order=key['target_order']),
                complete_source_predicates_field='/source_sealing'))]
    if excluded:
        supplemental.append(dict(id='finite57.legacy.strict_forced_entry_geometry',
            statement='At the exact verified original corner and proper native boundary, the missing supporting root is genuinely new, physical norm2; complete original-contact/end hypotheses give either a same-face short8 escape or the source-bound projection obstruction. No unrecorded contact or old omitted wall repairs entry.',
            facts=dict(complete_boundary_certificate_field='/source_sealing',
                       source_implication_field='/source_branch_conclusion')))
    else:
        supplemental += [
            dict(id='finite57.legacy.source_whole_component_and_branch_geometry',
                statement='Under NO same-face short8 escape, the exact source radical/saturated/protected proof seals precisely the complete paying factors and whole unmarked witnesses declared in EACH exhaustive exclusive physical branch. No indirect supporting-white tail has an unexamined first contact.',
                facts=dict(source_full_core_seal_field='/source_sealing',
                    exact_branch_conditions=[b['condition'] for b in bound.get('receipt_branches',[])],
                    complete_branch_receipts_field='/conditional_receipt_branches')),
            dict(id='finite57.legacy.current_scalar_canonical_roles',
                statement='The frozen complete source canonical images and actual ordered physical marks are the complete roles of the declared whole factors. CURRENT scalar values or fully validated exact actual terminal receipts apply to those SAME factors and distinct corners, with all current witness guards/debits retained.',
                facts=dict(current_literal=pin(LITERAL),canonical_manuscript=pin(TEX),
                    scalar_inputs_results_field='/consumer_literal_evaluations',
                    prior_actual_receipt_identities_field='/actual_terminal_receipts',
                    distinct_paying_corners_field='/current_paying_corners'))]
    acceptance = accepted_premises([p['id'] for p in supplemental])
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
        current_paying_corners=bound['paying_corners'],intrinsic_selected_current_lower=str(intrinsic),
        conditional_receipt_branches=bound.get('receipt_branches',[]),
        consumer_literal_evaluations=bound.get('consumer_literal_evaluations',[]),
        supplemental_source_premise_acceptance=acceptance,
        supplemental_source_premises=supplemental,
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
            implication=('For THIS exact full source/key and complete physical/native/contact/end/inverse predicates: NO same-face adjacent white-interior path<=8 excludes this actual completion.' if excluded else
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
