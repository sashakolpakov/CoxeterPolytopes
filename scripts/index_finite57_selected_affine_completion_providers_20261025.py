#!/usr/bin/env python3
"""SOURCE ONLY: selected affine provider keys from ten accepted source packets.

Exact stored full maps are joined to authoritative native high rows. Proper
exclusion implications remain source-computed/geometric premises; whole model
payments are recomputed with the DIRECT current receiver condition. This does
not assert that the resulting union is the complete G6 bridge.
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
SOURCE_PINS = 'results/finite57_selected_affine_completion_chain_source_pins_20261024.json'
SOURCE_PINS_SHA = '5f94aff40b8178158e8421c0cbfe1eb377469cd20fdc333aafe6e94e7fe611b5'
LOWER_PINS = 'results/finite57_g4_lower_cycle_source_inputs_20261025.json'
LOWER_PINS_SHA = '9f26c3f91c8c47126e8506e23d7a77497735974e012209a582da1344ddeb232b'
MAX_OUTPUT = 500 * 1024
STEMS = (
    'fresh_b14_b15_singleton_completion_transfer_20261017',
    'fresh_order15_affine_c_terminal_projection_transfer_20261017',
    'fresh_order15_affine_b_pair_reduction_20261017',
    'fresh_order15_affine_d_half_coupon_transfer_20261017',
    'fresh_order16_affine_b_singleton_projection_transfer_20261018',
    'fresh_order16_d_distance9_half_transfer_20261018',
    'fresh_order16_bd_sparse_unit_transfer_20261018',
    'fresh_order16_cd_sparse_sixteenth_transfer_20261018',
    'fresh_order17_affine_d_three_quarter_transfer_20261017',
    'fresh_order13_bc_sparse_terminal_transfer_20261018',
)
EXCLUDED = {
    'projection_exclusion', 'projection_exceeds_physical_norm_two',
    'projection_above_norm_two', 'saturated_entry_has_no_next_white',
    'forced_doubled_terminal_exceeds_norm2', 'excluded_forced_entry_projection',
    'excluded_multiple_core_contact_roots', 'excluded_projection_above_norm',
    'excluded_lorentz_index', 'excluded_forced_lorentz_index',
    'two_entry_exclusion_or_shortcut', 'forced_lorentz_index_exclusion',
    'physical_gauge_projection_exceeds_new_norm',
}


def need(ok, message):
    if not ok:
        raise ValueError(message)


def raw(relative):
    path = (ROOT / relative).resolve()
    need(path.is_relative_to(ROOT) and path.is_file(), 'bounded selected-affine source path')
    with path.open('rb') as stream:
        data = stream.read(2 * 1024**2 + 1)
    need(0 < len(data) <= 2 * 1024**2, 'bounded selected-affine source read')
    return data


def prepare():
    data = raw(LOADER)
    need(sha256(data).hexdigest() == LOADER_SHA, 'selected-affine loader identity')
    helper = types.ModuleType('finite57_selected_affine_loader')
    helper.__file__ = str(ROOT / LOADER)
    exec(compile(data, str(ROOT / LOADER), 'exec'), helper.__dict__)
    r, c = helper.prepare()
    data = raw(SOURCE_PINS)
    need(sha256(data).hexdigest() == SOURCE_PINS_SHA, 'accepted selected-affine source register identity')
    pins = json.loads(data, object_pairs_hook=c.unique)['sources']
    need(len(pins) == 68 and len({p['path'] for p in pins}) == 68, 'unique complete selected-affine source register')
    for p in pins:
        if p['path'] in c.PINS:
            need(c.PINS[p['path']] == p['sha256'], 'selected-affine/base pin collision')
        c.PINS[p['path']] = p['sha256']
    data = raw(LOWER_PINS)
    need(sha256(data).hexdigest() == LOWER_PINS_SHA, 'lower-order source register identity')
    lower = json.loads(data, object_pairs_hook=c.unique)['inputs']
    need(len(lower) == 64 and len({p['path'] for p in lower}) == 64, 'unique lower-order source register')
    for p in lower:
        if p['path'] in c.PINS:
            need(c.PINS[p['path']] == p['sha256'], 'lower-order/base pin collision')
        c.PINS[p['path']] = p['sha256']
    return helper, r, c


def dense_target(saved, canonical):
    n = len(canonical)
    need(type(saved) is list and len(saved) == n, 'saved canonical target order')
    for i, row in enumerate(saved):
        need(type(row) is list and len(row) == n, 'saved native target matrix schema')
        for j, label in enumerate(row):
            need(type(label) is int and label == canonical[i].get(j, 0), 'entire saved/canonical target identity')


def groups(stem, c):
    """Return explicit source case/map containers, never broad JSON name scans."""
    position = STEMS.index(stem)
    if position == 0:
        for k, rec in enumerate(c.field(stem, ('records',))):
            yield f'records/{k}', 'B', rec['order'], rec['target'], rec['cases'], 'maps'
    elif position in (6, 7, 9):
        section = 'rules' if position == 7 else 'records'
        order = 13 if position == 9 else 16
        for k, rec in enumerate(c.field(stem, (section,))):
            container = 'completions' if position == 7 and rec['kind'] == 'D' else 'rows'
            if position == 9:
                container = 'bc13'
            yield f'{section}/{k}', rec['kind'], order, rec['target'], rec['cases'], container
    else:
        kind, order, container = {
            1: ('C', 15, 'proof-maps'), 2: ('B', 15, 'finite-bound'),
            3: ('D', 15, 'proof-maps'), 4: ('B', 16, 'maps'),
            5: ('D', 16, 'completions'), 8: ('D', 17, 'proof-maps'),
        }[position]
        paths = [('cases',)]
        if position != 2:
            paths.append(('target',))
        else:
            paths.append(('removed_pairs',))
        doc = c.document(stem, *paths)
        saved = doc.get('target')
        yield 'root', kind, order, saved, doc['cases'], container


def model_receipts(c, encoded_case, path):
    result = []
    for location, gram, white, corners, gp, wp in c.complete_models(encoded_case, path):
        receipt = c.receipt('selected_affine_model', location, gram, white, corners, dict(gram=gp, white=wp))
        need(receipt['direct_current_receiver_condition'] and not receipt['legacy_pointwise_bounds_used'],
             'whole affine provider uses direct CURRENT receipt only')
        result.append(dict(path=location, context_sha256=receipt['context_sha256'],
            current_total=receipt['total'], gram=gram, white=white,
            current_corners=[dict(marks=q['marks'], kind=q['kind'], lower=q['weight']) for q in receipt['corners']]))
    need(result, 'whole selected-affine provider requires a complete current-funded model list')
    return result


def maps_and_body(p, container):
    if container == 'proof-maps':
        body = p['proof']
        return body['maps'], body
    if container == 'bc13':
        if 'inherited_full_core_seal' in p:
            old = p['inherited_full_core_seal']
            body = dict(gram=old['gram'], white=old['white'], corners=p['corners'])
            return None, body  # Source stores the exact orbit, not every full B map.
        return p['maps'], p
    return p[container], p


def physical_case(c, p, body, core):
    source = body.get('core_gram', body.get('gram', p.get('core_gram', p.get('gram'))))
    need(source is not None, 'source case original physical Gram')
    need(c.E.matrix(source) == c.E.matrix(core['gram']), 'source/high original physical Gram identity')


def whole_map_matches(c, r, model, case, kind, order, images):
    g = c.E.matrix(model['gram'])
    ids = sorted(set(model['white']) | set(case['preferred_corner']))
    parts = c.E.components(g, ids, c.FUEL)
    selected = next((p for p in parts if set(case['preferred_corner']) <= set(p)), None)
    if selected is None or len(selected) != order or not set(case['marked_component']) <= set(selected):
        return False
    target = c.W.classical_graph(kind, order, affine=True)
    labels = c.E.labels(g, selected)
    for full in r.embeddings(labels, target, c.FUEL, c.E):
        if [full[selected.index(i)] for i in case['marked_component']] == images:
            return True
    return False


def bc13_projection(c, payload, target, images, case, core):
    """Bind the unlabelled BC13 source rows by their exact physical entries."""
    base = case['marked_component']
    missing = sorted(set(range(13)) - set(images))
    required = sorted((new, old, base[i], target[new][old])
                      for new in missing for i, old in enumerate(images) if old in target[new])
    got = []
    for entry in payload['forced_entries']:
        new, old, node, label = (entry[k] for k in ('new_target_node', 'old_target_node', 'core_node', 'native_label'))
        need(label in (3, 4) and node in case['possible_contacts'], 'BC13 actual native entry/contact domain')
        square = F(1 if label == 3 else 2)
        gauge = F(core['gauge'][node])
        inverse = F(core['inverse'][node][node])
        need(F(entry['physical_native_contact_square']) == square
             and F(entry['old_gauge_squared']) == gauge
             and F(entry['inverse_diagonal']) == inverse, 'BC13 physical gauge/inverse associations')
        need(F(entry['projection_lower_bound']) == square * gauge * inverse > 2,
             'BC13 strict physical projection bound at this exact map contact')
        got.append((new, old, node, label))
    need(required and sorted(got) == required, 'complete BC13 required native entry set')


def run(number, c, r):
    stem = STEMS[number]
    script = f'scripts/{stem}.py'
    need(script in c.PINS and c.digest(script) == c.PINS[script], 'registered selected-affine source pin')
    c.READ.add(script)
    high, cores, partition = r.domain(c)
    result = []
    for group, kind, order, saved_target, cases, container in groups(stem, c):
        need(kind in ('B', 'C', 'D') and 13 <= order <= 17, 'explicit selected-affine target domain')
        target = c.W.classical_graph(kind, order, affine=True)
        if saved_target is not None:
            dense_target(saved_target, target)
        need(type(cases) is list and cases, 'nonempty complete registered case list')
        for ordinal, p in enumerate(cases):
            index = p.get('index', p.get('case'))
            need(type(index) is int and 0 <= index < len(high), 'authoritative selected-affine high index')
            case = high[index]
            core = cores[case['core_key']]
            if 'core_key' in p:
                need(p['core_key'] == case['core_key'], 'selected-affine original core key')
            base = case['marked_component']
            bound = case.get('component_order_bound')
            need(type(bound) is int and bound >= order and len(base) <= order,
                 'entire authoritative selected affine eligibility predicate')
            labels = c.E.labels(c.E.matrix(core['gram']), base)
            prefix = f'{group}/cases/{ordinal}'
            if container == 'finite-bound':
                physical_case(c, p, p, core)
                maximum = p['outcome']['maximum_order']
                need(type(maximum) is int and maximum < order, 'actual whole finite-component bound excludes affine target')
                need(p['outcome']['method'] in ('no_direct_contact_to_marked_component', 'zero_residual_endpoint_extension'),
                     'accepted exact B15 finite-component bound source')
                dropped = {tuple(q) for q in c.field(stem, ('removed_pairs',))}
                candidates = [dict(map=mp, outcome='actual_finite_component_order_below_affine_target')
                              for mp in r.embeddings(labels, target, c.FUEL, c.E)
                              if tuple(sorted(mp[base.index(m)] for m in case['preferred_corner'])) in dropped]
                body = p
            else:
                candidates, body = maps_and_body(p, container)
                physical_case(c, p, body, core)
                if container == 'bc13' and candidates is not None:
                    need(body['allowed_contacts'] == case['possible_contacts']
                         and body['end_sets'] == case['end_sets'], 'BC13 actual full contact/end domains')
                    need([F(x) for x in body['gauge_squared']] == [F(x) for x in core['gauge']],
                         'BC13 entire original physical gauge')
                    ends = body['end_sets']
                    need(len(ends) == 2 and not set(ends[0]) & set(ends[1])
                         and set(ends[0] + ends[1]) == set(case['possible_contacts']),
                         'BC13 disjoint complete contact partition')
                    need(core['inverse'] is not None and all(F(core['inverse'][i][j]) > 0
                         for end in ends for i in end for j in end), 'BC13 same-end inverse positivity')
                if candidates is None:
                    pair = {'B': [1, 10], 'C': [1, 11]}[kind]
                    candidates = [mp for mp in r.embeddings(labels, target, c.FUEL, c.E)
                                  if sorted(mp[base.index(m)] for m in case['preferred_corner']) == pair]
            for metadata in (p, body):
                for name in ('base_nodes', 'marked_component', 'selected_nodes'):
                    if name in metadata:
                        need(metadata[name] == base, 'source full-map original base-node ordering')
                for name in ('preferred_corner', 'selected_marks', 'marks'):
                    if name in metadata:
                        need(metadata[name] == case['preferred_corner'], 'source selected ordered mark association')
            need(type(candidates) is list and candidates, 'complete nonempty selected map-provider list')
            entries = []
            seen = set()
            whole = None
            model_catalog = {}
            for j, candidate in enumerate(candidates):
                if type(candidate) is list:
                    images, payload, outcome = candidate, body, 'whole_source_case'
                else:
                    need(type(candidate) is dict, 'selected map-provider record schema')
                    images = candidate.get('map', candidate.get('embedding', candidate.get('images')))
                    payload = candidate.get('proof', candidate)
                    outcome = payload.get('outcome')
                c.E.indices(images, order, len(base))
                need(all(labels[a][b] == target[images[a]].get(images[b], 0)
                         for a in range(len(base)) for b in range(a)), 'FULL induced native selected-affine map')
                key = dict(high_index=index, core_key=case['core_key'], base_nodes=base,
                    selected_marks=case['preferred_corner'], target_kind=kind, target_order=order,
                    affine=True, induced_map=images,
                    marked_image=sorted(images[base.index(m)] for m in case['preferred_corner']))
                encoded = c.canonical(key)
                need(encoded not in seen, 'duplicate full selected-affine map key in case')
                seen.add(encoded)
                if container == 'bc13' and outcome is None:
                    bc13_projection(c, payload, target, images, case, core)
                    outcome = 'physical_gauge_projection_exceeds_new_norm'
                if outcome in ('whole_source_case', 'sealed_completion'):
                    if outcome == 'whole_source_case':
                        if whole is None:
                            whole = model_receipts(c, body, prefix + '/proof')
                        funded = whole
                    else:
                        funded = model_receipts(c, payload, prefix + f'/maps/{j}')
                    complete_keys, matching = [], []
                    for q in funded:
                        model_key = sha256(c.canonical([q['gram'], q['white'],
                            [corner['marks'] for corner in q['current_corners']]]).encode()).hexdigest()
                        model_catalog.setdefault(model_key, dict(model_key=model_key,
                            path=q['path'], context_sha256=q['context_sha256'],
                            current_total=q['current_total'], current_corners=q['current_corners']))
                        if model_key not in complete_keys:
                            complete_keys.append(model_key)
                        if whole_map_matches(c, r, q, case, kind, order, images) and model_key not in matching:
                            matching.append(model_key)
                    need(matching, 'current-funded completed factor must realize THIS exact original-base map')
                    proof = dict(kind='CURRENT_FUNDED_WHOLE_AFFINE_VARIANTS',
                                 complete_funded_model_keys=complete_keys, matching_model_keys=matching,
                                 required_source_variant_exhaustiveness='Accepted source-specific full-variant theorem; not established by this graph search',
                                 source_variant_exhaustiveness_verified_by_this_index=False)
                else:
                    need(outcome in EXCLUDED or outcome == 'actual_finite_component_order_below_affine_target',
                         'explicit accepted selected-affine exclusion outcome')
                    proof = dict(kind='EXACT_PROPER_AFFINE_EXCLUSION', outcome=outcome,
                        computed_record_sha256=sha256(c.canonical(payload).encode()).hexdigest(),
                        premise='Accepted source-specific all-native contact/residual/inertia theorem at its exact saved computed identities, for THIS full map, modulo independent current S8.')
                entries.append(dict(key=key, source_map_path=prefix + f'/maps/{j}', provider=proof))
            result.append(dict(source=stem, group=group, source_case_ordinal=ordinal,
                               high_index=index, current_funded_model_catalog=list(model_catalog.values()),
                               complete_registered_map_providers=entries))
    return dict(source=stem, source_number=number, exact_registered_case_providers=result,
                scope='Only this explicit source packet. Global full-map coverage, older original/newer providers, S8/triangle/priority domains remain separate.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch', choices=tuple(f'affine-{i}' for i in range(len(STEMS))), required=True)
    parser.add_argument('--parent-lock', action='store_true', required=True)
    args = parser.parse_args()
    relative = str(Path(__file__).resolve().relative_to(ROOT))
    before = raw(relative)
    helper, r, c = prepare()
    out = run(int(args.batch.rsplit('-', 1)[1]), c, r)
    for path in c.READ:
        need(c.digest(path) == c.PINS[path], 'selected-affine source/input changed')
    need(raw(relative) == before and sha256(raw(LOADER)).hexdigest() == LOADER_SHA, 'affine index/loader changed')
    need(sha256(raw(SOURCE_PINS)).hexdigest() == SOURCE_PINS_SHA, 'selected-affine source register changed')
    need(sha256(raw(LOWER_PINS)).hexdigest() == LOWER_PINS_SHA, 'lower-order source register changed')
    need(sha256(helper.raw(helper.COLLECTOR)).hexdigest() == helper.COLLECTOR_SHA, 'direct collector changed')
    need(sha256(helper.raw(helper.ROUTES)).hexdigest() == helper.ROUTES_SHA, 'native route definitions changed')
    need(sha256(c.raw(c.MANIFEST)).hexdigest() == c.MANIFEST_SHA, 'base110 input register changed')
    payload = json.dumps(dict(status='SELECTED_AFFINE_PROVIDER_INDEX_ONLY', batch=args.batch,
        source_sha256=sha256(before).hexdigest(), current_terminal_bridge_complete=False,
        consumed_inputs={p: c.PINS[p] for p in sorted(c.READ)}, fuel_used=c.FUEL[0], result=out), sort_keys=True)
    need(len(payload.encode()) <= MAX_OUTPUT, 'affine provider output cap; refuse not truncate maps')
    print(payload)


if __name__ == '__main__':
    main()
