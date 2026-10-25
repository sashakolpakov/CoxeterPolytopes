#!/usr/bin/env python3
"""Bind complete actual provider exports to reviewed, conditional proof domains.

Only JSON/source/receipt associations are checked. No repository worker is
imported or executed, no native maps are enumerated, and no geometric or old
complete numerical premise is regenerated. Root alone launches this reader.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
import base64
from decimal import Decimal
import fcntl
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import re
import zlib

ROOT = Path(__file__).resolve().parents[1]
OWN = 'scripts/curate_finite57_provider_occurrences_20261025.py'
INPUTS = 'results/finite57_provider_occurrence_actual_inputs_20261025.json'
INPUT_SHA = '0ec63161c7ea3473b82673bf9c935a4dc6ad7f5b7e13cbe032a81a5054830f48'
REGISTRY = 'results/finite57_current_provider_proof_home_registry_candidate_20261025.json'
REGISTRY_SHA = '857ed9fe5363b6b4d40f97f7f51a0f3fb6eadfc6891c729fe8a33ab08ea22081'
HIGH = 'results/fresh_distance8_high_overlap_20261017.json'
RUNNER = 'scripts/run_bounded_research_check_20261018.py'
RUNNER_SHA = '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132'
POLICY = 'scripts/run_lean_batches_20261018.py'
AFFINE9 = 'scripts/index_finite57_selected_affine9_schema_safe_providers_20261025.py'
AFFINE9_SHA = '679f792f08adf3f51a6e63d5f0dc3a3f55074c73c382bf35c4513da0fdd3d4f5'
MAX_FILE = 40 * 1024**2
MAX_READ = 256 * 1024**2
MAX_RESULT = 8 * 1024**2
MAX_OUTPUT = 500 * 1024
MAX_VISITS = 2_000_000
FILES = {}
DOCUMENTS = {}
READ_BYTES = 0
VISITS = 0
SHA = re.compile(r'[0-9a-f]{64}\Z')
KEY_FIELDS = {'high_index','core_key','base_nodes','selected_marks','target_kind',
              'target_order','affine','induced_map','marked_image'}
CONTROLLERS = {RUNNER,POLICY,'scripts/run_proof_batch_20261024.py',
    'scripts/audit_existing_lean_module_20261024.py',
    'scripts/verify_manuscript_audited_20261024.py',
    'scripts/monitor_owned_proof_tree_20261024.py',
    'scripts/run_lean_batches_signal_safe_20261018.py',
    'scripts/build_illustrated_ra_manuscript_20261018.py',
    'scripts/build_compact_arithmetic_nine_field19_20261020.py','scripts/verify_manuscript.py'}
ORDINARY = {'actual-star','native-classification','complete-short-path-dag',
    'generic-extension-bounds','universal-receiver-geometry','bounded59-163',
    'unbounded-tail','global-incidence','current-rule-identity','literal35-retained'}


def need(ok, message):
    if not ok:
        raise ValueError(message)


def visit():
    global VISITS
    VISITS += 1
    need(VISITS <= MAX_VISITS, 'complete semantic traversal cap; refuse, never truncate')


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',', ':'),ensure_ascii=True)


def digest(value):
    need(type(value) is str and SHA.fullmatch(value), 'exact SHA256 string')
    return value


def unique(pairs):
    out = {}
    for key,value in pairs:
        need(key not in out, 'duplicate JSON field')
        out[key] = value
    return out


def loads(data, cap=MAX_FILE, exact=True):
    text = data.decode('utf-8') if type(data) is bytes else data
    need(type(text) is str and 0 < len(text.encode()) <= cap, 'bounded complete JSON')
    depth = 0; quoted = False; escaped = False
    for char in text:
        if quoted:
            if escaped: escaped = False
            elif char == '\\': escaped = True
            elif char == '"': quoted = False
        elif char == '"': quoted = True
        elif char in '[{':
            depth += 1; need(depth <= 48, 'bounded JSON nesting')
        elif char in ']}':
            depth -= 1; need(depth >= 0, 'JSON bracket balance')
    need(not quoted and depth == 0, 'complete JSON lexical structure')
    def source_float(token):
        value = float(token)
        need(value not in (float('inf'),float('-inf')), 'finite source JSON number')
        return value
    return json.loads(text,object_pairs_hook=unique,parse_float=Decimal if exact else source_float,
        parse_constant=lambda value: (_ for _ in ()).throw(ValueError('nonfinite JSON')))


def read(name, expected=None, cap=MAX_FILE):
    global READ_BYTES
    lexical = PurePosixPath(name)
    need(type(name) is str and 0 < len(name) <= 240 and not lexical.is_absolute()
         and '..' not in lexical.parts, 'bounded relative repository path')
    path = (ROOT/name).resolve()
    need(path.is_relative_to(ROOT) and path.is_file(), 'existing repository input')
    before = path.stat()
    identity = (before.st_dev,before.st_ino,before.st_size,before.st_mtime_ns,before.st_ctime_ns)
    need(0 < before.st_size <= cap, 'bounded complete input file')
    with path.open('rb') as stream:
        data = stream.read(cap+1)
    after = path.stat()
    need(identity == (after.st_dev,after.st_ino,after.st_size,after.st_mtime_ns,after.st_ctime_ns)
         and len(data) == before.st_size, 'stable complete input read')
    observed = sha256(data).hexdigest()
    if expected is not None:
        need(observed == digest(expected), 'source/input byte identity: '+name)
    if name in FILES:
        need(FILES[name] == (identity,observed), 'input changed or conflicting pin')
    else:
        READ_BYTES += len(data)
        need(READ_BYTES <= MAX_READ, 'total unique input cap')
        FILES[name] = identity,observed
    return data


def at(value, selector):
    for token in selector:
        visit()
        if type(value) is list:
            need(type(token) is int and 0 <= token < len(value), 'exact source list selector')
        else:
            need(type(value) is dict and type(token) is str and token in value,
                 'exact source object selector')
        value = value[token]
    return value


def pin(name, expected):
    if name in FILES:
        before = (ROOT/name).resolve().stat()
        identity = (before.st_dev,before.st_ino,before.st_size,before.st_mtime_ns,before.st_ctime_ns)
        need(FILES[name] == (identity,digest(expected)), 'consumed pin changed during curation')
    else:
        read(name,expected)


def source(stem, pins):
    name = 'results/'+stem+'.json'
    need(name in pins, 'original source must have been consumed by this actual producer')
    if name not in DOCUMENTS:
        DOCUMENTS[name] = loads(read(name,pins[name]),exact=False)
    else:
        need(FILES[name][1] == pins[name], 'shared original source pin agreement')
    return DOCUMENTS[name]


def source_path(stem, selector, pins):
    return at(source(stem,pins),selector)


def path_tokens(path):
    need(type(path) is str and path and not path.startswith('/'), 'literal source path')
    return [int(p) if p.isdecimal() else p for p in path.split('/')]


def full_key(key):
    need(type(key) is dict and set(key) == KEY_FIELDS, 'complete exact nine-field key')
    need(type(key['high_index']) is int and 0 <= key['high_index'] < 2614
         and type(key['core_key']) is str and 0 < len(key['core_key']) <= 256
         and type(key['affine']) is bool
         and key['target_kind'] in (('A','B','C','D') if key['affine'] else ('A','B','D'))
         and type(key['target_order']) is int and 10 <= key['target_order'] <= 17,
         'full key native target/high schema')
    nodes,marks,mp,pair = (key[n] for n in ('base_nodes','selected_marks','induced_map','marked_image'))
    for values in (nodes,marks,mp,pair):
        need(type(values) is list and all(type(v) is int and v >= 0 for v in values)
             and len(values) == len(set(values)), 'ordered integer key node/image lists')
    need(10 <= len(nodes) <= 17 and len(nodes) <= key['target_order'] and max(nodes) < 64
         and len(mp) == len(nodes)
         and set(mp) <= set(range(key['target_order'])) and len(marks) == len(pair) == 2
         and set(marks) <= set(nodes)
         and pair == sorted(mp[nodes.index(m)] for m in marks), 'complete induced map/ordered mark association')
    return sha256(canonical(key).encode()).hexdigest()


def funded(row, full=False):
    if full:
        need(row['direct_current_receiver_condition'] is True and
             row['legacy_pointwise_bounds_used'] is False, 'direct current receipt')
        corners = [(p['marks'],p['weight']) for p in row['corners']]
        total = row['total']
    else:
        corners = [(p['marks'],p['lower']) for p in row['current_corners']]
        total = row['current_total']
    need(1 <= len(corners) <= 64 and all(type(m) is list and len(m) == 2
         and all(type(v) is int and 0 <= v < 64 for v in m) and len(set(m)) == 2 for m,_ in corners)
         and len({tuple(sorted(m)) for m,_ in corners}) == len(corners),
         'distinct actual funded corner pairs')
    need(all(0 <= F(v) <= 1 for _,v in corners), 'exact current literal corner range')
    need(sum((F(v) for _,v in corners),F()) == F(total) >= 1,
         'complete exact current distinct-corner sum')
    digest(row['context_sha256'])


def models(rows):
    need(type(rows) is list, 'complete funded-model catalogue')
    out = {}
    for row in rows:
        visit(); name = digest(row['model_key'])
        need(name not in out, 'unique complete model keys')
        funded(row);out[name] = row
    return out


def percentage(resources, controllers):
    expected = dict(id='percentage-v1-20261025',source_sha256=controllers[POLICY],
        admission_disk_percent=2,admission_memory_percent=15,abort_disk_percent=1,
        abort_memory_percent=6,cleanup_disk_percent=2,rounding_digits=6)
    policy = resources['resource_policy']
    need(type(policy) is dict and policy.keys() == expected.keys()
         and all(type(policy[k]) is type(v) and policy[k] == v for k,v in expected.items()),
         'exact current resource policy identity')
    names = ('disk_free_bytes','disk_total_bytes','available_memory_bytes','physical_memory_bytes')
    need(all(type(resources[k]) is int for k in names), 'integer resource byte counts')
    disk,total,memory,physical = (resources[k] for k in names)
    need(0 <= disk <= total and total > 0 and 0 <= memory <= physical and physical > 0
         and 100*disk >= 2*total and 100*memory >= 15*physical, 'exact percentage admission ratios')
    for name,free,capacity in (('free_disk_percent',disk,total),('available_memory_percent',memory,physical)):
        value = resources[name]
        need(type(value) in (int,Decimal) and 0 <= value <= 100, 'finite resource percentage')
        value = F(value)
        need((value*10**6).denominator == 1 and abs(value-F(100*free,capacity)) <= F(1,2*10**6),
             'six-digit resource percentage agrees with bytes')


def actual(entry, review, controllers, expected_sha, expected_status):
    worker,batch = entry['worker'],entry['batch']
    need(review['source_sha256'] == expected_sha, 'fixed reviewed provider source')
    pin(worker,expected_sha);pin(review['path'],review['sha256'])
    pieces = []
    for field,hfield,sizefield,cap in (('path','sha256','bytes',2*1024**2),
            ('stdout_path','stdout_sha256','stdout_bytes',600*1024)):
        need(entry[field].startswith('results/'), 'actual evidence repository path')
        data = read(entry[field],entry[hfield],cap)
        need(type(entry[sizefield]) is int and len(data) == entry[sizefield], 'actual evidence exact bytes')
        pieces.append(data)
    outer = loads(pieces[0]);lines = pieces[1].decode().splitlines()
    need(outer['status'] == 'PASS_BOUNDED_RESEARCH_CHECK'
         and outer['command'] == [worker,'--batch',batch,'--parent-lock']
         and outer['checker_sha256'] == expected_sha and outer['runner_sha256'] == RUNNER_SHA,
         'actual named source/command PASS')
    percentage(outer['initial_resources'],controllers)
    need(type(outer['wall_limit_seconds']) is int and outer['wall_limit_seconds'] == 59
         and type(outer['rss_limit_bytes']) is int and outer['rss_limit_bytes'] == 1024**3
         and type(outer['threads']) is int and outer['threads'] == 1
         and outer['lock_owner'] == 'wrapper'
         and type(outer['elapsed_seconds']) in (int,Decimal) and 0 <= outer['elapsed_seconds'] < 59
         and type(outer['peak_process_tree_rss_bytes']) is int
         and 0 <= outer['peak_process_tree_rss_bytes'] < 1024**3, 'measured resource/lock contract')
    output = outer['stdout']
    need(type(output) is str and output.endswith('\n') and len(output.splitlines()) == 1
         and len(output.encode()) == outer['stdout_bytes'] < 512*1024
         and len(lines) == 3 and lines[0]+'\n' == output, 'complete public three-line stdout')
    need(loads(lines[1]) == {k:v for k,v in outer.items() if k != 'stdout'}, 'outer/public summary identity')
    bridge = loads(lines[2])
    need(bridge['controller_source_sha256'] == controllers and bridge['owned_tree_bridge_only'] is True
         and type(bridge['child_exit_code']) is int and bridge['child_exit_code'] == 0
         and bridge['cleanup']['signals'] == [] and bridge['cleanup']['verified_live_owned_processes_remaining'] == 0
         and bridge['assertions_enabled'] is True and type(bridge['python_optimize_flag']) is int
         and bridge['python_optimize_flag'] == 0 and bridge['child_pythonoptimize'] is None,
         'actual controller/ownership/assertion closure')
    inner = loads(output,exact=False)
    need(inner['status'] == entry['inner_status'] == expected_status and inner['batch'] == batch
         and inner['source_sha256'] == expected_sha
         and inner.get('current_terminal_bridge_complete',False) is False,
         'actual exact supported inner source/status/batch and local scope')
    pins = inner['consumed_inputs']
    need(type(pins) is dict and 'manuscripts/finite_volume57.tex' in pins, 'literal canonical source consumed')
    for name,expected_pin in pins.items():
        # Avoid rereading the same large physical source for every receipt.
        # All tracked files receive a complete final byte recheck below.
        pin(name,expected_pin)
    return inner,pins


def source_record(stem,selector,pins):
    value = source_path(stem,selector,pins)
    need(type(value) is dict, 'exact original source proof record')
    return value,dict(path='results/'+stem+'.json',selector=selector,
        source_sha256=pins['results/'+stem+'.json'],record_sha256=sha256(canonical(value).encode()).hexdigest())


def registration(worker_spec,batch,stem):
    rows = [r for r in worker_spec['rows'] if r['batch'] == batch and r['source'] == stem]
    need(len(rows) == 1, 'one exact registered batch/source proof domain')
    return rows[0]


def proof_homes(adapter,batch,reg):
    if adapter == 'affine':
        return ['sparse_BC13','S8'] if batch == 'affine-9' else ['selected_affine','S8']
    if adapter == 'old_affine':return ['original_affine','S8']
    if adapter == 'guard':return ['selected_guard','S8']
    names = ['S8']
    if adapter == 'newer':names.append('newer_separation')
    for name in ('finite_exchanges','late_receivers'):
        if name in reg['proof']:names.append(name)
    need(any(n in names for n in ('finite_exchanges','late_receivers')), 'concrete finite geometric proof home')
    return names


def affine_case(case, batch, pins):
    number = int(batch.rsplit('-',1)[1]);stem = case['source'];group = case['group']
    ordinal = case['source_case_ordinal']
    need(type(ordinal) is int and ordinal >= 0, 'exact affine original case ordinal')
    prefix = [] if group == 'root' else path_tokens(group)
    row,ref = source_record(stem,prefix+['cases',ordinal],pins)
    if number in (1,3,8):
        body = row['proof']; candidates = body['maps'];selector = prefix+['cases',ordinal,'proof','maps']
    elif number in (0,4):
        body = row;candidates = row['maps'];selector = prefix+['cases',ordinal,'maps']
    elif number == 5:
        body = row;candidates = row['completions'];selector = prefix+['cases',ordinal,'completions']
    elif number in (6,7):
        kind = source_path(stem,prefix,pins)['kind']
        container = 'completions' if number == 7 and kind == 'D' else 'rows'
        body = row;candidates = row[container];selector = prefix+['cases',ordinal,container]
    elif number == 2:
        body = row;candidates = None;selector = prefix+['cases',ordinal]
        need(body['outcome']['method'] in ('no_direct_contact_to_marked_component','zero_residual_endpoint_extension')
             and type(body['outcome']['maximum_order']) is int and body['outcome']['maximum_order'] < 15,
             'registered exact finite-order exclusion premise')
    else:
        need(number == 9 and type(row['proof']) is dict and type(row['maps']) is list,
             'exact nested sparse BC13 source case')
        body = row['proof'];candidates = row['maps'];selector = prefix+['cases',ordinal,'maps']
    return row,body,candidates,selector,ref


def origins(inner,spec,pins):
    """Yield every export in frozen435e order, without map enumeration."""
    result = inner['result'];adapter = spec['adapter'];batch = inner['batch']
    ordinal = 0
    if adapter == 'proper':
        for pi,packet in enumerate(result['packets']):
            reg = registration(spec,batch,packet['source'])
            for qi,proof in enumerate(packet['proper_map_providers']):
                for mi,mapping in enumerate(proof['maps']):
                    visit();selector = path_tokens(mapping['source_path'])
                    row,ref = source_record(packet['source'],selector,pins)
                    parent = source_path(packet['source'],selector[:-2],pins)
                    method = parent.get('method')
                    if packet['source_section'] == 'selected_proofs':
                        need(method in ('proper_completion_exclusion','impossible_proper_completion','saturated_completion_exclusion'),
                             'exact proper source filter')
                    need(ref['record_sha256'] == mapping['computed_record_sha256'] and mapping['forced_entries']
                         and mapping['geometric_premise'], 'proper complete record/contact association')
                    need(mapping['key']['high_index'] == proof['high_index'] == parent['index'], 'proper original high association')
                    need(parent['core_key'] == mapping['key']['core_key'] and
                         row.get('map',row.get('images')) == mapping['key']['induced_map'],
                         'proper exact original core/full-map association')
                    yield mapping['key'],dict(export_ordinal=ordinal,origin=['result','packets',pi,'proper_map_providers',qi,'maps',mi],
                        source_record=ref,method=method,source_outcome=mapping['outcome']),reg
                    ordinal += 1
    elif adapter in ('whole_schema_safe','newer'):
        if adapter == 'whole_schema_safe':
            need(result['shards'] == 8 and result['shard'] == int(batch[-1])
                 and result['complete_reference_interface_count'] == 564, 'complete whole reference scope')
            packets = result['packets']
            sequence = [(p['source'],proof,['result','packets',pi,'current_funded_whole_model_providers',qi],p)
                        for pi,p in enumerate(packets) for qi,proof in enumerate(p['current_funded_whole_model_providers'])]
        else:
            sequence = []
            for ri,row in enumerate(result['exact_newer_provider_rows']):
                funded(row['direct_current_receipt'],True)
                if row['domain'] != 'CURRENT_SELECTED_WHOLE_FINITE':
                    need(row['domain'] in ('UNIVERSAL_NAMED_CURRENT_OUTPUT','ORIGINAL_IDEAL17_CURRENT_UPDATE')
                         and row['high_row_id_assigned'] is False, 'nonselected newer domains stay separate')
                    continue
                sequence.append((row['source'],row['whole_provider'],['result','exact_newer_provider_rows',ri,'whole_provider'],row))
        for stem,proof,origin,container in sequence:
            reg = registration(spec,batch,stem)
            if adapter == 'newer':
                selector = path_tokens(container['source_path'])
            else:
                selector = path_tokens(proof['actual_source_row_path'])
            row,ref = source_record(stem,selector,pins)
            if adapter == 'newer':
                identifier = container['source_identifier']
                need(identifier == row.get('index',row.get('case',row.get('name','A14'))),
                     'exact newer source identifier at its actual selector')
                if type(identifier) is int:
                    need(identifier == proof['high_index'], 'newer current integer high identity')
                elif identifier in ('D15_1147','D15_1269'):
                    need(proof['high_index'] == {'D15_1147':1147,'D15_1269':1269}[identifier], 'two symbolic D15 high identities')
                else:
                    need(identifier == 'A14' and row['affected_contexts'] == [1639] and proof['high_index'] == 1639,
                         'sole current A14 singleton identity')
            else:
                need(row['index'] == proof['high_index'], 'whole actual source-row/high identity')
            catalogue = models(proof['all_current_funded_source_models'])
            need(proof['automatic_exhaustiveness_claim'] is False and proof['geometric_applicability_premise'],
                 'whole variant exhaustiveness stays an explicit source premise')
            for mi,mapping in enumerate(proof['full_completion_keys']):
                visit();keys = mapping['funded_model_keys']
                need(keys and len(keys) == len(set(keys)) and set(keys) <= set(catalogue), 'THIS full map funded-model subset')
                need(mapping['key']['high_index'] == proof['high_index'], 'whole proof/high identity')
                yield mapping['key'],dict(export_ordinal=ordinal,origin=origin+['full_completion_keys',mi],
                    source_record=ref,method=row.get('method'),funded_model_keys=keys,
                    finite_predicate_scope='Complete source variants, this native full-map realization, same physical source, DISTINCT current corners and whole marked/witness seals; model membership alone is insufficient.'),reg
                ordinal += 1
    elif adapter == 'affine':
        for ci,case in enumerate(result['exact_registered_case_providers']):
            reg = registration(spec,batch,case['source'])
            row,body,candidates,selector,ref = affine_case(case,batch,pins)
            need(row.get('index',row.get('case')) == case['high_index'], 'affine actual source case/high identity')
            catalogue = models(case['current_funded_model_catalog'])
            for mi,mapping in enumerate(case['complete_registered_map_providers']):
                visit();provider = mapping['provider']; key = mapping['key']
                need(key['high_index'] == case['high_index'], 'affine full-map/high identity')
                need(mapping['source_map_path'] == case['group']+'/cases/'+str(case['source_case_ordinal'])+'/maps/'+str(mi),
                     'exact normalized affine source map ordinal')
                if candidates is not None:
                    need(mi < len(candidates), 'complete source map ordinal')
                    candidate = candidates[mi]
                    images = candidate if type(candidate) is list else candidate.get('map',candidate.get('embedding',candidate.get('images')))
                    need(images == key['induced_map'], 'same original affine full map')
                if provider['kind'] == 'EXACT_PROPER_AFFINE_EXCLUSION':
                    if candidates is None:
                        payload = dict(map=key['induced_map'],outcome='actual_finite_component_order_below_affine_target')
                    else:
                        candidate = candidates[mi]
                        # BC13 outer maps are paired with nested physical-entry proofs.
                        if int(batch.rsplit('-',1)[1]) == 9 and type(candidate) is list:
                            nested = body['maps'][mi]
                            payload = nested.get('proof',nested)
                        else:
                            payload = body if type(candidate) is list else candidate.get('proof',candidate)
                    need(sha256(canonical(payload).encode()).hexdigest() == provider['computed_record_sha256']
                         and provider['premise'], 'proper affine exact computed record/premise')
                else:
                    need(provider['kind'] == 'CURRENT_FUNDED_WHOLE_AFFINE_VARIANTS'
                         and provider['source_variant_exhaustiveness_verified_by_this_index'] is False,
                         'explicit affine whole-variant premise')
                    complete,matching = provider['complete_funded_model_keys'],provider['matching_model_keys']
                    need(complete and matching and len(complete) == len(set(complete))
                         and len(matching) == len(set(matching)) and set(matching) <= set(complete) <= set(catalogue),
                         'THIS affine full map complete/matching funded models')
                yield key,dict(export_ordinal=ordinal,origin=['result','exact_registered_case_providers',ci,'complete_registered_map_providers',mi],
                    source_record=ref,source_map_selector=selector+[mi] if candidates is not None else selector,
                    method=body.get('method',body.get('outcome')),provider_kind=provider['kind'],
                    finite_predicate_scope='Exact all-native target/contact/gauge/inverse/inertia/completed-principal and physical-variant propositions; same map and DISTINCT current whole receivers, modulo independent S8.'),reg
                ordinal += 1
    elif adapter == 'old_affine':
        reg = next(r for r in spec['rows'] if r['batch'] == batch)
        for ri,row in enumerate(result['current_exception_providers']):
            funded(row['direct_current_receipt'],True)
            seal = row['sealing_witness']
            need(seal['method'] in ('positive_affine_kernel_and_strict_single_contact_projection','positive_affine_kernel_and_actual_zero_radical')
                 and len(seal['omitted_black']) == 2 and seal['positive_kernel'], 'original affine whole-star seal')
            section,items = ('exceptional_D15','cases') if batch == 'old-d16' else ('order17','replacement_cases')
            rows = source_path(row['old_source'],[section,items],pins)
            selected = [(i,p) for i,p in enumerate(rows) if p['case'] == row['old_case']]
            need(len(selected) == 1, 'one exact old source case')
            oi,_ = selected[0];_,ref = source_record(row['old_source'],[section,items,oi],pins)
            for mi,key in enumerate(row['complete_dropped_map_keys']):
                visit();need(key['high_index'] == row['new_case'], 'old/current full-map high identity')
                yield key,dict(export_ordinal=ordinal,origin=['result','current_exception_providers',ri,'complete_dropped_map_keys',mi],
                    source_record=ref,method=seal['method'],old_case=row['old_case'],
                    finite_predicate_scope='Exact complete old/new marked-core/white/corner/native-map correspondence, positive kernel, D16 strict inverse/signature or D17 actual-zero radical/opposite omitted signs, including indirect tails.'),reg
                ordinal += 1
    else:
        need(adapter == 'guard', 'supported exact provider adapter')
        stem = 'fresh_d13_three_rank15_context_transfer_20261018'
        source_rows = source_path(stem,['selected_proofs'],pins)
        reg = dict(batch=batch,source=stem,proof='selected_guard:Complete selected finite support and arbitrary extensions; fresh local-routing Section5; source review/root6ea')
        methods = {p['method']:p for p in spec['method_rules']}
        for pi,proof in enumerate(result['guarded_providers']):
            matches = [(i,p) for i,p in enumerate(source_rows) if p['index'] == proof['high_index']]
            need(len(matches) == 1 and proof['method'] in methods, 'one exact recognized selected guard source proof')
            si,row = matches[0];_,ref = source_record(stem,['selected_proofs',si],pins)
            need(row['method'] == proof['method'] and row['core_key'] == proof['core_key'], 'same selected guard method/core')
            for mi,mp in enumerate(proof['full_map_keys']):
                visit();yield dict(high_index=proof['high_index'],core_key=proof['core_key'],
                    base_nodes=row['selected_nodes'],selected_marks=row['selected_marks'],target_kind='D',
                    target_order=13,affine=False,induced_map=mp['images'],marked_image=mp['marked_images']),dict(
                    export_ordinal=None,origin=['result','guarded_providers',pi,'full_map_keys',mi],
                    guard_proof_origin=pi,guard_map_origin=mi,source_record=ref,method=proof['method'],
                    guard_hypothesis='TRUE: three DISTINCT ACTUAL whole rank15 A/D witnesses; FALSE has ordinary unit1.',
                    finite_predicate_scope=methods[proof['method']]['proof']+' Complete saved inverse/contact/radical and whole-small predicates; independent S8 escape remains explicit.'),reg


def encode(value):
    data = canonical(value).encode('ascii')
    need(0 < len(data) <= MAX_RESULT, 'complete curator result8MiB cap')
    compressed = zlib.compress(data,6)
    need(zlib.decompress(compressed) == data and loads(data,MAX_RESULT,False) == value,
         'complete curator byte/object round trip')
    return dict(encoding='zlib-base64-json-v1',decoded_bytes=len(data),decoded_sha256=sha256(data).hexdigest(),
        compressed_bytes=len(compressed),data=base64.b64encode(compressed).decode('ascii'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',choices=('catalogue',),required=True)
    parser.add_argument('--parent-lock',action='store_true',required=True)
    args = parser.parse_args()
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('r+') as stream:
        try:
            fcntl.flock(stream.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            pass
        else:
            fcntl.flock(stream.fileno(),fcntl.LOCK_UN)
            raise ValueError('shared wrapper proof lock must already be held')
    before = read(OWN)
    meta = loads(read(INPUTS,INPUT_SHA),exact=False)
    need(meta['schema'] == 'finite57-provider-occurrence-actual-inputs-v1', 'exact immutable curator input schema')
    need(meta['proof_home_registry'] == dict(path=REGISTRY,sha256=REGISTRY_SHA), 'exact reviewed proof-home guide')
    registry = loads(read(REGISTRY,REGISTRY_SHA),exact=False)
    for row in registry['input_registers']:read(row['path'],row['sha256'])
    controllers = meta['controller_source_sha256']
    need(type(controllers) is dict and set(controllers) == CONTROLLERS and controllers[RUNNER] == RUNNER_SHA,
         'complete current ten-controller closure')
    for name,pin in controllers.items():read(name,pin)
    accepted = meta['accepted_assumptions']
    need(type(accepted) is list and len(accepted) == 10 and {p['id'] for p in accepted} == ORDINARY
         and all(p['runtime_regenerated'] is False for p in accepted), 'exact ten ordinary retained premises; no assumed provider acceptance')
    for row in accepted:
        for field in ('source','review','root_acceptance'):read(row[field]['path'],row[field]['sha256'])
    for home in registry['proof_homes'].values():
        for field in ('source','review','root_acceptance','erratum'):
            if field in home:read(home[field]['path'],home[field]['sha256'])
    specs = {p['worker']:p for p in registry['selected_provider_workers']}
    affine = next(p for p in specs.values() if p['adapter'] == 'affine')
    specs[AFFINE9] = dict(affine,worker=AFFINE9,sha256=AFFINE9_SHA)
    expected = set()
    for worker,spec in specs.items():
        batches = spec.get('batches') or sorted({r['batch'] for r in spec['rows']})
        for batch in batches:
            if spec['adapter'] == 'affine' and ((worker == AFFINE9) != (batch == 'affine-9')):continue
            expected.add((worker,batch))
    entries = meta['actual_receipts'];reviews = {p['worker']:p for p in meta['source_review_bindings']}
    need(type(entries) is list and len(entries) == 36 and len(expected) == 36
         and {(p['worker'],p['batch']) for p in entries} == expected
         and len({p['path'] for p in entries}) == len({p['stdout_path'] for p in entries}) == 36,
         'exact complete chosen36 provider receipt set')
    need(len(reviews) == len(meta['source_review_bindings']) == 7 and set(reviews) == set(specs),
         'exact seven reviewed immutable provider sources')
    units = []
    allpins = {}
    for entry in entries:
        spec = specs[entry['worker']]
        inner,pins = actual(entry,reviews[entry['worker']],controllers,spec['sha256'],spec['inner_status'])
        for name,pin in pins.items():
            need(name not in allpins or allpins[name] == pin, 'shared consumed input identities')
            allpins[name] = pin
        units.append((entry,spec,inner,pins))
    need(HIGH in allpins, 'actual authoritative high source consumed')
    highdoc = loads(read(HIGH,allpins[HIGH]),exact=False)['isolation']
    high,cores = highdoc['cases'],highdoc['cores']
    need(len(high) == 2614 and [p['index'] for p in high] == list(range(2614)), 'complete authoritative high-index domain')
    lookup = {};contracts = {};accounting = []
    for entry,spec,inner,pins in units:
        counts = dict(exported=0,eligible=0,unused_out_of_selected_domain=0)
        for key,binding,reg in origins(inner,spec,pins):
            visit();counts['exported'] += 1;keyhash = full_key(key);p = high[key['high_index']]
            need(key['core_key'] == p['core_key'] and key['base_nodes'] == p['marked_component']
                 and key['selected_marks'] == p['preferred_corner'] and p['core_key'] in cores,
                 'same full physical high/core/ordered base/marks')
            h = p.get('component_order_bound')
            eligible = (p['method'] not in ('finite_leaf9_or_ideal_feature5','nine_context_adaptive_finite_payment')
                        and type(h) is int and len(key['base_nodes']) <= key['target_order'] <= h <= 17)
            if not eligible:
                counts['unused_out_of_selected_domain'] += 1
                continue
            counts['eligible'] += 1
            cid = sha256(canonical([spec['adapter'],reg]).encode()).hexdigest()
            if cid not in contracts:
                home_ids = proof_homes(spec['adapter'],entry['batch'],reg)
                contracts[cid] = dict(adapter=spec['adapter'],registry_row=reg,
                    proof_home_ids=home_ids,proof_homes={h:registry['proof_homes'][h] for h in home_ids},
                    ordinary_premise_ids=sorted(ORDINARY),
                    provider_implication='For THIS exact source occurrence/full map, at its complete source-specific finite/variant/whole-witness predicates, NO adjacent white-interior path<=8 in this SAME face implies exclusion or total CURRENT payment>=1.',
                    independent_S8_composition='Assume the same actual face unpaid. Independent conditional S8 forbids every such short path; the reviewed exhaustive no-S8 provider branch excludes its realized map or contradicts unpaidness. No observed path or intrinsic unit is asserted.',
                    numerical_and_geometric_premises_regenerated=False)
            binding['contract_id'] = cid
            identity = entry['worker']+'|'+entry['batch']+'|'+keyhash
            if identity not in lookup:
                lookup[identity] = dict(worker=entry['worker'],batch=entry['batch'],key_sha256=keyhash,
                    key=key,receipt=entry['path'],receipt_sha256=entry['sha256'],bindings=[])
            need(lookup[identity]['key'] == key, 'no key-hash identity collision')
            lookup[identity]['bindings'].append(binding)
        need(counts['exported'] == counts['eligible']+counts['unused_out_of_selected_domain'], 'complete export accounting')
        accounting.append(dict(worker=entry['worker'],batch=entry['batch'],receipt=entry['path'],**counts))
    need(lookup and all(p['bindings'] for p in lookup.values()), 'complete nonempty eligible occurrence table')
    result = dict(schema='finite57-provider-occurrence-bindings-v1',provider_entries=entries,
        complete_chosen_provider_set=True,complete_eligible_occurrence_associations=True,
        unknown_or_unresolved_eligible_occurrences=[],accounting=accounting,lookup=lookup,contracts=contracts,
        input_snapshot=dict(path=INPUTS,sha256=INPUT_SHA),proof_home_registry=dict(path=REGISTRY,sha256=REGISTRY_SHA),
        ordinary_assumptions=accepted,
        premise_scope='Source-bound applicability of the already reviewed implications at complete old finite/native/variant/whole-factor predicates; those predicates remain accepted inputs. This is not target-map enumeration, geometric formalization, an observed S8 escape, selected coverage or global verification.',
        per_source_numerical_premises_regenerated=False,selected_provider_union_complete=False,
        full_transitive_regeneration=False,geometric_formalization=False,global_finite57_verified=False)
    need(read(OWN) == before, 'curator source changed')
    for name,(_,pin) in list(FILES.items()):read(name,pin)
    payload = dict(status='PASS_PROVIDER_OCCURRENCE_ASSOCIATIONS_AT_EXPLICIT_PREMISES',batch=args.batch,
        source_sha256=sha256(before).hexdigest(),input_snapshot=dict(path=INPUTS,sha256=INPUT_SHA),result=result,
        consumed_inputs={name:pin for name,(_,pin) in sorted(FILES.items())},
        selected_provider_union_complete=False,full_transitive_regeneration=False,
        geometric_formalization=False,global_finite57_verified=False)
    output = json.dumps(payload,sort_keys=True)
    if len(output.encode())+1 > MAX_OUTPUT:
        packet = encode(payload)
        output = json.dumps(dict(status=payload['status'],batch=args.batch,source_sha256=payload['source_sha256'],
            complete_payload_encoding=packet),sort_keys=True)
    need(len(output.encode())+1 <= MAX_OUTPUT, 'complete curator output cap; never omit binding records')
    print(output)


if __name__ == '__main__':
    main()
