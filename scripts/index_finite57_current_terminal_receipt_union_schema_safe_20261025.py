#!/usr/bin/env python3
"""PROVISIONAL SOURCE: schema-safe whole-provider actual-key join; S8 stays open.

This reads preserved bounded receipts; it imports no mathematical worker and
does no enumeration. Named source-geometric/computed premises remain premises.
Unsupported priorities and every missing key stay OPEN. This is not a global
finite57, full-terminal or all64-shard verification claim.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
from decimal import Decimal
import fcntl
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import re

ROOT = Path(__file__).resolve().parents[1]
INPUTS = 'results/finite57_current_terminal_union_schema_safe_actual_inputs_20261025.json'
MAX_OUTPUT = 500 * 1024
MAX_TOTAL_READ = 256 * 1024**2
MAX_RECEIPTS = 192
MAX_VISITS = 2_000_000
VISITS = 0
READ_BYTES = 0
READ = {}
SHA = re.compile(r'[0-9a-f]{64}\Z')
KEY_FIELDS = {'high_index', 'core_key', 'base_nodes', 'selected_marks',
              'target_kind', 'target_order', 'affine', 'induced_map', 'marked_image'}
ROUTES = 'scripts/index_finite57_current_terminal_routes_syntax_safe_20261025.py'
ROUTES_SHA = '43a0c4e986cbdb92d112771bf0edc1d128ce2d0feced30fffd9d9825dda7c29d'
RUNNER = 'scripts/run_bounded_research_check_20261018.py'
RUNNER_SHA = '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132'
PRIMARY = 'manuscripts/finite_volume57.tex'
CONTROLLER_PATHS = {
    'scripts/run_proof_batch_20261024.py',
    'scripts/audit_existing_lean_module_20261024.py',
    'scripts/verify_manuscript_audited_20261024.py',
    'scripts/monitor_owned_proof_tree_20261024.py', RUNNER,
    'scripts/run_lean_batches_signal_safe_20261018.py',
    'scripts/run_lean_batches_20261018.py',
    'scripts/build_illustrated_ra_manuscript_20261018.py',
    'scripts/build_compact_arithmetic_nine_field19_20261020.py',
    'scripts/verify_manuscript.py',
}
REGISTERS = (
    ('results/finite57_terminal_adapter_inputs_small_affine_20261025.json',
     '8a0973d7b7ef08231e199bee9f04f45c6a8dc01ba66ca8de8acb2d8e513ee7b7', 'inputs', 110),
    ('results/finite57_selected_affine_completion_chain_source_pins_20261024.json',
     '5f94aff40b8178158e8421c0cbfe1eb377469cd20fdc333aafe6e94e7fe611b5', 'sources', 68),
    ('results/finite57_g4_lower_cycle_source_inputs_20261025.json',
     '9f26c3f91c8c47126e8506e23d7a77497735974e012209a582da1344ddeb232b', 'inputs', 64),
)
# (fixed worker SHA, exact inner status, supported batch names, adapter).
WORKERS = {
    ROUTES: (ROUTES_SHA, 'CURRENT_TERMINAL_OBLIGATION_INDEX_ONLY',
             tuple(f'selected-{i:02d}' for i in range(64)), 'obligations'),
    'scripts/index_finite57_all_proper_finite_completion_providers_20261025.py': (
        '0ef6647b6b9b702f076b264b730caef89a31e40ff2a502bf5a5e03ab1fe29ba3',
        'EXACT_PROPER_FINITE_PROVIDER_INDEX_ONLY', tuple(f'proper-{i}' for i in range(4)), 'proper'),
    'scripts/index_finite57_funded_selected_whole_providers_20261025.py': (
        '0a245deec21b283f0d47af0323cf9557fd6ba4465f6e25a2b7c47897f2a46e86',
        'CURRENT_FUNDED_WHOLE_PROVIDER_INDEX_ONLY', tuple(f'whole-{i}' for i in range(8)), 'disabled_whole_schema'),
    'scripts/index_finite57_funded_selected_whole_providers_schema_safe_20261025.py': (
        '6ace38a526ccf22d9913ffe0b814b7ee4177d74c3be478295a981fef5b7afb03',
        'CURRENT_FUNDED_WHOLE_PROVIDER_INDEX_ONLY', tuple(f'whole-{i}' for i in range(8)), 'whole_schema_safe'),
    'scripts/index_finite57_newer_selected_whole_providers_20261025.py': (
        '54f43fd07f84717d99a14545eb62edb397e74da648648d800ba2a292d8fb39ad',
        'NEWER_CURRENT_WHOLE_PROVIDER_INDEX_ONLY', tuple(f'newer-{i}' for i in range(4)), 'newer'),
    'scripts/index_finite57_selected_affine_completion_providers_20261025.py': (
        '63384dcd33cdfbfc0ae580e745fa339a65cbf8586d871fc7bd634ed437f9e98f',
        'SELECTED_AFFINE_PROVIDER_INDEX_ONLY', tuple(f'affine-{i}' for i in range(10)), 'affine'),
    'scripts/check_finite57_original_affine_seal_correspondences_20261025.py': (
        '5cedf8447f2d827f9c6f833aebdb783a12e69c0f22fe2b34fc2ce5404420db5d',
        'ORIGINAL_AFFINE_CURRENT_SEAL_BINDING_ONLY', ('old-d16', 'old-d17'), 'old_affine'),
    'scripts/check_finite57_selected_guard_bindings_syntax_safe_20261025.py': (
        'c293236c58d51aeb778fe1b40ef97671feaa478b18798a36390bb18d07f4d647',
        'PASS_NAMED_SELECTED_D13_GUARD_BINDING', tuple(f'guard-{i}' for i in range(8)), 'guard'),
}
PARTITION = {'radical': 1329, 'lorentz_residual': 823, 'relevant_component': 259,
             'marked_component_isolated_by_end_inverses': 120,
             'finite_leaf9_or_ideal_feature5': 74, 'nine_context_adaptive_finite_payment': 9}


def need(ok, message):
    if not ok:
        raise ValueError(message)


def charge():
    global VISITS
    VISITS += 1
    need(VISITS <= MAX_VISITS, 'union semantic traversal cap; refuse rather than truncate')


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True)


def unique(pairs):
    out = {}
    for name, value in pairs:
        need(name not in out, 'duplicate JSON field')
        out[name] = value
    return out


def loads(encoded):
    text = encoded.decode('utf-8') if type(encoded) is bytes else encoded
    need(type(text) is str and 0 < len(text.encode()) <= 2 * 1024**2, 'bounded JSON payload')
    depth = 0; quoted = False; escaped = False
    for char in text:
        if quoted:
            if escaped:
                escaped = False
            elif char == '\\':
                escaped = True
            elif char == '"':
                quoted = False
        elif char == '"':
            quoted = True
        elif char in '[{':
            depth += 1; need(depth <= 48, 'bounded JSON nesting')
        elif char in ']}':
            depth -= 1; need(depth >= 0, 'JSON bracket balance')
    need(not quoted and depth == 0, 'complete JSON lexical structure')
    return json.loads(text, object_pairs_hook=unique, parse_float=Decimal,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError('nonfinite JSON')))


def relative_path(name, prefix=None):
    need(type(name) is str and 0 < len(name) <= 240, 'bounded relative path')
    lexical = PurePosixPath(name)
    need(not lexical.is_absolute() and str(lexical) == name and '..' not in lexical.parts,
         'canonical repository path')
    if prefix:
        need(lexical.parts[0] == prefix, 'unexpected input path prefix')
    path = (ROOT / name).resolve()
    need(path.is_relative_to(ROOT) and path.is_file(), 'existing repository file required')
    return path


def raw(name, cap=2 * 1024**2, remember=True):
    global READ_BYTES
    path = relative_path(name)
    with path.open('rb') as stream:
        before = os.fstat(stream.fileno())
        data = stream.read(cap + 1)
        after = os.fstat(stream.fileno())
    identity = lambda st: (st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns)
    need(identity(before) == identity(after) == identity(path.stat()), 'input changed while reading')
    need(0 < len(data) <= cap, 'bounded input read; no truncation')
    READ_BYTES += len(data); need(READ_BYTES <= MAX_TOTAL_READ, 'total input-read cap')
    digest = sha256(data).hexdigest()
    if remember:
        need(name not in READ or READ[name] == digest, 'input hash changed')
        READ[name] = digest
    return data


def digest_string(value):
    need(type(value) is str and SHA.fullmatch(value), 'canonical SHA256 string')
    return value


def scalar(value):
    need(type(value) is str and re.fullmatch(r'-?[0-9]{1,96}(?:/[1-9][0-9]{0,95})?', value),
         'bounded exact rational string')
    return F(value)


def integers(value, bound, length=None):
    need(type(value) is list and (length is None or len(value) == length), 'integer-list schema')
    need(all(type(i) is int and 0 <= i < bound for i in value), 'integer-list domain')
    need(len(value) == len(set(value)), 'integer-list uniqueness')
    return value


def full_key(key):
    charge()
    need(type(key) is dict and set(key) == KEY_FIELDS, 'exact complete-map key fields')
    need(type(key['high_index']) is int and 0 <= key['high_index'] < 2614, 'high-row domain')
    need(type(key['core_key']) is str and 0 < len(key['core_key']) <= 256, 'physical-core key schema')
    base = integers(key['base_nodes'], 64)
    marks = integers(key['selected_marks'], 64, 2)
    need(10 <= len(base) <= 17 and set(marks) <= set(base), 'long selected base/mark domain')
    n = key['target_order']
    need(type(n) is int and len(base) <= n <= 17 and type(key['affine']) is bool, 'target-order/affine domain')
    need(key['target_kind'] in (('A', 'B', 'C', 'D') if key['affine'] else ('A', 'B', 'D')),
         'registered long target family')
    mp = integers(key['induced_map'], n, len(base))
    need(key['marked_image'] == sorted(mp[base.index(m)] for m in marks), 'full map/mark image association')
    return canonical(key)


def input_pins():
    pins = {}
    for path, expected, field, count in REGISTERS:
        data = raw(path, 128 * 1024)
        need(sha256(data).hexdigest() == expected, 'exact source-input register')
        rows = loads(data)[field]
        need(type(rows) is list and len(rows) == count, 'complete register size')
        names = set()
        for row in rows:
            name, digest = row['path'], digest_string(row['sha256'])
            need(name not in names and (name not in pins or pins[name] == digest), 'source-register collision')
            names.add(name); pins[name] = digest
    return pins


def actual(entry, pins, controllers):
    worker, batch = entry['worker'], entry['batch']
    need(worker in WORKERS, 'unregistered receipt worker')
    source_hash, status, batches, adapter = WORKERS[worker]
    need(batch in batches and adapter != 'disabled_whole_schema', 'unsupported/uncleared provider batch')
    need(sha256(raw(worker, 128 * 1024)).hexdigest() == source_hash, 'current worker source identity')
    for path_field, digest_field, size_field, cap in (
            ('path', 'sha256', 'bytes', 2 * 1024**2),
            ('stdout_path', 'stdout_sha256', 'stdout_bytes', 600 * 1024)):
        relative_path(entry[path_field], 'results')
        data = raw(entry[path_field], cap)
        need(type(entry[size_field]) is int and len(data) == entry[size_field]
             and sha256(data).hexdigest() == digest_string(entry[digest_field]), 'actual receipt/log byte binding')
        if path_field == 'path':
            outer_bytes = data
        else:
            log = data.decode('utf-8')
    outer = loads(outer_bytes)
    need(outer['status'] == 'PASS_BOUNDED_RESEARCH_CHECK', 'actual bounded PASS required')
    need(outer['command'] == [worker, '--batch', batch, '--parent-lock'], 'exact actual command')
    need(outer['checker_sha256'] == source_hash and outer['runner_sha256'] == RUNNER_SHA,
         'worker/unchanged bounded-runner source provenance')
    need(type(outer['wall_limit_seconds']) is int and outer['wall_limit_seconds'] == 59
         and type(outer['rss_limit_bytes']) is int and outer['rss_limit_bytes'] == 1024**3
         and type(outer['threads']) is int and outer['threads'] == 1
         and outer['lock_owner'] == 'wrapper', 'unchanged execution gates')
    need(type(outer['elapsed_seconds']) in (int, Decimal) and 0 <= outer['elapsed_seconds'] < 59
         and type(outer['peak_process_tree_rss_bytes']) is int
         and 0 <= outer['peak_process_tree_rss_bytes'] < 1024**3, 'actual wall/RSS bounds')
    need(outer['initial_resources']['free_disk_gib'] >= 40
         and outer['initial_resources']['available_memory_gib'] >= 6, 'proof admission floors')
    output = outer['stdout']
    need(type(output) is str and output.endswith('\n') and len(output.splitlines()) == 1
         and 0 < len(output.encode()) == outer['stdout_bytes'] < 524288, 'complete worker stdout')
    lines = log.splitlines()
    need(len(lines) == 3 and lines[0] + '\n' == output, 'exact public-entry stdout three-line schema')
    summary = loads(lines[1])
    need(summary == {k: v for k, v in outer.items() if k != 'stdout'}, 'outer ledger/stdout summary agreement')
    bridge = loads(lines[2])
    need(bridge['owned_tree_bridge_only'] is True and type(bridge['child_exit_code']) is int
         and bridge['child_exit_code'] == 0,
         'actual owned-tree bridge success')
    need(bridge['cleanup']['signals'] == [] and bridge['cleanup']['verified_live_owned_processes_remaining'] == 0,
         'no signalled or live owned stragglers')
    need(bridge['controller_source_sha256'] == controllers and bridge['assertions_enabled'] is True
         and type(bridge['python_optimize_flag']) is int and bridge['python_optimize_flag'] == 0
         and bridge['child_pythonoptimize'] is None,
         'exact controller/assertion/environment provenance')
    inner = loads(output)
    need(inner['status'] == status and inner['source_sha256'] == source_hash and inner['batch'] == batch,
         'exact actual worker payload status/source/batch')
    consumed = inner['consumed_inputs']
    need(type(consumed) is dict and PRIMARY in consumed
         and all(name in pins and pins[name] == digest for name, digest in consumed.items()),
         'consumed-input equality against accepted source registers')
    return inner, adapter


def funded(model, full=False):
    if full:
        need(model['direct_current_receiver_condition'] is True
             and model['legacy_pointwise_bounds_used'] is False, 'DIRECT current receiver only')
    corners = model['corners'] if full else model['current_corners']
    total = scalar(model['total'] if full else model['current_total'])
    need(type(corners) is list and 1 <= len(corners) <= 64, 'funded distinct-corner list')
    seen = set(); actual_total = F(0)
    for corner in corners:
        charge()
        marks = tuple(sorted(integers(corner['marks'], 64, 2)))
        need(marks not in seen, 'receiver corner counted twice')
        seen.add(marks)
        lower = scalar(corner['weight'] if full else corner['lower'])
        need(0 <= lower <= 1, 'current literal corner range')
        actual_total += lower
    need(total == actual_total and total >= 1, 'complete current receiver sum, never legacy surplus')
    digest_string(model['context_sha256'])


def catalog(models):
    need(type(models) is list and models, 'complete current-funded model catalog')
    result = {}
    for model in models:
        charge()
        funded(model)
        name = digest_string(model['model_key'])
        need(name not in result, 'unique funded model key')
        result[name] = model
    return result


def whole_entries(proof):
    models = catalog(proof['all_current_funded_source_models'])
    need(proof['automatic_exhaustiveness_claim'] is False, 'variant geometry must remain an explicit premise')
    for row in proof['full_completion_keys']:
        keys = row['funded_model_keys']
        need(type(keys) is list and keys and len(keys) == len(set(keys)) and set(keys) <= set(models),
             'exact whole key has a current-funded physical model')
        yield row['key'], 'WHOLE_CURRENT_RECEIVER_AT_SOURCE_EXHAUSTIVENESS_PREMISE'


def providers(inner, adapter, obligations):
    result = inner['result']
    if adapter == 'proper':
        for packet in result['packets']:
            for proof in packet['proper_map_providers']:
                for row in proof['maps']:
                    digest_string(row['computed_record_sha256'])
                    need(row['forced_entries'] and row['geometric_premise'], 'exact proper entry/premise evidence')
                    yield row['key'], 'PROPER_FINITE_EXCLUSION_AT_SOURCE_COMPUTED_PREMISES'
    elif adapter == 'whole_schema_safe':
        need(result['shards'] == 8 and type(result['shard']) is int
             and result['shard'] == int(inner['batch'][-1]), 'exact schema-safe whole shard')
        names = result['complete_packet_keys']
        need(type(names) is list and len(names) == 22 and len(set(names)) == 22
             and all(type(name) is str and name for name in names), 'complete twenty-two source packet keys')
        need(result['complete_reference_interface_count'] == 564,
             'exact source-bound whole interface count')
        packets = result['packets']
        need(type(packets) is list and [p['source'] for p in packets]
             == [name for i, name in enumerate(names) if i % 8 == result['shard']],
             'entire named whole packet shard')
        for packet in packets:
            rows = packet['current_funded_whole_model_providers']
            need(type(rows) is list, 'whole model provider row list')
            for proof in rows:
                need(type(proof['geometric_applicability_premise']) is str
                     and proof['geometric_applicability_premise'], 'whole applicability remains a premise')
                yield from whole_entries(proof)
    elif adapter == 'newer':
        for row in result['exact_newer_provider_rows']:
            funded(row['direct_current_receipt'], full=True)
            if row['domain'] == 'CURRENT_SELECTED_WHOLE_FINITE':
                yield from whole_entries(row['whole_provider'])
            else:
                need(row['domain'] in ('UNIVERSAL_NAMED_CURRENT_OUTPUT', 'ORIGINAL_IDEAL17_CURRENT_UPDATE'),
                     'recognized separated newer interface')
                # No fabricated selected high index or key is assigned here.
    elif adapter == 'affine':
        for case in result['exact_registered_case_providers']:
            rows = case['complete_registered_map_providers']
            models = catalog(case['current_funded_model_catalog']) if case['current_funded_model_catalog'] else {}
            for row in rows:
                proof = row['provider']
                if proof['kind'] == 'EXACT_PROPER_AFFINE_EXCLUSION':
                    digest_string(proof['computed_record_sha256'])
                    need(proof['premise'], 'proper affine geometric premise retained')
                    label = 'PROPER_AFFINE_EXCLUSION_AT_SOURCE_COMPUTED_PREMISES'
                else:
                    need(proof['kind'] == 'CURRENT_FUNDED_WHOLE_AFFINE_VARIANTS'
                         and proof['source_variant_exhaustiveness_verified_by_this_index'] is False,
                         'explicit whole affine variant boundary')
                    complete, matching = proof['complete_funded_model_keys'], proof['matching_model_keys']
                    need(complete and matching and len(complete) == len(set(complete))
                         and len(matching) == len(set(matching)) and set(matching) <= set(complete) <= set(models),
                         'THIS full affine map joins current-funded complete source variants')
                    label = 'WHOLE_AFFINE_CURRENT_RECEIVER_AT_SOURCE_EXHAUSTIVENESS_PREMISE'
                yield row['key'], label
    elif adapter == 'old_affine':
        for row in result['current_exception_providers']:
            funded(row['direct_current_receipt'], full=True)
            seal = row['sealing_witness']
            need(seal['method'] in ('positive_affine_kernel_and_strict_single_contact_projection',
                                    'positive_affine_kernel_and_actual_zero_radical'),
                 'exact visible affine-kernel sealing provider')
            need(len(seal['omitted_black']) == 2 and seal['positive_kernel'], 'visible full-kernel/two-omitted-wall scope')
            for key in row['complete_dropped_map_keys']:
                yield key, 'ORIGINAL_VISIBLE_AFFINE_KERNEL_SEAL_AND_DIRECT_CURRENT_RECEIVER'
    elif adapter == 'guard':
        need(result['shards'] == 8, 'exact eight-shard guard domain')
        for proof in result['guarded_providers']:
            for mp in proof['full_map_keys']:
                need(type(mp['images']) is list and mp['marked_images'] in ([0, 11], [0, 12]),
                     'exact changed D13 map role')
                for obligation in obligations:
                    charge()
                    key = obligation['key']
                    if (obligation['conditional_negative_D13_guard'] is True
                            and key['high_index'] == proof['high_index'] and key['core_key'] == proof['core_key']
                            and key['affine'] is False and key['target_kind'] == 'D' and key['target_order'] == 13
                            and key['induced_map'] == mp['images'] and key['marked_image'] == mp['marked_images']):
                        yield key, 'SELECTED_D13_GUARD_TRUE_EXCLUDED_GUARD_FALSE_UNIT'
    else:
        raise ValueError('unsupported provider extraction; no recursive key scavenging')


def parent_lock():
    with (ROOT / 'formal/.proof-batches/exclusive.lock').open('r+') as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return
        fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
    raise ValueError('shared parent lock must already be held')


def run(batch):
    manifest = loads(raw(INPUTS, 256 * 1024))
    need(manifest['schema'] == 'finite57-current-terminal-actual-receipts-v1', 'exact actual-input manifest schema')
    controllers = manifest['controller_source_sha256']
    need(type(controllers) is dict and set(controllers) == CONTROLLER_PATHS and controllers.get(RUNNER) == RUNNER_SHA,
         'complete ten-controller closure')
    for path, expected in controllers.items():
        relative_path(path, 'scripts')
        need(sha256(raw(path, 128 * 1024)).hexdigest() == digest_string(expected), 'current controller source identity')
    pins = input_pins()
    entries = manifest['actual_receipts']
    need(type(entries) is list and 1 <= len(entries) <= MAX_RECEIPTS, 'bounded actual-receipt inventory')
    names = [(e['worker'], e['batch']) for e in entries]
    need(len(names) == len(set(names)), 'one actual successful receipt per worker/batch')
    need(len({e['path'] for e in entries}) == len(entries)
         and len({e['stdout_path'] for e in entries}) == len(entries), 'unique receipt/log files')
    reviews = manifest['source_review_bindings']
    need(type(reviews) is list, 'source-review bindings are explicit external premises')
    accepted = {}
    for review in reviews:
        worker = review['worker']
        need(worker in WORKERS and worker not in accepted and review['source_sha256'] == WORKERS[worker][0],
             'unique exact reviewed worker source')
        relative_path(review['path'], 'results')
        data = raw(review['path'], 512 * 1024)
        need(sha256(data).hexdigest() == digest_string(review['sha256']), 'source-review document binding')
        accepted[worker] = review
    index_entries = [e for e in entries if e['worker'] == ROUTES and e['batch'] == batch]
    need(len(index_entries) == 1 and ROUTES in accepted, 'one actual reviewed selected-shard receipt required')
    inner, _ = actual(index_entries[0], pins, controllers)
    need(inner['collector_sha256'] == '930ffff5601e6096359d3b2f8846fa6d7d050f95fb2aabe449447e0d3d0bf9ca',
         'unchanged index collector provenance')
    result = inner['result']; shard = int(batch[-2:])
    need(result['shard'] == shard and result['shards'] == 64 and result['complete_authoritative_partition'] == PARTITION,
         'exact named shard of authoritative2614/six-method partition')
    obligations = result['unbound_current_completion_obligations']
    need(type(obligations) is list and len(obligations) <= 5000, 'bounded exact obligation list')
    wanted = {}
    for ordinal, row in enumerate(obligations):
        key = row['key']; token = full_key(key)
        need(key['high_index'] % 64 == shard and token not in wanted, 'disjoint complete-map shard keys')
        need(row['binding_status'] == 'UNBOUND' and 0 <= scalar(row['intrinsic_current_lower']) < 1
             and type(row['conditional_negative_D13_guard']) is bool, 'honest nonunit obligation schema')
        if row['conditional_negative_D13_guard']:
            need(not key['affine'] and (key['target_kind'], key['target_order']) == ('D', 13)
                 and key['marked_image'] in ([0, 11], [0, 12]) and scalar(row['intrinsic_current_lower']) == F(1, 2),
                 'exact conservative guard-TRUE scope')
        wanted[token] = dict(ordinal=ordinal, row=row, providers=[])
    paid = result['direct_whole_factor_unit_maps']
    need(type(paid) is int and paid >= 0 and result['all_induced_maps_visited'] == paid + len(wanted),
         'every visited full map is a direct unit or retained nonunit')
    digest_string(result['direct_map_key_sha256'])
    used = []; unmatched_provider_keys = 0
    for entry in entries:
        if entry['worker'] == ROUTES:
            continue  # Other shards are neither regenerated nor credited here.
        need(entry['worker'] in accepted, 'provider source requires bound different-author/root review')
        provider, adapter = actual(entry, pins, controllers)
        need(provider.get('current_terminal_bridge_complete', False) is False, 'provider is not global closure')
        seen = set(); matched = 0
        for key, nature in providers(provider, adapter, obligations):
            token = full_key(key)
            if (token, nature) in seen:
                continue  # Independent variants can yield the same exact key.
            seen.add((token, nature))
            if token not in wanted:
                unmatched_provider_keys += 1
                continue
            prerequisites = []
            if adapter in ('proper', 'whole_schema_safe', 'newer', 'affine'):
                prerequisites.append(dict(kind='independent_current_S8',
                    status='UNVERIFIED_EXTERNAL_TO_THIS_JOIN', discharged_by_key_match=False))
            wanted[token]['providers'].append(dict(receipt=entry['path'], adapter=adapter,
                nature=nature, application_prerequisites=prerequisites))
            matched += 1
        used.append(dict(receipt=entry['path'], worker=entry['worker'], batch=entry['batch'],
                         exact_matching_keys=matched, receipt_sha256=entry['sha256']))
    covered = sha256(); missing = []; joined = []; matched_count = 0
    for token, row in wanted.items():
        if row['providers']:
            matched_count += 1
            covered.update((token + '\n').encode())
            joined.append(dict(obligation_ordinal=row['ordinal'], full_key_sha256=sha256(token.encode()).hexdigest(),
                               providers=row['providers']))
        else:
            missing.append(dict(obligation_ordinal=row['ordinal'], key=row['row']['key'],
                                intrinsic_current_lower=row['row']['intrinsic_current_lower'],
                                conditional_negative_D13_guard=row['row']['conditional_negative_D13_guard']))
    prefix = result['prefix_obligations']
    need(type(prefix) is list and len(prefix) <= 41, 'exact per-shard prefix obligations')
    seen_prefix = set()
    for row in prefix:
        need(type(row['index']) is int and 0 <= row['index'] < 2614 and row['index'] % 64 == shard,
             'prefix belongs to this exact shard')
        need(row['index'] not in seen_prefix and row['method'] in (
            'finite_leaf9_or_ideal_feature5', 'nine_context_adaptive_finite_payment',
            'independent_S8', 'display_exceeds_proved_whole_bound'), 'unique recognized open prefix')
        seen_prefix.add(row['index'])
    return dict(shard=shard, total_shards_required=64, actual_selected_receipt=index_entries[0],
                direct_unit_map_count=paid, direct_unit_map_key_sha256=result['direct_map_key_sha256'],
                total_nonunit_full_keys=len(wanted), joined_nonunit_keys=matched_count,
                joined_nonunit_key_sha256=covered.hexdigest(), exact_joined_key_bindings=joined,
                uncovered_completion_keys=missing,
                exact_provider_receipts=used, provider_keys_outside_this_shard=unmatched_provider_keys,
                prefix_obligations_still_open=prefix,
                selected_nonunit_key_union_complete=not missing,
                selected_shard_all_routes_closed=False,
                current_terminal_bridge_complete=False, full_finite57_verified=False,
                source_geometric_and_computed_premises_proved_by_this_checker=False,
                unresolved_global_interfaces=['All64 actual selected shards and their exact provider unions',
                    'Current S8/protected74/adaptive9/extension-impossibility prefix applications',
                    'Universal/lower-order theorem-to-full-key applications; G4 requests are not evidence',
                    'Quiet affineB13 named64 consumer/provider joins and all other terminal collector domains',
                    'Complete bounded59--163 and unbounded-tail integration, separate from ambient58'],
                source_review_bindings=reviews)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch', choices=tuple(f'selected-{i:02d}' for i in range(64)), required=True)
    parser.add_argument('--parent-lock', action='store_true', required=True)
    args = parser.parse_args(); parent_lock()
    relative = str(Path(__file__).resolve().relative_to(ROOT))
    before = raw(relative, 128 * 1024)
    result = run(args.batch)
    for path, digest in list(READ.items()):
        need(sha256(raw(path, remember=False)).hexdigest() == digest, 'bound source/actual receipt changed before publication')
    need(raw(relative, 128 * 1024, remember=False) == before, 'union source changed')
    payload = json.dumps(dict(status='CURRENT_SELECTED_ACTUAL_RECEIPT_UNION_WITH_EXPLICIT_OPEN_GATES',
        batch=args.batch, source_sha256=sha256(before).hexdigest(), result=result,
        bound_input_sha256=READ, total_input_bytes_read=READ_BYTES, semantic_visits=VISITS,
        scope='Exact named-shard actual-receipt union at bound source proof/computed premises only. No all-terminal/global PASS, regeneration, new theorem or omitted-prefix discharge.'), sort_keys=True)
    need(len(payload.encode()) <= MAX_OUTPUT, 'union output cap: refuse rather than omit keys')
    print(payload)


if __name__ == '__main__':
    main()
