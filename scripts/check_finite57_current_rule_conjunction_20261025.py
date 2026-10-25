#!/usr/bin/env python3
"""SOURCE ONLY: current-rule conjunction at explicit accepted premises.

This successor preserves435e and historical evidence. It never asserts complete
transitive regeneration, a geometric formalization, or a global finite57 PASS.
Root review and measured bounded control/provider executions are prerequisites.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
from decimal import Decimal
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import types

ROOT = Path(__file__).resolve().parents[1]
UNION = 'scripts/index_finite57_current_terminal_receipt_union_schema_safe_20261025.py'
UNION_SHA = '435e2912d2e1bc82e97da8d33f499dc018802dbd10b2a00368a57e4f14d9fb71'
ROUTE = 'scripts/index_finite57_current_terminal_routes_syntax_safe_20261025.py'
ROUTE_SHA = '43a0c4e986cbdb92d112771bf0edc1d128ce2d0feced30fffd9d9825dda7c29d'
S8 = 'scripts/check_finite57_current_s8_applications_20261025.py'
UNIVERSAL = 'scripts/check_finite57_current_universal_receivers_20261025.py'
QUIET_CASE = 'scripts/check_finite57_final_quiet_cases_20261025.py'
QUIET_PARENT_SHA = '29ddbcf65cae2a7c7ff7703a1b4a81bad3ca1fdb0ade4bdbffda26fb9c8a1631'
LOSSLESS = 'scripts/check_finite57_large_d_receivers_compressed_20261025.py'
LOSSLESS_SHA = '6df8d216a1e2d8473cdfaea0c591fb9d74e7fb8e6d6aadcf6df8fc30ba17e5a0'
AFFINE9 = 'scripts/index_finite57_selected_affine9_schema_safe_providers_20261025.py'
AFFINE9_SHA = '679f792f08adf3f51a6e63d5f0dc3a3f55074c73c382bf35c4513da0fdd3d4f5'
CURATOR = 'scripts/curate_finite57_provider_occurrences_20261025.py'
CURATOR_SHA = 'aefd1a1db6a3edc52846316e969c189ecc12e22028d9151c2a264e0d09c3658c'
CURATOR_INPUT = dict(path='results/finite57_provider_occurrence_actual_inputs_20261025.json',
    sha256='0ec63161c7ea3473b82673bf9c935a4dc6ad7f5b7e13cbe032a81a5054830f48')
QUIET_CASE_INDICES = tuple(3+6*i+j for i in range(8) for j in range(3))
OWN = 'scripts/check_finite57_current_rule_conjunction_20261025.py'
LITERAL = 'scripts/finite57_literal_weights_audit_20261024.py'
LITERAL_SHA = 'f6eba7af77af50b821e9e99c589388c1352f9ba5a99052aaeee64393f5ada1e0'
SELECTED_INPUTS = 'results/finite57_current_rule_selected_inputs_20261025.json'
FINAL_INPUTS = 'results/finite57_current_rule_final_inputs_20261025.json'
POLICY_SOURCE = 'scripts/run_lean_batches_20261018.py'
LEGACY_REGISTRY = 'results/finite57_immutable_historical_nine_receipts_20261025.json'
LEGACY_REGISTRY_SHA = 'a00ae77d7680bc738656b78e1c366f9aa13fc48d428d5a9ff8e9bee3ce8fffe6'
LEGACY_RECORDS = None
MAX_OUTPUT = 500 * 1024
MAX_UNITS = 320
MAX_HASH_BYTES = 256 * 1024**2
HASH_BYTES = 0
FILES = {}
U = None
CODEC = None
LARGE_DECODED = {}
EXPANDED_BYTES = 0
MAX_EXPANDED_BYTES = 128 * 1024**2
SOURCE_DOCUMENTS = {}
# These are EXACT actual worker roles, not labels that establish theorems.
G4 = 'scripts/check_finite57_g4_lower_cycle_terminal_bindings_20261025.py'
FINITE = 'scripts/check_finite57_terminal_collectors_syntax_safe_20261025.py'
DIRECT = 'scripts/check_finite57_terminal_collectors_direct_receivers_20261025.py'
PROTECTED = 'scripts/index_finite57_current_protected74_order_safe_applications_20261025.py'
ADAPTIVE = 'scripts/check_finite57_current_adaptive9_copy_map_safe_20261025.py'
ZERO = 'scripts/check_finite57_current_d13_center_zero_outputs_20261025.py'
QUIET = 'scripts/check_finite57_quiet_b13_current_priority_bridge_20261025.py'
ONE_SIDE = 'scripts/check_finite57_final_guard_one_sided_current_receivers_20261025.py'
D12 = 'scripts/check_finite57_current_d12_terminal_bindings_guard_safe_20261025.py'
D14_CURRENT = 'scripts/check_finite57_d14_current24_witness_receivers_20261025.py'
SOURCE_APPLICATION_CONTROLS = 'scripts/test_finite57_current_rule_source_application_controls_20261025.py'
HIGH_PROJECTION_CONTROLS = 'scripts/test_finite57_source_application_high_projection_controls_20261025.py'
FINAL_INTEGRATION_CONTROLS = 'scripts/test_finite57_old_affine_adaptive_integration_controls_20261025.py'
LEGACY_APPLICATION = 'scripts/finite57_selected_legacy_receipt_applications_20261025.py'
FULL_D_APPLICATION = 'scripts/finite57_full_universal_affine_d_applications_20261025.py'
OLD_AFFINE_APPLICATION = 'scripts/finite57_selected_old_affine_source_applications_20261025.py'
DEEP_BATCHES = tuple(f'deep-{i:02d}' for i in range(16))
G4_BATCHES = ('order10-cycle','order10-b','order10-d','cycle11','order11-bcd',
    'order11-bc','order12-center','order12-near','order12-selected',
    'order13-sparse','cycle12','cycle13','cycle14','cycle16','cycle18-tree','large-d')
UNIVERSAL_BATCHES = ('order12-d','five-local','large-d-0','large-d-1',
                     'large-d-2','large-d-3')
CONTROL_WORKERS = {
 FINAL_INTEGRATION_CONTROLS,
 HIGH_PROJECTION_CONTROLS,
 SOURCE_APPLICATION_CONTROLS,
 'scripts/test_finite57_current_rule_conjunction_release_controls_20261025.py',
 'scripts/test_finite57_current_rule_conjunction_controls_20261025.py',
 'scripts/test_finite57_current_adaptive9_twenty_controls_20261025.py',
 'scripts/test_finite57_protected74_order_and_structure_controls_20261025.py',
 'scripts/test_finite57_proper_contact_associations_20261025.py',
 'scripts/test_finite57_whole_provider_schema_associations_20261025.py',
 'scripts/test_finite57_direct_receiver_condition_20261025.py',
 'scripts/test_finite57_selected_guard_and_partial_routes_syntax_safe_20261025.py',
}
# Whole affine target predicates have the separately accepted source meaning.
# They are applied only to exact full maps/whole target factors.
ROLE_THEOREMS = {
 ('C',10,0,9):'finite57.current.C10_endpoint_zero.whole_face.v1',
 **{('B',11,8,j):'finite57.current.BC11_zero.whole_face.v1' for j in (9,10)},
 **{('C',11,*p):'finite57.current.BC11_zero.whole_face.v1'
    for p in ((0,9),(0,10),(1,10))},
 ('B',12,0,9):'finite57.current.B12_center_terminal.whole_face.v1',
 **{('B',12,8,j):'finite57.current.BC12_near_terminal.whole_face.v1' for j in (10,11)},
 ('C',12,1,10):'finite57.current.C12_near_companion_zero.whole_face.v1',
 **{('D',12,i,j):'finite57.current.affine_order12_D_opposite_leaves.whole_face.v1'
    for i in (8,9) for j in (10,11)},
 **{('B',13,8,j):'finite57.current.BC13_inner_terminal.whole_face.v1' for j in (11,12)},
 ('C',13,2,10):'finite57.current.C13_inner_first_zero.whole_face.v1',
 ('D',13,0,8):'finite57.current.D13_center_zero.whole_face.v1',
 **{('D',13,*p):'finite57.current.D13_center_leaf.whole_face.v1'
    for p in ((0,11),(0,12),(8,9),(8,10))},
 **{('D',14,*p):'finite57.current.D14_near_fork.whole_face.v1'
    for p in ((1,12),(1,13),(8,10),(8,11))},
 **{('D',15,*p):'finite57.current.D15_basic_fractional.whole_face.v1'
    for p in ((2,13),(2,14),(8,11),(8,12))},
}


def need(ok, message):
    if not ok:
        raise ValueError(message)


def canonical(x):
    return json.dumps(x, sort_keys=True, separators=(',', ':'), ensure_ascii=True)


def source(name):
    path = (ROOT / name).resolve()
    need(path.is_relative_to(ROOT) and path.is_file(), 'existing repository source')
    with path.open('rb') as stream:
        data = stream.read(512 * 1024 + 1)
    need(0 < len(data) <= 512 * 1024, 'bounded source size')
    return data


def load(name, path, digest):
    data = source(path)
    need(sha256(data).hexdigest() == digest, 'exact reviewed helper source: ' + path)
    module = types.ModuleType(name)
    module.__file__ = str(ROOT / path)
    exec(compile(data, module.__file__, 'exec'), module.__dict__)
    return module


def codec():
    global CODEC
    if CODEC is None:
        hash_file(LOSSLESS, LOSSLESS_SHA)
        CODEC = load('finite57_complete_result_codec', LOSSLESS, LOSSLESS_SHA)
    return CODEC


def expanded_loads(encoded):
    """The frozen parser's exact JSON rules, with a bounded8MiB envelope."""
    text = encoded.decode('utf-8') if type(encoded) is bytes else encoded
    need(type(text) is str and 0 < len(text.encode()) <= 8*1024**2,
         'bounded expanded JSON payload')
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
    return json.loads(text, object_pairs_hook=U.unique, parse_float=Decimal,
        parse_constant=lambda value: (_ for _ in ()).throw(ValueError('nonfinite JSON')))


def expanded_packet(packet):
    global EXPANDED_BYTES
    decoded = codec().decode_result(packet)
    EXPANDED_BYTES += len(decoded)
    need(EXPANDED_BYTES <= MAX_EXPANDED_BYTES,
         'aggregate complete decoded payload cap; never truncate')
    return expanded_loads(decoded)


def hash_file(name, expected):
    global HASH_BYTES
    U.digest_string(expected)
    path = U.relative_path(name)
    before = path.stat()
    identity = (before.st_dev,before.st_ino,before.st_size,before.st_mtime_ns,before.st_ctime_ns)
    if name in FILES:
        need(FILES[name] == (identity,expected), 'source/input changed or conflicting pin')
        return
    need(0 < before.st_size <= 40 * 1024**2, 'bounded source/computed input')
    HASH_BYTES += before.st_size
    need(HASH_BYTES <= MAX_HASH_BYTES, 'total source recheck cap; never truncate')
    digest = sha256()
    with path.open('rb') as stream:
        while True:
            block = stream.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    after = path.stat()
    need(identity == (after.st_dev,after.st_ino,after.st_size,after.st_mtime_ns,after.st_ctime_ns)
         and digest.hexdigest() == expected, 'actual consumed source/input byte binding')
    FILES[name] = identity, expected


def application_source_document(name,expected):
    """Strict complete source JSON, read only after its exact byte binding.

    The original HIGH objects come from the pinned route domain, not this
    reader. Large arbitrary source files cannot silently enlarge the8MiB
    parse boundary. Cached objects remain subject to final identity checks.
    """
    hash_file(name,expected)
    identity=(name,expected)
    if identity not in SOURCE_DOCUMENTS:
        raw=U.raw(name,8*1024**2,remember=False)
        need(sha256(raw).hexdigest()==expected,'exact complete application-source bytes')
        SOURCE_DOCUMENTS[identity]=expanded_loads(raw)
    hash_file(name,expected)
    return SOURCE_DOCUMENTS[identity]


def source_application_helpers(reviews,cases,cores,c):
    """Load only independently reviewed definitions, with explicit premises."""
    out=[]
    for ordinal,path in enumerate((LEGACY_APPLICATION,FULL_D_APPLICATION)):
        need(path in reviews,'missing exact reviewed source-application helper')
        review=reviews[path]
        hash_file(path,review['source_sha256'])
        obj=load('finite57_source_application_'+str(ordinal),path,review['source_sha256'])
        obj.SOURCE_PINNER=hash_file
        obj.SOURCE_READER=application_source_document
        obj.APPLICATION_SOURCE_SHA=review['source_sha256']
        obj.SOURCE_ACCEPTANCE=dict(path=review['path'],sha256=review['sha256'])
        if path==LEGACY_APPLICATION:
            # Reuse the route's exact already-read HIGH projection. The pure
            # helper still pins the entire original file on every lookup.
            need(obj.HIGH=='results/'+c.HIGH+'.json' and
                 c.PINS[obj.HIGH]==obj.PINS[obj.HIGH],
                 'same original HIGH source for route and source applications')
            hash_file(obj.HIGH,obj.PINS[obj.HIGH])
            obj.CACHE[obj.HIGH]=dict(isolation=dict(cases=cases,cores=cores))
        out.append(obj)
    return out


def old_affine_application_helper(reviews,legacy):
    """Add reviewed old-affine transports without changing the two-helper API."""
    path=OLD_AFFINE_APPLICATION
    need(path in reviews,'missing exact reviewed old-affine source-application helper')
    review=reviews[path]
    hash_file(path,review['source_sha256'])
    obj=load('finite57_old_affine_source_applications',path,review['source_sha256'])
    obj.SOURCE_PINNER=hash_file
    obj.SOURCE_READER=application_source_document
    obj.APPLICATION_SOURCE_SHA=review['source_sha256']
    obj.SOURCE_ACCEPTANCE=dict(path=review['path'],sha256=review['sha256'])
    obj.LEGACY_HELPER=legacy
    return obj


def assumptions(meta):
    rows = meta['accepted_assumptions']
    need(type(rows) is list and rows, 'explicit accepted proof/input assumptions required')
    result = {}
    allowed = {'actual_geometry','classification','complete_finite_data',
               'source_geometric_implication','inherited_execution','rule_identity',
               'bounded_tail','unbounded_tail','incidence_contradiction'}
    for row in rows:
        need(type(row) is dict and row['kind'] in allowed and
             type(row['id']) is str and row['id'] not in result and
             row['acceptance_scope'] and row['statement'] and row['runtime_regenerated'] is False,
             'structured retained assumption; no fabricated regeneration')
        for field in ('source','review','root_acceptance'):
            pin = row[field]
            hash_file(pin['path'], pin['sha256'])
        result[row['id']] = row
    return result


def manifest(batch):
    name = FINAL_INPUTS if batch == 'final' else SELECTED_INPUTS
    raw = U.raw(name, 512 * 1024)
    meta = U.loads(raw)
    need(meta['schema'] == 'finite57-current-rule-actual-receipts-v1',
         'exact current-rule actual manifest schema')
    need(meta['controller_source_sha256'].keys() == U.CONTROLLER_PATHS,
         'complete current ten-controller closure')
    for path, digest in meta['controller_source_sha256'].items():
        hash_file(path,digest)
    rows = meta['actual_receipts']
    need(type(rows) is list and 1 <= len(rows) <= MAX_UNITS, 'bounded actual unit inventory')
    need(len({(r['worker'],r.get('batch')) for r in rows}) == len(rows) and
         len({r['path'] for r in rows}) == len(rows) and
         len({r['stdout_path'] for r in rows}) == len(rows), 'unique actual workers/ledgers/stdout')
    reviews = {}
    for row in meta['source_review_bindings']:
        worker = row['worker']
        need(worker not in reviews, 'one exact review binding per worker')
        U.digest_string(row['source_sha256'])
        hash_file(worker,row['source_sha256'])
        hash_file(row['path'],row['sha256'])
        reviews[worker] = row
    need(rows and all(r['worker'] in reviews for r in rows), 'reviewed actual sources only')
    need(OWN in reviews and reviews[OWN]['source_sha256']==sha256(source(OWN)).hexdigest(),
         'current conjunction source requires its own explicit accepted review binding')
    return meta, reviews, assumptions(meta), dict(path=name,sha256=sha256(raw).hexdigest())


def historical_receipts():
    global LEGACY_RECORDS
    if LEGACY_RECORDS is None:
        raw=U.raw(LEGACY_REGISTRY,128*1024)
        need(sha256(raw).hexdigest()==LEGACY_REGISTRY_SHA,'immutable historical nine-receipt registry')
        doc=U.loads(raw)
        need(doc['schema']=='finite57-immutable-historical-nine-receipts-v1' and
             len(doc['receipts'])==doc['receipt_count']==9,
             'exact historical nine-receipt inventory')
        rows=doc['receipts']
        need(len({canonical(r['entry']) for r in rows})==9,'unique historical receipt identities')
        for row in rows:
            need(row['controller_source_sha256'].keys()==U.CONTROLLER_PATHS and
                 row['runner_sha256']==U.RUNNER_SHA and
                 row['controller_source_sha256'][U.RUNNER]==row['runner_sha256'],
                 'exact historical controller and runner provenance')
            for digest in row['controller_source_sha256'].values():
                U.digest_string(digest)
        LEGACY_RECORDS=rows
    return LEGACY_RECORDS


def match_historical_receipt(entry,records):
    matches=[r for r in records if canonical(r['entry'])==canonical(entry)]
    need(len(matches)<=1,'unambiguous historical receipt identity')
    return matches[0] if matches else None


def percentage_policy(source_sha):
    U.digest_string(source_sha)
    return dict(id='percentage-v1-20261025',source_sha256=source_sha,
        admission_disk_percent=2,admission_memory_percent=15,
        abort_disk_percent=1,abort_memory_percent=6,cleanup_disk_percent=2,rounding_digits=6)


def percentage_resources(resources,source_sha):
    expected=percentage_policy(source_sha);policy=resources['resource_policy']
    need(type(policy) is dict and policy.keys()==expected.keys() and
         type(policy['id']) is str and policy['id']==expected['id'] and
         type(policy['source_sha256']) is str and policy['source_sha256']==source_sha and
         all(type(policy[k]) is int and policy[k]==v for k,v in expected.items()
             if k not in ('id','source_sha256')),'exact current percentage policy/source/thresholds')
    fields=('disk_free_bytes','disk_total_bytes','available_memory_bytes','physical_memory_bytes')
    need(all(type(resources[k]) is int for k in fields),'integer resource bytes, never booleans')
    disk,total,memory,physical=(resources[k] for k in fields)
    need(0<=disk<=total and total>0 and 0<=memory<=physical and physical>0,
         'positive capacities and bounded free/available bytes')
    need(100*disk>=2*total and 100*memory>=15*physical,
         'percentage admission uses exact byte ratios')
    for name,available,capacity in (('free_disk_percent',disk,total),
                                    ('available_memory_percent',memory,physical)):
        value=resources[name]
        need(type(value) in (int,Decimal) and 0<=value<=100,
             'finite numeric resource percentage')
        shown=F(value)
        need((shown*10**6).denominator==1 and
             abs(shown-F(100*available,capacity))<=F(1,2*10**6),
             'six-digit percentages agree with exact bytes')


def receipt_resources(entry,resources,controllers):
    if 'resource_policy' in resources:
        percentage_resources(resources,controllers[POLICY_SOURCE])
        return controllers,controllers[U.RUNNER]
    historical=match_historical_receipt(entry,historical_receipts())
    need(historical is not None,'unknown policy-less receipt is not historical evidence')
    need(type(resources['free_disk_gib']) in (int,Decimal) and
         type(resources['available_memory_gib']) in (int,Decimal) and
         resources['free_disk_gib']>=40 and resources['available_memory_gib']>=6,
         'historical actual40/6 admission remains unchanged')
    return historical['controller_source_sha256'],historical['runner_sha256']


def actual(entry, reviews, controllers):
    worker, batch = entry['worker'], entry.get('batch')
    review = reviews[worker]
    worker_sha = review['source_sha256']
    hash_file(worker,worker_sha)
    pieces = []
    for field, hfield, sizefield, cap in (
        ('path','sha256','bytes',2*1024**2),
        ('stdout_path','stdout_sha256','stdout_bytes',600*1024)):
        U.relative_path(entry[field],'results')
        data = U.raw(entry[field],cap)
        need(type(entry[sizefield]) is int and len(data) == entry[sizefield] and
             sha256(data).hexdigest() == U.digest_string(entry[hfield]), 'actual ledger/stdout byte binding')
        pieces.append(data)
    outer = U.loads(pieces[0]); lines = pieces[1].decode().splitlines()
    expected_args = [worker,'--parent-lock'] if batch is None else [worker,'--batch',batch,'--parent-lock']
    need(outer['status'] == 'PASS_BOUNDED_RESEARCH_CHECK' and
         outer['command'] == expected_args and outer['checker_sha256'] == worker_sha,
         'actual exact bounded command/source PASS')
    receipt_controllers,runner_sha=receipt_resources(entry,outer['initial_resources'],controllers)
    need(outer['runner_sha256']==runner_sha,'actual receipt-policy runner provenance')
    need(type(outer['wall_limit_seconds']) is int and outer['wall_limit_seconds'] == 59 and
         type(outer['rss_limit_bytes']) is int and outer['rss_limit_bytes'] == 1024**3 and
         type(outer['threads']) is int and outer['threads'] == 1 and outer['lock_owner'] == 'wrapper' and
         type(outer['elapsed_seconds']) in (int,Decimal) and 0 <= outer['elapsed_seconds'] < 59 and
         type(outer['peak_process_tree_rss_bytes']) is int and 0 <= outer['peak_process_tree_rss_bytes'] < 1024**3,
         'unchanged measured59s/1GiB/one-thread wrapper lock')
    output = outer['stdout']
    need(type(output) is str and output.endswith('\n') and len(output.splitlines()) == 1 and
         len(output.encode()) == outer['stdout_bytes'] < 512*1024 and
         len(lines) == 3 and lines[0]+'\n' == output, 'complete public three-line stdout')
    need(U.loads(lines[1]) == {k:v for k,v in outer.items() if k != 'stdout'}, 'outer/log equality')
    bridge = U.loads(lines[2])
    need(bridge['owned_tree_bridge_only'] is True and type(bridge['child_exit_code']) is int and bridge['child_exit_code'] == 0 and
         bridge['cleanup']['signals'] == [] and bridge['cleanup']['verified_live_owned_processes_remaining'] == 0 and
         bridge['controller_source_sha256'] == receipt_controllers and bridge['assertions_enabled'] is True and
         type(bridge['python_optimize_flag']) is int and bridge['python_optimize_flag'] == 0 and bridge['child_pythonoptimize'] is None,
         'actual owned process-tree/assertion/source closure')
    inner = U.loads(output)
    if worker in (OWN,CURATOR) and 'complete_payload_encoding' in inner:
        need(set(inner) == {'status','batch','source_sha256','complete_payload_encoding'},
             'exact complete-payload envelope')
        expanded = expanded_packet(inner['complete_payload_encoding'])
        need(type(expanded) is dict and all(expanded.get(k) == inner[k]
             for k in ('status','batch','source_sha256')),
             'complete payload agrees with actual envelope status/source/batch')
        inner = expanded
    # Two existing association drivers distinguish their own test hash from
    # the adapter source hash; neither field is silently relabelled.
    own_sha = inner.get('test_sha256',inner.get('source_sha256'))
    need(own_sha == worker_sha and inner.get('batch') == batch, 'actual inner source/batch')
    expected_status = entry['inner_status']
    need(type(expected_status) is str and inner['status'] == expected_status,
         'registered exact actual inner status, never status-prefix acceptance')
    if worker in U.WORKERS:
        need(worker_sha == U.WORKERS[worker][0] and expected_status == U.WORKERS[worker][1],
             'frozen435e provider contract')
    else:
        need(expected_status in KNOWN_STATUSES.get(worker,()), 'supported current-rule actual worker/status')
    if worker == AFFINE9:
        need(worker_sha == AFFINE9_SHA and batch == 'affine-9',
             'exact separate sparse BC13 successor, never old affine0..8')
    if worker == CURATOR:
        need(worker_sha == CURATOR_SHA and batch == 'catalogue',
             'exact complete provider-occurrence curator source and batch')
    consumed = inner.get('consumed_inputs',inner.get('read_inputs',inner.get('inputs',{})))
    need(type(consumed) is dict, 'explicit consumed source map')
    for path,digest in consumed.items():
        hash_file(path,digest)
    if worker not in CONTROL_WORKERS:
        need('manuscripts/finite_volume57.tex' in consumed or worker == OWN,
             'current literal canonical source must be consumed')
    if 'input_snapshot' in inner:
        hash_file(inner['input_snapshot']['path'],inner['input_snapshot']['sha256'])
    return inner


KNOWN_STATUSES = {
 FINAL_INTEGRATION_CONTROLS: ('PASS_OLD_AFFINE_ADAPTIVE_INTEGRATION_CONTROLS_ONLY',),
 HIGH_PROJECTION_CONTROLS: ('PASS_CURRENT_RULE_HIGH_PROJECTION_CONTROLS_ONLY',),
 SOURCE_APPLICATION_CONTROLS: ('PASS_CURRENT_RULE_SOURCE_APPLICATION_CONTROLS_ONLY',),
 D14_CURRENT: ('PASS_CURRENT_D14_ALL24_WITNESS_RECEIVERS_AT_EXPLICIT_PREMISES',),
 CURATOR: ('PASS_PROVIDER_OCCURRENCE_ASSOCIATIONS_AT_EXPLICIT_PREMISES',),
 LOSSLESS: ('PASS_COMPLETE_LARGE_D_RECEIVER_RESULT_LOSSLESSLY_ENCODED',),
 AFFINE9: ('SELECTED_AFFINE_PROVIDER_INDEX_ONLY',),
 DIRECT: ('PASS_NAMED_CURRENT_TERMINAL_ADAPTER',),
 FINITE: ('PASS_NAMED_CURRENT_TERMINAL_ADAPTER',),
 G4: ('PASS_G4_CURRENT_TERMINALS_DECLARED_SCOPE',),
 PROTECTED: ('CURRENT_PROTECTED74_APPLICATION_INDEX_WITH_EXPLICIT_OPEN_PRIORITIES',),
 ADAPTIVE: ('PASS_CURRENT_ADAPTIVE9_LITERAL_BINDING_WITH_EXPLICIT_S8_PREREQUISITE',),
 ZERO: ('PASS_CURRENT_D13_CENTER_ZERO_OUTPUT_BINDING_AT_EXPLICIT_PREMISES',),
 QUIET: ('PASS_NAMED_QUIET_B13_CURRENT_PRIORITY_BRIDGE',),
 ONE_SIDE: ('PASS_NAMED_FINAL_GUARD_EIGHT_ONE_SIDED_CURRENT_RECEIVERS',),
 D12: ('PASS_GUARD_SAFE_CURRENT_D12_TERMINAL_LITERAL_BINDING',),
 S8: ('PASS_CURRENT_S8_CONDITIONAL_APPLICATION_BINDINGS','PASS_CURRENT_S8_APPLICATION_ASSOCIATION_CONTROLS'),
 UNIVERSAL: ('PASS_NAMED_CURRENT_UNIVERSAL_RECEIVER_JOINS_AT_DECLARED_PREMISES',),
 QUIET_CASE: ('PASS_CURRENT_FINAL_QUIET_CASE_AT_DECLARED_PREMISES',),
 OWN: ('CURRENT_RULE_SELECTED_COMPLETE_AT_EXPLICIT_PREMISES',),
 'scripts/test_finite57_current_rule_conjunction_controls_20261025.py': ('PASS_CURRENT_RULE_CONJUNCTION_CONTROLS_ONLY',),
 'scripts/test_finite57_current_rule_conjunction_release_controls_20261025.py': ('PASS_CURRENT_RULE_CONJUNCTION_CONTROLS_ONLY',),
 'scripts/test_finite57_current_adaptive9_twenty_controls_20261025.py': ('PASS_ADAPTIVE9_ASSOCIATION_CONTROLS_ONLY',),
 'scripts/test_finite57_protected74_order_and_structure_controls_20261025.py': ('PASS_PROTECTED74_ORDER_AND_STRUCTURE_CONTROLS_ONLY',),
 'scripts/test_finite57_proper_contact_associations_20261025.py': ('PASS_TINY_PROPER_CONTACT_ASSOCIATIONS',),
 'scripts/test_finite57_whole_provider_schema_associations_20261025.py': ('PASS_TINY_WHOLE_PROVIDER_SCHEMA_ASSOCIATIONS',),
 'scripts/test_finite57_direct_receiver_condition_20261025.py': ('PASS_DIRECT_CURRENT_RECEIVER_REGRESSIONS',),
 'scripts/test_finite57_selected_guard_and_partial_routes_syntax_safe_20261025.py': ('PASS_TINY_PARTIAL_ROUTE_GUARD_REGRESSIONS',),
}


def inventory(meta,reviews):
    out = {}
    for entry in meta['actual_receipts']:
        out[entry['worker'],entry.get('batch')] = (entry,actual(entry,reviews,meta['controller_source_sha256']))
    return out


def required(units,worker,batch):
    key = worker,batch
    need(key in units, 'missing actual unit: '+worker+' '+str(batch))
    return units[key][1]


def large_d_receipt(units,shard):
    """One genuine four-shard mode: old complete JSON or its exact encoding."""
    need(type(shard) is int and 0 <= shard < 4, 'exact large-D shard')
    legacy = {(UNIVERSAL,f'large-d-{i}') for i in range(4)}
    encoded = {(LOSSLESS,f'large-d-{i}') for i in range(4)}
    need(not (legacy.intersection(units) and encoded.intersection(units)),
         'do not mix original and encoded large-D receipts')
    key = (LOSSLESS if encoded.intersection(units) else UNIVERSAL,f'large-d-{shard}')
    need(key in units, 'missing actual complete large-D shard')
    entry, inner = units[key]
    if key[0] == UNIVERSAL:
        return entry, inner
    need(inner['parent_source_sha256'] == QUIET_PARENT_SHA and
         inner['consumed_inputs'].get(UNIVERSAL) == QUIET_PARENT_SHA,
         'encoded large-D uses the exact unchanged parent')
    identity = (entry['path'],entry['sha256'])
    if identity not in LARGE_DECODED:
        result = expanded_packet(inner['result_encoding'])
        need(type(result) is dict, 'complete decoded large-D result object')
        LARGE_DECODED[identity] = dict(inner,result=result)
    return entry, LARGE_DECODED[identity]


def provider_adapter(worker):
    if worker == AFFINE9:
        return 'affine'
    return U.WORKERS[worker][3] if worker in U.WORKERS else None


def provider_curation(units,accepted):
    entry,doc = units.get((CURATOR,'catalogue'),(None,None))
    need(entry is not None, 'missing actual complete provider-occurrence curation')
    result = doc['result']
    need(result['schema'] == 'finite57-provider-occurrence-bindings-v1' and
         result['complete_chosen_provider_set'] is True and
         result['complete_eligible_occurrence_associations'] is True and
         result['unknown_or_unresolved_eligible_occurrences'] == [] and
         result['input_snapshot'] == CURATOR_INPUT and doc['input_snapshot'] == CURATOR_INPUT,
         'complete fixed source-occurrence catalogue at explicit premises')
    for field in ('per_source_numerical_premises_regenerated','selected_provider_union_complete',
                  'full_transitive_regeneration','geometric_formalization','global_finite57_verified'):
        need(result[field] is False, 'curation cannot claim selected coverage or geometry')
    actual = [e for (worker,batch),(e,inner) in units.items()
              if provider_adapter(worker) is not None and worker != ROUTE]
    recorded = result['provider_entries']
    need(type(recorded) is list and len(recorded) == len(actual) == 36 and
         len({canonical(e) for e in recorded}) == 36 and
         {canonical(e) for e in recorded} == {canonical(e) for e in actual},
         'every and only actual chosen provider entry has complete curation')
    premise = accepted['selected-provider-geometry']
    need(premise['actual_curator_receipt'] == entry and
         premise['provider_catalogue'] == CURATOR_INPUT,
         'accepted premise binds the exact actual curation and fixed catalogue')
    ordinary = result['ordinary_assumptions']
    need(type(ordinary) is list and len(ordinary) == 10 and
         len({p['id'] for p in ordinary}) == 10 and
         {p['id']:canonical(p) for p in ordinary} ==
         {k:canonical(p) for k,p in accepted.items() if k != 'selected-provider-geometry'},
         'curation and selected execution retain identical ordinary premises')
    need(type(result['lookup']) is dict and result['lookup'] and
         type(result['contracts']) is dict and result['contracts'],
         'complete nonempty occurrence/proof-contract tables')
    return result


def curated_occurrences(catalogue,entry,key):
    keyhash = sha256(U.full_key(key).encode()).hexdigest()
    identity = entry['worker']+'|'+entry['batch']+'|'+keyhash
    need(identity in catalogue['lookup'], 'uncurated exact full provider key')
    row = catalogue['lookup'][identity]
    need(row['worker'] == entry['worker'] and row['batch'] == entry['batch'] and
         row['receipt'] == entry['path'] and row['receipt_sha256'] == entry['sha256'] and
         row['key_sha256'] == keyhash and canonical(row['key']) == canonical(key),
         'same actual provider receipt and complete physical nine-field key')
    bindings = row['bindings']
    need(type(bindings) is list and bindings and
         len({canonical(b['origin']) for b in bindings}) == len(bindings),
         'all distinct stable source occurrences of this provider key')
    for binding in bindings:
        contract = catalogue['contracts'].get(binding['contract_id'])
        need(contract is not None and contract['adapter'] == provider_adapter(entry['worker']) and
             contract['proof_home_ids'] and contract['proof_homes'] and
             contract['numerical_and_geometric_premises_regenerated'] is False and
             binding['source_record']['source_sha256'] and binding['source_record']['record_sha256'],
             'exact reviewed occurrence theorem and retained finite predicates')
    return bindings


def proof_inputs(accepted):
    needed = {'actual-star','native-classification','complete-short-path-dag',
              'generic-extension-bounds','selected-provider-geometry','universal-receiver-geometry',
              'bounded59-163','unbounded-tail','global-incidence','current-rule-identity',
              'literal35-retained'}
    need(needed <= set(accepted), 'missing explicit accepted source/computed premise')
    return [accepted[k] for k in sorted(needed)]


def runtime_spine(units):
    # No protected/adaptive/largeD theorem establishes independentS8.
    required(units,S8,'association-controls')
    required_g4 = ('order10-cycle','cycle11','order10-b','order10-d',
                   'order11-bc','order11-bcd','cycle12','cycle14','cycle16','cycle18-tree')
    for batch in required_g4:
        required(units,G4,batch)
    zero = required(units,ZERO,'d13-center-zero')['universal_provider']
    need(zero['theorem'] == 'finite57.current.D13_center_zero.whole_face.v1' and
         zero['requires_independent_S8'] is False and
         zero['canonical_marked_images'] == [[0,8]], 'independent D13zero primitive')
    deep_keys=set(); ranks_seen=set()
    for residue,batch in enumerate(DEEP_BATCHES):
        doc=required(units,S8,batch)
        ranks=[n for n in range(9,59) if n%16==residue]
        need(doc['rank_partition']==dict(modulus=16,residue=residue,ranks=ranks),
             'complete exact deep rank partition')
        need(not ranks_seen.intersection(ranks),'disjoint deep rank domains')
        ranks_seen.update(ranks)
        for app in doc['deep_applications']:
            request=app['application_key']['original_deep_request']
            token=canonical(request)
            need(token not in deep_keys and request['factor_rank'] in ranks and
                 request['factor_kind'] in ('A','B','D') and request['marked_distance'] in (6,7,8)
                 and request['status']=='UNBOUND_G6_CURRENT_OUTPUT_BRIDGE' and
                 request['current_intrinsic_weight']=='0' and request['ambient_dimension']==58
                 and len(request['feature_depths'])==2 and min(request['feature_depths'])>=6,
                 'complete unique exact deep request and typed hypotheses')
            deep_keys.add(token)
    need(ranks_seen==set(range(9,59)),'all deep finite ranks9..58')
    tree=required(units,S8,'tree-immediate')
    branches=[canonical(p['application_key']['source_branch']) for p in tree['tree_applications']]
    profiles=('isolated','outer_end3','outer_end4','outer_next3','inner_end3',
              'inner_next3','both_ends3','inner_end4')
    expected={canonical(dict(P=p,c=a,d=b)) for p in (6,7,8) for a in profiles for b in profiles
              if 'inner_next3' in (a,b)}
    need(len(branches)==len(set(branches)) and set(branches)==expected,
         'every exact tree immediate branch, never a count-only inference')
    cfirst=required(units,UNIVERSAL,'five-local')['result']['theorem_applications']
    matches=[p for p in cfirst if p['theorem']=='finite57.current.C13_inner_first_zero.whole_face.v1'
             and p['application_key']==['C',13,2,10]]
    need(len(matches)==1 and matches[0]['independent_S8_dependency'] is False,
         'exact independent C-FIRST application; not sparseB13')
    # The C-FIRST clause is a distinct source theorem, never sparseB13.
    return dict(theorem='finite57.current.S8.whole_face.v1',
        status='RUNTIME_APPLICATION_INTERFACES_AT_ACCEPTED_COMPLETE_SOURCE_PREMISES',
        dag=['S5','F6','S6','F7','S7','F8','S8'],
        actual_units=[[G4,b] for b in required_g4]+[[ZERO,'d13-center-zero']]+
                     [[S8,b] for b in DEEP_BATCHES+('tree-immediate',)],
        protected74_used_as_premise=False,adaptive9_used_as_premise=False,
        largeD_half_used_as_premise=False,order12_D_used_as_premise=False)


def receiver_sum(row):
    corners=row['current_receivers']
    marks=[tuple(sorted(p['current']['marks'])) for p in corners]
    need(corners and len(marks)==len(set(marks)),'distinct original receiver corners')
    total=sum((F(p['current']['weight']) for p in corners),F())
    need(total==F(row['current_total']) and total>=1,'actual complete current receiver sum')
    U.digest_string(row['physical_context_sha256'])


def old_largeD_priority(mapping,W,units):
    U.full_key(mapping['key'])
    theorem=mapping['conditional_whole_face_priority']
    if theorem is None:
        return None
    application=whole_application(mapping['key'],W,units)
    need(application['theorem']==theorem,
         'exact OLD largeD full-map conditional priority theorem')
    return dict(requested_theorem=theorem,application=application)


def retain_selected_assumptions(accepted,rows):
    need(type(rows) is list and rows,'selected receipt must retain its ordinary premises')
    seen=set()
    for row in rows:
        need(type(row) is dict and type(row.get('id')) is str and
             row['id'] not in seen and row.get('runtime_regenerated') is False,
             'unique explicit selected premise records without regeneration')
        name=row['id'];seen.add(name)
        need(name in accepted,'selected premise missing from final accepted premises')
        need(canonical(row)==canonical(accepted[name]),
             'selected/final accepted premise conflict')
    return sorted(seen)


def quiet_applications(units):
    """One actual full parent OR all24 actual case receipts, never a mixture."""
    full_key = UNIVERSAL, 'final-quiet'
    case_keys = [(QUIET_CASE, f'quiet-case-{index:02d}') for index in QUIET_CASE_INDICES]
    observed_cases = {key for key in units if key[0] == QUIET_CASE}
    need(observed_cases <= set(case_keys), 'only the exact24 quiet case batch names')
    any_case = bool(observed_cases)
    need(not (full_key in units and any_case), 'do not mix full and partial quiet receipts')
    witnesses = [(a,b,c) for a in ('A','D') for b in ('A','D') for c in ('A','D')]
    expected = {(index,*kinds) for index in QUIET_CASE_INDICES for kinds in witnesses}
    evidence = []
    if full_key in units:
        doc = required(units,*full_key)
        need(doc['source_sha256'] == QUIET_PARENT_SHA and
             doc['result']['quiet24_complete_at_declared_premises'] is True,
             'genuine unchanged full quiet parent receipt')
        rows = doc['result']['theorem_applications']
        evidence.append(dict(worker=full_key[0],batch=full_key[1],receipt=units[full_key][0]['path']))
        mode = 'ONE_ACTUAL_FULL_PARENT_RECEIPT'
    else:
        rows = []
        consumed = None
        for index,key in zip(QUIET_CASE_INDICES,case_keys,strict=True):
            doc = required(units,*key)
            part = doc['result']
            need(doc['base_universal_source_sha256'] == QUIET_PARENT_SHA and
                 type(part['case_index']) is int and part['case_index'] == index and
                 part['original24_case_domain'] == list(QUIET_CASE_INDICES) and
                 part['case_complete_at_declared_premises'] is True and
                 part['quiet24_complete_at_declared_premises'] is False,
                 'THIS exact quiet case and its partial scope')
            if consumed is None:
                consumed = doc['consumed_inputs']
            need(doc['consumed_inputs'] == consumed,
                 'same exact rule/source/premise input pins across quiet cases')
            apps = part['theorem_applications']
            keys = [tuple(p['application_key']) for p in apps]
            need(len(keys) == len(set(keys)) == 8 and
                 set(keys) == {(index,*kinds) for kinds in witnesses},
                 'every ordered witness key of THIS quiet case')
            rows.extend(apps)
            evidence.append(dict(worker=key[0],batch=key[1],receipt=units[key][0]['path']))
        mode = 'TWENTY_FOUR_ACTUAL_CASE_RECEIPTS'
    keys = [tuple(p['application_key']) for p in rows]
    need(len(keys) == len(set(keys)) == 192 and set(keys) == expected,
         'complete exact192 quiet product, not count-only coverage')
    return rows,dict(mode=mode,actual_receipts=evidence)


def universal_domains(units,W):
    order12=required(units,UNIVERSAL,'order12-d')['result']
    pairs=[(a,b) for a in (8,9) for b in (10,11)]
    expected32={canonical([list(p),side,branch]) for p in pairs for side in ('left','right')
                for branch in ('isolated4','mate3_1','D5','A1_A3')}
    expected_separate={canonical([list(p),a,b,small]) for p in pairs for a in (8,9)
                       for b in (10,11) for small in (12,13)}
    for field,expected in (('current32',expected32),('separate32',expected_separate)):
        rows=order12[field];keys=[canonical(p['application_key']) for p in rows]
        need(len(keys)==len(set(keys)) and set(keys)==expected,'exact physical '+field+' receiver keys')
        for row in rows:
            receiver_sum(row)
    complete=None;seen=set();old_priority_applications=[]
    for shard in range(4):
        data=large_d_receipt(units,shard)[1]['result']
        keys=data['complete_old40plus16_keys']
        tokens=[canonical(k) for k in keys]
        need(len(tokens)==len(set(tokens)),'unique source-bound OLD40+16 physical keys')
        if complete is None:
            complete=tokens
        need(tokens==complete and data['shard']==shard and data['shards']==4 and
             data['old_split_is_not_generic_current_protected74'] is True,
             'same complete acyclic old split in every shard')
        apps=data['theorem_applications']
        expected={canonical(k) for k in keys if k['high_index']%4==shard}
        got=[canonical(p['application_key']) for p in apps]
        need(len(got)==len(set(got)) and set(got)==expected and not seen.intersection(got),
             'entire disjoint actual OLD40+16 application shard')
        seen.update(got)
        for app in apps:
            need(app['theorem']=='finite57.current.largeD_leaf9_distinct_half.whole_face.v1'
                 and app['generic_current_protected74_used'] is False,
                 'largeD theorem cannot depend on genericprotected74')
            for option in app['protected_options']:
                for mapping in option['complete_small_native_maps']:
                    resolution=old_largeD_priority(mapping,W,units)
                    if resolution is not None:
                        old_priority_applications.append(dict(
                            old_largeD_consumer_key=app['application_key'],
                            actual_old_protected_corner=option['actual_corner'],
                            physical_context_sha256=option['physical_context_sha256'],
                            **resolution))
    need(seen==set(complete),'all exact old-split applications, not merely their cardinality')
    quiet,quiet_evidence=quiet_applications(units)
    keys=[canonical(p['application_key']) for p in quiet]
    need(len(keys)==len(set(keys)),'distinct exact quiet physical witness-copy contexts')
    priority=required(units,QUIET,'quiet-b13')['exact_quiet_consumer_joins']
    providers={canonical(p['consumer_key']):p for p in priority}
    used=set()
    for app in quiet:
        need(app['same_three_witness_copies_at_both_corners'] is True and
             app['intrinsic_zero_credited_as_unit'] is False,
             'same actual rank15 witness copies, no intrinsic zero promotion')
        if app['distinct_current_total'] is None:
            token=canonical(app['application_key'])
            need(token in providers and providers[token]['intrinsic_current_weight']=='0',
                 'exact quietB13 consumer/provider key, not a family-name substitute')
            provider=providers[token]
            need(app['original_source_context_sha256']==provider['physical_context_sha256'] and
                 app['original_source_selector'][1]=='nonextreme_zero_cases' and
                 app['original_source_selector'][2]==provider['provider_key']['quiet_source_index'] and
                 app['paying_receiver']['actual_corner']==provider['original_paying_marks'] and
                 app['paying_receiver']['current']['marked_images']==provider['provider_key']['marked_images'] and
                 app['paying_receiver']['current']['weight']=='0',
                 'same ORIGINAL source context/index/paying marks/complete B13 images')
            used.add(token)
        else:
            need(F(app['distinct_current_total'])>=1,'quiet distinct current receiver sum')
    need(used==set(providers),'every quietB13 consumer joins its actual quiet context')
    return dict(order12_current32_and_separate32='exact physical keys and current sums checked',
                old_largeD_application_keys=complete,
                old_largeD_exact_priority_applications=old_priority_applications,quiet_context_keys=keys,
                quiet_receiver_evidence=quiet_evidence,
                quiet_B13_consumer_keys=sorted(used))


def theorem_evidence(theorem,units):
    if theorem == 'finite57.current.largeD_leaf9_distinct_half.whole_face.v1':
        entries = [large_d_receipt(units,i)[0] for i in range(4)]
        return [dict(worker=e['worker'],batch=e['batch'],receipt=e['path']) for e in entries]
    table = {
      'finite57.current.C10_endpoint_zero.whole_face.v1':[(G4,'order10-b'),(G4,'order11-bc')],
      'finite57.current.BC11_zero.whole_face.v1':[(G4,'order11-bc')],
      'finite57.current.B12_center_terminal.whole_face.v1':[(G4,'order12-center')],
      'finite57.current.BC12_near_terminal.whole_face.v1':[(G4,'order12-near'),(UNIVERSAL,'five-local')],
      'finite57.current.C12_near_companion_zero.whole_face.v1':[(UNIVERSAL,'five-local')],
      'finite57.current.affine_order12_D_opposite_leaves.whole_face.v1':[(UNIVERSAL,'order12-d')],
      'finite57.current.BC13_inner_terminal.whole_face.v1':[(UNIVERSAL,'five-local'),(FINITE,'universal-bc13')],
      'finite57.current.C13_inner_first_zero.whole_face.v1':[(UNIVERSAL,'five-local')],
      'finite57.current.D13_center_zero.whole_face.v1':[(ZERO,'d13-center-zero')],
      'finite57.current.D13_center_leaf.whole_face.v1':[(UNIVERSAL,'five-local'),(FINITE,'universal-d13')],
      'finite57.current.D14_near_fork.whole_face.v1':[(UNIVERSAL,'five-local'),(FINITE,'universal-d14-near')],
      'finite57.current.D15_basic_fractional.whole_face.v1':[(UNIVERSAL,'five-local'),(FINITE,'universal-d15-basic')],
      'finite57.current.largeD_leaf9_distinct_half.whole_face.v1':
          [(UNIVERSAL,f'large-d-{i}') for i in range(4)],
    }
    need(theorem in table, 'unsupported source-specific whole-face theorem')
    for worker,batch in table[theorem]:
        required(units,worker,batch)
    return [dict(worker=w,batch=b,receipt=units[w,b][0]['path']) for w,b in table[theorem]]


def whole_application(key,W,units):
    U.full_key(key)
    need(key['affine'] is True, 'universal whole target is affine')
    role = (key['target_kind'],key['target_order'],*key['marked_image'])
    theorem = ROLE_THEOREMS.get(role)
    need(theorem is not None, 'no supported universal exact target-role application')
    if theorem in {'finite57.current.BC12_near_terminal.whole_face.v1',
            'finite57.current.C12_near_companion_zero.whole_face.v1',
            'finite57.current.BC13_inner_terminal.whole_face.v1',
            'finite57.current.C13_inner_first_zero.whole_face.v1',
            'finite57.current.D13_center_leaf.whole_face.v1',
            'finite57.current.D14_near_fork.whole_face.v1',
            'finite57.current.D15_basic_fractional.whole_face.v1',
            'finite57.current.affine_order12_D_opposite_leaves.whole_face.v1'}:
        batch='order12-d' if role[1]==12 and role[0]=='D' else 'five-local'
        apps=required(units,UNIVERSAL,batch)['result']['theorem_applications']
        matching=[p for p in apps if p['theorem']==theorem and
                  p['application_key']==list(role)]
        need(len(matching)==1 and matching[0]['status']=='BOUND_CONDITIONAL_WHOLE_FACE_APPLICATION',
             'THIS exact whole target role joins its actual universal theorem application')
    graph = W.classical_graph(role[0],role[1],affine=True)
    pair = key['marked_image']
    intrinsic = W.affine_weight(role[0],role[1],*pair)
    return dict(status='BOUND_CONDITIONAL_WHOLE_FACE_APPLICATION',
        theorem=theorem,selected_full_key=key,
        typed_hypotheses=dict(whole_actual_affine_factor=True,affine_kind=role[0],
            affine_order=role[1],affine_rank=role[1]-1,canonical_marks=pair,
            marked_distance=W._distances(graph,pair[0])[pair[1]],
            current_intrinsic_weight=str(intrinsic)),
        actual_theorem_interface_receipts=theorem_evidence(theorem,units),
        assumptions=['actual-star','native-classification','universal-receiver-geometry'],
        conclusion=dict(quantity='total_current_charge_of_same_actual_face',lower='1'),
        intrinsic_corner_unit=False,geometric_premises_proved_here=False)


def prefix_identity(p,core,row,c):
    consumer=dict(high_index=p['index'],method=p['method'],core_key=p['core_key'],
                  selected_marks=p['preferred_corner'])
    need(all(row['consumer_key'].get(k)==v for k,v in consumer.items()),
         'same high/method/core/ordered marks')
    expected=sha256(c.canonical([core['gram'],p['marked_component'],p['preferred_corner']]).encode()).hexdigest()
    need(row['physical_context_sha256']==expected,'same actual prefix physical core')
    return consumer,expected


def adaptive_prefix_identity(p, core, row, c):
    """Bind an ALREADY F.actual-validated adaptive row to its actual corner.

    Returns identity only. Whole-support geometry and independent same-face
    S8 must remain separate accepted inputs in the enclosing application.
    """
    method = 'nine_context_adaptive_finite_payment'
    key = row['consumer_key']
    fields = {'high_index', 'method', 'core_key', 'local_case', 'selected_marks'}
    need(type(key) is dict and set(key) == fields,
         'complete adaptive consumer key, including local source case')
    need(p['method'] == method and type(p['index']) is int and
         0 <= p['index'] < 2614 and type(p['core_key']) is str and
         type(p['local_case']) is int and 0 <= p['local_case'] < 9,
         'actual adaptive high/core/local-case identity')
    gram = core['gram']
    need(type(gram) is list and 1 <= len(gram) <= 64 and
         all(type(r) is list and len(r) == len(gram) for r in gram),
         'bounded complete physical core Gram')
    marks = key['selected_marks']
    need(type(marks) is list and len(marks) == 2 and
         all(type(v) is int and 0 <= v < len(gram) for v in marks) and
         len(set(marks)) == 2 and type(key['high_index']) is int and
         type(key['local_case']) is int,
         'ordered distinct adaptive actual marks and strict integer key')
    consumer = dict(high_index=p['index'], method=method,
                    core_key=p['core_key'], local_case=p['local_case'],
                    selected_marks=marks)
    need(key == consumer and type(p['payment']) is dict and
         p['payment']['corner'] == marks,
         'same adaptive high/method/core/local-case and actual payment marks')
    corners = p['all_corners']
    need(type(corners) is list and len(corners) == 4 and
         all(type(q) is dict for q in corners),
         'four original actual corner records')
    matches = [(i, q) for i, q in enumerate(corners)
               if q.get('preferred_corner') == marks]
    need(len(matches) == 1, 'unique actual adaptive paying-corner selector')
    ordinal, chosen = matches[0]
    nodes = chosen['marked_component']
    need(type(nodes) is list and len(nodes) == len(set(nodes)) and
         all(type(v) is int and 0 <= v < len(gram) for v in nodes) and
         set(marks) <= set(nodes), 'complete selected adaptive component nodes')
    expected = sha256(c.canonical([gram, nodes, marks]).encode()).hexdigest()
    need(row['physical_context_sha256'] == expected,
         'same adaptive actual core/component/ordered-mark digest')
    local = p['local_case']
    branch_method = ('sealed_marked_leaf' if local in (1, 3, 5, 8) else
                     'adaptive_untouched_marked_leaf' if local in (0, 4, 7) else
                     'adaptive_leaf_or_bounded_relevant_tail')
    need(row['method'] == branch_method and
         type(row['current_alternative_branches']) is list and
         row['current_alternative_branches'] and
         row['required_independent_theorem'] == 'finite57.current.S8.whole_face.v1' and
         row['independent_S8_verified'] is False,
         'same accepted adaptive branch; independent S8 remains a prerequisite')
    path = row['retained_path']
    need(type(path) is list and len(path) == 10 and len(set(path)) == 10 and
         all(type(v) is int and v in nodes for v in path) and
         [path[0], path[-1]] == marks,
         'same ordered ten-node retained source path')
    selector = dict(high_row_selector=['isolation', 'cases', p['index']],
                    actual_corner_selector=['all_corners', ordinal],
                    component_nodes=nodes, selected_marks=marks,
                    local_case=local, identity_only=True,
                    whole_face_payment_proved_here=False,
                    intrinsic_corner_unit=False)
    return consumer, expected, selector



def prefix_rows(route_result,cases,cores,r,c,s8,units):
    applications=[]
    protected={}
    for i in range(4):
        doc=required(units,PROTECTED,f'protected-{i}')
        for row in doc['current_prefix_applications']:
            k=row['consumer_key']['high_index']
            need(k not in protected, 'duplicate protected physical high key')
            protected[k]=row
    adaptive={}
    for row in required(units,ADAPTIVE,'adaptive9')['current_prefix_providers']:
        k=row['consumer_key']['high_index']
        need(k not in adaptive, 'duplicate adaptive physical high key')
        adaptive[k]=row
    for request in route_result['prefix_obligations']:
        p=cases[request['index']];core=cores[p['core_key']]
        method=request['method']
        if method=='independent_S8':
            application=s8.apply_physical_s8(p,core,c.E,c.FUEL)
        elif method=='display_exceeds_proved_whole_bound':
            consumer=dict(high_index=p['index'],method=p['method'],core_key=p['core_key'],
                          selected_marks=p['preferred_corner'])
            h=p['component_order_bound'];nodes=p['marked_component']
            need(type(h) is int and 0<=h<len(nodes) and p['method'] in
                 ('radical','lorentz_residual','relevant_component',
                  'marked_component_isolated_by_end_inverses'),
                 'exact displayed component violates its source-bound extension certificate')
            application=dict(status='BOUND_CONDITIONAL_COMPONENT_IMPOSSIBILITY',
                consumer_key=consumer,displayed_nodes=nodes,whole_order_bound=h,
                exact_bound_certificate=p,
                assumptions=['generic-extension-bounds','actual-star'],
                conclusion='The displayed nodes cannot lie in the claimed actual whole component.',
                payment_claimed=False)
        elif method in ('finite_leaf9_or_ideal_feature5','nine_context_adaptive_finite_payment'):
            rows=protected if method=='finite_leaf9_or_ideal_feature5' else adaptive
            need(p['index'] in rows, 'missing actual physical prefix provider')
            row=rows[p['index']]
            actual_corner_identity=None
            if method=='nine_context_adaptive_finite_payment':
                consumer,expected,actual_corner_identity=adaptive_prefix_identity(p,core,row,c)
            else:
                consumer,expected=prefix_identity(p,core,row,c)
            priority=[]
            if method=='finite_leaf9_or_ideal_feature5':
                for item in row['exact_small_current_priority_requests']:
                    app=whole_application(item['key'],c.W,units)
                    need(app['theorem']==item['required_current_whole_face_provider'],
                         'exact protected whole-face priority role')
                    priority.append(app)
                theorem_evidence('finite57.current.largeD_leaf9_distinct_half.whole_face.v1',units)
                required(units,UNIVERSAL,'order12-d')
            application=dict(status='BOUND_CONDITIONAL_PHYSICAL_PREFIX_APPLICATION',
                consumer_key=consumer,physical_context_sha256=expected,
                actual_corner_identity=actual_corner_identity,
                actual_prefix_provider=row,priority_applications=priority,
                independent_S8_runtime=runtime_spine(units),
                assumptions=['actual-star','native-classification','generic-extension-bounds',
                             'universal-receiver-geometry'],
                conclusion=dict(quantity='total_current_charge_of_same_actual_face',lower='1'),
                intrinsic_corner_unit=False)
        else:
            raise ValueError('unsupported original prefix route')
        applications.append(dict(original_prefix=request,application=application))
    return applications


def reused_route_result(units,batch,collector_sha,partition):
    """Consume one actual index; never reconstruct a missing enumeration."""
    need(type(batch) is str and batch in {f'selected-{i:02d}' for i in range(64)},
         'exact named selected route batch')
    original=required(units,ROUTE,batch)
    need(original['batch']==batch and original['source_sha256']==ROUTE_SHA and
         original['collector_sha256']==collector_sha,
         'preserved actual route-index source/batch/collector identity')
    result=original['result']
    need(type(result['shard']) is int and result['shard']==int(batch[-2:]) and
         type(result['shards']) is int and result['shards']==64 and
         result['complete_authoritative_partition']==U.PARTITION==partition,
         'exact authoritative reused shard/partition')
    return result,dict(kind='REUSED_ACTUAL_INDEX_WITHOUT_REENUMERATION',
                       receipt=units[ROUTE,batch][0])


def selected(batch,meta,reviews,accepted,snapshot):
    units=inventory(meta,reviews)
    spine=runtime_spine(units);proof_inputs(accepted)
    curated = provider_curation(units,accepted)
    s8=load('finite57_final_s8',S8,reviews[S8]['source_sha256'])
    r=load('finite57_final_routes',ROUTE,ROUTE_SHA);c=r.prepare()
    cases,cores,partition=r.domain(c)
    legacy_application,full_d_application=source_application_helpers(reviews,cases,cores,c)
    old_affine_application=old_affine_application_helper(reviews,legacy_application)
    shard=int(batch[-2:])
    # All64 named indexes already ran successfully. Reuse their complete
    # actual outputs; a missing receipt is not permission to repeat a search.
    result,origin=reused_route_result(units,batch,r.COLLECTOR_SHA,partition)
    obligations=result['unbound_current_completion_obligations'];wanted={}
    for ordinal,row in enumerate(obligations):
        token=U.full_key(row['key'])
        need(token not in wanted and row['key']['high_index']%64==shard and
             row['binding_status']=='UNBOUND' and 0<=U.scalar(row['intrinsic_current_lower'])<1,
             'exact unique conservative nonunit full-map obligation')
        wanted[token]=dict(ordinal=ordinal,row=row,providers=[])
    for (worker,provider_batch),(entry,inner) in units.items():
        adapter=provider_adapter(worker)
        if adapter is None or worker==ROUTE:
            continue
        need(adapter!='disabled_whole_schema','old malformed whole0a remains disabled')
        for provider_ordinal,(key,nature) in enumerate(U.providers(inner,adapter,obligations)):
            token=U.full_key(key)
            if token not in wanted:
                continue
            occurrence_bindings=curated_occurrences(curated,entry,key)
            application=None
            if adapter in ('proper','whole_schema_safe','newer','affine','old_affine','guard'):
                application=s8.independent_s8_application(key,c.W)
                premise=application['provider_composition_premise']
                need(premise['status']=='REQUIRES_EXACT_MATCHED_ACCEPTED_PROVIDER_IMPLICATION' and
                     premise['same_full_key']==key,
                     'exact accepted provider implication, never a requiresS8 flag alone')
                premise['matched_provider']=dict(worker=worker,batch=provider_batch,
                    receipt=entry['path'],receipt_sha256=entry['sha256'],
                    exact_extraction_ordinal=provider_ordinal,source_review=reviews[worker],
                    actual_curated_source_occurrences=occurrence_bindings,
                    accepted_assumption_id='selected-provider-geometry')
            wanted[token]['providers'].append(dict(receipt=entry['path'],adapter=adapter,nature=nature,
                exact_extraction_ordinal=provider_ordinal,source_review=reviews[worker],
                actual_curated_source_occurrences=occurrence_bindings,
                conditional_S8_escape=application,
                proof_logic=None if application is None else dict(
                    accepted_exhaustiveness_premise='selected-provider-geometry',
                    split='The exact accepted provider proof excludes this full map or pays this SAME face in the no-S8 branch; the only alternative is a same-face short-path escape.',
                    unpaid_face_argument='Assume total current charge<1. Conditional universal S8 then rules out every adjacent white-interior path<=8 in this face; the exhaustive no-S8 provider branch excludes the map or contradicts its unpaid charge.',
                    actual_short_path_asserted=False)))
    joined=[]
    for token,item in wanted.items():
        if not item['providers']:
            key=item['row']['key']
            graph=c.W.classical_graph(key['target_kind'],key['target_order'],affine=key['affine'])
            distance=c.W._distances(graph,key['marked_image'][0])[key['marked_image'][1]]
            if distance<=8:
                # The accepted independent theorem is about actual paths,
                # including shorter affine-cycle arcs, not only finite deep
                # zeros. This is a whole-face payment, never a corner unit.
                item['providers'].append(dict(adapter='direct_actual_short_path_S8',
                    application=s8.apply_s8(key,c.W)))
        if not item['providers']:
            # Exact universal WHOLE target implication is distinct from an
            # intrinsic unit or matching a family name in a historical list.
            key=item['row']['key']
            role=(key['target_kind'],key['target_order'],*key['marked_image'])
            if key['affine'] and role in ROLE_THEOREMS:
                item['providers'].append(dict(adapter='typed_universal_whole_face',
                    application=whole_application(key,c.W,units)))
        if not item['providers']:
            key=item['row']['key']
            application=full_d_application.apply_selected(key,c.W,units,
                dict(same_full_key=key,whole_actual_affine_factor=True,actual_runtime_S8=spine))
            if application is not None:
                item['providers'].append(dict(adapter='full_universal_affine_D_source_application',
                    application=application))
        if not item['providers']:
            application=legacy_application.apply_selected(item['row']['key'],c.W,units)
            if application is not None:
                item['providers'].append(dict(adapter='legacy_exact_source_application',
                    application=application))
        if not item['providers']:
            application=old_affine_application.apply_selected(item['row']['key'],c.W,units)
            if application is not None:
                item['providers'].append(dict(adapter='old_affine_exact_source_application',
                    application=application))
        need(item['providers'],('uncovered exact full map',item['row']['key']))
        joined.append(dict(obligation_ordinal=item['ordinal'],key=item['row']['key'],
            intrinsic_current_lower=item['row']['intrinsic_current_lower'],
            conditional_negative_D13_guard=item['row']['conditional_negative_D13_guard'],
            providers=item['providers']))
    paid=result['direct_whole_factor_unit_maps']
    need(type(paid)is int and paid>=0 and result['all_induced_maps_visited']==paid+len(joined),
         'every enumerated map direct-unit or explicit joined nonunit')
    U.digest_string(result['direct_map_key_sha256'])
    prefixes=prefix_rows(result,cases,cores,r,c,s8,units)
    for path in c.READ:
        hash_file(path,c.PINS[path])
    return dict(status='CURRENT_RULE_SELECTED_COMPLETE_AT_EXPLICIT_PREMISES',
        batch=batch,source_sha256=sha256(source(OWN)).hexdigest(),input_snapshot=snapshot,
        result=dict(shard=shard,shards=64,complete_authoritative_partition=partition,
            authoritative_route_result=result,route_origin=origin,exact_nonunit_applications=joined,
            exact_prefix_applications=prefixes,independent_S8=spine,
            source_application_contracts=full_d_application.application_contracts(),
            selected_routes_complete_at_explicit_premises=True,
            full_transitive_regeneration=False,geometric_formalization=False,
            global_finite57_verified=False),
        accepted_assumptions=list(accepted.values()),
        consumed_inputs={p:c.PINS[p] for p in sorted(c.READ)})


def original_domain(units):
    rows=[];domain=None
    for shard in range(4):
        doc=required(units,DIRECT,f'original-{shard}')['result']
        keys=doc['complete_domain_keys']
        need(type(keys)is list and len(keys)==len({canonical(k) for k in keys})==216,'exact complete original domain keys')
        if domain is None:
            domain=keys
        need(keys==domain,'same complete original physical source/case domain')
        expected=domain[shard::4]
        current=doc['rows']
        need(len(current)==len(expected),'complete original named shard')
        actual_keys=[]
        for row in current:
            U.funded(row,full=True)
            actual_keys.append([row['source'],canonical(row['index'])])
            need(row['final_collector_joins'],'exact original final-corner joins retained')
        need(actual_keys==expected,'exact source/case original physical receiver shard')
        rows.extend(current)
    need(len(rows)==216,'all exact original receivers occur once')
    return dict(complete_domain_keys=domain,physical_receivers=216,
                receipts=[units[DIRECT,f'original-{i}'][0]['path'] for i in range(4)])


def final(meta,reviews,accepted,snapshot):
    units=inventory(meta,reviews);proof_inputs(accepted);spine=runtime_spine(units)
    provider_curation(units,accepted)
    original=original_domain(units)
    for batch in G4_BATCHES:
        required(units,G4,batch)
    for batch in UNIVERSAL_BATCHES:
        if batch.startswith('large-d-'):
            large_d_receipt(units,int(batch[-1]))
        else:
            required(units,UNIVERSAL,batch)
    hash_file(LITERAL,LITERAL_SHA)
    literal=load('finite57_final_old_priority_literal',LITERAL,LITERAL_SHA)
    universal=universal_domains(units,literal)
    for i in range(4):
        required(units,PROTECTED,f'protected-{i}')
    required(units,ADAPTIVE,'adaptive9')
    required(units,QUIET,'quiet-b13');required(units,ONE_SIDE,'one-sided-eight')
    for batch in ('d12-first','d12-rest'):
        required(units,D12,batch)
    for batch in ('universal-bc13','universal-d13','universal-d14-near','universal-d15-basic'):
        required(units,FINITE,batch)
    for worker in CONTROL_WORKERS:
        required(units,worker,None)
    selected=[];all_prefixes=set();all_maps=set();application_contracts={}
    for shard in range(64):
        doc=required(units,OWN,f'selected-{shard:02d}')
        premise_ids=retain_selected_assumptions(accepted,doc['accepted_assumptions'])
        result=doc['result'];route=result['authoritative_route_result']
        need(result['shard']==shard and result['shards']==64 and
             result['selected_routes_complete_at_explicit_premises'] is True and
             result['complete_authoritative_partition']==U.PARTITION and
             result['full_transitive_regeneration'] is False and
             result['geometric_formalization'] is False and result['global_finite57_verified'] is False,
             'same exact conditional selected scope; no historical false flag discarded')
        joined=result['exact_nonunit_applications'];prefix=result['exact_prefix_applications']
        contracts=result['source_application_contracts']
        need(type(contracts) is dict,'expanded selected source-application contracts required')
        for cid,contract in contracts.items():
            need(type(cid) is str and type(contract) is dict and
                 contract['contract_id']==cid and contract['assumptions'] and
                 contract['geometric_premises_proved_by_code'] is False and
                 contract['full_transitive_regeneration'] is False and
                 contract['global_finite57_verified'] is False,
                 'complete explicit conditional source-application contract')
            if cid in application_contracts:
                need(canonical(application_contracts[cid])==canonical(contract),
                     'same shared contract means the same expanded source and premise tables')
            application_contracts[cid]=contract
        need([j['key'] for j in joined]==[p['key'] for p in route['unbound_current_completion_obligations']]
             and [p['original_prefix'] for p in prefix]==route['prefix_obligations'],
             'every original full map/prefix retains its exact actual application')
        for row in joined:
            token=U.full_key(row['key'])
            need(token not in all_maps and row['key']['high_index']%64==shard and row['providers'],
                 'disjoint complete actual nonunit applications')
            all_maps.add(token)
            for provider in row['providers']:
                if provider['adapter']=='full_universal_affine_D_source_application':
                    app=provider['application']
                    need(app['selected_full_key']==row['key'] and app['contract_id'] in contracts,
                         'every full-D application retains its exact key and expanded contract')
                elif provider['adapter'] in ('legacy_exact_source_application',
                                              'old_affine_exact_source_application'):
                    app=provider['application']
                    need(app['application_key']['selected_full_key']==row['key'] and
                         app['supplemental_source_premises'] and
                         app['supplemental_source_premise_acceptance'],
                         'every legacy application retains its exact key and supplemental premises')
        for row in prefix:
            key=row['original_prefix']['index']
            need(key not in all_prefixes and key%64==shard,'disjoint actual prefix applications')
            all_prefixes.add(key)
        selected.append(dict(batch=doc['batch'],receipt=units[OWN,doc['batch']][0]['path'],
            retained_accepted_assumption_ids=premise_ids,
            direct_unit_full_key_sha256=route['direct_map_key_sha256'],
            direct_unit_maps=route['direct_whole_factor_unit_maps'],
            joined_nonunit_maps=len(joined),physical_prefix_applications=len(prefix)))
    # Inherited numerical/profile/tail propositions remain explicit accepted
    # premises. This conjunction never rewrites them as fresh executions.
    return dict(status='CURRENT_RULE_RUNTIME_CONJUNCTION_COMPLETE_AT_EXPLICIT_PREMISES',
        batch='final',source_sha256=sha256(source(OWN)).hexdigest(),input_snapshot=snapshot,
        actual_original_domain=original,actual_selected_shards=selected,independent_S8=spine,
        actual_universal_domain_joins=universal,
        source_application_contracts=application_contracts,
        accepted_assumptions=list(accepted.values()),actual_unit_receipts=[e for e,_ in units.values()],
        current_rule_application_runtime_conjunction_complete=True,
        full_transitive_regeneration=False,geometric_formalization=False,
        global_finite57_verified=False,new_bound_or_construction=False,
        scope='Current-rule application/runtime conjunction ONLY at the structured accepted geometry, complete finite-data, inherited tail and incidence premises. Every required actual source is separately bound; no full scientific regeneration or geometric formalization.')


def main():
    global U
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',choices=('final',)+tuple(f'selected-{i:02d}' for i in range(64)),required=True)
    parser.add_argument('--parent-lock',action='store_true',required=True)
    args=parser.parse_args()
    U=load('finite57_final_receipt_parser',UNION,UNION_SHA)
    U.parent_lock()
    hash_file(UNION,UNION_SHA)
    before=source(OWN)
    codec()  # Bind the exact lossless representation before final source rechecks.
    meta,reviews,accepted,snapshot=manifest(args.batch)
    payload=final(meta,reviews,accepted,snapshot) if args.batch=='final' else selected(args.batch,meta,reviews,accepted,snapshot)
    need(source(OWN)==before,'final conjunction source changed')
    for name,(identity,digest) in FILES.items():
        st=U.relative_path(name).stat()
        need(identity==(st.st_dev,st.st_ino,st.st_size,st.st_mtime_ns,st.st_ctime_ns),
             'consumed source/input changed before publication')
    for name,digest in list(U.READ.items()):
        need(sha256(U.raw(name,2*1024**2,remember=False)).hexdigest()==digest,
             'actual manifest/receipt/log changed before publication')
    payload['bound_source_sha256']={p:d for p,(_,d) in sorted(FILES.items())}
    encoded=json.dumps(payload,sort_keys=True)
    if len(encoded.encode())+1 > MAX_OUTPUT:
        packet = codec().encode_result(payload)
        need(expanded_loads(codec().decode_result(packet)) == payload,
             'complete conjunction output byte/object round trip')
        encoded=json.dumps(dict(status=payload['status'],batch=payload['batch'],
            source_sha256=payload['source_sha256'],complete_payload_encoding=packet),sort_keys=True)
    need(len(encoded.encode())+1<=MAX_OUTPUT,'complete current-rule output cap; no omitted application keys')
    print(encoded)


if __name__=='__main__':
    main()

