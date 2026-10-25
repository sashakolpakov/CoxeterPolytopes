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
OWN = 'scripts/check_finite57_current_rule_conjunction_20261025.py'
SELECTED_INPUTS = 'results/finite57_current_rule_selected_inputs_20261025.json'
FINAL_INPUTS = 'results/finite57_current_rule_final_inputs_20261025.json'
MAX_OUTPUT = 500 * 1024
MAX_UNITS = 256
MAX_HASH_BYTES = 256 * 1024**2
HASH_BYTES = 0
FILES = {}
U = None
# These are EXACT actual worker roles, not labels that establish theorems.
G4 = 'scripts/check_finite57_g4_lower_cycle_terminal_bindings_20261025.py'
FINITE = 'scripts/check_finite57_terminal_collectors_syntax_safe_20261025.py'
DIRECT = 'scripts/check_finite57_terminal_collectors_direct_receivers_20261025.py'
PROTECTED = 'scripts/index_finite57_current_protected74_order_safe_applications_20261025.py'
ADAPTIVE = 'scripts/check_finite57_current_adaptive9_path_binding_20261025.py'
ZERO = 'scripts/check_finite57_current_d13_center_zero_outputs_20261025.py'
QUIET = 'scripts/check_finite57_quiet_b13_current_priority_bridge_20261025.py'
ONE_SIDE = 'scripts/check_finite57_final_guard_one_sided_current_receivers_20261025.py'
D12 = 'scripts/check_finite57_current_d12_terminal_bindings_guard_safe_20261025.py'
DEEP_BATCHES = tuple(f'deep-{i:02d}' for i in range(16))
G4_BATCHES = ('order10-cycle','order10-b','order10-d','cycle11','order11-bcd',
    'order11-bc','order12-center','order12-near','order12-selected',
    'order13-sparse','cycle12','cycle13','cycle14','cycle16','cycle18-tree','large-d')
UNIVERSAL_BATCHES = ('order12-d','five-local','large-d-0','large-d-1',
                     'large-d-2','large-d-3','final-quiet')
CONTROL_WORKERS = {
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
         'complete unchanged ten-controller closure')
    for path, digest in meta['controller_source_sha256'].items():
        hash_file(path,digest)
    need(meta['controller_source_sha256'][U.RUNNER] == U.RUNNER_SHA,
         'unchanged bounded runner')
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
         outer['command'] == expected_args and outer['checker_sha256'] == worker_sha and
         outer['runner_sha256'] == U.RUNNER_SHA, 'actual exact bounded command/source PASS')
    need(type(outer['wall_limit_seconds']) is int and outer['wall_limit_seconds'] == 59 and
         type(outer['rss_limit_bytes']) is int and outer['rss_limit_bytes'] == 1024**3 and
         type(outer['threads']) is int and outer['threads'] == 1 and outer['lock_owner'] == 'wrapper' and
         type(outer['elapsed_seconds']) in (int,Decimal) and 0 <= outer['elapsed_seconds'] < 59 and
         type(outer['peak_process_tree_rss_bytes']) is int and 0 <= outer['peak_process_tree_rss_bytes'] < 1024**3,
         'unchanged measured59s/1GiB/one-thread wrapper lock')
    need(outer['initial_resources']['free_disk_gib'] >= 40 and
         outer['initial_resources']['available_memory_gib'] >= 6, 'actual40/6 admission')
    output = outer['stdout']
    need(type(output) is str and output.endswith('\n') and len(output.splitlines()) == 1 and
         len(output.encode()) == outer['stdout_bytes'] < 512*1024 and
         len(lines) == 3 and lines[0]+'\n' == output, 'complete public three-line stdout')
    need(U.loads(lines[1]) == {k:v for k,v in outer.items() if k != 'stdout'}, 'outer/log equality')
    bridge = U.loads(lines[2])
    need(bridge['owned_tree_bridge_only'] is True and type(bridge['child_exit_code']) is int and bridge['child_exit_code'] == 0 and
         bridge['cleanup']['signals'] == [] and bridge['cleanup']['verified_live_owned_processes_remaining'] == 0 and
         bridge['controller_source_sha256'] == controllers and bridge['assertions_enabled'] is True and
         type(bridge['python_optimize_flag']) is int and bridge['python_optimize_flag'] == 0 and bridge['child_pythonoptimize'] is None,
         'actual owned process-tree/assertion/source closure')
    inner = U.loads(output)
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
 OWN: ('CURRENT_RULE_SELECTED_COMPLETE_AT_EXPLICIT_PREMISES',),
 'scripts/test_finite57_current_rule_conjunction_controls_20261025.py': ('PASS_CURRENT_RULE_CONJUNCTION_CONTROLS_ONLY',),
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


def universal_domains(units):
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
    complete=None;seen=set()
    for shard in range(4):
        data=required(units,UNIVERSAL,f'large-d-{shard}')['result']
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
                    U.full_key(mapping['key'])
    need(seen==set(complete),'all exact old-split applications, not merely their cardinality')
    quiet=required(units,UNIVERSAL,'final-quiet')['result']['theorem_applications']
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
                old_largeD_application_keys=complete,quiet_context_keys=keys,
                quiet_B13_consumer_keys=sorted(used))


def theorem_evidence(theorem,units):
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
        consumer=dict(high_index=p['index'],method=p['method'],core_key=p['core_key'],
                      selected_marks=p['preferred_corner'])
        method=request['method']
        if method=='independent_S8':
            application=s8.apply_physical_s8(p,core,c.E,c.FUEL)
        elif method=='display_exceeds_proved_whole_bound':
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


def selected(batch,meta,reviews,accepted,snapshot):
    units=inventory(meta,reviews)
    spine=runtime_spine(units);proof_inputs(accepted)
    s8=load('finite57_final_s8',S8,reviews[S8]['source_sha256'])
    r=load('finite57_final_routes',ROUTE,ROUTE_SHA);c=r.prepare()
    cases,cores,partition=r.domain(c)
    shard=int(batch[-2:])
    if shard==0:
        original=required(units,ROUTE,'selected-00')
        result=original['result']
        need(original['collector_sha256']==r.COLLECTOR_SHA, 'preserved selected00 collector identity')
        origin=dict(kind='REUSED_ACTUAL_INDEX_WITHOUT_REENUMERATION',receipt=units[ROUTE,'selected-00'][0])
    else:
        result=r.selected(shard,c)
        origin=dict(kind='AUTHORITATIVE_SHARD_ENUMERATED_ONCE_INSIDE_THIS_APPLICATION_JOB')
    need(result['shard']==shard and result['shards']==64 and
         result['complete_authoritative_partition']==U.PARTITION==partition,
         'exact authoritative shard/partition')
    obligations=result['unbound_current_completion_obligations'];wanted={}
    for ordinal,row in enumerate(obligations):
        token=U.full_key(row['key'])
        need(token not in wanted and row['key']['high_index']%64==shard and
             row['binding_status']=='UNBOUND' and 0<=U.scalar(row['intrinsic_current_lower'])<1,
             'exact unique conservative nonunit full-map obligation')
        wanted[token]=dict(ordinal=ordinal,row=row,providers=[])
    for (worker,provider_batch),(entry,inner) in units.items():
        if worker not in U.WORKERS or worker==ROUTE:
            continue
        adapter=U.WORKERS[worker][3]
        need(adapter!='disabled_whole_schema','old malformed whole0a remains disabled')
        for provider_ordinal,(key,nature) in enumerate(U.providers(inner,adapter,obligations)):
            token=U.full_key(key)
            if token not in wanted:
                continue
            application=None
            if adapter in ('proper','whole_schema_safe','newer','affine','old_affine'):
                application=s8.independent_s8_application(key,c.W)
                premise=application['provider_composition_premise']
                need(premise['status']=='REQUIRES_EXACT_MATCHED_ACCEPTED_PROVIDER_IMPLICATION' and
                     premise['same_full_key']==key,
                     'exact accepted provider implication, never a requiresS8 flag alone')
                premise['matched_provider']=dict(worker=worker,batch=provider_batch,
                    receipt=entry['path'],receipt_sha256=entry['sha256'],
                    exact_extraction_ordinal=provider_ordinal,source_review=reviews[worker],
                    accepted_assumption_id='selected-provider-geometry')
            wanted[token]['providers'].append(dict(receipt=entry['path'],adapter=adapter,nature=nature,
                exact_extraction_ordinal=provider_ordinal,source_review=reviews[worker],
                conditional_S8_escape=application,
                proof_logic=None if application is None else dict(
                    accepted_exhaustiveness_premise='selected-provider-geometry',
                    split='The exact accepted provider proof excludes this full map or pays this SAME face in the no-S8 branch; the only alternative is a same-face short-path escape.',
                    unpaid_face_argument='Assume total current charge<1. Conditional universal S8 then rules out every adjacent white-interior path<=8 in this face; the exhaustive no-S8 provider branch excludes the map or contradicts its unpaid charge.',
                    actual_short_path_asserted=False)))
    joined=[]
    for token,item in wanted.items():
        if not item['providers']:
            # Exact universal WHOLE target implication is distinct from an
            # intrinsic unit or matching a family name in a historical list.
            item['providers'].append(dict(adapter='typed_universal_whole_face',
                application=whole_application(item['row']['key'],c.W,units)))
        need(item['providers'],'uncovered exact full map')
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
    original=original_domain(units)
    for batch in G4_BATCHES:
        required(units,G4,batch)
    for batch in UNIVERSAL_BATCHES:
        required(units,UNIVERSAL,batch)
    universal=universal_domains(units)
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
    selected=[];all_prefixes=set();all_maps=set()
    for shard in range(64):
        doc=required(units,OWN,f'selected-{shard:02d}')
        result=doc['result'];route=result['authoritative_route_result']
        need(result['shard']==shard and result['shards']==64 and
             result['selected_routes_complete_at_explicit_premises'] is True and
             result['complete_authoritative_partition']==U.PARTITION and
             result['full_transitive_regeneration'] is False and
             result['geometric_formalization'] is False and result['global_finite57_verified'] is False,
             'same exact conditional selected scope; no historical false flag discarded')
        joined=result['exact_nonunit_applications'];prefix=result['exact_prefix_applications']
        need([j['key'] for j in joined]==[p['key'] for p in route['unbound_current_completion_obligations']]
             and [p['original_prefix'] for p in prefix]==route['prefix_obligations'],
             'every original full map/prefix retains its exact actual application')
        for row in joined:
            token=U.full_key(row['key'])
            need(token not in all_maps and row['key']['high_index']%64==shard and row['providers'],
                 'disjoint complete actual nonunit applications')
            all_maps.add(token)
        for row in prefix:
            key=row['original_prefix']['index']
            need(key not in all_prefixes and key%64==shard,'disjoint actual prefix applications')
            all_prefixes.add(key)
        selected.append(dict(batch=doc['batch'],receipt=units[OWN,doc['batch']][0]['path'],
            direct_unit_full_key_sha256=route['direct_map_key_sha256'],
            direct_unit_maps=route['direct_whole_factor_unit_maps'],
            joined_nonunit_maps=len(joined),physical_prefix_applications=len(prefix)))
    # Inherited numerical/profile/tail propositions remain explicit accepted
    # premises. This conjunction never rewrites them as fresh executions.
    return dict(status='CURRENT_RULE_RUNTIME_CONJUNCTION_COMPLETE_AT_EXPLICIT_PREMISES',
        batch='final',source_sha256=sha256(source(OWN)).hexdigest(),input_snapshot=snapshot,
        actual_original_domain=original,actual_selected_shards=selected,independent_S8=spine,
        actual_universal_domain_joins=universal,
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
    need(len(encoded.encode())+1<=MAX_OUTPUT,'complete current-rule output cap; no omitted application keys')
    print(encoded)


if __name__=='__main__':
    main()

