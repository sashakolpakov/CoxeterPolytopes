#!/usr/bin/env python3
"""SOURCE ONLY: twenty-three synthetic receipt/parser/whole-schema control regressions.

Memory fixtures are deliberately NOT actual executions. They exercise the
pinned union reader's parser and logic, never produce a scientific receipt,
and load no mathematical worker. Filesystem/monitor operation is not tested.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
import fcntl
from hashlib import sha256
import io
import json
from pathlib import Path
import types
import unittest

ROOT = Path(__file__).resolve().parents[1]
UNION = 'scripts/index_finite57_current_terminal_receipt_union_schema_safe_20261025.py'
UNION_SHA = '435e2912d2e1bc82e97da8d33f499dc018802dbd10b2a00368a57e4f14d9fb71'
PROPER = 'scripts/index_finite57_all_proper_finite_completion_providers_20261025.py'
WHOLE = 'scripts/index_finite57_funded_selected_whole_providers_schema_safe_20261025.py'
SOURCE_SHA = {
    WHOLE: '6ace38a526ccf22d9913ffe0b814b7ee4177d74c3be478295a981fef5b7afb03',
    UNION: UNION_SHA,
    'scripts/index_finite57_current_terminal_routes_syntax_safe_20261025.py':
        '43a0c4e986cbdb92d112771bf0edc1d128ce2d0feced30fffd9d9825dda7c29d',
    PROPER: '0ef6647b6b9b702f076b264b730caef89a31e40ff2a502bf5a5e03ab1fe29ba3',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
}
PRIMARY_SHA = 'c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0'
U = None
SOURCES = {}


def need(ok, message):
    if not ok:
        raise ValueError(message)


def raw(relative):
    with (ROOT / relative).open('rb') as stream:
        data = stream.read(128 * 1024 + 1)
    need(0 < len(data) <= 128 * 1024, 'bounded union-test source read')
    return data


def parent_lock():
    with (ROOT / 'formal/.proof-batches/exclusive.lock').open('r+') as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return
        fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
    raise ValueError('shared parent lock must already be held')


def key(reverse=False):
    images = list(range(10))
    if reverse:
        images.reverse()
    return dict(high_index=64, core_key='SYNTHETIC_ASSOCIATION_ONLY', base_nodes=list(range(10)),
                selected_marks=[0, 9], target_kind='A', target_order=12, affine=False,
                induced_map=images, marked_image=[0, 9])


def selected(keys=None, prefix=None):
    keys = [key()] if keys is None else keys
    return dict(shard=0, shards=64, complete_authoritative_partition=dict(U.PARTITION),
        all_induced_maps_visited=len(keys), direct_whole_factor_unit_maps=0,
        direct_map_key_sha256=sha256(b'').hexdigest(), prefix_obligations=prefix or [],
        unbound_current_completion_obligations=[dict(key=k, intrinsic_current_lower='0',
            conditional_negative_D13_guard=False, required_link='synthetic provider', binding_status='UNBOUND') for k in keys])


def proper(k):
    return dict(packets=[dict(proper_map_providers=[dict(maps=[dict(key=k,
        outcome='projection_exceeds_two', computed_record_sha256='a' * 64,
        forced_entries=[[10, 9, 9, 3]],
        geometric_premise='SYNTHETIC fixture: source exclusion only modulo independent current S8.')])])])


class Memory:
    """Patch only reader I/O; ledger contents remain subject to real parsers."""
    def __init__(self):
        self.mem = dict(SOURCES)
        self.controllers = {}
        for name in U.CONTROLLER_PATHS:
            data = SOURCES[name] if name in SOURCES else ('# SYNTHETIC controller ' + name + '\n').encode()
            self.mem[name] = data
            self.controllers[name] = sha256(data).hexdigest()
        self.entries = []; self.reviews = []

    def add(self, worker, batch, result, outer_change=None, line_change=None):
        number = len(self.entries)
        source_hash, status, batches, adapter = U.WORKERS[worker]
        inner = dict(status=status, source_sha256=source_hash, batch=batch, result=result,
                     consumed_inputs={U.PRIMARY: PRIMARY_SHA}, current_terminal_bridge_complete=False)
        if worker == U.ROUTES:
            inner['collector_sha256'] = '930ffff5601e6096359d3b2f8846fa6d7d050f95fb2aabe449447e0d3d0bf9ca'
        output = json.dumps(inner, sort_keys=True) + '\n'
        outer = dict(status='PASS_BOUNDED_RESEARCH_CHECK', command=[worker, '--batch', batch, '--parent-lock'],
            elapsed_seconds=.01, peak_process_tree_rss_bytes=1,
            initial_resources=dict(free_disk_gib=40.5, available_memory_gib=6.5),
            wall_limit_seconds=59, rss_limit_bytes=1024**3, threads=1, stdout_bytes=len(output.encode()),
            stdout=output, runner_sha256=U.RUNNER_SHA, checker_sha256=source_hash, lock_owner='wrapper')
        if outer_change:
            outer_change(outer)
        bridge = dict(owned_tree_bridge_only=True, child_exit_code=0,
            cleanup=dict(signals=[], verified_live_owned_processes_remaining=0),
            controller_source_sha256=dict(self.controllers), assertions_enabled=True,
            python_optimize_flag=0, child_pythonoptimize=None)
        lines = [output.rstrip('\n'), json.dumps({k: v for k, v in outer.items() if k != 'stdout'}), json.dumps(bridge)]
        if line_change:
            line_change(lines)
        ledger, log = json.dumps(outer).encode(), ('\n'.join(lines) + '\n').encode()
        path, stdout_path = f'results/SYNTHETIC_union_{number}.json', f'results/SYNTHETIC_union_{number}.stdout'
        self.mem[path], self.mem[stdout_path] = ledger, log
        entry = dict(worker=worker, batch=batch, path=path, sha256=sha256(ledger).hexdigest(), bytes=len(ledger),
                     stdout_path=stdout_path, stdout_sha256=sha256(log).hexdigest(), stdout_bytes=len(log))
        self.entries.append(entry)
        if not any(p['worker'] == worker for p in self.reviews):
            review_path = f'results/SYNTHETIC_union_review_{number}.md'
            review = b'Synthetic external-review premise only. This is not an acceptance record.\n'
            self.mem[review_path] = review
            self.reviews.append(dict(worker=worker, source_sha256=source_hash,
                path=review_path, sha256=sha256(review).hexdigest()))
        return entry

    def __enter__(self):
        self.saved = U.raw, U.relative_path, U.input_pins, U.READ, U.READ_BYTES, U.VISITS
        U.raw = self.read
        U.relative_path = self.path
        U.input_pins = lambda: {U.PRIMARY: PRIMARY_SHA}
        U.READ, U.READ_BYTES, U.VISITS = {}, 0, 0
        return self

    def __exit__(self, *exc):
        U.raw, U.relative_path, U.input_pins, U.READ, U.READ_BYTES, U.VISITS = self.saved

    def path(self, name, prefix=None):
        if name not in self.mem:
            raise FileNotFoundError(name)
        need(prefix is None or name.startswith(prefix + '/'), 'synthetic path prefix')
        return ROOT / name

    def read(self, name, cap=2 * 1024**2, remember=True):
        self.path(name)
        data = self.mem[name]
        need(0 < len(data) <= cap, 'synthetic input cap')
        if remember:
            U.READ[name] = sha256(data).hexdigest()
        return data

    def run(self):
        self.mem[U.INPUTS] = json.dumps(dict(schema='finite57-current-terminal-actual-receipts-v1',
            controller_source_sha256=self.controllers, actual_receipts=self.entries,
            source_review_bindings=self.reviews)).encode()
        with self:
            return U.run('selected-00')


def actual_fixture(outer_change=None, line_change=None):
    memory = Memory()
    entry = memory.add(U.ROUTES, 'selected-00', selected(), outer_change, line_change)
    return memory, entry


def receiver(values):
    pairs = [[0, 1], [1, 2], [2, 3], [0, 3]]
    from fractions import Fraction
    return dict(corners=[dict(marks=m, weight=v) for m, v in zip(pairs, values, strict=True)],
        total=str(sum(map(Fraction, values))), context_sha256='a' * 64,
        direct_current_receiver_condition=True, legacy_pointwise_bounds_used=False,
        legacy_pointwise_shortfalls=[dict(saved_legacy_lower='1/4', status='UNUSED_IN_DIRECT_RECEIVER_PROOF')])


def whole(k):
    model = dict(model_key='b' * 64, context_sha256='a' * 64, current_total='1',
        current_corners=[dict(marks=p, lower='1/4') for p in [[0,1],[1,2],[2,3],[0,3]]])
    proof = dict(all_current_funded_source_models=[model], automatic_exhaustiveness_claim=False,
        full_completion_keys=[dict(key=k, funded_model_keys=['b' * 64])],
        geometric_applicability_premise='Synthetic whole applicability modulo independent current S8.')
    names = ['synthetic_source_' + str(i) for i in range(22)]
    return dict(shard=6, shards=8, complete_packet_keys=names, complete_reference_interface_count=564,
        packets=[dict(source=names[i], current_funded_whole_model_providers=[proof] if i == 6 else [])
                 for i in range(22) if i % 8 == 6])


class UnionControls(unittest.TestCase):
    def test_synthetic_three_line_actual_schema_accepted(self):
        m, e = actual_fixture()
        with m:
            inner, adapter = U.actual(e, {U.PRIMARY: PRIMARY_SHA}, m.controllers)
        self.assertEqual(adapter, 'obligations')
        self.assertEqual(inner['batch'], 'selected-00')

    def test_missing_ledger_refused(self):
        m, e = actual_fixture(); del m.mem[e['path']]
        with m, self.assertRaises(FileNotFoundError):
            U.actual(e, {U.PRIMARY: PRIMARY_SHA}, m.controllers)

    def test_failed_outer_even_with_consistent_hash_refused(self):
        m, e = actual_fixture(outer_change=lambda x: x.update(status='FAIL_BOUNDED_RESEARCH_CHECK'))
        with m, self.assertRaisesRegex(ValueError, 'actual bounded PASS required'):
            U.actual(e, {U.PRIMARY: PRIMARY_SHA}, m.controllers)

    def test_changed_ledger_bytes_refused(self):
        m, e = actual_fixture(); m.mem[e['path']] = b' ' + m.mem[e['path']]
        with m, self.assertRaisesRegex(ValueError, 'receipt/log byte binding'):
            U.actual(e, {U.PRIMARY: PRIMARY_SHA}, m.controllers)

    def test_public_worker_line_cannot_differ_from_embedded_stdout(self):
        m, e = actual_fixture(line_change=lambda lines: lines.__setitem__(0, lines[0] + ' '))
        with m, self.assertRaisesRegex(ValueError, 'stdout three-line schema'):
            U.actual(e, {U.PRIMARY: PRIMARY_SHA}, m.controllers)

    def test_missing_owned_tree_receipt_refused(self):
        m, e = actual_fixture(line_change=lambda lines: lines.pop())
        with m, self.assertRaisesRegex(ValueError, 'stdout three-line schema'):
            U.actual(e, {U.PRIMARY: PRIMARY_SHA}, m.controllers)

    def test_controller_hash_reassociation_refused(self):
        def mutate(lines):
            doc = json.loads(lines[2]); doc['controller_source_sha256'][U.RUNNER] = '0' * 64
            lines[2] = json.dumps(doc)
        m, e = actual_fixture(line_change=mutate)
        with m, self.assertRaisesRegex(ValueError, 'controller/assertion/environment provenance'):
            U.actual(e, {U.PRIMARY: PRIMARY_SHA}, m.controllers)

    def test_summary_cannot_change_the_recorded_thread_limit(self):
        def mutate(lines):
            doc = json.loads(lines[1]); doc['threads'] = 2; lines[1] = json.dumps(doc)
        m, e = actual_fixture(line_change=mutate)
        with m, self.assertRaisesRegex(ValueError, 'ledger/stdout summary agreement'):
            U.actual(e, {U.PRIMARY: PRIMARY_SHA}, m.controllers)

    def test_duplicate_json_fields_refused(self):
        with self.assertRaisesRegex(ValueError, 'duplicate JSON field'):
            U.loads('{"same":0,"same":1}')

    def test_same_marked_image_different_full_map_does_not_pay(self):
        m = Memory(); m.add(U.ROUTES, 'selected-00', selected()); m.add(PROPER, 'proper-0', proper(key(True)))
        result = m.run()
        self.assertEqual(result['joined_nonunit_keys'], 0)
        self.assertEqual(result['provider_keys_outside_this_shard'], 1)
        self.assertEqual(result['uncovered_completion_keys'][0]['key'], key())

    def test_duplicate_full_obligation_cannot_hide_in_a_count(self):
        m = Memory(); m.add(U.ROUTES, 'selected-00', selected([key(), key()]))
        with self.assertRaisesRegex(ValueError, 'disjoint complete-map shard keys'):
            m.run()

    def test_direct_unit_is_funded_despite_unused_legacy_quarter(self):
        U.funded(receiver(['0', '1/2', '0', '1/2']), full=True)

    def test_current_three_quarters_refused_despite_unused_old_surplus(self):
        fixture = receiver(['0', '1/2', '0', '1/4'])
        fixture['unused_old_total'] = '5/4'
        with self.assertRaisesRegex(ValueError, 'complete current receiver sum'):
            U.funded(fixture, full=True)

    def test_unbound_g4_request_is_not_a_provider(self):
        request = dict(result=dict(key=key(), status='UNBOUND_G6_CURRENT_OUTPUT_BRIDGE'))
        with self.assertRaisesRegex(ValueError, 'unsupported provider extraction'):
            list(U.providers(request, 'G4', []))

    def test_semantic_cap_refuses_instead_of_truncating(self):
        before = U.VISITS
        try:
            U.VISITS = U.MAX_VISITS
            with self.assertRaisesRegex(ValueError, 'semantic traversal cap'):
                U.full_key(key())
        finally:
            U.VISITS = before

    def test_local_key_coverage_never_discards_prefix_or_global_open_flags(self):
        prefix = [dict(index=1536, method='finite_leaf9_or_ideal_feature5', required_provider='synthetic open prefix')]
        m = Memory(); m.add(U.ROUTES, 'selected-00', selected(prefix=prefix)); m.add(PROPER, 'proper-0', proper(key()))
        result = m.run()
        self.assertTrue(result['selected_nonunit_key_union_complete'])
        self.assertEqual(result['prefix_obligations_still_open'], prefix)
        self.assertFalse(result['selected_shard_all_routes_closed'])
        self.assertFalse(result['current_terminal_bridge_complete'])
        self.assertFalse(result['full_finite57_verified'])
        self.assertTrue(any('Current S8' in item for item in result['unresolved_global_interfaces']))


    def test_schema_safe_whole_exact_key_is_indexed_with_open_s8(self):
        m = Memory(); m.add(U.ROUTES, 'selected-00', selected()); m.add(WHOLE, 'whole-6', whole(key()))
        result = m.run()
        self.assertTrue(result['selected_nonunit_key_union_complete'])
        provider = result['exact_joined_key_bindings'][0]['providers'][0]
        self.assertEqual(provider['adapter'], 'whole_schema_safe')
        self.assertEqual(provider['application_prerequisites'], [dict(kind='independent_current_S8',
            status='UNVERIFIED_EXTERNAL_TO_THIS_JOIN', discharged_by_key_match=False)])
        self.assertFalse(result['selected_shard_all_routes_closed'])
        self.assertFalse(result['current_terminal_bridge_complete'])
        self.assertFalse(result['full_finite57_verified'])

    def test_schema_safe_whole_reversed_full_map_is_not_a_match(self):
        m = Memory(); m.add(U.ROUTES, 'selected-00', selected()); m.add(WHOLE, 'whole-6', whole(key(True)))
        result = m.run()
        self.assertFalse(result['selected_nonunit_key_union_complete'])
        self.assertEqual(len(result['uncovered_completion_keys']), 1)

    def test_schema_safe_whole_unfunded_model_is_rejected(self):
        result = whole(key()); result['packets'][0]['current_funded_whole_model_providers'][0][
            'all_current_funded_source_models'][0]['current_total'] = '0'
        m = Memory(); m.add(U.ROUTES, 'selected-00', selected()); m.add(WHOLE, 'whole-6', result)
        with self.assertRaisesRegex(ValueError, 'complete current receiver sum'):
            m.run()

    def test_schema_safe_whole_exhaustiveness_is_never_inferred(self):
        result = whole(key()); result['packets'][0]['current_funded_whole_model_providers'][0][
            'automatic_exhaustiveness_claim'] = True
        m = Memory(); m.add(U.ROUTES, 'selected-00', selected()); m.add(WHOLE, 'whole-6', result)
        with self.assertRaisesRegex(ValueError, 'variant geometry'):
            m.run()

    def test_schema_safe_whole_wrong_packet_shard_is_rejected(self):
        result = whole(key()); result['packets'][0]['source'] = 'synthetic_source_7'
        m = Memory(); m.add(U.ROUTES, 'selected-00', selected()); m.add(WHOLE, 'whole-6', result)
        with self.assertRaisesRegex(ValueError, 'entire named whole packet shard'):
            m.run()

    def test_schema_safe_whole_empty_applicability_premise_is_rejected(self):
        result = whole(key()); result['packets'][0]['current_funded_whole_model_providers'][0][
            'geometric_applicability_premise'] = ''
        m = Memory(); m.add(U.ROUTES, 'selected-00', selected()); m.add(WHOLE, 'whole-6', result)
        with self.assertRaisesRegex(ValueError, 'whole applicability'):
            m.run()

    def test_proper_exact_key_still_requires_independent_current_s8(self):
        m = Memory(); m.add(U.ROUTES, 'selected-00', selected()); m.add(PROPER, 'proper-0', proper(key()))
        result = m.run(); self.assertTrue(result['selected_nonunit_key_union_complete'])
        premise = result['exact_joined_key_bindings'][0]['providers'][0]['application_prerequisites'][0]
        self.assertEqual(premise['kind'], 'independent_current_S8')
        self.assertEqual(premise['status'], 'UNVERIFIED_EXTERNAL_TO_THIS_JOIN')
        self.assertFalse(premise['discharged_by_key_match'])
        self.assertFalse(result['source_geometric_and_computed_premises_proved_by_this_checker'])


def main():
    global U, SOURCES
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent-lock', action='store_true', required=True)
    parser.parse_args(); parent_lock()
    relative = str(Path(__file__).resolve().relative_to(ROOT)); before = raw(relative)
    SOURCES = {p: raw(p) for p in SOURCE_SHA}
    need({p: sha256(b).hexdigest() for p, b in SOURCES.items()} == SOURCE_SHA, 'exact union-test source closure')
    U = types.ModuleType('finite57_actual_receipt_union_control_tests'); U.__file__ = str(ROOT / UNION)
    exec(compile(SOURCES[UNION], str(ROOT / UNION), 'exec'), U.__dict__)
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(UnionControls)
    need(suite.countTestCases() == 23, 'exact twenty-three union/whole-schema control cases')
    stream = io.StringIO(); result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    need(raw(relative) == before and all(raw(p) == b for p, b in SOURCES.items()), 'union-test source changed')
    payload = json.dumps(dict(status='PASS_TINY_ACTUAL_RECEIPT_UNION_CONTROLS' if result.wasSuccessful()
        else 'FAIL_TINY_ACTUAL_RECEIPT_UNION_CONTROLS', source_sha256=sha256(before).hexdigest(),
        sources=SOURCE_SHA, tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        output=stream.getvalue(),
        scope='Twenty-three in-memory synthetic parser/key/receiver-control cases only. Faked receipts are not actual runs. No scientific worker is imported; source bytes are read only. Filesystem/monitor correctness, actual provider execution, geometric premises, all64 coverage and finite57 are not verified.'), sort_keys=True)
    need(len(payload.encode()) <= 500 * 1024, 'tiny union-control output cap')
    print(payload)
    if not result.wasSuccessful():
        raise SystemExit(1)


if __name__ == '__main__':
    main()
