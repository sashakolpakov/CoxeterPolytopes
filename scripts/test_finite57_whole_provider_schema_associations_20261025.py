#!/usr/bin/env python3
"""SOURCE ONLY: fifteen tiny whole-provider schema/path-association tests.

Synthetic metadata and context hashes test selectors and associations only.
No physical Gram, saved scientific JSON, completion search or worker runs.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
from copy import deepcopy
import fcntl
from hashlib import sha256
import io
import json
from pathlib import Path
import types
import unittest

ROOT = Path(__file__).resolve().parents[1]
PROVIDER = 'scripts/index_finite57_funded_selected_whole_providers_schema_safe_20261025.py'
PROVIDER_SHA = '6ace38a526ccf22d9913ffe0b814b7ee4177d74c3be478295a981fef5b7afb03'
P = None


def need(ok, message):
    if not ok:
        raise ValueError(message)


def raw(relative):
    with (ROOT / relative).open('rb') as stream:
        data = stream.read(128 * 1024 + 1)
    need(0 < len(data) <= 128 * 1024, 'bounded whole-schema test source read')
    return data


def parent_lock():
    with (ROOT / 'formal/.proof-batches/exclusive.lock').open('r+') as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return
        fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
    raise ValueError('shared parent lock must already be held')


def alternate():
    indices = [1035, 1049, 1905, 1964, 2298, 2357, 2421, 2429]
    return dict(proofs=[dict(index=i, old_seal=dict(index=i)) for i in indices],
        complete_eligible_patterns=394, validated_saved_maps=578,
        affected_contexts=[dict(index=i, maps=[{'fixture': 'first'}, {'fixture': 'second'}]) for i in indices],
        affected_maps=16, all_selected_components_whole_D16=True,
        proper_displayed_subcomponent_cases=0, geometry_census_repeated=False)


def interfaces():
    path = 'selected_proofs/proofs/0/corners/2'
    required = {path: dict(marks=[3, 7], context_sha256='a' * 64)}
    wanted = {path: dict(path=path, marks=[3, 7])}
    covered = {path: dict(marks=[3, 7], context_sha256='a' * 64, model_key='b' * 64)}
    return path, required, wanted, covered


class WholeSchema(unittest.TestCase):
    def test_ordinary_list_keeps_its_actual_selector(self):
        rows = [dict(index=1664)]
        selector, got = P.indexed_rows('ordinary_fixture', 'selected_proofs', rows)
        self.assertEqual(selector, ['selected_proofs'])
        self.assertIs(got, rows)

    def test_actual_dictionary_wrapper_keeps_proofs_selector(self):
        container = alternate()
        selector, rows = P.indexed_rows(P.ALT, 'selected_proofs', container)
        self.assertEqual(selector, ['selected_proofs', 'proofs'])
        self.assertIs(rows, container['proofs'])
        self.assertEqual([p['index'] for p in rows], P.ALT_INDICES)

    def test_dictionary_at_an_unregistered_source_refused(self):
        with self.assertRaisesRegex(ValueError, 'ordinary indexed proof section must be a list'):
            P.indexed_rows('ordinary_fixture', 'selected_proofs', alternate())

    def test_bare_list_at_the_actual_dictionary_source_refused(self):
        with self.assertRaisesRegex(ValueError, 'exact D16 alternate dictionary wrapper schema'):
            P.indexed_rows(P.ALT, 'selected_proofs', alternate()['proofs'])

    def test_missing_proofs_field_refused(self):
        row = alternate(); del row['proofs']
        with self.assertRaisesRegex(ValueError, 'exact D16 alternate dictionary wrapper schema'):
            P.indexed_rows(P.ALT, 'selected_proofs', row)

    def test_unknown_wrapper_field_refused(self):
        row = alternate(); row['another_proof_list'] = row['proofs']
        with self.assertRaisesRegex(ValueError, 'exact D16 alternate dictionary wrapper schema'):
            P.indexed_rows(P.ALT, 'selected_proofs', row)

    def test_equal_counts_do_not_excuse_affected_case_reassociation(self):
        row = alternate()
        row['affected_contexts'][0], row['affected_contexts'][1] = row['affected_contexts'][1], row['affected_contexts'][0]
        with self.assertRaisesRegex(ValueError, 'proof/affected-context exact association'):
            P.indexed_rows(P.ALT, 'selected_proofs', row)

    def test_inherited_seal_must_belong_to_same_high_row(self):
        rows = [dict(index=1664, old_seal=dict(index=1665))]
        with self.assertRaisesRegex(ValueError, 'inherited seal high-row association'):
            P.indexed_rows('ordinary_fixture', 'selected_proofs', rows)

    def test_boolean_is_not_a_high_index(self):
        with self.assertRaisesRegex(ValueError, 'typed authoritative high ID'):
            P.indexed_rows('ordinary_fixture', 'selected_proofs', [dict(index=True)])

    def test_actual_dictionary_path_is_funded_without_rewriting(self):
        path, required, wanted, covered = interfaces()
        self.assertEqual(P.bind_interfaces(required, wanted, covered), {path: covered[path]})

    def test_artificial_slash_zero_alias_does_not_match_actual_path(self):
        path, required, wanted, covered = interfaces()
        wanted['selected_proofs/0/corners/2'] = wanted.pop(path)
        with self.assertRaisesRegex(ValueError, 'actual/reference interface path equality'):
            P.bind_interfaces(required, wanted, covered)

    def test_unfunded_required_occurrence_refused(self):
        path, required, wanted, covered = interfaces(); covered.clear()
        with self.assertRaisesRegex(ValueError, 'lacks a full DIRECT current receipt'):
            P.bind_interfaces(required, wanted, covered)

    def test_same_unordered_marks_do_not_excuse_ordered_mismatch(self):
        path, required, wanted, covered = interfaces(); covered[path]['marks'] = [7, 3]
        with self.assertRaisesRegex(ValueError, 'exact ordered marks'):
            P.bind_interfaces(required, wanted, covered)

    def test_right_path_wrong_physical_context_refused(self):
        path, required, wanted, covered = interfaces(); covered[path]['context_sha256'] = 'c' * 64
        with self.assertRaisesRegex(ValueError, 'exact original or extended Gram/white'):
            P.bind_interfaces(required, wanted, covered)

    def test_extra_unlisted_required_path_refused(self):
        path, required, wanted, covered = interfaces()
        required[path + '/extra'] = deepcopy(required[path])
        covered[path + '/extra'] = deepcopy(covered[path])
        with self.assertRaisesRegex(ValueError, 'actual/reference interface path equality'):
            P.bind_interfaces(required, wanted, covered)


def main():
    global P
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent-lock', action='store_true', required=True)
    parser.parse_args(); parent_lock()
    relative = str(Path(__file__).resolve().relative_to(ROOT))
    before, encoded = raw(relative), raw(PROVIDER)
    need(sha256(encoded).hexdigest() == PROVIDER_SHA, 'exact schema-safe whole-provider source')
    P = types.ModuleType('finite57_schema_safe_whole_provider_tests')
    P.__file__ = str(ROOT / PROVIDER)
    exec(compile(encoded, str(ROOT / PROVIDER), 'exec'), P.__dict__)
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(WholeSchema)
    need(suite.countTestCases() == 15, 'exact fifteen new schema-association cases')
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    need(raw(relative) == before and raw(PROVIDER) == encoded, 'whole-schema test/provider source changed')
    payload = json.dumps(dict(
        status='PASS_TINY_WHOLE_PROVIDER_SCHEMA_ASSOCIATIONS' if result.wasSuccessful()
        else 'FAIL_TINY_WHOLE_PROVIDER_SCHEMA_ASSOCIATIONS',
        source_sha256=sha256(before).hexdigest(), sources={PROVIDER: PROVIDER_SHA},
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors), output=stream.getvalue(),
        scope='Fifteen synthetic selector/path/mark/context association cases only. No actual Gram, inverse, enumeration, saved scientific JSON, whole22 replay, Lean/TeX or final finite57 proof is verified.'),
        sort_keys=True)
    need(len(payload.encode()) <= 500 * 1024, 'tiny whole-schema output cap')
    print(payload)
    if not result.wasSuccessful():
        raise SystemExit(1)


if __name__ == '__main__':
    main()
