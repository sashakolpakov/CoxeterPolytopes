#!/usr/bin/env python3
"""SOURCE ONLY: thirteen tiny proper-provider contact-association regressions.

The inverse entries below are synthetic association fixtures. They are not
claimed Gram inverses, geometric examples or a replay of a saved certificate.
Only the pinned provider's bind_contacts helper is called at test runtime.
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
PROVIDER = 'scripts/index_finite57_all_proper_finite_completion_providers_20261025.py'
PROVIDER_SHA = '0ef6647b6b9b702f076b264b730caef89a31e40ff2a502bf5a5e03ab1fe29ba3'
P = None


def need(ok, message):
    if not ok:
        raise ValueError(message)


def raw(relative):
    with (ROOT / relative).open('rb') as stream:
        data = stream.read(128 * 1024 + 1)
    need(0 < len(data) <= 128 * 1024, 'bounded proper-contact source read')
    return data


def parent_lock():
    with (ROOT / 'formal/.proof-batches/exclusive.lock').open('r+') as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return
        fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
    raise ValueError('shared parent lock must already be held')


def fixture(saturated=False):
    base, images = [7, 8], [0, 1]
    inverse = [['3' if i == j else '0' for j in range(9)] for i in range(9)]
    if saturated:
        target = [{1: 3}, {0: 3, 2: 3, 3: 3}, {1: 3}, {1: 3}]
        inverse[8][8] = '2'
        row = dict(missing_target_nodes=[2, 3],
                   required_contacts={'2': [[1, 8, 3]], '3': [[1, 8, 3]]},
                   forced_node=8, outcome='two_saturated_short_leaves_coincide',
                   projection_square='2')
    else:
        target = [{1: 3}, {0: 3, 2: 3}, {1: 3}]
        row = dict(missing_target_nodes=[2], required_contacts=[[1, 8, 3]],
                   outcome='forced_new_white_projection_exceeds_norm2',
                   projection_lower_bound='3')
    return row, target, images, base, dict(inverse=inverse)


def bind(parts, method='proper_completion_exclusion'):
    row, target, images, base, core = parts
    return P.bind_contacts(None, row, target, images, base, core, method)


class ProperContacts(unittest.TestCase):
    def test_strict_single_missing_list_accepted(self):
        self.assertEqual(bind(fixture()), [(2, 1, 8, 3)])

    def test_saturated_two_missing_dictionary_accepted(self):
        self.assertEqual(bind(fixture(True)), [(2, 1, 8, 3), (3, 1, 8, 3)])

    def test_direct_strict_forced_contact_branch_accepted(self):
        parts = fixture()
        row = parts[0]
        del row['required_contacts']
        row.update(forced_old_target_node=1, core_contact=8, label=3)
        self.assertEqual(bind(parts, 'impossible_proper_completion'), [(2, 1, 8, 3)])

    def test_wrong_missing_target_rejected(self):
        parts = fixture(); parts[0]['missing_target_nodes'] = [1]
        with self.assertRaisesRegex(ValueError, 'exact full-map missing target nodes'):
            bind(parts)

    def test_equal_inverse_at_wrong_contact_still_rejected(self):
        parts = fixture(); parts[0]['required_contacts'] = [[1, 7, 3]]
        self.assertEqual(parts[4]['inverse'][7][7], parts[4]['inverse'][8][8])
        with self.assertRaisesRegex(ValueError, 'complete exact rank14 required-contact list'):
            bind(parts)

    def test_saved_inverse_wrong_at_actual_contact_rejected(self):
        parts = fixture(); parts[0]['projection_lower_bound'] = '4'
        with self.assertRaisesRegex(ValueError, 'strict projection bound at its actual forced contact'):
            bind(parts)

    def test_equal_two_is_not_a_strict_exclusion(self):
        parts = fixture()
        parts[0]['projection_lower_bound'] = '2'; parts[4]['inverse'][8][8] = '2'
        with self.assertRaisesRegex(ValueError, 'strict projection bound at its actual forced contact'):
            bind(parts)

    def test_saturated_branch_requires_exact_two(self):
        parts = fixture(True)
        parts[0]['projection_square'] = '3'; parts[4]['inverse'][8][8] = '3'
        with self.assertRaisesRegex(ValueError, 'saturated fork contact at its actual inverse entry'):
            bind(parts)

    def test_missing_dictionary_entry_rejected(self):
        parts = fixture(True); del parts[0]['required_contacts']['3']
        with self.assertRaisesRegex(ValueError, 'complete exact rank14 required-contact list'):
            bind(parts)

    def test_noncanonical_dictionary_node_rejected(self):
        parts = fixture(True)
        parts[0]['required_contacts']['02'] = parts[0]['required_contacts'].pop('2')
        with self.assertRaisesRegex(ValueError, 'canonical target node key'):
            bind(parts)

    def test_wrong_native_contact_label_rejected(self):
        parts = fixture(); parts[0]['required_contacts'][0][2] = 4
        with self.assertRaisesRegex(ValueError, 'complete exact rank14 required-contact list'):
            bind(parts)

    def test_single_missing_list_cannot_hide_a_second_missing_node(self):
        parts = fixture(True); parts[0]['required_contacts'] = [[1, 8, 3]]
        with self.assertRaisesRegex(ValueError, 'single-missing-node contact-list schema'):
            bind(parts)

    def test_changed_full_map_cannot_reuse_saved_contacts(self):
        parts = fixture(); parts[2][:] = [0, 2]; parts[0]['missing_target_nodes'] = [1]
        with self.assertRaisesRegex(ValueError, 'complete exact rank14 required-contact list'):
            bind(parts)


def main():
    global P
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent-lock', action='store_true', required=True)
    parser.parse_args()
    parent_lock()
    relative = str(Path(__file__).resolve().relative_to(ROOT))
    before, encoded = raw(relative), raw(PROVIDER)
    need(sha256(encoded).hexdigest() == PROVIDER_SHA, 'exact proper-provider source identity')
    P = types.ModuleType('finite57_proper_contact_test_provider')
    P.__file__ = str(ROOT / PROVIDER)
    exec(compile(encoded, str(ROOT / PROVIDER), 'exec'), P.__dict__)
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(ProperContacts)
    need(suite.countTestCases() == 13, 'exact thirteen new association regressions')
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    need(raw(relative) == before and raw(PROVIDER) == encoded, 'contact-test/provider source changed')
    payload = json.dumps(dict(
        status='PASS_TINY_PROPER_CONTACT_ASSOCIATIONS' if result.wasSuccessful()
        else 'FAIL_TINY_PROPER_CONTACT_ASSOCIATIONS',
        source_sha256=sha256(before).hexdigest(), sources={PROVIDER: PROVIDER_SHA},
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        output=stream.getvalue(),
        scope='Thirteen synthetic bind_contacts association tests only. No saved inverse, geometric theorem, full-map enumeration, scientific JSON, Lean/TeX or final-provider coverage is verified.'),
        sort_keys=True)
    need(len(payload.encode()) <= 500 * 1024, 'tiny contact-test output cap')
    print(payload)
    if not result.wasSuccessful():
        raise SystemExit(1)


if __name__ == '__main__':
    main()
