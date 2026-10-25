#!/usr/bin/env python3
"""SOURCE ONLY: two direct-current receiver control regressions.

Corner values are explicit stubs. These are receiver-logic tests, not
geometric examples or a replay of any scientific certificate.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')
import argparse
from hashlib import sha256
import io
import json
from pathlib import Path
import types
import unittest

ROOT = Path(__file__).resolve().parents[1]
COLLECTOR = 'scripts/check_finite57_terminal_collectors_direct_receivers_20261025.py'
COLLECTOR_SHA = 'a2a9d61a7637e45d5a4d7af1c3283e120ca7685d07e9f686ae7b6ee585520eca'
ENGINE = 'scripts/finite57_terminal_literal_engine_small_affine_20261025.py'
ENGINE_SHA = '39dd55d2b2abef3db045e646edcc95bc1c8363f483d9f853d8c1d4e1b65d0e1f'
C = E = None


def need(ok, message):
    if not ok:
        raise ValueError(message)


def raw(relative):
    with (ROOT / relative).open('rb') as stream:
        data = stream.read(128 * 1024 + 1)
    need(0 < len(data) <= 128 * 1024, 'bounded receiver-test source read')
    return data


def module(relative, digest, name):
    data = raw(relative)
    need(sha256(data).hexdigest() == digest, 'receiver-test source identity')
    obj = types.ModuleType(name)
    obj.__file__ = str(ROOT / relative)
    exec(compile(data, str(ROOT / relative), 'exec'), obj.__dict__)
    return obj


def receipt(current):
    marks = [[0, 1], [1, 2], [2, 3], [0, 3]]
    lookup = {tuple(m): value for m, value in zip(marks, current, strict=True)}
    old = ['0', '1/2', '1/4', '1/2']
    C.value = lambda gram, white, pair, *args, **kwargs: dict(
        weight=lookup[tuple(pair)], marks=pair, kind='finite', forest=[],
        sealed_witness_components=[], ignored_unsealed_components=[])
    C.FUEL = [0]
    return C.receipt('synthetic_receiver', 26, [], [],
                     [dict(marks=m) for m in marks], 'explicit stub values',
                     [dict(marks=m, new_lower=value) for m, value in zip(marks, old, strict=True)])


class DirectReceiver(unittest.TestCase):
    def test_current_unit_accepted_and_unused_quarter_reported(self):
        result = receipt(['0', '1/2', '0', '1/2'])
        self.assertEqual(result['total'], '1')
        self.assertTrue(result['direct_current_receiver_condition'])
        self.assertFalse(result['legacy_pointwise_bounds_used'])
        shortfalls = result['legacy_pointwise_shortfalls']
        self.assertEqual(len(shortfalls), 1)
        self.assertEqual(shortfalls[0]['marks'], [2, 3])
        self.assertEqual(shortfalls[0]['saved_legacy_lower'], '1/4')
        self.assertEqual(shortfalls[0]['conservative_current_lower'], '0')
        self.assertEqual(shortfalls[0]['status'], 'UNUSED_IN_DIRECT_RECEIVER_PROOF')

    def test_current_three_quarters_refused_despite_old_surplus(self):
        with self.assertRaisesRegex(ValueError, 'current distinct-corner receipt below one'):
            receipt(['0', '1/2', '0', '1/4'])


def main():
    global C, E
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent-lock', action='store_true', required=True)
    parser.parse_args()
    relative = str(Path(__file__).resolve().relative_to(ROOT))
    before = raw(relative)
    E = module(ENGINE, ENGINE_SHA, 'direct_receiver_test_engine')
    C = module(COLLECTOR, COLLECTOR_SHA, 'direct_receiver_test_collector')
    C.E = E
    C.parent_lock()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(DirectReceiver)
    need(suite.countTestCases() == 2, 'exact two new receiver regressions')
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    need(raw(relative) == before and sha256(raw(COLLECTOR)).hexdigest() == COLLECTOR_SHA
         and sha256(raw(ENGINE)).hexdigest() == ENGINE_SHA, 'receiver-test sources changed')
    payload = json.dumps(dict(
        status='PASS_DIRECT_CURRENT_RECEIVER_REGRESSIONS' if result.wasSuccessful()
        else 'FAIL_DIRECT_CURRENT_RECEIVER_REGRESSIONS',
        source_sha256=sha256(before).hexdigest(),
        sources={COLLECTOR: COLLECTOR_SHA, ENGINE: ENGINE_SHA},
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        output=stream.getvalue(),
        scope='Two new explicit-value receiver-logic tests only; no geometry, scientific JSON, worker, Lean/TeX or full finite57 check.'), sort_keys=True)
    need(len(payload.encode()) <= 500 * 1024, 'receiver-test output cap')
    print(payload)
    if not result.wasSuccessful():
        raise SystemExit(1)


if __name__ == '__main__':
    main()
