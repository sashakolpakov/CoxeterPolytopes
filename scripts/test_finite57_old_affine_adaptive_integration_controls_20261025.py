#!/usr/bin/env python3
"""Retain56 controls; add five adaptive identities and one helper-wiring case.

Only pinned definitions and synthetic objects are loaded. No actual HIGH,
old-affine application helper or scientific data is loaded. No mathematical
worker is prepared or run.
Root must freeze FINAL_SHA and accept the source before bounded execution.
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
PREVIOUS = 'scripts/test_finite57_source_application_high_projection_controls_20261025.py'
PREVIOUS_SHA = 'b2630300393ef13c1594275de8a6cee755818630f565c9216db080320dd3bc99'
OLD55_SHA = 'b9d1d4f56654f34fb8cb9155b9afc5b57da6283e73343e1c569a4e906114d7ea'
OLD47_SHA = 'e51d172fb858b904999bcf8c580e98e0e8492a40ed28853abf9d675bac6c4d9a'
OLD23_SHA = '7f9801491d5b9f21ba669155f46d9d5469b4328880770f39fbec4ef013be4640'
ADAPTIVE = 'results/finite57_adaptive_prefix_identity_proposal_20261025.py'
ADAPTIVE_SHA = 'c16d95611578fae45e08327c6fccfb421c5ba216ce91672aa9d3bcd87e1afc06'
FINAL = 'scripts/check_finite57_current_rule_conjunction_20261025.py'
FINAL_SHA = '3b9843ba745e8e760007dc5647cdcddfbb78354121109401f35f75b5c3c7240c'
F = None


def raw(name):
    with (ROOT / name).open('rb') as stream:
        data = stream.read(128 * 1024 + 1)
    if not 0 < len(data) <= 128 * 1024:
        raise ValueError('bounded integration-control source')
    return data


def load(name, path, digest):
    data = raw(path)
    if sha256(data).hexdigest() != digest:
        raise ValueError('exact integration-control source identity')
    obj = types.ModuleType(name)
    obj.__file__ = str(ROOT / path)
    exec(compile(data, obj.__file__, 'exec'), obj.__dict__)
    return obj


class OldAffineHelperWiringControls(unittest.TestCase):
    def test_exact_pin_before_load_explicit_injections_and_missing_review_refusal(self):
        path = F.OLD_AFFINE_APPLICATION
        digest = 'a' * 64
        review = dict(source_sha256=digest,
                      path='results/SYNTHETIC_ACCEPTANCE_ONLY.json', sha256='b' * 64)
        legacy = types.SimpleNamespace(synthetic_no_realization=True)
        helper = types.SimpleNamespace(synthetic_no_realization=True)
        events = []
        saved = F.load, F.hash_file

        def pin(actual_path, actual_digest):
            if (actual_path, actual_digest) != (path, digest):
                raise ValueError('synthetic exact old-affine source pin')
            events.append(('pin', actual_path, actual_digest))

        def fake_load(name, actual_path, actual_digest):
            self.assertEqual((actual_path, actual_digest), (path, digest))
            self.assertEqual(events, [('pin', path, digest)])
            events.append(('load', name, actual_path, actual_digest))
            return helper

        try:
            F.load = fake_load
            F.hash_file = pin
            got = F.old_affine_application_helper({path: review}, legacy)
            self.assertIs(got, helper)
            self.assertIs(got.SOURCE_PINNER, pin)
            self.assertIs(got.SOURCE_READER, F.application_source_document)
            self.assertEqual(got.APPLICATION_SOURCE_SHA, digest)
            self.assertEqual(got.SOURCE_ACCEPTANCE,
                             dict(path=review['path'], sha256=review['sha256']))
            self.assertIs(got.LEGACY_HELPER, legacy)
            self.assertEqual(len(events), 2)
            self.assertEqual(events[0], ('pin', path, digest))
            self.assertEqual(events[1][0], 'load')
            self.assertEqual(events[1][2:], (path, digest))
            before = list(events)
            with self.assertRaisesRegex(ValueError, 'missing exact reviewed old-affine'):
                F.old_affine_application_helper({}, legacy)
            self.assertEqual(events, before)
        finally:
            F.load, F.hash_file = saved


def main():
    global F
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent-lock', action='store_true', required=True)
    parser.parse_args()
    own = str(Path(__file__).resolve().relative_to(ROOT))
    before = raw(own)
    W = load('finite57_unchanged56_definitions', PREVIOUS, PREVIOUS_SHA)
    if W.PREVIOUS_SHA != OLD55_SHA:
        raise ValueError('unchanged55 source identity')
    V = W.load('finite57_unchanged55_definitions', W.PREVIOUS, OLD55_SHA)
    if V.PREVIOUS_SHA != OLD47_SHA:
        raise ValueError('unchanged47 source identity')
    P = V.load('finite57_unchanged47_definitions', V.PREVIOUS, OLD47_SHA)
    if P.OLD_SHA != OLD23_SHA:
        raise ValueError('unchanged23 source identity')
    T = P.load('finite57_unchanged23_definitions', P.OLD, OLD23_SHA)
    T.parent_lock()
    T.SOURCES = {path: raw(path) for path in T.SOURCE_SHA}
    if {p: sha256(b).hexdigest() for p, b in T.SOURCES.items()} != T.SOURCE_SHA:
        raise ValueError('unchanged predecessor fixture source closure')
    T.U = P.load('finite57_control_union435e', T.UNION, T.UNION_SHA)
    F = load('finite57_old_affine_adaptive_conjunction', FINAL, FINAL_SHA)
    F.U = T.U
    P.F = F
    P.T = T
    P.C = P.load('finite57_release_codec_definitions', P.CODEC, P.CODEC_SHA)
    S = V.load('finite57_unchanged_s8_definitions', V.S8, V.S8_SHA)
    V.F, V.P, V.S = F, P, S
    W.F, W.V = F, V
    A = load('finite57_adaptive_identity_fixture_definitions', ADAPTIVE, ADAPTIVE_SHA)
    A.IDENTITY = F.adaptive_prefix_identity
    if sha256(raw(F.LEGACY_REGISTRY)).hexdigest() != F.LEGACY_REGISTRY_SHA:
        raise ValueError('unchanged historical fixture identity')
    classes = [T.UnionControls, P.FinalControls, P.CodecControls,
               P.ExpandedJSONControls, P.ConsumerAdapterControls,
               P.CurationTrustControls, V.RouteReuseControls, V.DirectS8Controls,
               V.SourceReaderControls, W.HighProjectionControls,
               A.AdaptivePrefixIdentityControls, OldAffineHelperWiringControls]
    suites = [unittest.defaultTestLoader.loadTestsFromTestCase(cls) for cls in classes]
    if [s.countTestCases() for s in suites] != [23, 15, 4, 1, 2, 2, 4, 3, 1, 1, 5, 1]:
        raise ValueError('exact56 retained definitions plus five adaptive and one wiring case')
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(unittest.TestSuite(suites))
    sources = {**T.SOURCE_SHA, PREVIOUS: PREVIOUS_SHA, W.PREVIOUS: OLD55_SHA,
               V.PREVIOUS: OLD47_SHA, P.OLD: OLD23_SHA,
               FINAL: FINAL_SHA, ADAPTIVE: ADAPTIVE_SHA, V.S8: V.S8_SHA,
               P.CODEC: P.CODEC_SHA, F.LEGACY_REGISTRY: F.LEGACY_REGISTRY_SHA}
    if raw(own) != before or any(sha256(raw(p)).hexdigest() != h for p, h in sources.items()):
        raise ValueError('source changed during integration controls')
    output = dict(
        status=('PASS' if result.wasSuccessful() else 'FAIL') +
               '_OLD_AFFINE_ADAPTIVE_INTEGRATION_CONTROLS_ONLY',
        source_sha256=sha256(before).hexdigest(), sources=sources,
        tests_run=result.testsRun, prepared_case_count=62,
        retained_high_projection_cases=56, new_adaptive_identity_cases=5,
        new_old_affine_wiring_cases=1, failures=len(result.failures),
        errors=len(result.errors), output=stream.getvalue(),
        scientific_data_read=False, geometric_theorems_tested=False,
        application_helpers_executed=False, actual_receipt_execution_fixtures=False,
        scope='Exactly56 frozen synthetic definitions against the new reader, five '
              'adaptive actual-corner identity fixtures and one mocked old-affine '
              'loader-wiring fixture. No actual HIGH, application helper, native '
              'enumeration, numerical geometry or global S8 is executed.')
    encoded = json.dumps(output, sort_keys=True)
    if len(encoded.encode()) + 1 > 500 * 1024:
        raise ValueError('complete integration-control output cap')
    print(encoded)
    if not result.wasSuccessful():
        raise SystemExit(1)


if __name__ == '__main__':
    main()
