#!/usr/bin/env python3
"""Synthetic source-association controls for the immutable provider curator.

Five test methods use only tiny in-memory schema fixtures. These are not
diagrams, realized examples, geometric arguments, actual receipt coverage or
proof-premise acceptance. The curator's main function is never called.
Root alone launches this script after source review under the shared lock.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
from contextlib import contextmanager
from copy import deepcopy
import fcntl
from hashlib import sha256
import io
import json
from pathlib import Path, PurePosixPath
import types
import unittest

ROOT = Path(__file__).resolve().parents[1]
OWN = 'scripts/test_finite57_provider_occurrence_curation_controls_20261025.py'
SOURCE = 'scripts/curate_finite57_provider_occurrences_20261025.py'
SOURCE_SHA = 'aefd1a1db6a3edc52846316e969c189ecc12e22028d9151c2a264e0d09c3658c'
SOURCE_CAP = 128 * 1024
OUTPUT_CAP = 500 * 1024
C = None


def need(ok, message):
    if not ok:
        raise ValueError(message)


def raw(name):
    lexical = PurePosixPath(name)
    need(type(name) is str and not lexical.is_absolute() and '..' not in lexical.parts,
         'bounded control-source path')
    path = (ROOT/name).resolve()
    need(path.is_relative_to(ROOT) and path.is_file(), 'existing control source')
    before = path.stat()
    identity = (before.st_dev,before.st_ino,before.st_size,before.st_mtime_ns,before.st_ctime_ns)
    with path.open('rb') as stream:
        data = stream.read(SOURCE_CAP+1)
    after = path.stat()
    need(0 < len(data) <= SOURCE_CAP and len(data) == before.st_size
         and identity == (after.st_dev,after.st_ino,after.st_size,after.st_mtime_ns,after.st_ctime_ns),
         'stable complete bounded control source')
    return data


def parent_lock():
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('r+') as stream:
        try:
            fcntl.flock(stream.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            pass
        else:
            fcntl.flock(stream.fileno(),fcntl.LOCK_UN)
            raise ValueError('shared wrapper proof lock must already be held')


def key(affine=False):
    # A nine-field representation fixture, with no asserted graph realization.
    return dict(high_index=7,core_key='SYNTHETIC_CORE_NO_REALIZATION',
        base_nodes=[0,1,2,3,4,5,6,7,8,9],selected_marks=[0,9],
        target_kind='C' if affine else 'D',target_order=13 if affine else 12,
        affine=affine,induced_map=[0,1,2,3,4,5,6,7,8,9],marked_image=[0,9])


def pin_for(stem):
    return {'results/'+stem+'.json':'e'*64}


def spec(adapter,batch,stem):
    return dict(adapter=adapter,rows=[dict(batch=batch,source=stem,
        proof=dict(finite_exchanges='SYNTHETIC_NOT_AN_ACCEPTED_THEOREM'))])


@contextmanager
def synthetic_sources(documents):
    """Exercise selectors against in-memory records; forbid all source I/O."""
    saved = C.source,C.read,C.pin
    def fixture_source(stem,pins):
        name = 'results/'+stem+'.json'
        C.need(stem in documents and pins.get(name) == 'e'*64,
               'exact synthetic source/pin association')
        return documents[stem]
    def no_read(*args,**kwargs):
        raise AssertionError('scientific source I/O is forbidden in synthetic controls')
    C.source,C.read,C.pin = fixture_source,no_read,no_read
    try:
        yield
    finally:
        C.source,C.read,C.pin = saved


def exported(inner,worker_spec,pins,documents):
    with synthetic_sources(documents):
        rows = list(C.origins(inner,worker_spec,pins))
        for complete_key,_,_ in rows:
            C.full_key(complete_key)
        return rows


def proper_fixture():
    stem = 'SYNTHETIC_PROPER_SOURCE_NOT_A_FILE'
    complete_key = key()
    record = dict(map=deepcopy(complete_key['induced_map']),
                  proof=dict(kind='SYNTHETIC_CONTACT_RECORD'))
    original = dict(index=7,core_key=complete_key['core_key'],
                    method='proper_completion_exclusion',maps=[record])
    mapping = dict(key=complete_key,source_path='selected_proofs/0/maps/0',
        computed_record_sha256=sha256(C.canonical(record).encode()).hexdigest(),
        forced_entries=[dict(synthetic_contact=True)],
        geometric_premise='SYNTHETIC_REQUIRED_PREMISE_NOT_ACCEPTANCE',
        outcome='SYNTHETIC_EXCLUSION_NOT_A_GEOMETRIC_THEOREM')
    packet = dict(source=stem,source_section='selected_proofs',
        proper_map_providers=[dict(high_index=7,maps=[mapping])])
    inner = dict(batch='proper-0',result=dict(packets=[packet]))
    return inner,spec('proper','proper-0',stem),pin_for(stem),{stem:dict(selected_proofs=[original])}


def whole_fixture():
    stem = 'SYNTHETIC_WHOLE_SOURCE_NOT_A_FILE'
    model = dict(model_key='a'*64,context_sha256='b'*64,current_total='1',
        current_corners=[dict(marks=[0,1],lower='1/2'),dict(marks=[1,2],lower='1/2')])
    proof = dict(high_index=7,actual_source_row_path='selected_proofs/0',
        all_current_funded_source_models=[model],automatic_exhaustiveness_claim=False,
        geometric_applicability_premise='SYNTHETIC_UNREGENERATED_VARIANT_PREMISE',
        full_completion_keys=[dict(key=key(),funded_model_keys=['a'*64])])
    inner = dict(batch='whole-0',result=dict(shards=8,shard=0,
        complete_reference_interface_count=564,packets=[dict(source=stem,
            current_funded_whole_model_providers=[proof])]))
    source_row = dict(index=7,method='SYNTHETIC_WHOLE_SOURCE_METHOD')
    return inner,spec('whole_schema_safe','whole-0',stem),pin_for(stem),{stem:dict(selected_proofs=[source_row])}


def affine9_fixture():
    stem = 'SYNTHETIC_AFFINE9_SOURCE_NOT_A_FILE'
    complete_key = key(affine=True)
    payload = dict(method='SYNTHETIC_NESTED_CONTACT',physical_entries='SYNTHETIC_NOT_A_GRAM')
    row = dict(index=7,maps=[deepcopy(complete_key['induced_map'])],
        proof=dict(method='SYNTHETIC_BC13_NESTED_METHOD',maps=[dict(proof=payload)]))
    other = deepcopy(row);other['index'] = 8
    provider = dict(kind='EXACT_PROPER_AFFINE_EXCLUSION',
        computed_record_sha256=sha256(C.canonical(payload).encode()).hexdigest(),
        premise='SYNTHETIC_REQUIRED_PREMISE_NOT_ACCEPTANCE')
    case = dict(source=stem,group='BC13',source_case_ordinal=0,high_index=7,
        current_funded_model_catalog=[],complete_registered_map_providers=[
            dict(key=complete_key,provider=provider,source_map_path='BC13/cases/0/maps/0')])
    inner = dict(batch='affine-9',result=dict(exact_registered_case_providers=[case]))
    return inner,spec('affine','affine-9',stem),pin_for(stem),{stem:dict(BC13=dict(cases=[row,other]))}


def guard_fixture():
    # This literal adapter stem resolves only in synthetic_sources, never on disk.
    stem = 'fresh_d13_three_rank15_context_transfer_20261018'
    complete_key = key();method = 'SYNTHETIC_GUARD_METHOD'
    source_row = dict(index=7,method=method,core_key=complete_key['core_key'],
        selected_nodes=deepcopy(complete_key['base_nodes']),
        selected_marks=deepcopy(complete_key['selected_marks']))
    unrelated = deepcopy(source_row);unrelated['index'] = 8
    first = dict(images=deepcopy(complete_key['induced_map']),marked_images=[0,9])
    second = deepcopy(first);second['images'][1:3] = [2,1]
    proof = dict(high_index=7,method=method,core_key=complete_key['core_key'],
        full_map_keys=[first,second])
    inner = dict(batch='guard-0',result=dict(guarded_providers=[proof]))
    worker_spec = dict(adapter='guard',method_rules=[dict(method=method,
        proof='SYNTHETIC_SOURCE_IMPLICATION_NOT_A_THEOREM')])
    return inner,worker_spec,pin_for(stem),{stem:dict(selected_proofs=[unrelated,source_row])}


class CurationControls(unittest.TestCase):
    def setUp(self):
        C.VISITS = 0

    def test_full_key_positive_and_finite_C_or_large_physical_node_refusals(self):
        valid = key()
        self.assertEqual(len(C.full_key(valid)),64)
        self.assertEqual(len(C.full_key(key(affine=True))),64)
        finite_c = deepcopy(valid);finite_c['target_kind'] = 'C'
        with self.assertRaisesRegex(ValueError,'full key native target/high schema'):
            C.full_key(finite_c)
        large = deepcopy(valid)
        large['base_nodes'][-1] = 64;large['selected_marks'][-1] = 64
        with self.assertRaisesRegex(ValueError,'complete induced map/ordered mark association'):
            C.full_key(large)

    def test_proper_exact_source_record_and_full_map_associations(self):
        fixture = proper_fixture()
        rows = exported(*fixture)
        self.assertEqual(len(rows),1)
        complete_key,binding,registration = rows[0]
        self.assertEqual(complete_key,key())
        self.assertEqual(binding['export_ordinal'],0)
        self.assertEqual(binding['origin'],['result','packets',0,'proper_map_providers',0,'maps',0])
        self.assertEqual(binding['source_record']['selector'],['selected_proofs',0,'maps',0])
        self.assertEqual(registration,fixture[1]['rows'][0])
        broken = deepcopy(fixture)
        broken[0]['result']['packets'][0]['proper_map_providers'][0]['maps'][0]['computed_record_sha256'] = '0'*64
        with self.assertRaisesRegex(ValueError,'proper complete record/contact association'):
            exported(*broken)
        wrong_map = deepcopy(fixture)
        wrong_map[0]['result']['packets'][0]['proper_map_providers'][0]['maps'][0]['key']['induced_map'][1:3] = [2,1]
        with self.assertRaisesRegex(ValueError,'proper exact original core/full-map association'):
            exported(*wrong_map)

    def test_whole_funded_model_subset_and_original_source_high_associations(self):
        fixture = whole_fixture()
        rows = exported(*fixture)
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0][1]['funded_model_keys'],['a'*64])
        self.assertEqual(rows[0][1]['source_record']['selector'],['selected_proofs',0])
        missing_model = deepcopy(fixture)
        missing_model[0]['result']['packets'][0]['current_funded_whole_model_providers'][0]['full_completion_keys'][0]['funded_model_keys'] = ['c'*64]
        with self.assertRaisesRegex(ValueError,'THIS full map funded-model subset'):
            exported(*missing_model)
        wrong_high = deepcopy(fixture)
        next(iter(wrong_high[3].values()))['selected_proofs'][0]['index'] = 8
        with self.assertRaisesRegex(ValueError,'whole actual source-row/high identity'):
            exported(*wrong_high)

    def test_affine9_nested_record_and_original_case_map_ordinals(self):
        fixture = affine9_fixture()
        rows = exported(*fixture)
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0][1]['source_map_selector'],['BC13','cases',0,'maps',0])
        self.assertEqual(rows[0][1]['source_record']['selector'],['BC13','cases',0])
        self.assertEqual(rows[0][1]['method'],'SYNTHETIC_BC13_NESTED_METHOD')
        self.assertEqual(rows[0][1]['provider_kind'],'EXACT_PROPER_AFFINE_EXCLUSION')
        bad_record = deepcopy(fixture)
        bad_record[0]['result']['exact_registered_case_providers'][0]['complete_registered_map_providers'][0]['provider']['computed_record_sha256'] = '0'*64
        with self.assertRaisesRegex(ValueError,'proper affine exact computed record/premise'):
            exported(*bad_record)
        wrong_ordinal = deepcopy(fixture)
        wrong_ordinal[0]['result']['exact_registered_case_providers'][0]['source_case_ordinal'] = 1
        with self.assertRaisesRegex(ValueError,'affine actual source case/high identity'):
            exported(*wrong_ordinal)
        wrong_map_path = deepcopy(fixture)
        wrong_map_path[0]['result']['exact_registered_case_providers'][0]['complete_registered_map_providers'][0]['source_map_path'] = 'BC13/cases/0/maps/1'
        with self.assertRaisesRegex(ValueError,'exact normalized affine source map ordinal'):
            exported(*wrong_map_path)

    def test_guard_source_proof_and_stable_map_origins(self):
        fixture = guard_fixture()
        rows = exported(*fixture)
        self.assertEqual(len(rows),2)
        for ordinal,(complete_key,binding,_) in enumerate(rows):
            self.assertEqual(binding['export_ordinal'],None)
            self.assertEqual(binding['source_record']['selector'],['selected_proofs',1])
            self.assertEqual(binding['guard_proof_origin'],0)
            self.assertEqual(binding['guard_map_origin'],ordinal)
            self.assertEqual(binding['origin'],['result','guarded_providers',0,'full_map_keys',ordinal])
            self.assertEqual(complete_key['induced_map'],
                fixture[0]['result']['guarded_providers'][0]['full_map_keys'][ordinal]['images'])
        self.assertNotEqual(rows[0][0],rows[1][0])
        wrong_source = deepcopy(fixture)
        next(iter(wrong_source[3].values()))['selected_proofs'][1]['method'] = 'OTHER_SYNTHETIC_METHOD'
        with self.assertRaisesRegex(ValueError,'same selected guard method/core'):
            exported(*wrong_source)
        ambiguous_source = deepcopy(fixture)
        doc = next(iter(ambiguous_source[3].values()))
        doc['selected_proofs'][0] = deepcopy(doc['selected_proofs'][1])
        with self.assertRaisesRegex(ValueError,'one exact recognized selected guard source proof'):
            exported(*ambiguous_source)
        wrong_mark_image = deepcopy(fixture)
        wrong_mark_image[0]['result']['guarded_providers'][0]['full_map_keys'][1]['marked_images'] = [9,0]
        with self.assertRaisesRegex(ValueError,'complete induced map/ordered mark association'):
            exported(*wrong_mark_image)


def main():
    global C
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent-lock',action='store_true',required=True)
    parser.parse_args()
    parent_lock()
    own_before = raw(OWN);source_before = raw(SOURCE)
    need(sha256(source_before).hexdigest() == SOURCE_SHA,'exact reviewed curator source')
    C = types.ModuleType('finite57_synthetic_provider_curator_definitions')
    C.__file__ = str(ROOT/SOURCE)
    exec(compile(source_before,C.__file__,'exec'),C.__dict__)
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(CurationControls)
    need(suite.countTestCases() == 5,'five focused synthetic association methods')
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    need(raw(OWN) == own_before and raw(SOURCE) == source_before,
         'control or curator source changed during synthetic controls')
    output = dict(status='PASS_PROVIDER_OCCURRENCE_CURATION_CONTROLS_ONLY' if result.wasSuccessful()
        else 'FAIL_PROVIDER_OCCURRENCE_CURATION_CONTROLS_ONLY',
        source_sha256=sha256(own_before).hexdigest(),sources={SOURCE:SOURCE_SHA},
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        output=stream.getvalue(),fixture_kind='SYNTHETIC_IN_MEMORY_SCHEMA_ASSOCIATIONS_ONLY',
        curator_main_called=False,scientific_data_read=False,actual_receipts_read=False,
        native_maps_enumerated=False,geometric_theorems_tested=False,
        provider_occurrence_acceptance_established=False,selected_provider_union_complete=False,
        current_terminal_bridge_complete=False,full_transitive_regeneration=False,
        geometric_formalization=False,global_finite57_verified=False,
        scope='Five positive/negative synthetic source-association methods for key schema, proper record/map, whole model/high, nested affine9 ordinals and guard proof/map origins. No actual provider coverage, finite predicates, geometric applicability or mathematical premise is verified.')
    encoded = json.dumps(output,sort_keys=True)
    need(len(encoded.encode())+1 <= OUTPUT_CAP,'complete synthetic control output cap')
    print(encoded)
    if not result.wasSuccessful():
        raise SystemExit(1)


if __name__ == '__main__':
    main()
