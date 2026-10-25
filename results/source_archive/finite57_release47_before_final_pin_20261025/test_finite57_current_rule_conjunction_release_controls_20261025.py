#!/usr/bin/env python3
"""Release controls: retained38 tests and nine codec/JSON/consumer properties.

The38 predecessor/semantic tests remain unchanged. Codec checks use only
synthetic bytes and nested JSON; its pinned definitions are loaded in main.
No scientific data or mathematical parent preparation is read or executed.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
import binascii
from contextlib import contextmanager
from copy import deepcopy
from hashlib import sha256
import io
import json
from pathlib import Path
import types
import unittest

ROOT=Path(__file__).resolve().parents[1]
FINAL='scripts/check_finite57_current_rule_conjunction_20261025.py'
FINAL_SHA='fcf224c092a9cc76bec5c22dc6e9da7ef522acbd9b7ddfc80dab1a3bdfaf5f29'
OLD='scripts/test_finite57_actual_receipt_union_schema_safe_controls_20261025.py'
OLD_SHA='7f9801491d5b9f21ba669155f46d9d5469b4328880770f39fbec4ef013be4640'
CODEC='scripts/check_finite57_large_d_receivers_compressed_20261025.py'
CODEC_SHA='6df8d216a1e2d8473cdfaea0c591fb9d74e7fb8e6d6aadcf6df8fc30ba17e5a0'
F=None
T=None
C=None


def raw(name):
    with (ROOT/name).open('rb') as stream:
        data=stream.read(128*1024+1)
    if not 0<len(data)<=128*1024:
        raise ValueError('bounded synthetic control source')
    return data


def load(name,path,digest):
    data=raw(path)
    if sha256(data).hexdigest()!=digest:
        raise ValueError('exact synthetic control/helper source')
    obj=types.ModuleType(name);obj.__file__=str(ROOT/path)
    exec(compile(data,obj.__file__,'exec'),obj.__dict__)
    return obj


@contextmanager
def no_filesystem_hashes():
    old=F.hash_file
    # Actual ledger/stdout/source and semantic fields are still checked by
    # F.actual. This fixture does not test the filesystem hashing routine.
    F.hash_file=lambda path,digest:None
    try:
        yield
    finally:
        F.hash_file=old


def receiver(index):
    return dict(source='SYNTHETIC_NO_REALIZATION',index=index,context_sha256='a'*64,
        direct_current_receiver_condition=True,legacy_pointwise_bounds_used=False,
        total='1',corners=[dict(marks=[0,1],weight='1/2'),dict(marks=[1,2],weight='1/2')],
        final_collector_joins=[dict(synthetic_identity_only=True)])


def originals():
    domain=[['SYNTHETIC_NO_REALIZATION',str(i)] for i in range(216)]
    return {(F.DIRECT,f'original-{s}'):(dict(path=f'NO_FILE_{s}'),
        dict(result=dict(complete_domain_keys=deepcopy(domain),
                         rows=[receiver(i) for i in range(s,216,4)]))) for s in range(4)}


def quiet_fixtures():
    """Synthetic keys only; no geometry, actual receipt or filesystem evidence."""
    witnesses=[(a,b,c) for a in ('A','D') for b in ('A','D') for c in ('A','D')]
    rows=[dict(application_key=[index,*kinds])
          for index in F.QUIET_CASE_INDICES for kinds in witnesses]
    full={(F.UNIVERSAL,'final-quiet'):(dict(path='NO_FILE_FULL_QUIET'),
        dict(source_sha256=F.QUIET_PARENT_SHA,
             result=dict(quiet24_complete_at_declared_premises=True,theorem_applications=deepcopy(rows))))}
    micro={}
    for index in F.QUIET_CASE_INDICES:
        batch=f'quiet-case-{index:02d}'
        micro[F.QUIET_CASE,batch]=(dict(path='NO_FILE_'+batch),
            dict(base_universal_source_sha256=F.QUIET_PARENT_SHA,
                 consumed_inputs={'SYNTHETIC_SOURCE_ONLY':'a'*64},
                 result=dict(case_index=index,original24_case_domain=list(F.QUIET_CASE_INDICES),
                     case_complete_at_declared_premises=True,quiet24_complete_at_declared_premises=False,
                     theorem_applications=[deepcopy(row) for row in rows if row['application_key'][0]==index])))
    return full,micro


def percentage_sample(controllers):
    gib=1024**3
    return dict(free_disk_gib=3,available_memory_gib=4,
        disk_free_bytes=3*gib,disk_total_bytes=100*gib,
        available_memory_bytes=4*gib,physical_memory_bytes=20*gib,
        free_disk_percent=3,available_memory_percent=20,
        resource_policy=F.percentage_policy(controllers[F.POLICY_SOURCE]))


def percentage_fixture(outer_change=None):
    memory=T.Memory()
    def change(outer):
        outer['initial_resources']=percentage_sample(memory.controllers)
        if outer_change:
            outer_change(outer)
    entry=memory.add(T.U.ROUTES,'selected-00',T.selected(),outer_change=change)
    return memory,entry


class FinalControls(unittest.TestCase):
    def test_new_actual_reader_accepts_exact_three_line_fixture(self):
        memory,entry=percentage_fixture()
        entry['inner_status']=T.U.WORKERS[T.U.ROUTES][1]
        with memory,no_filesystem_hashes():
            F.U=T.U
            doc=F.actual(entry,{p['worker']:p for p in memory.reviews},memory.controllers)
        self.assertEqual(doc['batch'],'selected-00')

    def test_new_actual_reader_refuses_invented_registered_status(self):
        memory,entry=percentage_fixture()
        entry['inner_status']='INVENTED_PASS'
        with memory,no_filesystem_hashes(),self.assertRaisesRegex(ValueError,'exact actual inner status'):
            F.U=T.U
            F.actual(entry,{p['worker']:p for p in memory.reviews},memory.controllers)

    def test_new_actual_reader_refuses_boolean_thread_count(self):
        memory,entry=percentage_fixture(outer_change=lambda d:d.update(threads=True))
        entry['inner_status']=T.U.WORKERS[T.U.ROUTES][1]
        with memory,no_filesystem_hashes(),self.assertRaisesRegex(ValueError,'wrapper lock'):
            F.U=T.U
            F.actual(entry,{p['worker']:p for p in memory.reviews},memory.controllers)

    def test_exact_original_source_case_shards_accepted(self):
        F.U=T.U
        result=F.original_domain(originals())
        self.assertEqual(result['physical_receivers'],216)

    def test_original_equal_count_wrong_source_case_refused(self):
        F.U=T.U
        data=originals()
        data[F.DIRECT,'original-0'][1]['result']['rows'][0]['index']=1
        with self.assertRaisesRegex(ValueError,'exact source/case'):
            F.original_domain(data)

    def test_original_legacy_surplus_does_not_fund_current_shortfall(self):
        F.U=T.U
        data=originals();row=data[F.DIRECT,'original-0'][1]['result']['rows'][0]
        row['corners'][1]['weight']='1/4';row['total']='3/4'
        row['saved_legacy_total']='2'
        with self.assertRaisesRegex(ValueError,'complete current receiver sum'):
            F.original_domain(data)

    def test_physical_prefix_same_count_different_core_refused(self):
        p=dict(index=64,method='finite_leaf9_or_ideal_feature5',core_key='CORE_A',
               marked_component=[0,1],preferred_corner=[0,1])
        core=dict(gram=[['2','-1'],['-1','2']])
        digest=sha256(F.canonical([core['gram'],p['marked_component'],p['preferred_corner']]).encode()).hexdigest()
        row=dict(consumer_key=dict(high_index=64,method=p['method'],core_key='CORE_B',selected_marks=[0,1]),
                 physical_context_sha256=digest)
        with self.assertRaisesRegex(ValueError,'same high/method/core/ordered marks'):
            F.prefix_identity(p,core,row,types.SimpleNamespace(canonical=F.canonical))

    def test_physical_prefix_reversed_ordered_marks_refused(self):
        p=dict(index=64,method='finite_leaf9_or_ideal_feature5',core_key='CORE_A',
               marked_component=[0,1],preferred_corner=[0,1])
        row=dict(consumer_key=dict(high_index=64,method=p['method'],core_key='CORE_A',selected_marks=[1,0]),
                 physical_context_sha256='a'*64)
        with self.assertRaisesRegex(ValueError,'same high/method/core/ordered marks'):
            F.prefix_identity(p,dict(gram=[]),row,types.SimpleNamespace(canonical=F.canonical))

    def test_universal_wrong_full_role_refused_even_with_batch_present(self):
        F.U=T.U
        key=T.key();key.update(target_kind='D',target_order=12,affine=True,marked_image=[8,10])
        # Keep the full map/marks consistent; only its claimed actual theorem
        # application belongs to the other opposite-leaf pair.
        key['induced_map']=[8,0,1,2,3,4,5,6,7,10]
        theorem='finite57.current.affine_order12_D_opposite_leaves.whole_face.v1'
        application=dict(theorem=theorem,application_key=['D',12,8,11])
        payload=dict(result=dict(theorem_applications=[application]))
        units={(F.UNIVERSAL,'order12-d'):(dict(path='NO_FILE'),payload)}
        with self.assertRaisesRegex(ValueError,'THIS exact whole target role'):
            F.whole_application(key,None,units)

    def test_old_largeD_priority_requires_exact_theorem_and_actual_evidence(self):
        F.U=T.U
        key=T.key();key.update(target_kind='D',target_order=12,affine=True,marked_image=[8,10])
        key['induced_map']=[8,0,1,2,3,4,5,6,7,10]
        theorem='finite57.current.affine_order12_D_opposite_leaves.whole_face.v1'
        app=dict(theorem=theorem,application_key=['D',12,8,10],
                 status='BOUND_CONDITIONAL_WHOLE_FACE_APPLICATION')
        payload=dict(result=dict(theorem_applications=[app]))
        units={(F.UNIVERSAL,'order12-d'):(dict(path='NO_FILE_SYNTHETIC'),payload)}
        class LiteralFixture:
            def classical_graph(self,kind,order,affine=False):return {}
            def _distances(self,graph,mark):return {8:0,10:9}
            def affine_weight(self,kind,order,*marks):return 0
        mapping=dict(key=key,conditional_whole_face_priority=theorem)
        bound=F.old_largeD_priority(mapping,LiteralFixture(),units)
        self.assertEqual(bound['application']['selected_full_key'],key)
        self.assertEqual(bound['application']['actual_theorem_interface_receipts'][0]['receipt'],
                         'NO_FILE_SYNTHETIC')
        wrong=dict(mapping,conditional_whole_face_priority='finite57.current.BC11_zero.whole_face.v1')
        with self.assertRaisesRegex(ValueError,'exact OLD largeD full-map conditional priority theorem'):
            F.old_largeD_priority(wrong,LiteralFixture(),units)
        with self.assertRaisesRegex(ValueError,'missing actual unit'):
            F.old_largeD_priority(mapping,LiteralFixture(),{})

    def test_selected_premise_cannot_drift_or_disappear_in_final(self):
        pin=dict(path='NO_FILE_SYNTHETIC',sha256='a'*64)
        row=dict(id='source-premise',kind='source_geometric_implication',
                 statement='SYNTHETIC conditional statement',acceptance_scope='CONTROL_ONLY',
                 source=deepcopy(pin),review=deepcopy(pin),root_acceptance=deepcopy(pin),
                 runtime_regenerated=False)
        extra=dict(deepcopy(row),id='final-extra-premise')
        accepted={row['id']:deepcopy(row),extra['id']:extra}
        self.assertEqual(F.retain_selected_assumptions(accepted,[deepcopy(row)]),[row['id']])
        drift=deepcopy(row);drift['source']['sha256']='b'*64
        with self.assertRaisesRegex(ValueError,'selected/final accepted premise conflict'):
            F.retain_selected_assumptions(accepted,[drift])
        with self.assertRaisesRegex(ValueError,'selected premise missing from final'):
            F.retain_selected_assumptions({extra['id']:extra},[row])

    def test_quiet_full_or_all24_cases_preserve_exact192_keys_and_provenance(self):
        full,micro=quiet_fixtures()
        old_rows,old_evidence=F.quiet_applications(full)
        new_rows,new_evidence=F.quiet_applications(micro)
        self.assertEqual(old_rows,new_rows)
        self.assertEqual(len(new_rows),192)
        self.assertEqual(old_evidence['mode'],'ONE_ACTUAL_FULL_PARENT_RECEIPT')
        self.assertEqual(new_evidence['mode'],'TWENTY_FOUR_ACTUAL_CASE_RECEIPTS')
        self.assertEqual(len(old_evidence['actual_receipts']),1)
        self.assertEqual(len(new_evidence['actual_receipts']),24)
        self.assertTrue(all(r['worker']==F.QUIET_CASE for r in new_evidence['actual_receipts']))

    def test_quiet_incomplete_mixed_misassociated_or_drifted_cases_refuse(self):
        full,micro=quiet_fixtures()
        first=(F.QUIET_CASE,f'quiet-case-{F.QUIET_CASE_INDICES[0]:02d}')
        missing=deepcopy(micro);del missing[first]
        wrong_case=deepcopy(micro);wrong_case[first][1]['result']['case_index']=F.QUIET_CASE_INDICES[1]
        wrong_key=deepcopy(micro);wrong_key[first][1]['result']['theorem_applications'][0]['application_key'][0]=F.QUIET_CASE_INDICES[1]
        drift=deepcopy(micro);drift[first][1]['consumed_inputs']['SYNTHETIC_SOURCE_ONLY']='b'*64
        unknown=deepcopy(micro);unknown[F.QUIET_CASE,'quiet-case-00']=deepcopy(micro[first])
        for name,data,message in [('missing',missing,'missing actual unit'),
                ('mixed',{**full,**micro},'do not mix'),
                ('wrong case',wrong_case,'THIS exact quiet case'),
                ('wrong witness key',wrong_key,'every ordered witness key'),
                ('source drift',drift,'same exact rule/source/premise input pins'),
                ('unknown batch',unknown,'only the exact24')]:
            with self.subTest(name=name),self.assertRaisesRegex(ValueError,message):
                F.quiet_applications(data)

    def test_percentage_policy_uses_raw_bytes_and_rejects_inconsistent_metadata(self):
        F.U=T.U
        resources=percentage_sample({F.POLICY_SOURCE:'a'*64})
        F.percentage_resources(resources,'a'*64)
        self.assertLess(resources['free_disk_gib'],40)
        self.assertLess(resources['available_memory_gib'],6)
        below=deepcopy(resources);below['disk_free_bytes']=2*1024**3-1
        below['free_disk_percent']=2
        display=deepcopy(resources);display['free_disk_percent']=4
        threshold=deepcopy(resources);threshold['resource_policy']['abort_memory_percent']=5
        boolean=deepcopy(resources);boolean['disk_free_bytes']=True
        wrong_source=deepcopy(resources);wrong_source['resource_policy']['source_sha256']='b'*64
        for name,value,message in [('rounded below admission',below,'exact byte ratios'),
                ('wrong displayed percentage',display,'six-digit percentages'),
                ('changed abort threshold',threshold,'exact current percentage policy'),
                ('boolean byte count',boolean,'integer resource bytes'),
                ('wrong policy source',wrong_source,'exact current percentage policy')]:
            with self.subTest(name=name),self.assertRaisesRegex(ValueError,message):
                F.percentage_resources(value,'a'*64)

    def test_only_exact_historical_nine_keep_their_recorded_controller_maps(self):
        F.U=T.U
        records=F.historical_receipts()
        self.assertEqual(len(records),9)
        d12=[];later=[]
        for row in records:
            matched=F.match_historical_receipt(deepcopy(row['entry']),records)
            self.assertEqual(matched['controller_source_sha256'],row['controller_source_sha256'])
            monitor=matched['controller_source_sha256']['scripts/monitor_owned_proof_tree_20261024.py']
            (d12 if row['entry']['worker']==F.D12 else later).append(monitor)
        self.assertEqual(d12,['517c0c78347c0aa81cd16489afcf7f11b7bce064490919b9cf5ff94c88e59cda']*2)
        self.assertEqual(later,['b07f84fbcc5b35e484696f0515b8e233ea975f13c2da6ff75c0d3b54b99ed5db']*7)
        unknown=deepcopy(records[0]['entry']);unknown['sha256']='0'*64
        self.assertIsNone(F.match_historical_receipt(unknown,records))
        with self.assertRaisesRegex(ValueError,'unknown policy-less receipt'):
            F.receipt_resources(unknown,dict(free_disk_gib=80,available_memory_gib=16),
                                {F.POLICY_SOURCE:'a'*64,T.U.RUNNER:'b'*64})


def codec_fixture():
    """Representation data only; its complete fields must survive decoding."""
    return dict(scope='SYNTHETIC_CODEC_ONLY_NO_GEOMETRIC_REALIZATION',
        source_keys=[['source',0],['source',1]],
        theorem_applications=[dict(application_key=['case','A','D'],
            physical_context_sha256='a'*64,
            receivers=[dict(marks=[3,7],current='1/2'),dict(marks=[2,8],current='1/2')],
            assumptions=[dict(kind='synthetic_whole_factor_premise',runtime_verified=False)])],
        optional=None,flags=[False,True],nested=dict(empty=[],rational='17/32'))


class CodecControls(unittest.TestCase):
    def test_complete_nested_result_byte_and_object_round_trip(self):
        value=codec_fixture()
        expected=json.dumps(value,sort_keys=True,separators=(',',':')).encode('ascii')
        packet=C.encode_result(value)
        decoded=C.decode_result(packet)
        self.assertEqual(decoded,expected)
        self.assertEqual(json.loads(decoded),value)
        self.assertEqual(packet['decoded_bytes'],len(expected))
        self.assertEqual(packet['decoded_sha256'],sha256(expected).hexdigest())
        self.assertEqual(packet['compressed_bytes'],
                         len(C.base64.b64decode(packet['data'],validate=True)))

    def test_digest_exact_lengths_types_and_decoded_bound_refused(self):
        packet=C.encode_result(codec_fixture())
        digest=packet['decoded_sha256']
        wrong_digest=('0' if digest[0]!='0' else '1')+digest[1:]
        changes=[('wrong digest','decoded_sha256',wrong_digest,'complete bounded decoded'),
            ('decoded length smaller','decoded_bytes',packet['decoded_bytes']-1,'complete bounded decoded'),
            ('decoded length larger','decoded_bytes',packet['decoded_bytes']+1,'complete bounded decoded'),
            ('compressed length smaller','compressed_bytes',packet['compressed_bytes']-1,'exact compressed length'),
            ('compressed length larger','compressed_bytes',packet['compressed_bytes']+1,'exact compressed length'),
            ('boolean decoded length','decoded_bytes',True,'bounded compressed-result lengths'),
            ('boolean compressed length','compressed_bytes',True,'bounded compressed-result lengths'),
            ('advertised decoded bytes exceed8MiB','decoded_bytes',8*1024**2+1,'bounded compressed-result lengths')]
        for name,field,value,message in changes:
            changed=dict(packet);changed[field]=value
            with self.subTest(name=name),self.assertRaisesRegex(ValueError,message):
                C.decode_result(changed)

    def test_malformed_base64_refused(self):
        packet=C.encode_result(codec_fixture())
        packet['data']='!'+packet['data'][1:]
        with self.assertRaises((ValueError,binascii.Error)):
            C.decode_result(packet)

    def test_truncated_or_trailing_stream_refused_without_partial_result(self):
        packet=C.encode_result(codec_fixture())
        compressed=C.base64.b64decode(packet['data'],validate=True)
        variants=[('truncated checksum',compressed[:-1]),
                  ('trailing byte',compressed+b'\x00'),
                  ('second complete stream',compressed+compressed)]
        for name,data in variants:
            changed=dict(packet,compressed_bytes=len(data),
                         data=C.base64.b64encode(data).decode('ascii'))
            with self.subTest(name=name),self.assertRaisesRegex(ValueError,'complete bounded decoded'):
                C.decode_result(changed)


class ExpandedJSONControls(unittest.TestCase):
    def test_expanded_size_exact_numbers_duplicate_nonfinite_and_depth_contract(self):
        F.U=T.U
        payload=b'{"payload":"'+b'['*(2*1024**2)+b'"}'
        self.assertGreater(len(payload),2*1024**2)
        self.assertLess(len(payload),8*1024**2)
        decoded=F.expanded_loads(payload)
        self.assertEqual(set(decoded),{'payload'})
        self.assertEqual(decoded['payload'],'['*(2*1024**2))
        number=F.expanded_loads(b'{"value":0.123456789012345678901}')['value']
        self.assertIsInstance(number,F.Decimal)
        self.assertEqual(str(number),'0.123456789012345678901')
        nested=F.expanded_loads(b'['*48+b'0'+b']'*48)
        for _ in range(48):
            self.assertIsInstance(nested,list)
            self.assertEqual(len(nested),1)
            nested=nested[0]
        self.assertEqual(nested,0)
        for name,data in [('duplicate fields',b'{"same":1,"same":2}'),
                ('NaN',b'{"value":NaN}'),('positive infinity',b'{"value":Infinity}'),
                ('negative infinity',b'{"value":-Infinity}'),
                ('depth49',b'['*49+b'0'+b']'*49)]:
            with self.subTest(name=name),self.assertRaises(ValueError):
                F.expanded_loads(data)


@contextmanager
def isolated_codec_consumer():
    saved=(F.codec,F.CODEC,F.LARGE_DECODED,F.EXPANDED_BYTES,F.U)
    F.codec=lambda:C  # Only this synthetic fixture bypasses codec source loading.
    F.LARGE_DECODED={}
    F.EXPANDED_BYTES=0
    F.U=T.U
    try:
        yield
    finally:
        F.codec,F.CODEC,F.LARGE_DECODED,F.EXPANDED_BYTES,F.U=saved


class ConsumerAdapterControls(unittest.TestCase):
    def test_large_d_complete_result_preserves_receipt_and_refuses_wrong_mode_or_parent(self):
        value=codec_fixture()
        entry=dict(path='NO_FILE_SYNTHETIC_ENCODED_LARGE_D',sha256='a'*64,
                   worker=F.LOSSLESS,batch='large-d-0')
        inner=dict(batch='large-d-0',source_sha256=CODEC_SHA,
            parent_source_sha256=F.QUIET_PARENT_SHA,
            consumed_inputs={F.UNIVERSAL:F.QUIET_PARENT_SHA},result_encoding=C.encode_result(value))
        units={(F.LOSSLESS,'large-d-0'):(entry,inner)}
        with isolated_codec_consumer():
            actual,decoded=F.large_d_receipt(units,0)
            self.assertIs(actual,entry)
            self.assertEqual(actual['worker'],F.LOSSLESS)
            self.assertEqual(decoded['result'],value)
            self.assertEqual({k:decoded[k] for k in inner},inner)
            charged=F.EXPANDED_BYTES
            self.assertGreater(charged,0)
            cached_entry,cached=F.large_d_receipt(units,0)
            self.assertIs(cached_entry,entry)
            self.assertIs(cached,decoded)
            self.assertEqual(F.EXPANDED_BYTES,charged)
            old_entry=dict(path='NO_FILE_SYNTHETIC_ORIGINAL_LARGE_D',sha256='b'*64,
                           worker=F.UNIVERSAL,batch='large-d-0')
            old_inner=dict(result=deepcopy(value))
            legacy={(F.UNIVERSAL,'large-d-0'):(old_entry,old_inner)}
            self.assertEqual(F.large_d_receipt(legacy,0),(old_entry,old_inner))
            wrong_parent=deepcopy(units)
            wrong_parent[F.LOSSLESS,'large-d-0'][1]['parent_source_sha256']='b'*64
            wrong_consumed=deepcopy(units)
            wrong_consumed[F.LOSSLESS,'large-d-0'][1]['consumed_inputs'][F.UNIVERSAL]='b'*64
            for name,data,shard,message in [
                    ('mixed mode',{**units,**legacy},0,'do not mix original and encoded'),
                    ('missing requested shard',units,1,'missing actual complete large-D shard'),
                    ('wrong parent identity',wrong_parent,0,'exact unchanged parent'),
                    ('wrong consumed parent',wrong_consumed,0,'exact unchanged parent'),
                    ('boolean shard',units,True,'exact large-D shard'),
                    ('shard outside0..3',units,4,'exact large-D shard')]:
                with self.subTest(name=name),self.assertRaisesRegex(ValueError,message):
                    F.large_d_receipt(data,shard)

    def test_affine9_and_original_affine_alias_leave_frozen_registry_unchanged(self):
        original='scripts/index_finite57_selected_affine_completion_providers_20261025.py'
        registry=T.U.WORKERS
        before=deepcopy(registry)
        with isolated_codec_consumer():
            self.assertEqual(F.provider_adapter(F.AFFINE9),'affine')
            self.assertEqual(F.provider_adapter(original),'affine')
            self.assertNotIn(F.AFFINE9,registry)
            self.assertIsNone(F.provider_adapter('NO_SUCH_SYNTHETIC_WORKER'))
        self.assertIs(T.U.WORKERS,registry)
        self.assertEqual(registry,before)


@contextmanager
def isolated_curation_consumer():
    saved_module,saved_visits=F.U,T.U.VISITS
    F.U=T.U
    try:
        yield
    finally:
        T.U.VISITS=saved_visits
        F.U=saved_module


def synthetic_curation():
    """Identity fixtures only; these36 names are not registered worker batches."""
    worker='scripts/index_finite57_all_proper_finite_completion_providers_20261025.py'
    entries=[dict(worker=worker,batch=f'SYNTHETIC_BATCH_ONLY_{i:02d}',
        path=f'NO_FILE_PROVIDER_{i:02d}',sha256=f'{i+1:064x}') for i in range(36)]
    ordinary=[dict(id=f'SYNTHETIC_ORDINARY_{i:02d}',kind='complete_finite_data',
        statement='SYNTHETIC retained premise only',runtime_regenerated=False) for i in range(10)]
    entry=dict(worker=F.CURATOR,batch='catalogue',path='NO_FILE_CURATOR',sha256='a'*64)
    result=dict(schema='finite57-provider-occurrence-bindings-v1',
        complete_chosen_provider_set=True,complete_eligible_occurrence_associations=True,
        unknown_or_unresolved_eligible_occurrences=[],input_snapshot=deepcopy(F.CURATOR_INPUT),
        per_source_numerical_premises_regenerated=False,selected_provider_union_complete=False,
        full_transitive_regeneration=False,geometric_formalization=False,global_finite57_verified=False,
        provider_entries=deepcopy(entries),ordinary_assumptions=deepcopy(ordinary),
        lookup={'SYNTHETIC_KEY_ONLY':{}},contracts={'SYNTHETIC_CONTRACT_ONLY':{}})
    doc=dict(input_snapshot=deepcopy(F.CURATOR_INPUT),result=result)
    units={(p['worker'],p['batch']):(p,{}) for p in entries}
    units[F.CURATOR,'catalogue']=(entry,doc)
    accepted={p['id']:deepcopy(p) for p in ordinary}
    accepted['selected-provider-geometry']=dict(actual_curator_receipt=deepcopy(entry),
        provider_catalogue=deepcopy(F.CURATOR_INPUT))
    return units,accepted,entries[0]


class CurationTrustControls(unittest.TestCase):
    def test_curation_binds_exact36_inventory_receipt_catalogue_and_ordinary_premises(self):
        with isolated_curation_consumer():
            units,accepted,first=synthetic_curation()
            curator=(F.CURATOR,'catalogue')
            self.assertIs(F.provider_curation(units,accepted),units[curator][1]['result'])
            missing=dict(units);del missing[curator]
            wrong_entry=deepcopy(units)
            wrong_entry[first['worker'],first['batch']][0]['sha256']='f'*64
            missing_provider=deepcopy(units);del missing_provider[first['worker'],first['batch']]
            drift=deepcopy(units)
            drift[curator][1]['result']['input_snapshot']['sha256']='f'*64
            promoted=deepcopy(units)
            promoted[curator][1]['result']['global_finite57_verified']=True
            ordinary=deepcopy(accepted)
            ordinary['SYNTHETIC_ORDINARY_00']['statement']='SYNTHETIC DRIFT'
            wrong_receipt=deepcopy(accepted)
            wrong_receipt['selected-provider-geometry']['actual_curator_receipt']['sha256']='f'*64
            for name,data,premises,message in [
                    ('missing curator',missing,accepted,'missing actual complete provider-occurrence'),
                    ('same-count provider identity drift',wrong_entry,accepted,'every and only actual chosen provider'),
                    ('missing provider',missing_provider,accepted,'every and only actual chosen provider'),
                    ('fixed catalogue drift',drift,accepted,'complete fixed source-occurrence catalogue'),
                    ('global claim promotion',promoted,accepted,'cannot claim selected coverage or geometry'),
                    ('ordinary premise drift',units,ordinary,'identical ordinary premises'),
                    ('curator receipt identity drift',units,wrong_receipt,'exact actual curation and fixed catalogue')]:
                with self.subTest(name=name),self.assertRaisesRegex(ValueError,message):
                    F.provider_curation(data,premises)

    def test_curated_lookup_binds_full_key_actual_receipt_origin_and_contract(self):
        with isolated_curation_consumer():
            key=T.key()
            worker='scripts/index_finite57_all_proper_finite_completion_providers_20261025.py'
            entry=dict(worker=worker,batch='SYNTHETIC_BATCH_ONLY',
                       path='NO_FILE_EXACT_PROVIDER',sha256='a'*64)
            keyhash=sha256(T.U.full_key(key).encode()).hexdigest()
            identity=worker+'|'+entry['batch']+'|'+keyhash
            binding=dict(origin=dict(source='SYNTHETIC_NO_REALIZATION',selector=['cases',0,'maps',0]),
                contract_id='SYNTHETIC_PROPER_CONTRACT',
                source_record=dict(source_sha256='b'*64,record_sha256='c'*64))
            contract=dict(adapter='proper',proof_home_ids=['SYNTHETIC_HOME'],
                proof_homes=[dict(id='SYNTHETIC_HOME',path='NO_FILE_PROOF_HOME')],
                numerical_and_geometric_premises_regenerated=False)
            row=dict(worker=worker,batch=entry['batch'],receipt=entry['path'],
                receipt_sha256=entry['sha256'],key_sha256=keyhash,key=deepcopy(key),
                bindings=[deepcopy(binding)])
            catalogue=dict(lookup={identity:row},contracts={binding['contract_id']:contract})
            bound=F.curated_occurrences(catalogue,entry,key)
            self.assertIs(bound,row['bindings'])
            self.assertEqual(bound,[binding])
            self.assertEqual(catalogue['contracts'][bound[0]['contract_id']],contract)
            changed_key=deepcopy(catalogue);changed_key['lookup'][identity]['key']=T.key(reverse=True)
            unknown=deepcopy(catalogue);unknown['lookup'][identity]['bindings'][0]['contract_id']='UNKNOWN'
            adapter=deepcopy(catalogue);adapter['contracts'][binding['contract_id']]['adapter']='affine'
            duplicate=deepcopy(catalogue);duplicate['lookup'][identity]['bindings'].append(deepcopy(binding))
            missing_hash=deepcopy(catalogue)
            missing_hash['lookup'][identity]['bindings'][0]['source_record']['record_sha256']=''
            for name,data,actual,selected,message in [
                    ('different full map',catalogue,entry,T.key(reverse=True),'uncurated exact full provider key'),
                    ('row full-key drift',changed_key,entry,key,'complete physical nine-field key'),
                    ('wrong receipt path',catalogue,dict(entry,path='NO_FILE_OTHER_PROVIDER'),key,'same actual provider receipt'),
                    ('wrong receipt hash',catalogue,dict(entry,sha256='f'*64),key,'same actual provider receipt'),
                    ('unknown proof binding',unknown,entry,key,'exact reviewed occurrence theorem'),
                    ('wrong adapter',adapter,entry,key,'exact reviewed occurrence theorem'),
                    ('duplicate source origin',duplicate,entry,key,'distinct stable source occurrences'),
                    ('missing source-record identity',missing_hash,entry,key,'exact reviewed occurrence theorem')]:
                with self.subTest(name=name),self.assertRaisesRegex(ValueError,message):
                    F.curated_occurrences(data,actual,selected)


def main():
    global F,T,C
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--parent-lock',action='store_true',required=True);p.parse_args()
    own=str(Path(__file__).resolve().relative_to(ROOT));before=raw(own)
    T=load('finite57_incorporated_old23_controls',OLD,OLD_SHA)
    T.parent_lock()
    T.SOURCES={path:raw(path) for path in T.SOURCE_SHA}
    if {p:sha256(b).hexdigest() for p,b in T.SOURCES.items()}!=T.SOURCE_SHA:
        raise ValueError('unchanged predecessor fixture source closure')
    T.U=load('finite57_control_union435e',T.UNION,T.UNION_SHA)
    F=load('finite57_control_final_conjunction',FINAL,FINAL_SHA);F.U=T.U
    C=load('finite57_release_codec_definitions',CODEC,CODEC_SHA)
    historical=raw(F.LEGACY_REGISTRY)
    if sha256(historical).hexdigest()!=F.LEGACY_REGISTRY_SHA:
        raise ValueError('immutable historical metadata fixture')
    predecessor=unittest.defaultTestLoader.loadTestsFromTestCase(T.UnionControls)
    new=unittest.defaultTestLoader.loadTestsFromTestCase(FinalControls)
    codec=unittest.defaultTestLoader.loadTestsFromTestCase(CodecControls)
    expanded=unittest.defaultTestLoader.loadTestsFromTestCase(ExpandedJSONControls)
    consumer=unittest.defaultTestLoader.loadTestsFromTestCase(ConsumerAdapterControls)
    curation=unittest.defaultTestLoader.loadTestsFromTestCase(CurationTrustControls)
    if (predecessor.countTestCases()!=23 or new.countTestCases()!=15 or
            codec.countTestCases()!=4 or expanded.countTestCases()!=1 or
            consumer.countTestCases()!=2 or curation.countTestCases()!=2):
        raise ValueError('exact retained38 plus nine codec/JSON/consumer/curation properties')
    suite=unittest.TestSuite([predecessor,new,codec,expanded,consumer,curation])
    stream=io.StringIO();result=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    if (raw(own)!=before or sha256(raw(FINAL)).hexdigest()!=FINAL_SHA or
            sha256(raw(OLD)).hexdigest()!=OLD_SHA or raw(F.LEGACY_REGISTRY)!=historical or
            sha256(raw(CODEC)).hexdigest()!=CODEC_SHA):
        raise ValueError('synthetic control source changed')
    output=dict(status='PASS_CURRENT_RULE_CONJUNCTION_CONTROLS_ONLY' if result.wasSuccessful()
        else 'FAIL_CURRENT_RULE_CONJUNCTION_CONTROLS_ONLY',
        source_sha256=sha256(before).hexdigest(),
        sources={FINAL:FINAL_SHA,OLD:OLD_SHA,F.LEGACY_REGISTRY:F.LEGACY_REGISTRY_SHA,CODEC:CODEC_SHA},
        tests_run=result.testsRun,incorporated_predecessor_cases=23,new_semantic_cases=15,
        codec_properties=4,expanded_json_properties=1,consumer_properties=2,
        curation_trust_properties=2,codec_fixture_only=True,
        failures=len(result.failures),errors=len(result.errors),output=stream.getvalue(),
        actual_execution_fixtures=False,scientific_data_read=False,
        filesystem_hashing_and_monitoring_tested=False,geometric_theorems_tested=False,
        historical_metadata_fixture_read=True,
        scope='Retains all38 predecessor/semantic controls and adds four codec, one expanded-JSON, two receipt/adapter and two curation-trust properties on synthetic data only. No scientific input, actual large-D producer or curator execution, provider execution or geometry theorem is tested.')
    encoded=json.dumps(output,sort_keys=True)
    if len(encoded.encode())+1>500*1024:
        raise ValueError('complete control output cap')
    print(encoded)
    if not result.wasSuccessful():
        raise SystemExit(1)


if __name__=='__main__':
    main()

