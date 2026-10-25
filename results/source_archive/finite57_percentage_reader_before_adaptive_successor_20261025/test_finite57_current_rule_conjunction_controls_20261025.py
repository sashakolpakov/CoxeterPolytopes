#!/usr/bin/env python3
"""Synthetic controls for final conditional conjunction; source preparation only.

The23 never-executed predecessor fixtures are incorporated once. New controls
exercise exact receiver identities, physical prefix identity and universal role
association. No scientific data or mathematical producer is imported here.
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
from hashlib import sha256
import io
import json
from pathlib import Path
import types
import unittest

ROOT=Path(__file__).resolve().parents[1]
FINAL='scripts/check_finite57_current_rule_conjunction_20261025.py'
FINAL_SHA='a5cfc91010933a85b34c4875a9d19e2db8262b93c084154042926b35b44c90ab'
OLD='scripts/test_finite57_actual_receipt_union_schema_safe_controls_20261025.py'
OLD_SHA='7f9801491d5b9f21ba669155f46d9d5469b4328880770f39fbec4ef013be4640'
F=None
T=None


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


def main():
    global F,T
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
    historical=raw(F.LEGACY_REGISTRY)
    if sha256(historical).hexdigest()!=F.LEGACY_REGISTRY_SHA:
        raise ValueError('immutable historical metadata fixture')
    predecessor=unittest.defaultTestLoader.loadTestsFromTestCase(T.UnionControls)
    new=unittest.defaultTestLoader.loadTestsFromTestCase(FinalControls)
    if predecessor.countTestCases()!=23 or new.countTestCases()!=15:
        raise ValueError('exact incorporated23 plus15 new semantic controls')
    suite=unittest.TestSuite([predecessor,new])
    stream=io.StringIO();result=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    if (raw(own)!=before or sha256(raw(FINAL)).hexdigest()!=FINAL_SHA or
            sha256(raw(OLD)).hexdigest()!=OLD_SHA or raw(F.LEGACY_REGISTRY)!=historical):
        raise ValueError('synthetic control source changed')
    output=dict(status='PASS_CURRENT_RULE_CONJUNCTION_CONTROLS_ONLY' if result.wasSuccessful()
        else 'FAIL_CURRENT_RULE_CONJUNCTION_CONTROLS_ONLY',
        source_sha256=sha256(before).hexdigest(),
        sources={FINAL:FINAL_SHA,OLD:OLD_SHA,F.LEGACY_REGISTRY:F.LEGACY_REGISTRY_SHA},
        tests_run=result.testsRun,incorporated_predecessor_cases=23,new_semantic_cases=15,
        failures=len(result.failures),errors=len(result.errors),output=stream.getvalue(),
        actual_execution_fixtures=False,scientific_data_read=False,
        filesystem_hashing_and_monitoring_tested=False,geometric_theorems_tested=False,
        historical_metadata_fixture_read=True,
        scope='One successor suite incorporates23 never-executed memory controls plus15 new exact-identity/source-case/current-sum/universal-role/premise-retention/quiet-case-partition/resource-policy controls. The historical fixture contains provenance metadata only; synthetic fixtures are not provider executions or geometry evidence.')
    encoded=json.dumps(output,sort_keys=True)
    if len(encoded.encode())+1>500*1024:
        raise ValueError('complete control output cap')
    print(encoded)
    if not result.wasSuccessful():
        raise SystemExit(1)


if __name__=='__main__':
    main()

