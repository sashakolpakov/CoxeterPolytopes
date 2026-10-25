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
FINAL_SHA='26ad0b1455add72e719c0d1df2e11bd0b9d7ccd2679979e293ac4dd2f99c57dd'
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


class FinalControls(unittest.TestCase):
    def test_new_actual_reader_accepts_exact_three_line_fixture(self):
        memory,entry=T.actual_fixture()
        entry['inner_status']=T.U.WORKERS[T.U.ROUTES][1]
        with memory,no_filesystem_hashes():
            F.U=T.U
            doc=F.actual(entry,{p['worker']:p for p in memory.reviews},memory.controllers)
        self.assertEqual(doc['batch'],'selected-00')

    def test_new_actual_reader_refuses_invented_registered_status(self):
        memory,entry=T.actual_fixture()
        entry['inner_status']='INVENTED_PASS'
        with memory,no_filesystem_hashes(),self.assertRaisesRegex(ValueError,'exact actual inner status'):
            F.U=T.U
            F.actual(entry,{p['worker']:p for p in memory.reviews},memory.controllers)

    def test_new_actual_reader_refuses_boolean_thread_count(self):
        memory,entry=T.actual_fixture(outer_change=lambda d:d.update(threads=True))
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
    predecessor=unittest.defaultTestLoader.loadTestsFromTestCase(T.UnionControls)
    new=unittest.defaultTestLoader.loadTestsFromTestCase(FinalControls)
    if predecessor.countTestCases()!=23 or new.countTestCases()!=9:
        raise ValueError('exact incorporated23 plus9 new semantic controls')
    suite=unittest.TestSuite([predecessor,new])
    stream=io.StringIO();result=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    if raw(own)!=before or sha256(raw(FINAL)).hexdigest()!=FINAL_SHA or sha256(raw(OLD)).hexdigest()!=OLD_SHA:
        raise ValueError('synthetic control source changed')
    output=dict(status='PASS_CURRENT_RULE_CONJUNCTION_CONTROLS_ONLY' if result.wasSuccessful()
        else 'FAIL_CURRENT_RULE_CONJUNCTION_CONTROLS_ONLY',
        source_sha256=sha256(before).hexdigest(),sources={FINAL:FINAL_SHA,OLD:OLD_SHA},
        tests_run=result.testsRun,incorporated_predecessor_cases=23,new_semantic_cases=9,
        failures=len(result.failures),errors=len(result.errors),output=stream.getvalue(),
        actual_execution_fixtures=False,scientific_data_read=False,
        filesystem_hashing_and_monitoring_tested=False,geometric_theorems_tested=False,
        scope='One successor suite incorporates23 never-executed memory controls plus9 new exact-identity/source-case/current-sum/universal-role controls. Synthetic fixtures are not provider executions or geometry evidence.')
    encoded=json.dumps(output,sort_keys=True)
    if len(encoded.encode())+1>500*1024:
        raise ValueError('complete control output cap')
    print(encoded)
    if not result.wasSuccessful():
        raise SystemExit(1)


if __name__=='__main__':
    main()

