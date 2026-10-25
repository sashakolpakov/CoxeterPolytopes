#!/usr/bin/env python3
"""Source-only successor: retain47 controls and check new consumer boundaries.

Only synthetic fixtures and pinned definitions are loaded. No mathematical
worker preparation, scientific data, geometry or native-map search runs.
The root must freeze FINAL_SHA and review this source before execution.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
import io
import json
from pathlib import Path
import types
import unittest

ROOT=Path(__file__).resolve().parents[1]
PREVIOUS='scripts/test_finite57_current_rule_conjunction_release_controls_20261025.py'
PREVIOUS_SHA='e51d172fb858b904999bcf8c580e98e0e8492a40ed28853abf9d675bac6c4d9a'
FINAL='scripts/check_finite57_current_rule_conjunction_20261025.py'
FINAL_SHA='0000000000000000000000000000000000000000000000000000000000000000'
S8='scripts/check_finite57_current_s8_applications_20261025.py'
S8_SHA='e9a8074d79b6861f04678d10fd9abe20259a36d7c7e04ed21f815836b81c5755'
F=None
P=None
S=None


def raw(name):
    with (ROOT/name).open('rb') as stream:data=stream.read(128*1024+1)
    if not 0<len(data)<=128*1024:raise ValueError('bounded control source')
    return data


def load(name,path,digest):
    data=raw(path)
    if sha256(data).hexdigest()!=digest:raise ValueError('exact control/helper source')
    obj=types.ModuleType(name);obj.__file__=str(ROOT/path)
    exec(compile(data,obj.__file__,'exec'),obj.__dict__)
    return obj


def route_fixture(shard=31):
    batch=f'selected-{shard:02d}'
    result=dict(shard=shard,shards=64,complete_authoritative_partition=dict(F.U.PARTITION),
        all_induced_maps_visited=0,direct_whole_factor_unit_maps=0,
        direct_map_key_sha256=sha256(b'').hexdigest(),prefix_obligations=[],
        unbound_current_completion_obligations=[])
    entry=dict(path='NO_FILE_SYNTHETIC_ACTUAL_RECEIPT_ONLY',sha256='a'*64)
    doc=dict(batch=batch,source_sha256=F.ROUTE_SHA,collector_sha256='b'*64,result=result)
    return batch,{(F.ROUTE,batch):(entry,doc)},result


class RouteReuseControls(unittest.TestCase):
    def test_nonzero_shard_reuses_complete_original_without_search(self):
        batch,units,result=route_fixture()
        got,origin=F.reused_route_result(units,batch,'b'*64,dict(F.U.PARTITION))
        self.assertIs(got,result)
        self.assertIs(origin['receipt'],units[F.ROUTE,batch][0])
        self.assertEqual(origin['kind'],'REUSED_ACTUAL_INDEX_WITHOUT_REENUMERATION')

    def test_missing_receipt_refuses_instead_of_enumerating(self):
        with self.assertRaisesRegex(ValueError,'missing actual unit'):
            F.reused_route_result({},'selected-31','b'*64,dict(F.U.PARTITION))

    def test_same_count_source_batch_collector_shard_and_partition_drift_refused(self):
        batch,units,_=route_fixture()
        mutations=[('source',lambda d:d.__setitem__('source_sha256','f'*64)),
            ('batch',lambda d:d.__setitem__('batch','selected-32')),
            ('collector',lambda d:d.__setitem__('collector_sha256','f'*64)),
            ('shard',lambda d:d['result'].__setitem__('shard',32)),
            ('boolean shard',lambda d:d['result'].__setitem__('shard',True)),
            ('shards',lambda d:d['result'].__setitem__('shards',63)),
            ('partition',lambda d:d['result']['complete_authoritative_partition'].__setitem__('radical',1330))]
        for name,mutation in mutations:
            changed=deepcopy(units);mutation(changed[F.ROUTE,batch][1])
            with self.subTest(name=name),self.assertRaises(ValueError):
                F.reused_route_result(changed,batch,'b'*64,dict(F.U.PARTITION))

    def test_out_of_domain_batch_and_caller_partition_refused(self):
        batch,units,_=route_fixture()
        for value in ('selected-64','selected-1',True):
            with self.subTest(batch=value),self.assertRaises(ValueError):
                F.reused_route_result(units,value,'b'*64,dict(F.U.PARTITION))
        wrong=dict(F.U.PARTITION);wrong['radical']+=1
        with self.assertRaisesRegex(ValueError,'reused shard/partition'):
            F.reused_route_result(units,batch,'b'*64,wrong)


class SyntheticCycle:
    """Small graph fixture; no scientific rule or actual polytope assertion."""
    @staticmethod
    def classical_graph(kind,order,affine=False):
        if kind!='A' or not affine or order!=17:raise ValueError('synthetic cycle only')
        return [{(i-1)%order:3,(i+1)%order:3} for i in range(order)]

    @staticmethod
    def _distances(graph,start):
        distances={start:0};queue=[start]
        for i in queue:
            for j in graph[i]:
                if j not in distances:distances[j]=distances[i]+1;queue.append(j)
        return distances

    @staticmethod
    def affine_weight(kind,order,a,b):return Fraction(1,2)


def cycle_key():
    # Displayed path has distance9; its WHOLE cycle has a shorter arc8.
    return dict(high_index=31,core_key='SYNTHETIC_NO_REALIZATION',base_nodes=list(range(10)),
        selected_marks=[0,9],target_kind='A',target_order=17,affine=True,
        induced_map=list(range(10)),marked_image=[0,9])


class DirectS8Controls(unittest.TestCase):
    def test_affine_cycle_shorter_arc_uses_whole_face_not_intrinsic_unit(self):
        key=cycle_key();got=S.apply_s8(key,SyntheticCycle)
        self.assertEqual(got['application_key']['marked_distance'],8)
        self.assertEqual(got['application_key']['selected_full_key'],key)
        self.assertEqual(got['application_key']['feature_depths'],[None,None])
        self.assertEqual(got['conclusion']['lower'],'1')
        self.assertFalse(got['intrinsic_corner_unit'])
        self.assertFalse(got['theorem_execution'])

    def test_full_map_mark_mismatch_and_noninjection_refused(self):
        badmark=cycle_key();badmark['marked_image']=[0,8]
        noninjective=cycle_key();noninjective['induced_map'][1]=0
        for key in (badmark,noninjective):
            with self.assertRaises(ValueError):S.apply_s8(key,SyntheticCycle)

    def test_long_marked_path_cannot_be_called_observed_s8(self):
        class LongFixture(SyntheticCycle):
            @staticmethod
            def _distances(graph,start):return {i:abs(i-start) for i in range(len(graph))}
        with self.assertRaisesRegex(ValueError,'actual marked white-interior path<=8'):
            S.apply_s8(cycle_key(),LongFixture)


class SourceReaderControls(unittest.TestCase):
    def test_cached_source_rechecks_identity_and_keeps_strict_complete_parser(self):
        name='results/SYNTHETIC_SOURCE_NO_REALIZATION.json'
        data=b'{"exact":0.123456789012345678901}'
        state=dict(data=data,reads=0,pins=0)
        saved=(F.U,F.hash_file,F.SOURCE_DOCUMENTS)
        def pin(path,digest):
            state['pins']+=1
            if path!=name or sha256(state['data']).hexdigest()!=digest:
                raise ValueError('synthetic byte pin drift')
        def read(path,cap,remember):
            self.assertEqual(path,name);self.assertEqual(cap,8*1024**2)
            self.assertFalse(remember);state['reads']+=1
            return state['data']
        try:
            F.U=types.SimpleNamespace(raw=read,unique=saved[0].unique)
            F.hash_file=pin;F.SOURCE_DOCUMENTS={}
            digest=sha256(data).hexdigest()
            first=F.application_source_document(name,digest)
            self.assertEqual(str(first['exact']),'0.123456789012345678901')
            self.assertIsInstance(first['exact'],F.Decimal)
            self.assertIs(F.application_source_document(name,digest),first)
            self.assertEqual(state['reads'],1);self.assertEqual(state['pins'],4)
            state['data']=b'{"changed":true}'
            with self.assertRaisesRegex(ValueError,'synthetic byte pin drift'):
                F.application_source_document(name,digest)
            self.assertEqual(state['reads'],1)
            for invalid in (b'{"x":1,"x":2}',b'{"x":NaN}',b'['*49+b'0'+b']'*49):
                state['data']=invalid
                with self.assertRaises(ValueError):
                    F.application_source_document(name,sha256(invalid).hexdigest())
        finally:
            F.U,F.hash_file,F.SOURCE_DOCUMENTS=saved


def main():
    global F,P,S
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent-lock',action='store_true',required=True);parser.parse_args()
    own=str(Path(__file__).resolve().relative_to(ROOT));before=raw(own)
    P=load('finite57_unchanged_release47_definitions',PREVIOUS,PREVIOUS_SHA)
    T=P.load('finite57_unchanged_old23_definitions',P.OLD,P.OLD_SHA)
    T.parent_lock()
    T.SOURCES={path:raw(path) for path in T.SOURCE_SHA}
    if {p:sha256(b).hexdigest() for p,b in T.SOURCES.items()}!=T.SOURCE_SHA:
        raise ValueError('unchanged predecessor fixture source closure')
    T.U=P.load('finite57_control_union435e',T.UNION,T.UNION_SHA)
    F=load('finite57_new_source_application_conjunction',FINAL,FINAL_SHA);F.U=T.U
    P.F=F;P.T=T;P.C=P.load('finite57_release_codec_definitions',P.CODEC,P.CODEC_SHA)
    S=load('finite57_unchanged_s8_definitions',S8,S8_SHA)
    historical=raw(F.LEGACY_REGISTRY)
    if sha256(historical).hexdigest()!=F.LEGACY_REGISTRY_SHA:raise ValueError('historical fixture')
    classes=[T.UnionControls,P.FinalControls,P.CodecControls,P.ExpandedJSONControls,
        P.ConsumerAdapterControls,P.CurationTrustControls,RouteReuseControls,DirectS8Controls,
        SourceReaderControls]
    suites=[unittest.defaultTestLoader.loadTestsFromTestCase(cls) for cls in classes]
    counts=[s.countTestCases() for s in suites]
    if counts!=[23,15,4,1,2,2,4,3,1]:raise ValueError('exact47 retained plus8 new boundary controls')
    stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream,verbosity=2).run(unittest.TestSuite(suites))
    sources={PREVIOUS:PREVIOUS_SHA,FINAL:FINAL_SHA,S8:S8_SHA,P.OLD:P.OLD_SHA,
        P.CODEC:P.CODEC_SHA,F.LEGACY_REGISTRY:F.LEGACY_REGISTRY_SHA}
    if raw(own)!=before or any(sha256(raw(p)).hexdigest()!=h for p,h in sources.items()):
        raise ValueError('source changed during controls')
    output=dict(status=('PASS' if result.wasSuccessful() else 'FAIL')+'_CURRENT_RULE_SOURCE_APPLICATION_CONTROLS_ONLY',
        source_sha256=sha256(before).hexdigest(),sources=sources,tests_run=result.testsRun,
        retained_release_cases=47,new_route_reuse_cases=4,new_direct_s8_cases=3,new_source_reader_cases=1,
        failures=len(result.failures),errors=len(result.errors),output=stream.getvalue(),
        scientific_data_read=False,geometric_theorems_tested=False,actual_receipt_execution_fixtures=False,
        scope='Synthetic consumer-boundary tests. Reuses unchanged47 definitions against the new final reader; tests all-shard actual-index reuse and affine-cycle path/full-map/S8 boundaries. No native-map search or geometry is tested.')
    encoded=json.dumps(output,sort_keys=True)
    if len(encoded.encode())+1>500*1024:raise ValueError('complete control output cap')
    print(encoded)
    if not result.wasSuccessful():raise SystemExit(1)


if __name__=='__main__':main()
