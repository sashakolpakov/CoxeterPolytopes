#!/usr/bin/env python3
"""Retain the55 controls and test the real HIGH stem/full-path interface.

Only pinned definitions and synthetic helper objects are used. The actual
HIGH source, application helpers and mathematical workers are not loaded.
Root must set FINAL_SHA and review the frozen source before execution.
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

ROOT=Path(__file__).resolve().parents[1]
PREVIOUS='scripts/test_finite57_current_rule_source_application_controls_20261025.py'
PREVIOUS_SHA='b9d1d4f56654f34fb8cb9155b9afc5b57da6283e73343e1c569a4e906114d7ea'
RELEASE_SHA='e51d172fb858b904999bcf8c580e98e0e8492a40ed28853abf9d675bac6c4d9a'
OLD_SHA='7f9801491d5b9f21ba669155f46d9d5469b4328880770f39fbec4ef013be4640'
FINAL='scripts/check_finite57_current_rule_conjunction_20261025.py'
FINAL_SHA='5276975ed175ae921386b0b075c2ba0d2b4740ba8e2bc187b438ca1d00d8d066'
F=None
V=None


def raw(name):
    with (ROOT/name).open('rb') as stream:data=stream.read(128*1024+1)
    if not 0<len(data)<=128*1024:raise ValueError('bounded synthetic control source')
    return data


def load(name,path,digest):
    data=raw(path)
    if sha256(data).hexdigest()!=digest:raise ValueError('exact source/control identity')
    obj=types.ModuleType(name);obj.__file__=str(ROOT/path)
    exec(compile(data,obj.__file__,'exec'),obj.__dict__)
    return obj


class HighProjectionControls(unittest.TestCase):
    def test_real_loader_stem_full_path_pin_and_projection_identity(self):
        stem='fresh_distance8_high_overlap_20261017'
        full='results/'+stem+'.json'
        cases=[dict(index=0,synthetic_no_realization=True)]
        cores={'SYNTHETIC_CORE':dict(synthetic_no_realization=True)}
        source_sha='d'*64
        reviews={F.LEGACY_APPLICATION:dict(source_sha256='a'*64,
                     path='results/SYNTHETIC_ACCEPTANCE_ONLY.json',sha256='c'*64),
                 F.FULL_D_APPLICATION:dict(source_sha256='b'*64,
                     path='results/SYNTHETIC_ACCEPTANCE_ONLY.json',sha256='c'*64)}
        route=types.SimpleNamespace(HIGH=stem,PINS={full:source_sha})
        state=dict(helper_high=full,pins=[],loads=[])
        saved=F.load,F.hash_file

        def pin(path,digest):
            allowed={full:source_sha,**{p:r['source_sha256'] for p,r in reviews.items()}}
            if allowed.get(path)!=digest:raise ValueError('synthetic exact byte-pin drift')
            state['pins'].append((path,digest))

        def helper(name,path,digest):
            if path not in reviews or digest!=reviews[path]['source_sha256']:
                raise ValueError('synthetic helper source association')
            state['loads'].append((name,path,digest))
            return types.SimpleNamespace(HIGH=state['helper_high'],
                PINS={state['helper_high']:source_sha},CACHE={})

        try:
            F.load=helper;F.hash_file=pin
            self.assertNotIn(stem,route.PINS)
            legacy,full_d=F.source_application_helpers(reviews,cases,cores,route)
            projection=legacy.CACHE[full]
            self.assertEqual(set(projection),{'isolation'})
            self.assertEqual(set(projection['isolation']),{'cases','cores'})
            self.assertIs(projection['isolation']['cases'],cases)
            self.assertIs(projection['isolation']['cores'],cores)
            self.assertEqual(full_d.CACHE,{})
            for path,obj in ((F.LEGACY_APPLICATION,legacy),(F.FULL_D_APPLICATION,full_d)):
                self.assertIs(obj.SOURCE_PINNER,pin)
                self.assertIs(obj.SOURCE_READER,F.application_source_document)
                self.assertEqual(obj.APPLICATION_SOURCE_SHA,reviews[path]['source_sha256'])
                self.assertEqual(obj.SOURCE_ACCEPTANCE,
                    dict(path=reviews[path]['path'],sha256=reviews[path]['sha256']))
            self.assertEqual([p for _,p,_ in state['loads']],
                             [F.LEGACY_APPLICATION,F.FULL_D_APPLICATION])
            self.assertIn((full,source_sha),state['pins'])
            self.assertNotIn(stem,[p for p,_ in state['pins']])

            wrong_pin=types.SimpleNamespace(HIGH=stem,PINS={full:'e'*64})
            with self.assertRaisesRegex(ValueError,'same original HIGH source'):
                F.source_application_helpers(reviews,cases,cores,wrong_pin)
            state['helper_high']='results/WRONG_SYNTHETIC_HIGH.json'
            with self.assertRaisesRegex(ValueError,'same original HIGH source'):
                F.source_application_helpers(reviews,cases,cores,route)
        finally:
            F.load,F.hash_file=saved


def main():
    global F,V
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent-lock',action='store_true',required=True);parser.parse_args()
    own=str(Path(__file__).resolve().relative_to(ROOT));before=raw(own)
    V=load('finite57_unchanged55_definitions',PREVIOUS,PREVIOUS_SHA)
    if V.PREVIOUS_SHA!=RELEASE_SHA:raise ValueError('unchanged47 source identity')
    P=V.load('finite57_unchanged47_definitions',V.PREVIOUS,RELEASE_SHA)
    if P.OLD_SHA!=OLD_SHA:raise ValueError('unchanged23 source identity')
    T=P.load('finite57_unchanged23_definitions',P.OLD,OLD_SHA)
    T.parent_lock()
    T.SOURCES={path:raw(path) for path in T.SOURCE_SHA}
    if {p:sha256(b).hexdigest() for p,b in T.SOURCES.items()}!=T.SOURCE_SHA:
        raise ValueError('unchanged predecessor fixture source closure')
    T.U=P.load('finite57_control_union435e',T.UNION,T.UNION_SHA)
    F=load('finite57_high_projection_conjunction',FINAL,FINAL_SHA);F.U=T.U
    P.F=F;P.T=T;P.C=P.load('finite57_release_codec_definitions',P.CODEC,P.CODEC_SHA)
    S=V.load('finite57_unchanged_s8_definitions',V.S8,V.S8_SHA)
    V.F=F;V.P=P;V.S=S
    if sha256(raw(F.LEGACY_REGISTRY)).hexdigest()!=F.LEGACY_REGISTRY_SHA:
        raise ValueError('unchanged historical fixture identity')
    classes=[T.UnionControls,P.FinalControls,P.CodecControls,P.ExpandedJSONControls,
        P.ConsumerAdapterControls,P.CurationTrustControls,V.RouteReuseControls,V.DirectS8Controls,
        V.SourceReaderControls,HighProjectionControls]
    suites=[unittest.defaultTestLoader.loadTestsFromTestCase(cls) for cls in classes]
    if [s.countTestCases() for s in suites]!=[23,15,4,1,2,2,4,3,1,1]:
        raise ValueError('exact55 retained definitions plus one HIGH integration regression')
    stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream,verbosity=2).run(unittest.TestSuite(suites))
    sources={PREVIOUS:PREVIOUS_SHA,V.PREVIOUS:RELEASE_SHA,FINAL:FINAL_SHA,
        V.S8:V.S8_SHA,P.OLD:OLD_SHA,P.CODEC:P.CODEC_SHA,
        F.LEGACY_REGISTRY:F.LEGACY_REGISTRY_SHA}
    if raw(own)!=before or any(sha256(raw(p)).hexdigest()!=h for p,h in sources.items()):
        raise ValueError('source changed during controls')
    output=dict(status=('PASS' if result.wasSuccessful() else 'FAIL')+'_CURRENT_RULE_HIGH_PROJECTION_CONTROLS_ONLY',
        source_sha256=sha256(before).hexdigest(),sources=sources,tests_run=result.testsRun,
        retained_source_application_cases=55,new_high_projection_cases=1,
        failures=len(result.failures),errors=len(result.errors),output=stream.getvalue(),
        scientific_data_read=False,geometric_theorems_tested=False,
        application_helpers_executed=False,actual_receipt_execution_fixtures=False,
        scope='Retains55 synthetic definitions against the new final reader and tests its actual stem/full-path HIGH adapter, exact source-pin refusal and original projection object identity. No actual HIGH source or mathematical application helper is loaded.')
    encoded=json.dumps(output,sort_keys=True)
    if len(encoded.encode())+1>500*1024:raise ValueError('complete control output cap')
    print(encoded)
    if not result.wasSuccessful():raise SystemExit(1)


if __name__=='__main__':main()
