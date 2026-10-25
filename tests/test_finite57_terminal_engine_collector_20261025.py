#!/usr/bin/env python3
"""SOURCE ONLY: focused bounded synthetic terminal-engine/collector regressions.

No scientific certificate is loaded and no native worker is imported. Root
alone may run this source after review, under the ordinary proof-batch gates.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE'in os.environ:
    raise SystemExit('Optimization environment is forbidden')
import argparse
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import io
import json
from pathlib import Path
import types
import unittest

ROOT=Path(__file__).resolve().parents[1]
SOURCES={
 'scripts/finite57_literal_weights_audit_20261024.py':'f6eba7af77af50b821e9e99c589388c1352f9ba5a99052aaeee64393f5ada1e0',
 'scripts/finite57_terminal_literal_engine_small_affine_20261025.py':'39dd55d2b2abef3db045e646edcc95bc1c8363f483d9f853d8c1d4e1b65d0e1f',
 'scripts/check_finite57_terminal_collectors_coverage_safe_20261025.py':'68f7a441f4e9701090e30d9a049430a23a3e5cae660e86d1f248a63a26a87dee',
}
MAX_OUTPUT=500*1024
W=E=None


def need(ok,message):
    if not ok:raise ValueError(message)


def raw(path):
    with (ROOT/path).open('rb')as stream:data=stream.read(256*1024+1)
    need(0<len(data)<=256*1024,'bounded synthetic-test source read')
    return data


def module(path,name):
    data=raw(path);need(sha256(data).hexdigest()==SOURCES[path],'synthetic-test exact source pin')
    obj=types.ModuleType(name);obj.__file__=str(ROOT/path)
    exec(compile(data,str(ROOT/path),'exec'),obj.__dict__)
    return obj


def gram(n,edges=()):
    out=[['2'if i==j else'0'for j in range(n)]for i in range(n)]
    for i,j,*magnitude in edges:
        q=str(magnitude[0])if magnitude else'1'
        out[i][j]=out[j][i]='-('+q+')'
    return out


def classical(kind,n,extras=()):
    graph=W.classical_graph(kind,n)
    return gram(n+sum(extras),[(i,j,'sqrt(2)'if lab==4 else'1')for i,row in enumerate(graph)for j,lab in row.items()if i<j])


def coverage_fixture(unfunded=False):
    """Exactly564 occurrence keys in22 tiny fake source packets.

    Literal values are stubbed here: this fixture tests coverage control and
    context identity, independently of the separate real literal-engine tests.
    """
    stems=[f'synthetic_{i:02d}'for i in range(22)]
    rows={s:[]for s in stems};expected={s:[]for s in stems};b13=[];d13=[]
    corners=[dict(marks=m,components=[])for m in([0,1],[1,2],[2,3],[0,3])]
    g=gram(4,[(0,1),(1,2),(2,3)])
    for i in range(141):
        stem=stems[i%22];offset=len(rows[stem]);name='loose_outputs'if unfunded and i==0 else'corners'
        rows[stem].append(dict(gram=g,white=[],**{name:deepcopy(corners)}))
        for j,c in enumerate(corners):
            path=f'selected_proofs/{offset}/{name}/{j}'
            expected[stem].append(dict(path=path,marks=c['marks']))
            entry=dict(source=stem,path=path,marks=c['marks'],images=[])
            b13.append(deepcopy(entry));d13.append(deepcopy(entry))
    reference=[dict(source=s,corners=expected[s])for s in stems]
    return stems,rows,reference,b13,d13


class TerminalRegression(unittest.TestCase):
    def setUp(self):
        self.c=module('scripts/check_finite57_terminal_collectors_coverage_safe_20261025.py','terminal_regression_collector')
        self.c.W=W;self.c.E=E;self.c.WITNESS_POOL={};self.c.FUEL=[0]

    def test_positive_gauge_invariance(self):
        g=gram(3,[(0,1),(1,2)]);scale=[2,3,5]
        gauged=[[str(E.number(x).rational()*scale[i]*scale[j])for j,x in enumerate(row)]for i,row in enumerate(g)]
        a=E.corner(g,[1],[0,2],W,[0]);b=E.corner(gauged,[1],[0,2],W,[0])
        self.assertEqual(a,b);self.assertEqual(a['weight'],'1')

    def test_wrong_edge_and_nonedge_refused(self):
        target=W.classical_graph('A',4)
        wrong=E.labels(E.matrix(gram(4,[(0,1),(1,2),(2,3),(0,3)])),list(range(4)))
        self.assertEqual(E.embeddings(wrong,target,[0,3],[0]),[])
        wrong=E.labels(E.matrix(gram(4,[(0,1),(1,2,'sqrt(2)'),(2,3)])),list(range(4)))
        self.assertEqual(E.embeddings(wrong,target,[0,3],[0]),[])

    def test_complete_canonical_images(self):
        g=E.matrix(gram(4,[(0,1),(1,2),(2,3)]));lab=E.labels(g,list(range(4)));target=W.classical_graph('A',4)
        self.assertEqual(E.embeddings(lab,target,[0],[0]),[(0,),(3,)])
        self.assertEqual(E.embeddings(lab,target,[0,1],[0]),[(0,1),(2,3)])

    def test_shared_fuel_refuses_not_partial(self):
        lab=E.labels(E.matrix(gram(4,[(0,1),(1,2),(2,3)])),list(range(4)))
        with self.assertRaisesRegex(ValueError,'fuel exhausted'):
            E.embeddings(lab,W.classical_graph('A',4),[0],[E.MAX_FUEL])

    def test_unsealed_positive_witness_is_ignored(self):
        g=classical('A',13,(1,));white=list(range(1,12))+[13]
        ordinary=E.corner(g,white,[0,12],W,[0])
        witnessed=E.corner(g,white,[0,12],W,[0],sealed_witnesses=[[13]])
        self.assertEqual(ordinary['weight'],'0');self.assertEqual(ordinary['ignored_unsealed_components'],[[13]])
        self.assertEqual(witnessed['weight'],'1');self.assertEqual(witnessed['sealed_witness_components'],[[13]])

    def test_partial_D13_guard_stays_conservative(self):
        g=classical('D',13,(2,));g[13][14]=g[14][13]='-1'
        result=E.corner(g,[i for i in range(15)if i not in(0,11)],[0,11],W,[0])
        self.assertEqual(result['weight'],'1/2');self.assertTrue(result['guard_possible'])
        self.assertEqual(result['ignored_unsealed_components'],[[13,14]])
        with self.assertRaisesRegex(ValueError,'complete augmentation'):
            E.corner(g,[i for i in range(15)if i not in(0,11)],[0,11],W,[0],complete=True)

    def test_incomplete_or_marked_witness_refused(self):
        g=classical('A',13,(2,));g[13][14]=g[14][13]='-1';white=list(range(1,12))+[13,14]
        for false_witness in([[13]],[[0,1]]):
            with self.assertRaisesRegex(ValueError,'distinct unmarked displayed component'):
                E.corner(g,white,[0,12],W,[0],sealed_witnesses=false_witness)

    def test_duplicate_corner_refused(self):
        g=gram(3,[(0,1),(1,2)]);corners=[dict(marks=[0,2]),dict(marks=[2,0])]
        with self.assertRaisesRegex(ValueError,'corner counted twice'):
            self.c.receipt('synthetic',0,g,[1],corners,'synthetic distinct-corner test')

    def test_exceptional_affine_E8_recognized(self):
        # Independent explicit arms1,2,5 around node0; order9.
        edges=[(0,1),(0,2),(2,3),(0,4),(4,5),(5,6),(6,7),(7,8)]
        got=E.corner(gram(9,edges),list(range(1,8)),[0,8],W,[0])
        self.assertEqual(got['kind'],'ideal');self.assertEqual(got['forest'],[['E',9]])
        self.assertEqual(got['weight'],'1')

    def test_exact_golden_label_and_unsafe_AST(self):
        g=E.matrix(gram(2,[(0,1,'(1+sqrt(5))/2')]))
        self.assertEqual(E.labels(g,[0,1]),[[0,5],[5,0]])
        for expression in("__import__('os')",'sqrt(2,3)','1/0'):
            with self.assertRaises((ValueError,ZeroDivisionError)):
                E.number(expression)

    def test_extended_model_context_is_preserved(self):
        old=gram(3,[(0,1),(1,2)]);new=gram(4,[(0,1),(1,2),(2,3)])
        corners=[dict(marks=m)for m in([0,1],[1,2],[2,3],[0,3])]
        fixture=dict(gram=old,white=[],extended_gram=new,extended_white=[1],after=corners)
        found=list(self.c.complete_models(fixture,'root'))
        self.assertEqual(len(found),1);path,g,w,cs,gp,wp=found[0]
        self.assertEqual((path,g,w,gp,wp),('root/after',new,[1],'root/extended_gram','root/extended_white'))

    def install_coverage(self,unfunded):
        c=self.c;stems,rows,reference,b13,d13=coverage_fixture(unfunded)
        c.LATER={s:('selected_proofs',)for s in stems}
        def field(stem,path):
            if stem==c.COMPOSE and path==('later_finite_output_interfaces',):return reference
            if stem==c.B13 and path==('current_interfaces','later_corners'):return b13
            if stem==c.FINAL and path==('current_interfaces','later_corners'):return d13
            raise ValueError('unexpected synthetic field route')
        c.field=field;c.document=lambda stem,*paths:dict(selected_proofs=rows[stem])
        c.value=lambda g,w,m,*a,**k:dict(weight='1',marks=m,kind='finite')
        c.current_collector_join=lambda g,w,m,b,d:dict(synthetic_join=True)
        return c

    def test_every_required_path_has_funded_model(self):
        c=self.install_coverage(False);out=c.later_batch(0)
        self.assertEqual(len(out['packets']),6)
        for packet in out['packets']:
            self.assertEqual(set(packet['interfaces']),set(packet['required_path_to_funded_model']))
            paid={r['model_key']for r in packet['receipts']}
            self.assertTrue(all(r['model_key']in paid for r in packet['required_path_to_funded_model'].values()))

    def test_required_path_outside_full_receipt_refused(self):
        c=self.install_coverage(True)
        with self.assertRaisesRegex(ValueError,'lacks a full current-funded polygon list'):
            c.later_batch(0)

    def test_new_corners_are_a_complete_receipt_list(self):
        g=gram(4,[(0,1),(1,2),(2,3)])
        cs=[dict(marks=m)for m in([0,1],[1,2],[2,3],[0,3])]
        found=list(self.c.complete_models(dict(old_seal=dict(gram=g,white=[]),new_corners=cs),'root'))
        self.assertEqual(len(found),1);self.assertEqual(found[0][0],'root/new_corners')
        self.assertEqual(found[0][1:4],(g,[],cs))

    def test_pending_provider_cannot_emit_plain_named_PASS(self):
        self.assertEqual(self.c.pending_provider_occurrences(dict(a=[dict(status='UNBOUND_G4_G6_PROVIDER')],b=dict(status='PASS'))),1)
        self.assertEqual(self.c.pending_provider_occurrences(dict(a=[dict(status='PASS')])),0)


def main():
    global W,E
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent-lock',action='store_true',required=True)
    parser.parse_args();relative=str(Path(__file__).resolve().relative_to(ROOT));before=raw(relative)
    W=module('scripts/finite57_literal_weights_audit_20261024.py','terminal_regression_literal')
    E=module('scripts/finite57_terminal_literal_engine_small_affine_20261025.py','terminal_regression_engine')
    gate=module('scripts/check_finite57_terminal_collectors_coverage_safe_20261025.py','terminal_regression_gate');gate.parent_lock()
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(TerminalRegression)
    need(suite.countTestCases()==15,'exact intended15 synthetic regressions')
    stream=io.StringIO();result=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    after={p:sha256(raw(p)).hexdigest()for p in SOURCES}
    need(after==SOURCES and raw(relative)==before,'synthetic test/source changed during batch')
    payload=json.dumps(dict(status='PASS_TERMINAL_SYNTHETIC_REGRESSIONS'if result.wasSuccessful()else'FAIL_TERMINAL_SYNTHETIC_REGRESSIONS',
        source_sha256=sha256(before).hexdigest(),sources_before=SOURCES,sources_after=after,
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),output=stream.getvalue(),
        scope='15 tiny synthetic representation/coverage tests only; zero scientific certificate/worker/Lean/TeX inputs; not a terminal-payment or global finite57 verification'),sort_keys=True)
    need(len(payload.encode())<=MAX_OUTPUT,'synthetic receipt cap')
    print(payload)
    if not result.wasSuccessful():raise SystemExit(1)


if __name__=='__main__':main()
