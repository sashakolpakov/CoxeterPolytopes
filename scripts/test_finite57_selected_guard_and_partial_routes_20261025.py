#!/usr/bin/env python3
"""SOURCE ONLY: tiny partial-injection/G5 schema-control regressions.

Synthetic inverse entries are fixtures for association checks, not claimed
Gram inverses or geometric examples. No scientific certificate is loaded.
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
 'scripts/index_finite57_current_terminal_routes_coverage_safe_20261025.py':'65c7a7e018800fa7bd259a984ec24249656410569aee578f668a6e04f5fde5d1',
 'scripts/check_finite57_selected_guard_bindings_contact_safe_20261025.py':'480d613a3e44f2038add91a09f68ad0625a45ddf6b3f36cce82304499183ab5d',
}
W=E=C=R=G=None


def need(ok,message):
    if not ok:raise ValueError(message)


def raw(relative):
    with(ROOT/relative).open('rb')as stream:data=stream.read(128*1024+1)
    need(0<len(data)<=128*1024,'tiny regression bounded source read')
    return data


def module(relative,name):
    data=raw(relative);need(sha256(data).hexdigest()==SOURCES[relative],'tiny regression source identity')
    obj=types.ModuleType(name);obj.__file__=str(ROOT/relative)
    exec(compile(data,str(ROOT/relative),'exec'),obj.__dict__)
    return obj


def fixture(proper=False):
    n=16;base=list(range(12 if proper else 13));marks=[0,11]
    target=W.classical_graph('A'if proper else'D',len(base))
    gram=[['2'if i==j else'0'for j in range(n)]for i in range(n)]
    for i,row in enumerate(target):
        for j,label in row.items():gram[i][j]='-1'
    inverse=[['2'if i==j else'0'for j in range(n)]for i in range(n)]
    inverse[13][13]='4/3';inverse[10][10]='3';inverse[9][9]='3'
    high=dict(index=1,core_key='synthetic-only',method='lorentz_residual',marked_component=base,
              preferred_corner=marks,left=[0,[0,14]],right=[0,[11,15]],right_copy_map=list(range(n)),
              possible_contacts=[10]if proper else[13])
    core=dict(gram=gram,inverse=inverse,radical=[])
    p=dict(index=1,core_key='synthetic-only',gram=gram,white=[i for i in range(n)if i not in(0,11,14,15)],
           selected_nodes=base,selected_marks=marks)
    maps=[dict(images=list(base),marked_images=marks)]
    if proper:
        p.update(method='genuine_new_fork_leaf_projection_exclusion',
                 maps=[dict(images=list(base),marked_images=marks,missing=[12],forced_contact=10,projection_lower='3')])
    else:
        p.update(method='one_contact_inverse_rank_bound',old_small_factor=dict(nodes=[13],type='A1'),
                 maximum_small_rank=4,remaining_inverse_diagonals=[[13,'4/3']],possible_contacts=[13],
                 contact=13,inverse_diagonal='4/3',maximum_added_rank=3,first_forbidden_rank_energy='16/15')
    return p,high,core,maps


class PartialAndGuard(unittest.TestCase):
    def setUp(self):C.FUEL=[0]

    def path(self):return [[0,3,0],[3,0,3],[0,3,0]]

    def test_partial_three_path_has_four_full_maps(self):
        maps=list(R.embeddings(self.path(),W.classical_graph('A',4),[0],E))
        self.assertEqual(sorted(maps),[[0,1,2],[1,2,3],[2,1,0],[3,2,1]])
        self.assertEqual(len(maps),4)
        self.assertEqual(len({tuple(sorted((m[0],m[2])))for m in maps}),2)

    def test_required_nonedge_mutation_rejected(self):
        pattern=self.path();pattern[0][2]=pattern[2][0]=3
        self.assertEqual(list(R.embeddings(pattern,W.classical_graph('A',4),[0],E)),[])

    def test_required_label_mutation_rejected(self):
        pattern=self.path();pattern[0][1]=pattern[1][0]=4
        self.assertEqual(list(R.embeddings(pattern,W.classical_graph('A',4),[0],E)),[])

    def test_fuel_after_first_map_refuses_complete_list(self):
        with self.assertRaisesRegex(ValueError,'fuel exhausted'):
            list(R.embeddings(self.path(),W.classical_graph('A',4),[E.MAX_FUEL-4],E))

    def test_actual_old_small_factor_schema(self):
        p,h,core,maps=fixture()
        self.assertNotIn('whole_small_factor',p)
        out=G.proof(C,p,h,core,maps)
        self.assertEqual(out['maximum_whole_small_rank'],4)
        self.assertEqual(out['whole_small_nodes'],[13])

    def test_wrong_one_contact_node_rejected(self):
        p,h,core,maps=fixture();p['contact']=12
        with self.assertRaisesRegex(ValueError,'unique actual small-factor contact'):
            G.proof(C,p,h,core,maps)

    def test_arbitrary_possible_subset_rejected(self):
        p,h,core,maps=fixture();p['possible_contacts']=[]
        with self.assertRaisesRegex(ValueError,'exact remaining feasible-contact subset'):
            G.proof(C,p,h,core,maps)

    def test_wrong_one_contact_inverse_rejected(self):
        p,h,core,maps=fixture();p['inverse_diagonal']='3/2'
        with self.assertRaisesRegex(ValueError,'rho joined to actual core inverse'):
            G.proof(C,p,h,core,maps)

    def test_exact_proper_missing_leaf_fixture(self):
        p,h,core,maps=fixture(True)
        out=G.proof(C,p,h,core,maps)
        self.assertEqual(out['method'],'genuine_new_fork_leaf_projection_exclusion')

    def test_wrong_missing_target_node_rejected(self):
        p,h,core,maps=fixture(True);p['maps'][0]['missing']=[13]
        with self.assertRaisesRegex(ValueError,'exact required missing target node'):
            G.proof(C,p,h,core,maps)

    def test_right_number_wrong_map_contact_rejected(self):
        p,h,core,maps=fixture(True);p['maps'][0]['forced_contact']=9;h['possible_contacts']=[9,10]
        self.assertEqual(core['inverse'][9][9],p['maps'][0]['projection_lower'])
        with self.assertRaisesRegex(ValueError,'actual full-map target neighbor'):
            G.proof(C,p,h,core,maps)

    def test_right_contact_wrong_inverse_entry_rejected(self):
        p,h,core,maps=fixture(True);p['maps'][0]['projection_lower']='4'
        with self.assertRaisesRegex(ValueError,'actual core inverse entry'):
            G.proof(C,p,h,core,maps)


def main():
    global W,E,C,R,G
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent-lock',action='store_true',required=True)
    parser.parse_args();relative=str(Path(__file__).resolve().relative_to(ROOT));before=raw(relative)
    W=module('scripts/finite57_literal_weights_audit_20261024.py','partial_guard_literal')
    E=module('scripts/finite57_terminal_literal_engine_small_affine_20261025.py','partial_guard_engine')
    C=module('scripts/check_finite57_terminal_collectors_coverage_safe_20261025.py','partial_guard_collector');C.W=W;C.E=E;C.parent_lock()
    R=module('scripts/index_finite57_current_terminal_routes_coverage_safe_20261025.py','partial_guard_routes')
    G=module('scripts/check_finite57_selected_guard_bindings_contact_safe_20261025.py','partial_guard_guard')
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(PartialAndGuard)
    need(suite.countTestCases()==12,'exact intended12 tiny regressions')
    stream=io.StringIO();result=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    after={p:sha256(raw(p)).hexdigest()for p in SOURCES}
    need(after==SOURCES and raw(relative)==before,'tiny regression source changed')
    payload=json.dumps(dict(status='PASS_TINY_PARTIAL_ROUTE_GUARD_REGRESSIONS'if result.wasSuccessful()else'FAIL_TINY_PARTIAL_ROUTE_GUARD_REGRESSIONS',
        source_sha256=sha256(before).hexdigest(),sources_before=SOURCES,sources_after=after,
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),output=stream.getvalue(),
        scope='12 synthetic control/association tests only. Synthetic inverse entries are not geometric certificates. Zero scientific JSON/worker/Lean/TeX inputs; no G5/G6 numerical-domain replay or global finite57 PASS.'),sort_keys=True)
    need(len(payload.encode())<=500*1024,'tiny regression output cap')
    print(payload)
    if not result.wasSuccessful():raise SystemExit(1)


if __name__=='__main__':main()
