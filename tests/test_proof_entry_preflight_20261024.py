#!/usr/bin/env python3
"""Infrastructure tests only: no mathematical worker, Lean or TeX is run.

The live cleanup test owns two tiny sleeping Python processes in new sessions.
All other subprocesses either refuse optimized startup or produce a plan.
"""
import contextlib
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import selectors
import signal
import subprocess
import sys
import time
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import run_proof_batch_20261024 as entry
import audit_existing_lean_module_20261024 as audit
from monitor_owned_proof_tree_20261024 import OwnedTree, Process, ProcessReader, execution_sources

ENTRY = ROOT / 'scripts/run_proof_batch_20261024.py'


class EntryTests(unittest.TestCase):
    def test_optimization_refused_before_any_downstream_job(self):
        for options, setting in [(['-O'], None), (['-OO'], None), ([], '1')]:
            with self.subTest(options=options, setting=setting):
                environment = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
                environment.pop('PYTHONOPTIMIZE', None)
                if setting is not None:
                    environment['PYTHONOPTIMIZE'] = setting
                result = subprocess.run([sys.executable, '-B', *options, str(ENTRY),
                    'manuscript', 'finite57', '--plan'], cwd=ROOT,
                    env=environment, capture_output=True, text=True, timeout=10)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('ASSERTIONS_REQUIRED', result.stderr)
                self.assertEqual(result.stdout, '')

    def test_check_defaults_to_shared_lock(self):
        command = entry.command_for('check', ['--',
            'scripts/fresh_d13_context_dimension57_20261018.py', '--verify'])
        self.assertIn('--parent-lock', command)
        self.assertTrue(command[0].endswith('run_bounded_research_check_20261018.py'))

    def test_self_locking_allowlist_matches_existing_catalogue(self):
        catalogue = json.loads((ROOT / 'manuscripts/catalogue.json').read_text())
        checks = catalogue['checks']
        if isinstance(checks, dict):
            checks = checks.values()
        declared = {}
        for check in checks:
            if check.get('owns_lock'):
                argv = check['argv']
                if Path(argv[0]).name in ('python', 'python3'):
                    argv = argv[1:]
                if argv[:1] == ['-B']:
                    argv = argv[1:]
                declared[Path(argv[0]).name] = tuple(argv[1:])
        self.assertEqual(entry.SELF_LOCKING_CHECKERS, set(declared))
        self.assertEqual(entry.SELF_LOCKING_ARGUMENTS, declared)

    def test_self_locking_checker_cannot_double_acquire(self):
        worker = 'scripts/check_compact_arithmetic20_composition_20261018.py'
        self.assertNotIn('--parent-lock', entry.command_for('check', ['--', worker, '--verify']))
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            entry.command_for('check', ['--parent-lock', '--', worker, '--verify'])

    def test_internal_worker_cannot_inherit_filename_lock_exemption(self):
        command = entry.command_for('check', ['--',
            'scripts/fresh_compact21_self_hub_exact40_reconstruction_20261018.py', '--worker'])
        self.assertIn('--parent-lock', command)

    def test_lean_routes_keep_outer_monitor_and_inner_lock(self):
        for mode in ('lean-compile', 'lean-audit'):
            command = entry.command_for(mode, ['CoxeterProofs.CompactDimensionSeventeenComposition'])
            self.assertTrue(command[0].endswith('run_bounded_research_check_20261018.py'))
            self.assertNotIn('--parent-lock', command)
        command = entry.command_for('lean-audit', ['CoxeterProofs.CompactDimensionSeventeenComposition'])
        self.assertIn('scripts/audit_existing_lean_module_20261024.py', command)
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            entry.command_for('lean-compile', ['CoxeterProofs'])

    def test_plan_has_separate_phases_without_running_proofs(self):
        environment = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
        environment.pop('PYTHONOPTIMIZE', None)
        for mode, expected in [('full', (136, 9, 5, 1)), ('core', (1, 9, 5, 1))]:
            result = subprocess.run([sys.executable, '-B', str(ENTRY),
                'manuscript', 'finite57', '--mode', mode, '--pdf', '--plan'],
                cwd=ROOT, env=environment, capture_output=True, text=True, timeout=15)
            self.assertEqual(result.returncode, 0, result.stderr)
            plan = json.loads(result.stdout)
            self.assertEqual(plan['status'], 'PLAN_ONLY_NOT_VERIFIED')
            jobs = plan['jobs']
            counts = (sum(j['kind'] == 'script' for j in jobs),
                      sum(j.get('phase') == 'compile' for j in jobs),
                      sum(j.get('phase') == 'separate_actual_axiom_audit' for j in jobs),
                      sum(j['kind'] == 'typesetting' for j in jobs))
            self.assertEqual(counts, expected)
            ids = [j['id'] for j in jobs]
            for job in jobs:
                if job.get('phase') == 'separate_actual_axiom_audit':
                    self.assertLess(ids.index(job['id'].removesuffix('-audit')+'-compile'),
                                    ids.index(job['id']))


class AxiomOutputTests(unittest.TestCase):
    def test_fully_qualified_and_namespace_local_queries(self):
        commands = ['#print axioms local_theorem', '#print axioms CoxeterProofs.T.empty']
        output = ("'CoxeterProofs.T.local_theorem' depends on axioms: [propext,\n"
                  " Classical.choice,\n Quot.sound]\n"
                  "'CoxeterProofs.T.empty' does not depend on any axioms\n").encode()
        parsed = audit.validate_axiom_output(output, commands)
        self.assertEqual(len(parsed), 2)
        self.assertEqual(parsed['CoxeterProofs.T.empty'], [])

    def test_missing_duplicate_malformed_or_unapproved_results_fail(self):
        good = "'CoxeterProofs.T.a' depends on axioms: [propext]\n"
        cases = ['', good+good,
                 "'CoxeterProofs.T.a' depends on axioms: [sorryAx]\n",
                 "'CoxeterProofs.T.a' depends on axioms: [Other.axiom]\n",
                 "'CoxeterProofs.T.a' depends on axioms: [,,]\n",
                 "'CoxeterProofs.T.a' depends on axioms: [propext,,Quot.sound]\n",
                 "'CoxeterProofs.T.a' depends on axioms: [propext,propext]\n",
                 good+"'CoxeterProofs.T.b' depends on axioms: [broken\n",
                 good+"'CoxeterProofs.T.b' does not depend on any axiomsGARBAGE\n"]
        for text in cases:
            with self.subTest(text=text), self.assertRaises(RuntimeError):
                audit.validate_axiom_output(text.encode(), ['#print axioms CoxeterProofs.T.a'])

    def test_ambiguous_namespace_suffix_fails(self):
        output = ("'CoxeterProofs.X.a' depends on axioms: []\n"
                  "'CoxeterProofs.Y.a' depends on axioms: []\n").encode()
        with self.assertRaises(RuntimeError):
            audit.validate_axiom_output(output, ['#print axioms a'])

    def test_real_preserved_stdout_accepts_short_queries(self):
        raw = (ROOT / 'results/ResidualBellmanCertificate_axiom_stdout_20261020.txt').read_bytes()
        source = (ROOT / 'formal/CoxeterProofs/ResidualBellmanCertificate.lean').read_text()
        commands = [line.strip() for line in source.splitlines()
                    if line.lstrip().startswith('#print axioms ')]
        parsed = audit.validate_axiom_output(raw, commands)
        self.assertEqual(len(parsed), len(commands))


class FakeChild:
    pid = 700001
    def wait(self, timeout=None):
        return 0


class FakeReader:
    def __init__(self):
        self.rows = {
            700001: Process(700001, os.getpid(), 700001, (10, 1), False),
            700002: Process(700002, 700001, 700002, (10, 2), False),
            700003: Process(700003, 700002, 700003, (10, 3), False),
        }
    def process(self, pid):
        return self.rows.get(pid)
    def table(self):
        return {pid: (item.ppid, item.pgid) for pid, item in self.rows.items()}


class OwnedTreeTests(unittest.TestCase):
    def test_reparented_owned_groups_cleaned_deepest_first(self):
        reader = FakeReader()
        tree = OwnedTree(FakeChild(), reader)
        tree.sample()
        reader.rows[700003] = Process(700003, 1, 700003, (10, 3), False)
        signals = []
        def terminate(group, signum):
            signals.append((group, signum))
            reader.rows = {pid: item for pid, item in reader.rows.items()
                           if item.pgid != group}
        with patch('os.killpg', side_effect=terminate):
            receipt = tree.stop(grace_seconds=0)
        self.assertEqual([group for group, _ in signals], [700003, 700002, 700001])
        self.assertEqual(receipt['verified_live_owned_processes_remaining'], 0)

    def test_recycled_pid_does_not_inherit_ownership(self):
        reader = FakeReader()
        reader.rows = {700001: reader.rows[700001]}
        tree = OwnedTree(FakeChild(), reader)
        reader.rows[700001] = Process(700001, 1, 700001, (99, 1), False)
        with patch('os.killpg') as kill:
            receipt = tree.stop(grace_seconds=0)
        kill.assert_not_called()
        self.assertEqual(receipt['reused_ids_skipped'], [700001])

    def test_live_owned_nested_sessions_leave_no_child(self):
        reader = ProcessReader()
        reader.probe()  # Permission/ABI failure occurs before either child exists.
        leaf = "import signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); print('READY',flush=True); time.sleep(20)"
        parent = ("import subprocess,sys,signal,time,json,os; "
                  "signal.signal(signal.SIGTERM,signal.SIG_IGN); "
                  f"c=subprocess.Popen([sys.executable,'-B','-c',{leaf!r}],start_new_session=True,stdout=subprocess.PIPE,text=True); "
                  "c.stdout.readline(); print(json.dumps({'leaf':c.pid}),flush=True); time.sleep(20)")
        read_fd, write_fd = os.pipe()
        gated = (f"import os\nready=os.read({read_fd},1)\nos.close({read_fd})\n"
                 "if ready!=b'1': raise SystemExit(125)\n" + parent)
        child = subprocess.Popen([sys.executable, '-B', '-c', gated],
                                 start_new_session=True, stdout=subprocess.PIPE,
                                 pass_fds=(read_fd,))
        os.close(read_fd)
        tree = None
        try:
            tree = OwnedTree(child, reader)
            os.write(write_fd, b'1')
            os.close(write_fd)
            write_fd = None
            readiness = bytearray()
            deadline = time.monotonic() + 5
            with selectors.DefaultSelector() as selector:
                selector.register(child.stdout, selectors.EVENT_READ)
                while b'\n' not in readiness:
                    remaining = deadline - time.monotonic()
                    if remaining <= 0 or not selector.select(remaining):
                        raise RuntimeError('Owned stand-in readiness timed out')
                    chunk = os.read(child.stdout.fileno(), 4096)
                    if not chunk or len(readiness) + len(chunk) > 4096:
                        raise RuntimeError('Invalid owned stand-in readiness')
                    readiness.extend(chunk)
            info = json.loads(readiness.split(b'\n', 1)[0])
            tree.sample()
            receipt = tree.stop(grace_seconds=.1)
            self.assertEqual(receipt['verified_live_owned_processes_remaining'], 0)
            groups = [item['pgid'] for item in receipt['signals']]
            self.assertLess(groups.index(info['leaf']), groups.index(child.pid))
            self.assertIsNotNone(child.poll())
        finally:
            if write_fd is not None:
                os.close(write_fd)
            if tree is not None:
                tree.stop(grace_seconds=.1)
            else:
                # The unreaped direct child is still behind the gate; no leaf
                # exists and it has not installed its TERM-ignore handler.
                child.kill()
                child.wait(timeout=2)
            child.stdout.close()


if __name__ == '__main__':
    started = time.monotonic()
    controllers = execution_sources()
    with Path(__file__).open('rb') as stream:
        own_hash = hashlib.file_digest(stream, 'sha256').hexdigest()
    print(json.dumps(dict(infrastructure_preflight_only=True,
        created_utc=datetime.now(timezone.utc).isoformat(),
        test_source_sha256=own_hash, controller_source_sha256=controllers)), flush=True)
    program = unittest.main(verbosity=2, exit=False)
    with Path(__file__).open('rb') as stream:
        unchanged = hashlib.file_digest(stream, 'sha256').hexdigest() == own_hash
    unchanged = unchanged and execution_sources() == controllers
    success = program.result.wasSuccessful() and unchanged
    print(json.dumps(dict(status='PASS_INFRASTRUCTURE_PREFLIGHT' if success else 'FAILED',
        tests_run=program.result.testsRun, seconds=round(time.monotonic()-started, 3),
        control_and_test_sources_unchanged=unchanged,
        test_source_sha256=own_hash, controller_source_sha256=controllers,
        mathematical_workers_launched=0, lean_or_tex_launched=0)), flush=True)
    sys.exit(0 if success else 1)
