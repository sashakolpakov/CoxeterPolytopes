#!/usr/bin/env python3
"""Assertion-safe public entry; retain the original resource-gated runners.

This file performs no mathematics.  Planning stays planning.  Actual checks,
Lean and TeX still pass through the existing serial monitors and admission
gates.  Frozen workers, runners and their source-bound receipts are unchanged.
"""
import os
import sys

# These must be real conditionals, before any repository module is imported.
if not __debug__ or sys.flags.optimize:
    raise SystemExit('ASSERTIONS_REQUIRED: Python -O/-OO cannot verify proofs')
if os.environ.get('PYTHONOPTIMIZE', '') not in ('', '0'):
    raise SystemExit('ASSERTIONS_REQUIRED: remove PYTHONOPTIMIZE before verification')

import argparse
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / 'scripts'
SELF_LOCKING_CHECKERS = {
    'fresh_compact21_self_hub_exact40_reconstruction_20261018.py',
    'check_compact21_self_hub_primal_independent_20261018.py',
    'check_compact_arithmetic20_composition_20261018.py',
    'check_ideal_ra4_octahedral_path_roots_20261018.py',
    'fresh_ideal_ra6_cusp_transport_20261018.py',
    'fresh_ideal_ra6_sharp_reciprocal_20261018.py',
    'fresh_compact21_supported_forest_triangle_coupling_20261018.py',
}
SELF_LOCKING_ARGUMENTS = {
    'fresh_compact21_self_hub_exact40_reconstruction_20261018.py': ('--verify',),
    'check_compact21_self_hub_primal_independent_20261018.py': (
        '--candidate', 'results/fresh_compact21_self_hub_exact40_reconstruction_20261018.json'),
    'check_compact_arithmetic20_composition_20261018.py': ('--verify',),
    'check_ideal_ra4_octahedral_path_roots_20261018.py': (),
    'fresh_ideal_ra6_cusp_transport_20261018.py': ('--verify',),
    'fresh_ideal_ra6_sharp_reciprocal_20261018.py': ('--verify',),
    'fresh_compact21_supported_forest_triangle_coupling_20261018.py': ('--verify',),
}


def command_for(mode, arguments):
    """Return an argv, never a shell command or a weakened replacement runner."""
    if mode == 'check':
        if '--' not in arguments:
            raise ValueError('Use check [--ledger PATH] -- scripts/CHECKER.py [arguments]')
        boundary = arguments.index('--')
        options, worker = arguments[:boundary], arguments[boundary+1:]
        parser = argparse.ArgumentParser(prog=mode)
        parser.add_argument('--ledger')
        parser.add_argument('--parent-lock', action='store_true')
        args = parser.parse_args(options)
        if not worker:
            parser.error('Missing checker')
        source = (ROOT / worker[0]).resolve()
        if source.parent != SCRIPTS or source.suffix != '.py' or not source.is_file():
            parser.error('Choose an existing repository checker')
        # A filename alone is insufficient: one frozen source has an internal
        # --worker branch that runs before its public main acquires the lock.
        owns_lock = (source.name in SELF_LOCKING_CHECKERS
                     and tuple(worker[1:]) == SELF_LOCKING_ARGUMENTS[source.name])
        if owns_lock and args.parent_lock:
            parser.error('This checker owns the shared lock; omit --parent-lock')
        command = [str(SCRIPTS / 'run_bounded_research_check_20261018.py')]
        if args.ledger:
            command += ['--ledger', args.ledger]
        if not owns_lock:
            command.append('--parent-lock')
        return [*command, '--', *worker]
    if mode == 'manuscript':
        return [str(SCRIPTS / 'verify_manuscript_audited_20261024.py'), *arguments]
    if mode in ('lean-compile', 'lean-audit'):
        parser = argparse.ArgumentParser(prog=mode)
        parser.add_argument('module')
        parser.add_argument('--ledger')
        args = parser.parse_args(arguments)
        if not re.fullmatch(r'CoxeterProofs(?:\.[A-Za-z0-9_]+)+', args.module):
            parser.error('Choose one named module, not the whole-project root')
        target = ('run_lean_batches_signal_safe_20261018.py'
                  if mode == 'lean-compile'
                  else 'audit_existing_lean_module_20261024.py')
        command = [str(SCRIPTS / 'run_bounded_research_check_20261018.py')]
        if args.ledger:
            command += ['--ledger', args.ledger]
        # The inner Lean compiler/auditor owns exclusive.lock.  Taking it here
        # with --parent-lock would cause a second acquisition, not extra safety.
        return command + ['--', 'scripts/' + target, args.module]
    if mode == 'tex':
        parser = argparse.ArgumentParser(prog=mode, add_help=False)
        parser.add_argument('--builder')
        args, rest = parser.parse_known_args(arguments)
        builder = args.builder or 'scripts/build_illustrated_ra_manuscript_20261018.py'
        allowed = {
            'scripts/build_illustrated_ra_manuscript_20261018.py',
            'scripts/build_compact_arithmetic_nine_field19_20261020.py',
        }
        if builder not in allowed:
            parser.error('Unknown bounded TeX builder')
        return [str(ROOT / builder), *rest]
    raise ValueError('Unknown mode: ' + mode)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('check', 'manuscript', 'lean-compile',
                                         'lean-audit', 'tex'))
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = command_for(args.mode, args.arguments)
    environment = dict(os.environ)
    environment.pop('PYTHONOPTIMIZE', None)
    environment['PYTHONDONTWRITEBYTECODE'] = '1'
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        environment[name] = '1'
    # exec retains the owner's process relationship.  -E also prevents Python
    # environment options from optimizing the selected runner; its descendants
    # inherit the explicitly cleaned environment.
    os.chdir(ROOT)
    selected = [sys.executable, '-E', '-B', *command]
    if args.mode != 'manuscript':
        selected = [sys.executable, '-E', '-B',
                    str(SCRIPTS / 'monitor_owned_proof_tree_20261024.py'),
                    '--', *selected]
    os.execve(sys.executable, selected, environment)


if __name__ == '__main__':
    main()
