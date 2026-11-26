#!/usr/bin/env python3
"""Test one SDPA tolerance on the unchanged, gated common-margin input.

This is a numerical diagnostic. Native status never replaces checks in the
original exact model. The historical launcher and preflight-bound sources stay
unchanged.
"""
import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from research.global_formulation import checkpoint_solver, conic, feasibility, sdpa_cli, solver_gate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native-input', type=Path, required=True)
    parser.add_argument('--preflight', type=Path, required=True)
    parser.add_argument('--initial', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--iterations', type=int, default=6)
    parser.add_argument('--seconds', type=int, default=780)
    parser.add_argument('--epsilon-dash', default='1e-8')
    parser.add_argument('--baseline', type=Path,
                        help='Require the identical restart and all other parameter bytes of this prior run')
    args = parser.parse_args()
    if not 1 <= args.iterations <= 1000 or not 1 <= args.seconds <= 800:
        raise ValueError('Invalid numerical diagnostic limits')
    epsilon = conic.F(args.epsilon_dash)
    if not conic.F('1e-12') <= epsilon <= conic.F('1e-6'):
        raise ValueError('Unexpected feasibility tolerance')
    started = time.monotonic()
    gate = solver_gate.verify_native(args.native_input, args.preflight)
    args.output.mkdir()
    paths = [Path(__file__), args.native_input, args.preflight, args.initial,
             sdpa_cli.EXECUTABLE, *[ROOT/'research/global_formulation'/f'{n}.py'
             for n in ['checkpoint_solver', 'conic', 'feasibility', 'sdpa_cli', 'solver_gate']]]
    if args.baseline:
        paths.extend([args.baseline,args.baseline.parent/'param.sdpa'])
    def bindings():
        return {str(p): {'bytes': p.stat().st_size, 'sha256': solver_gate.file_digest(p)}
                for p in paths}
    before = bindings()
    conic.write_json(args.output/'source_bindings.json', before)
    variables, sizes = checkpoint_solver.shape(args.native_input)
    initial_check = checkpoint_solver.checkpoint(args.initial, variables, sizes, require_interior=True)
    conic.write_json(args.output/'initial_checkpoint_check.json', initial_check)
    parameters = args.output/'param.sdpa'
    sdpa_cli.parameters(parameters, args.iterations)
    original_parameters = parameters.read_text()
    lines = original_parameters.splitlines()
    if len(lines) != 14 or lines[9] != '1e-10':
        raise ValueError('Unexpected base parameter layout')
    lines[9] = args.epsilon_dash
    parameters.write_text('\n'.join(lines)+'\n')
    changes = [{'zero_based_line': i, 'old': x, 'new': y}
               for i, (x, y) in enumerate(zip(original_parameters.splitlines(), lines, strict=True)) if x != y]
    if any(p['zero_based_line'] != 9 for p in changes):
        raise ValueError('Changed a parameter other than epsilonDash')
    matched = None
    if args.baseline:
        baseline = json.loads(args.baseline.read_text())
        if baseline['native']['input_sha256'] != gate['native_input_sha256']:
            raise ValueError('Baseline used different problem data')
        if baseline['native']['initial_source_sha256'] != solver_gate.file_digest(args.initial):
            raise ValueError('Baseline used a different initial state')
        baseline_parameters = (args.baseline.parent/'param.sdpa').read_text().splitlines()
        if len(baseline_parameters) != len(lines):
            raise ValueError('Baseline parameter layout differs')
        baseline_changes = [{'zero_based_line': i,'old':x,'new':y}
            for i,(x,y) in enumerate(zip(baseline_parameters,lines,strict=True)) if x!=y]
        if any(v['zero_based_line'] != 9 for v in baseline_changes):
            raise ValueError('Matched comparison changed more than epsilonDash')
        matched = {'baseline':str(args.baseline),'identical_ordered_input':True,
                   'identical_initial_source_bytes':True,
                   'parameter_file_changes':baseline_changes,
                   'baseline_numerical_error_messages':baseline.get('numerical_error_messages',[])}
    initial_path = args.output/'native_initial.dat'
    initial_path.write_text('\n'.join(feasibility.native_array(args.initial, field)
                                     for field in ['xVec', 'xMat', 'yMat'])+'\n')
    command = [str(sdpa_cli.EXECUTABLE), '-ds', str(args.native_input.resolve()),
               '-o', str((args.output/'native.result').resolve()), '-p', str(parameters.resolve()),
               '-numThreads', '1', '-id', str(initial_path.resolve())]
    metadata = {'record_kind':'direct SDPA parameter-probe invocation',
                'command':command,'input_sha256':gate['native_input_sha256'],
                'preflight_check':gate,'initial_source_path':str(args.initial),
                'initial_source_sha256':solver_gate.file_digest(args.initial),
                'initial_sha256':solver_gate.file_digest(initial_path),
                'parameters_sha256':solver_gate.file_digest(parameters),
                'executable_sha256':solver_gate.file_digest(sdpa_cli.EXECUTABLE),
                'one_thread':True,'matched_comparison':matched}
    conic.write_json(args.output/'native_metadata.json',metadata)
    timed_out = False
    with (args.output/'solver.stdout').open('x') as stream:
        try:
            code = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                                  timeout=args.seconds, check=False).returncode
        except subprocess.TimeoutExpired:
            code = None
            timed_out = True
    text = (args.output/'solver.stdout').read_text()
    record = {'status': 'PARAMETER_PROBE_NO_FEASIBILITY_DECISION',
              'verified_preflight': gate, 'command': command, 'returncode': code,
              'timed_out': timed_out, 'only_parameter_changes': changes,
              'iteration_budget': args.iterations, 'native_seconds': args.seconds,
              'geometric_constraints_changed': False, 'accepted_bounds_changed': False,
              'source_bound_geometric_review_complete': False,
              'parameter_sha256': solver_gate.file_digest(parameters),
              'initial_array_file_sha256': solver_gate.file_digest(initial_path),
              'matched_comparison':matched,
              'source_bindings': before}
    errors = [s for s in text.splitlines() if re.search(r'cannot move:|Strange behavior|Cholesky.*failed|\bERROR\b', s)]
    record['numerical_error_messages'] = errors
    if not timed_out and code == 0:
        native = args.output/'native.result'
        with native.open() as stream:
            header = stream.read(32768)
        formats = dict(re.findall(r'^\s*(xPrint|XPrint|YPrint|infPrint)\s*=\s*(\S+)', header, re.M))
        if formats != {n: '%+26.17e' for n in ['xPrint', 'XPrint', 'YPrint', 'infPrint']}:
            raise ValueError('Unexpected native output formats')
        record['actual_output_formats'] = formats
        record['checkpoint'] = checkpoint_solver.checkpoint(native, variables, sizes)
        match = re.search(r'Iteration\s*=\s*(\d+)', text)
        record['reported_iterations'] = int(match[1]) if match else None
        record['iteration_limit_reached'] = 'maxIteration is reached' in text
        if record['checkpoint']['native_phase'] in checkpoint_solver.TERMINAL_PHASES:
            record['status'] = 'NATIVE_TERMINAL_STATUS_REQUIRES_ORIGINAL_CHECK'
        elif errors:
            record['status'] = 'NATIVE_SOLVER_NUMERICAL_STOP'
        elif record['iteration_limit_reached']:
            record['status'] = 'PARAMETER_PROBE_CHECKPOINT_REQUIRES_CONTINUATION'
    if bindings() != before:
        raise ValueError('Probe input or source changed')
    record['elapsed_seconds'] = time.monotonic()-started
    conic.write_json(args.output/'result.json', record)
    print(json.dumps({k:v for k,v in record.items() if k not in ['source_bindings','command']}), flush=True)


if __name__ == '__main__':
    main()
