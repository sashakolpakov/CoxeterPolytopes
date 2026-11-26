#!/usr/bin/env python3
"""Prepare a bound one-iteration GMP job; running it needs a longer wall policy.

The current 900-second watchdog cannot authorize the longer native call.
The plan changes no mathematical input, precision, thread or RSS settings.
"""
import argparse
import json
import os
from pathlib import Path
from types import SimpleNamespace
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from research.global_formulation import conic,solver_gate
from scripts import run_balanced_dual_gmp_20261126 as gmp
from scripts import sdpa_gmp_sparse_graph_20261126 as graph_backend


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['plan','run'])
    for name in ('build','prepared','verified','plan'):
        parser.add_argument('--'+name,type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if not os.environ.get('COXETER_CHECK_OWNER'):raise RuntimeError('Use an authorized scientific watchdog')
    args.output=args.output.resolve()
    if args.output.exists() or not args.output.is_relative_to(ROOT/'results'):raise ValueError('Use fresh output')
    if args.action=='plan':
        verified=gmp.read(args.verified/'result.json')
        if (verified.get('status')!='EXACT_DUAL_GMP_ALL_NATIVE_COEFFICIENTS_REPLAYED'
                or verified.get('solver_ready') is not True
                or verified['prepared_result_sha256']!=solver_gate.file_digest(args.prepared/'result.json')
                or verified['native_input_sha256']!=solver_gate.file_digest(args.prepared/'input.dat-s')):
            raise ValueError('Missing complete GMP input replay')
        gmp.check_bindings(verified['source_bindings'])
        build=gmp.read(args.build/'result.json')
        if build['status']!='PRIVATE_GMP_GRAPH_BUILD_AND_COMPILED_EQUIVALENCE_CHECKS_PASSED':
            raise ValueError('Incomplete private solver build')
        gmp.check_bindings(build['source_bindings'])
        if solver_gate.file_digest(build['executable'])!=build['executable_sha256']:
            raise ValueError('Changed private binary')
        args.output.mkdir();gmp.parameters(args.output/'one_iteration.sdpa',256,1)
        report=dict(status='GMP_CHECKPOINT_INPUT_READY_LONGER_WALL_AUTHORIZATION_REQUIRED',
            native_execution_performed=False,build=str(args.build),prepared=str(args.prepared),verified=str(args.verified),
            native_input_sha256=verified['native_input_sha256'],executable=build['executable'],
            executable_sha256=build['executable_sha256'],precision_bits=256,native_iteration_limit=1,
            requested_outer_wall_seconds=3600,requested_native_wall_seconds=3300,
            threads=1,rss_limit_bytes=10_000_000_000,shared_lock_and_percentage_gates_required=True,
            source_of_current_wall_limit='AGENTS.md:118: future research batches: 900 seconds',
            explicit_user_wall_time_override_received=False,
            purpose='Complete one Newton iteration and write all native arrays before a restart.',
            checkpoint_protocol_fixture='results/global_gmp_checkpoint_checks02/result.json',
            parameters_sha256=solver_gate.file_digest(args.output/'one_iteration.sdpa'),
            source_bindings=gmp.bind([Path(__file__),Path(gmp.__file__),Path(graph_backend.__file__),
                args.build/'result.json',args.prepared/'result.json',args.verified/'result.json']))
        conic.write_json(args.output/'result.json',report)
    else:
        if int(os.environ.get('COXETER_EXPLICIT_GMP_WALL_SECONDS','0'))<3600:
            raise RuntimeError('This run requires an explicitly authorized longer monitored wall policy')
        plan=gmp.read(args.plan/'result.json');gmp.check_bindings(plan['source_bindings'])
        if plan['status']!='GMP_CHECKPOINT_INPUT_READY_LONGER_WALL_AUTHORIZATION_REQUIRED':
            raise ValueError('Wrong pending job plan')
        args.output.mkdir()
        graph_backend.backend(SimpleNamespace(action='solve',build=Path(plan['build']),
            prepared=Path(plan['prepared']),verified=Path(plan['verified']),output=args.output,
            precision=256,iterations=1,seconds=3300))
        result=gmp.read(args.output/'result.json')
        result['longer_wall_plan']=str(args.plan/'result.json')
        result['longer_wall_plan_sha256']=solver_gate.file_digest(args.plan/'result.json')
        result['authorized_wall_seconds_from_monitor']=int(os.environ['COXETER_EXPLICIT_GMP_WALL_SECONDS'])
        (args.output/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(dict(status=gmp.read(args.output/'result.json')['status'],output=str(args.output))),flush=True)


if __name__=='__main__':main()
