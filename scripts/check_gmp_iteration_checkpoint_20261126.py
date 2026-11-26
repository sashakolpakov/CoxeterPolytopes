#!/usr/bin/env python3
"""Check actual one-iteration GMP checkpoints and a native warm restart."""
from fractions import Fraction as F
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from research.global_formulation import conic,feasibility,sdpa_cli,solver_gate


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',type=Path,required=True)
    parser.add_argument('--fixture',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if not os.environ.get('COXETER_CHECK_OWNER'):raise RuntimeError('Use the serial scientific watchdog')
    if args.output.exists() or not args.output.resolve().is_relative_to(ROOT/'results'):
        raise ValueError('Use a fresh results folder')
    build=json.loads((args.build/'result.json').read_text())
    if build['status']!='PRIVATE_GMP_GRAPH_BUILD_AND_COMPILED_EQUIVALENCE_CHECKS_PASSED':
        raise ValueError('The private GMP build is incomplete')
    executable=Path(build['executable'])
    if solver_gate.file_digest(executable)!=build['executable_sha256']:raise ValueError('Changed native binary')
    args.output.mkdir();input_path=args.output/'input.dat-s'
    shutil.copyfile(args.fixture,input_path)
    sdpa_cli.GMP_EXECUTABLE=executable
    first=args.output/'first';resumed=args.output/'resumed'
    first.mkdir();resumed.mkdir()
    a=sdpa_cli.call(input_path,first,1,40,precision=256,fixture=True)
    b=sdpa_cli.call(input_path,resumed,1,40,initial=first/'native.result',precision=256,fixture=True)
    for folder in (first,resumed):
        for field in ('xVec','xMat','yMat'):feasibility.native_array(folder/'native.result',field)
    initial=(resumed/'gmp_initial.dat').read_text()
    expected='\n'.join(feasibility.native_array(first/'native.result',field)
        for field in ('xVec','xMat','yMat'))+'\n'
    if initial!=expected or not initial.startswith('{'):
        raise ValueError('Restart did not preserve exactly the three native arrays')
    old_mu=re.search(r'^\s*mu\s*=\s*(\S+)',(first/'native.result').read_text(),re.M)
    new_mu=re.search(r'^\s*0\s+(\S+)',(resumed/'solver.stdout').read_text(),re.M)
    if not old_mu or not new_mu or abs(F(old_mu[1])-F(new_mu[1]))>F(6,100)*abs(F(old_mu[1])):
        raise ValueError('Restart complementarity differs from the saved checkpoint')
    if abs(F(new_mu[1])-100)<F(1):raise ValueError('The warm restart used the default initialization')
    report=dict(status='GMP_ONE_ITERATION_CHECKPOINT_AND_WARM_RESTART_PASSED',
        fixture_only=True,research_feasibility_decision=False,
        precision_bits=256,iteration_limit_per_native_call=1,
        all_three_arrays_written_at_iteration_limit=True,
        restart_contains_only_the_three_arrays=True,
        native_report_header_excluded_from_initial_state=True,
        saved_complementarity=str(F(old_mu[1])),restarted_initial_complementarity=str(F(new_mu[1])),
        cold=a,warm=b,source_bindings={str(p):dict(sha256=solver_gate.file_digest(p)) for p in
            [Path(__file__),ROOT/'research/global_formulation/sdpa_cli.py',
             ROOT/'research/global_formulation/feasibility.py',args.build/'result.json']})
    conic.write_json(args.output/'result.json',report)
    print(json.dumps(dict(status=report['status'],output=str(args.output))),flush=True)


if __name__=='__main__':main()
