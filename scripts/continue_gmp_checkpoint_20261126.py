#!/usr/bin/env python3
"""Restart the verified full G(u) SDP from native GMP arrays without binary64.

The exact-export and private-build gates remain in force. Cone positivity of a
barrier state authorizes a restart, not a feasibility claim for G(u).
"""
import argparse
from fractions import Fraction as F
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from research.global_formulation import conic,feasibility,field_rounding,solver_gate
from scripts import run_balanced_dual_gmp_20261126 as gmp
RUN_DEADLINE=None


def private_solver(build_folder):
    b=gmp.read(build_folder/'result.json')
    if b['status']!='PRIVATE_GMP_GRAPH_BUILD_AND_COMPILED_EQUIVALENCE_CHECKS_PASSED':
        raise ValueError('Incomplete private solver build')
    gmp.check_bindings(b['source_bindings'])
    if solver_gate.file_digest(b['executable'])!=b['executable_sha256']:raise ValueError('Changed solver')
    gmp.EXECUTABLE=Path(b['executable']);gmp.EXECUTABLE_SHA=b['executable_sha256']
    gmp.INSTALLATION=Path(b['installation'])
    return b


def verified_plan(folder):
    plan=gmp.read(folder/'result.json');gmp.check_bindings(plan['source_bindings'])
    if plan['status']!='GMP_CHECKPOINT_INPUT_READY_LONGER_WALL_AUTHORIZATION_REQUIRED':
        raise ValueError('Wrong bound GMP plan')
    prepared=Path(plan['prepared']);verified=Path(plan['verified'])
    gate=gmp.read(verified/'result.json')
    if (gate.get('status')!='EXACT_DUAL_GMP_ALL_NATIVE_COEFFICIENTS_REPLAYED'
            or gate.get('solver_ready') is not True
            or gate['prepared_result_sha256']!=solver_gate.file_digest(prepared/'result.json')
            or gate['native_input_sha256']!=solver_gate.file_digest(prepared/'input.dat-s')
            or gate['native_input_sha256']!=plan['native_input_sha256']):
        raise ValueError('No matching complete exact matrix export replay')
    gmp.check_bindings(gate['source_bindings'])
    private_solver(Path(plan['build']))
    return plan,gmp.read(prepared/'result.json')


def state(path,m,sizes):
    fields={k:feasibility.native_array(path,k) for k in ('xVec','xMat','yMat')}
    vector=[F(x) for x in conic.NUMBER.findall(fields['xVec'])]
    if len(vector)!=m:raise ValueError('Wrong checkpoint free-vector dimension')
    degree=sum(abs(n) for n in sizes)
    field=field_rounding.Field(0);zero=field.zero
    def values(body):
        return [F(x) for x in conic.NUMBER.findall(body)]
    matrices={k:values(fields[k]) for k in ('xMat','yMat')}
    expected=sum(n*n if n>0 else -n for n in sizes)
    if any(len(v)!=expected for v in matrices.values()):raise ValueError('Wrong checkpoint cone dimensions')
    max_relative_asymmetry=F(0)
    for name,flat in matrices.items():
        position=0
        for n in sizes:
            if n<0:
                if any(v<=0 for v in flat[position:position-n]):
                    raise ValueError('Nonpositive checkpoint scalar cone in '+name)
                position-=n;continue
            block=flat[position:position+n*n];position+=n*n
            scale=max(abs(v) for v in block)
            asym=max(abs(block[i*n+j]-block[j*n+i]) for i in range(n) for j in range(n))
            if asym and (not scale or asym>scale/F(10**70)):
                raise ValueError('Checkpoint matrix has material asymmetry')
            if scale:max_relative_asymmetry=max(max_relative_asymmetry,asym/scale)
            # Check the native symmetric triangle exactly; preserve both printed
            # triangles unaltered in the restart file.
            M=[[field.element(block[min(i,j)*n+max(i,j)]) for j in range(n)] for i in range(n)]
            pivots=field_rounding.psd(M,field)
            if any(field.sign(p)<=0 for p in pivots):
                raise ValueError('Checkpoint matrix is not strictly positive definite')
    dot=F(0);position=0
    for n in sizes:
        if n<0:
            dot+=sum((matrices['xMat'][position+i]*matrices['yMat'][position+i] for i in range(-n)),F(0))
            position-=n;continue
        for i in range(n):
            for j in range(n):dot+=matrices['xMat'][position+i*n+j]*matrices['yMat'][position+j*n+i]
        position+=n*n
    return fields,dict(free_variables=m,cone_state_strictly_positive=True,
        exact_rational_upper_triangle_LDL=True,
        maximum_relative_printed_asymmetry=str(max_relative_asymmetry),
        checkpoint_barrier_mean_complementarity=gmp.decimal(dot/degree),
        exact_barrier_mean_complementarity=str(dot/degree),
        affine_feasibility_not_inferred_from_barrier_state=True)


def native(input_path,out,seconds,iterations,initial=None,shape=None):
    solver=gmp.executable_check();gmp.parameters(out/'param.sdpa',256,iterations)
    command=[str(gmp.EXECUTABLE),'-ds',str(input_path),'-o',str(out/'native.result'),'-p',str(out/'param.sdpa')]
    initial_info=None
    if initial is not None:
        fields,initial_info=state(initial,*shape)
        checkpoint=out/'native_initial.dat'
        checkpoint.write_text('\n'.join(fields[k] for k in ('xVec','xMat','yMat'))+'\n')
        del fields
        if not checkpoint.read_text().startswith('{'):raise ValueError('Report header leaked into initial state')
        command+=['-id',str(checkpoint)]
        initial_info.update(source=str(initial),source_sha256=solver_gate.file_digest(initial),
            initialization_sha256=solver_gate.file_digest(checkpoint),only_three_arrays=True)
        conic.write_json(out/'initial_state.json',initial_info)
    conic.write_json(out/'native_launch.json',dict(command=command,solver=solver,
        precision_bits=256,iterations=iterations,warm_start=initial is not None,
        native_input_sha256=solver_gate.file_digest(input_path)))
    started=time.monotonic();stop=None;peak=0
    if RUN_DEADLINE is not None:seconds=min(seconds,max(1,RUN_DEADLINE-started-180))
    with (out/'solver.stdout').open('x') as stream:
        worker=subprocess.Popen(command,stdout=stream,stderr=subprocess.STDOUT)
        while worker.poll() is None:
            current=gmp.process_rss(os.getpid());peak=max(peak,current)
            if current>8_800_000_000:stop='OWNED_NATIVE_MEMORY_LIMIT'
            if time.monotonic()-started>seconds:stop='NATIVE_BATCH_TIME_LIMIT'
            if stop:
                worker.terminate()
                try:worker.wait(timeout=8)
                except subprocess.TimeoutExpired:worker.kill();worker.wait()
                break
            time.sleep(.5)
    path=out/'native.result'
    text=path.read_text() if path.exists() and path.stat().st_size<128*1024**2 else ''
    stdout=(out/'solver.stdout').read_text()
    phase=re.search(r'phase\.value\s*=\s*(\S+)',text)
    precision=re.search(r'^\s*precision\s*=\s*(\d+)',text,re.M)
    formats=dict(re.findall(r'^\s*(xPrint|XPrint|YPrint|infPrint)\s*=\s*(\S+)',text,re.M))
    if worker.returncode==0 and not stop and phase:
        if not precision or int(precision[1])!=256 or formats!={k:'%+100.90Fe' for k in ('xPrint','XPrint','YPrint','infPrint')}:
            raise ValueError('Wrong native arithmetic precision or output formats')
    metrics={k:v for k,v in re.findall(r'^\s*(mu|objValPrimal|objValDual|p\.feas\.error|d\.feas\.error|relative eps|gap)\s*=\s*(\S+)',text,re.M)}
    return dict(status=stop or ('NATIVE_TERMINAL_STATUS' if phase else 'NATIVE_ERROR_NO_TERMINAL_PHASE'),
        phase=phase[1] if phase else None,returncode=worker.returncode,command=command,solver=solver,
        precision_bits=256,reported_precision_bits=int(precision[1]) if precision else None,
        full_precision_output_formats=formats,metrics=metrics,
        solution_arrays_present={k:k in text for k in ('xVec','xMat','yMat')},
        native_result_sha256=solver_gate.file_digest(path) if path.exists() else None,
        observed_owned_process_tree_rss_bytes=peak,elapsed_seconds=time.monotonic()-started,
        actual_native_wall_limit_seconds=seconds,
        initialization=initial_info,iteration_limit=iterations,
        completed_checkpoint_is_not_a_feasibility_certificate=True)


def check(args):
    private_solver(args.build)
    inp=args.output/'input.dat-s';shutil.copyfile(args.fixture,inp)
    fixture=solver_gate.verify_native(inp,None,fixture=True)
    first=args.output/'first';warm=args.output/'warm';first.mkdir();warm.mkdir()
    a=native(inp,first,40,1);b=native(inp,warm,40,1,first/'native.result',(1,[2,-1]))
    s1=state(first/'native.result',1,[2,-1])[1];s2=state(warm/'native.result',1,[2,-1])[1]
    if b['initialization']['exact_barrier_mean_complementarity']!=s1['exact_barrier_mean_complementarity']:
        raise ValueError('Transport changed checkpoint complementarity')
    saved=F(s1['exact_barrier_mean_complementarity'])
    mu=re.search(r'^\s*0\s+(\S+)',(warm/'solver.stdout').read_text(),re.M)
    if not mu or abs(F(mu[1])-saved)>F(6,100)*abs(saved):raise ValueError('Native restart ignored its saved state')
    conic.write_json(args.output/'result.json',dict(status='GMP_CONTINUATION_ARRAYS_PD_AND_NATIVE_RESTART_FIXTURE_PASSED',
        fixture_only=True,research_feasibility_decision=False,fixture_gate=fixture,
        cold=a,warm=b,cold_checkpoint=s1,warm_checkpoint=s2,
        source_bindings=gmp.bind([Path(__file__),Path(gmp.__file__)])))


def run(args):
    if int(os.environ.get('COXETER_EXPLICIT_GMP_WALL_SECONDS','0'))<3600:
        raise RuntimeError('Use the authorized one-hour monitor')
    plan,prepared=verified_plan(args.plan)
    previous=gmp.read(args.resume/'result.json');answer=args.resume/'native.result'
    if (Path(previous['prepared']).resolve()!=Path(plan['prepared']).resolve()
            or previous['native']['returncode']!=0
            or not all(previous['native']['solution_arrays_present'].values())
            or previous['native']['native_result_sha256']!=solver_gate.file_digest(answer)):
        raise ValueError('No intact completed checkpoint for this exact problem')
    info=native(Path(plan['prepared'])/'input.dat-s',args.output,3300,args.iterations,answer,
        (prepared['variables'],prepared['native_block_sizes']))
    validation=None
    if info['returncode']==0 and all(info['solution_arrays_present'].values()):
        validation=state(args.output/'native.result',prepared['variables'],prepared['native_block_sizes'])[1]
    conic.write_json(args.output/'result.json',dict(status='GMP_FULL_MODEL_CONTINUATION_BATCH_COMPLETED',
        prepared=plan['prepared'],separate_input_replay=str(Path(plan['verified'])/'result.json'),
        native=info,checkpoint_state=validation,resumed_from=str(args.resume),
        iterations_requested_in_batch=args.iterations,problem_role='Explicit normalized dual G(u), trace objective',
        original_feasibility_status='UNDECIDED',exact_certificate_obtained=False,
        source_bindings=gmp.bind([Path(__file__),Path(gmp.__file__),args.resume/'result.json',
            args.plan/'result.json',Path(plan['prepared'])/'result.json',Path(plan['verified'])/'result.json'])))


def main():
    global RUN_DEADLINE
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['check','run'])
    for name in ('build','fixture','plan','resume','output'):parser.add_argument('--'+name,type=Path,required=name=='output')
    parser.add_argument('--iterations',type=int,default=1)
    args=parser.parse_args()
    if not os.environ.get('COXETER_CHECK_OWNER'):raise RuntimeError('Use a scientific monitor')
    args.output=args.output.resolve()
    if args.output.exists() or not args.output.is_relative_to(ROOT/'results'):raise ValueError('Use a fresh results folder')
    if not 1<=args.iterations<=10:raise ValueError('Use short native iteration checkpoints')
    if args.action=='run':RUN_DEADLINE=time.monotonic()+3500
    args.output.mkdir();globals()[args.action](args)
    print(json.dumps(dict(status=gmp.read(args.output/'result.json')['status'],output=str(args.output))),flush=True)


if __name__=='__main__':main()
