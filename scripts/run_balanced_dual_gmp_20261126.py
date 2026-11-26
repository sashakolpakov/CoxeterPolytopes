#!/usr/bin/env python3
"""Export, replay and run SDPA-GMP on the replayed full explicit dual.

The constraint set is unchanged. A positive multiple of trace(G) supplies
a bounded-below objective and the exactly feasible strict dual point h I.
Coefficients come from Q, never from the binary64 conic encoding.
"""
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as F
from functools import lru_cache
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from research.global_formulation import conic,solver_gate

EXECUTABLE=Path.home()/'.local/opt/sdpa-gmp/bin/sdpa_gmp'
EXECUTABLE_SHA='430440d3810d78b9b6a24f657aa84f1797ca45b2072777e92265614287f6461f'
INSTALLATION=EXECUTABLE.parent.parent/'share/sdpa-gmp/installation.json'
DIGITS=100


def read(path):
    return json.loads(Path(path).read_text())


def bind(paths):
    return {str(p):dict(bytes=Path(p).stat().st_size,sha256=solver_gate.file_digest(p)) for p in paths}


def check_bindings(rows):
    for p,row in rows.items():
        if solver_gate.file_digest(p)!=row['sha256']:
            raise ValueError('Changed input or source: '+p)


def executable_check():
    installation=read(INSTALLATION)
    if (solver_gate.file_digest(EXECUTABLE)!=EXECUTABLE_SHA
            or installation['binary_sha256']!=EXECUTABLE_SHA):
        raise ValueError('Changed SDPA-GMP executable')
    return dict(executable=str(EXECUTABLE),sha256=EXECUTABLE_SHA,
        upstream=installation['upstream'],upstream_commit=installation['commit'],
        version=installation['version'],installation=bind([INSTALLATION]))


@lru_cache(maxsize=200000)
def rational(text,exponent=0):
    q=F(text)
    return q*(1<<exponent) if exponent>=0 else q/F(1<<(-exponent))


@lru_cache(maxsize=200000)
def decimal(q):
    with localcontext() as context:
        context.prec=DIGITS
        return format(Decimal(q.numerator)/Decimal(q.denominator),'+.99E')


def check_decimal(token,q):
    # A separate exact rational bound on each distinct emitted number.
    a=F(token)
    if (not q and a) or (q and (not a or abs(a-q)*10**98>abs(q))):
        raise ValueError('The native coefficient exceeds its 1e-98 relative input error bound')


def source_model(args):
    from scripts.search_balanced_dual_field_20261126 import checked_model,checked_cone_units
    from scripts.solve_balanced_certificate_20261126 import verified_balanced_model
    from research.global_formulation import conic_backend as backend
    verified_balanced_model(args.balanced,args.native_input,args.preflight)
    replay=checked_model(args.model,args.replay)
    data,layout=backend.load_encoding(args.model)
    _,exponents=checked_cone_units(args.units,data,layout,args.model)
    coordinates=read(args.model/'coordinates.json')
    with gzip.open(args.model/'blocks.json.gz','rt') as stream:blocks=json.load(stream)
    sources=[Path(__file__),ROOT/'research/global_formulation/conic.py',
        ROOT/'research/global_formulation/solver_gate.py',
        ROOT/'scripts/search_balanced_dual_field_20261126.py',
        ROOT/'scripts/solve_balanced_certificate_20261126.py',
        ROOT/'scripts/rescale_balanced_dual_20261126.py',
        args.model/'result.json',args.model/'forms.exact.jsonl.gz',args.model/'coordinates.json',
        args.model/'blocks.json.gz',args.replay,args.units/'result.json',
        args.units/'congruence.exact.json',args.units/'units.npz',args.preflight]
    return coordinates,blocks,list(map(int,exponents)),replay,bind(sources)


def slots(blocks):
    full=[b for b in blocks if b['size']>1]
    full_ids={b['start']:k+1 for k,b in enumerate(full)}
    scalar=0
    for block in blocks:
        n=block['size']
        if n==1:
            scalar+=1
            yield len(full)+1,scalar,scalar,True,block['name']
        else:
            for i in range(n):
                for j in range(i,n):
                    yield full_ids[block['start']],i+1,j+1,i==j,block['name']


def forms(path):
    with gzip.open(path,'rt') as stream:
        for line in stream:
            yield json.loads(line)


def header_and_objective(meta,blocks,exponents,path):
    m=len(meta['variables'])
    objective=[F(0) for _ in range(m)]
    constant=F(0)
    diagonal_count=sum(b['size'] for b in blocks)
    unit=F(1,1<<diagonal_count.bit_length())
    count=0
    for row,(_,_,_,diagonal,_),e in zip(forms(path),slots(blocks),exponents,strict=True):
        count+=1
        if diagonal:
            for k,text in row:
                q=rational(text,e)*unit
                if k==-1:constant+=q
                else:objective[k]+=q
    if count!=meta['forms']:raise ValueError('Wrong matrix entry count')
    sizes=[b['size'] for b in blocks if b['size']>1]
    scalars=sum(b['size']==1 for b in blocks)
    if scalars:sizes.append(-scalars)
    header=['"Full explicit dual, Q coefficients, 100 significant decimal digits',
        str(m),str(len(sizes)),' '.join(map(str,sizes)),' '.join(decimal(q) for q in objective)]
    return header,objective,dict(positive_trace_unit=str(unit),trace_constant=str(constant),
        feasibility_constraints_unchanged=True,
        objective='h * trace(G(u)); the constant is omitted from the native linear objective',
        exact_strict_dual_point='h I on every retained cone',
        objective_lower_bound=str(-constant),
        objective_is_not_an_additional_geometric_condition=True),sizes


def prepare(args):
    meta,blocks,exponents,replay,bindings=source_model(args)
    header,objective,objective_meta,sizes=header_and_objective(meta,blocks,exponents,args.model/'forms.exact.jsonl.gz')
    count=0;largest_support=0;PSD_support=[];current=set();last=None
    # Only support sizes are needed for the native sparse/dense selection audit.
    with (args.output/'input.dat-s').open('x') as stream:
        stream.write('\n'.join(header)+'\n')
        for row,(block,i,j,_,_),e in zip(forms(args.model/'forms.exact.jsonl.gz'),slots(blocks),exponents,strict=True):
            if block!=last and last is not None and last<=len(sizes)-(sizes[-1]<0):
                PSD_support.append(len(current));current=set()
            last=block
            support={k for k,text in row if k>=0 and rational(text,e)}
            if block<=len(sizes)-(sizes[-1]<0):current.update(support)
            else:largest_support=max(largest_support,len(support))
            for k,text in row:
                q=rational(text,e)
                if not q:continue
                matrix=k+1
                if k==-1:q=-q
                stream.write(f'{matrix} {block} {i} {j} {decimal(q)}\n');count+=1
    if current:PSD_support.append(len(current))
    conic.write_json(args.output/'objective.exact.json',dict(coefficients=list(map(str,objective)),**objective_meta))
    conic.write_json(args.output/'source_bindings.json',bindings)
    conic.write_json(args.output/'result.json',dict(status='EXACT_DUAL_GMP_EXPORT_REQUIRES_SEPARATE_REPLAY',
        solver_ready=False,model=str(args.model),separate_model_replay=str(args.replay),
        cone_units=str(args.units),variables=len(objective),native_block_sizes=sizes,
        scalar_cones=meta['scalar_cones'],PSD_blocks=len(meta['PSD_sizes']),
        sparse_nonzeros=count,decimal_significant_digits=DIGITS,arithmetic_source='Q, actual matrix entries',
        input_binary64_conversion_used=False,input_relative_decimal_error_bound='1e-98',
        largest_scalar_variable_support=largest_support,largest_PSD_variable_support=max(PSD_support,default=0),
        source_model_exact_sha256=replay['exact_model_sha256'],objective=objective_meta,
        native_input_sha256=solver_gate.file_digest(args.output/'input.dat-s'),
        source_bindings=bindings,original_feasibility_status='UNDECIDED'))


def verify(args):
    prepared=read(args.prepared/'result.json')
    if prepared['status']!='EXACT_DUAL_GMP_EXPORT_REQUIRES_SEPARATE_REPLAY':raise ValueError('Wrong export status')
    check_bindings(prepared['source_bindings'])
    if solver_gate.file_digest(args.prepared/'input.dat-s')!=prepared['native_input_sha256']:
        raise ValueError('Changed GMP input')
    meta,blocks,exponents,_,bindings=source_model(args)
    header,objective,objective_meta,sizes=header_and_objective(meta,blocks,exponents,args.model/'forms.exact.jsonl.gz')
    count=0;checked=set()
    with (args.prepared/'input.dat-s').open() as stream:
        for expected in header:
            if stream.readline().rstrip('\n')!=expected:raise ValueError('Changed GMP header/objective')
        for token,q in zip(header[-1].split(),objective,strict=True):check_decimal(token,q)
        for row,(block,i,j,_,_),e in zip(forms(args.model/'forms.exact.jsonl.gz'),slots(blocks),exponents,strict=True):
            for k,text in row:
                q=rational(text,e)
                if not q:continue
                if k==-1:q=-q
                entry=stream.readline().split()
                if (len(entry)!=5 or list(map(int,entry[:4]))!=[k+1,block,i,j]):
                    raise ValueError('Changed coefficient index or omitted matrix entry')
                if (entry[4],q) not in checked:
                    check_decimal(entry[4],q);checked.add((entry[4],q))
                count+=1
        if stream.read().strip():raise ValueError('Extra native coefficients')
    if count!=prepared['sparse_nonzeros']:raise ValueError('Wrong coefficient count')
    if read(args.prepared/'objective.exact.json')!=dict(coefficients=list(map(str,objective)),**objective_meta):
        raise ValueError('Changed exact trace objective')
    conic.write_json(args.output/'result.json',dict(status='EXACT_DUAL_GMP_ALL_NATIVE_COEFFICIENTS_REPLAYED',
        solver_ready=True,prepared=str(args.prepared),prepared_result_sha256=solver_gate.file_digest(args.prepared/'result.json'),
        native_input_sha256=prepared['native_input_sha256'],native_coefficients_checked=count,
        exact_trace_objective_checked=True,decimal_error_bound_checked_over_Q='1e-98 relative',
        distinct_decimal_coefficients_checked=len(checked),source_bindings=bindings,
        variables=len(objective),native_block_sizes=sizes,original_feasibility_status='UNDECIDED'))


def parameters(path,precision,iterations):
    lines=[str(iterations),'1e-45','10','10000','-1e100','1e100','0.1','0.3','0.9','1e-45',str(precision)]
    lines+=['%+100.90Fe = '+key for key in ('xPrint','XPrint','YPrint','infPrint')]
    path.write_text('\n'.join(lines)+'\n')


def process_rss(pid):
    rows=subprocess.check_output(['ps','-axo','pid=,ppid=,rss='],text=True)
    table={}
    for line in rows.splitlines():
        fields=line.split()
        if len(fields)==3:table[int(fields[0])]=(int(fields[1]),int(fields[2])*1024)
    owned={pid};changed=True
    while changed:
        newer={p for p,(parent,_) in table.items() if parent in owned}
        changed=bool(newer-owned);owned|=newer
    return sum(table.get(p,(0,0))[1] for p in owned)


def native(input_path,out,precision,iterations,seconds):
    solver=executable_check()
    parameters(out/'param.sdpa',precision,iterations)
    command=[str(EXECUTABLE),'-ds',str(input_path),'-o',str(out/'native.result'),'-p',str(out/'param.sdpa')]
    started=time.monotonic();stop=None;peak=0
    with (out/'solver.stdout').open('x') as stream:
        worker=subprocess.Popen(command,stdout=stream,stderr=subprocess.STDOUT)
        while worker.poll() is None:
            current=process_rss(os.getpid());peak=max(peak,current)
            # Leave room for the outer watchdog, stop only this owned solver.
            if current>8_800_000_000:stop='OWNED_NATIVE_MEMORY_LIMIT'
            if time.monotonic()-started>seconds:stop='NATIVE_BATCH_TIME_LIMIT'
            if stop:
                worker.terminate()
                try:worker.wait(timeout=8)
                except subprocess.TimeoutExpired:worker.kill();worker.wait()
                break
            time.sleep(.5)
    result_path=out/'native.result'
    text=result_path.read_text() if result_path.is_file() and result_path.stat().st_size<128*1024**2 else ''
    stdout=(out/'solver.stdout').read_text()
    phase=re.search(r'phase\.value\s*=\s*(\S+)',text)
    reported=re.search(r'^\s*precision\s*=\s*(\d+)',text,re.M)
    actual_formats=dict(re.findall(r'^\s*(xPrint|XPrint|YPrint|infPrint)\s*=\s*(\S+)',text,re.M))
    iters=re.findall(r'^\s*(\d+)\s+[+-]?(?:\d+\.)?\d+[eE][+-]\d+',stdout,re.M)
    arrays={name:(name in text) for name in ('xVec','xMat','yMat')}
    expected={key:'%+100.90Fe' for key in ('xPrint','XPrint','YPrint','infPrint')}
    if worker.returncode==0 and not stop and phase:
        if not reported or int(reported[1])!=precision or actual_formats!=expected:
            raise ValueError('Native precision or output formats do not match request')
    return dict(status=stop or ('NATIVE_TERMINAL_STATUS' if phase else 'NATIVE_ERROR_NO_TERMINAL_PHASE'),
        command=command,solver=solver,returncode=worker.returncode,phase=None if phase is None else phase[1],
        requested_precision_bits=precision,reported_precision_bits=None if reported is None else int(reported[1]),
        decimal_input_significant_digits=DIGITS,feasibility_tolerance='1e-45',gap_tolerance='1e-45',
        full_precision_output_formats=actual_formats,solution_arrays_present=arrays,
        last_printed_iteration=int(iters[-1]) if iters else None,elapsed_seconds=time.monotonic()-started,
        observed_owned_process_tree_rss_bytes=peak,scientific_watchdog_memory_limit_bytes=10_000_000_000,
        inner_owned_memory_stop_bytes=8_800_000_000,
        native_matrix_mode='SPARSE' if 'SPARSE computation' in stdout+text else ('DENSE' if 'DENSE computation' in stdout+text else None),
        native_result_sha256=solver_gate.file_digest(result_path) if result_path.is_file() else None,
        numerical_status_is_not_exact_certificate=True)


def solve(args):
    prepared=read(args.prepared/'result.json');verified=read(args.verified/'result.json')
    if (verified.get('status')!='EXACT_DUAL_GMP_ALL_NATIVE_COEFFICIENTS_REPLAYED'
            or verified.get('solver_ready') is not True
            or verified['prepared_result_sha256']!=solver_gate.file_digest(args.prepared/'result.json')
            or verified['native_input_sha256']!=solver_gate.file_digest(args.prepared/'input.dat-s')):
        raise ValueError('No matching complete multiprecision export replay')
    check_bindings(verified['source_bindings'])
    info=native(args.prepared/'input.dat-s',args.output,args.precision,args.iterations,args.seconds)
    conic.write_json(args.output/'result.json',dict(status='SDPA_GMP_FULL_EXPLICIT_DUAL_ATTEMPT_COMPLETED',
        prepared=str(args.prepared),separate_input_replay=str(args.verified/'result.json'),
        native=info,problem_role='Normalized dual certificate feasibility for reduced certificate search C',
        objective=prepared['objective'],original_feasibility_status='UNDECIDED',
        exact_certificate_obtained=False,geometric_constraints_changed=False,accepted_bounds_changed=False,
        source_bindings=bind([Path(__file__),args.prepared/'result.json',args.verified/'result.json'])))


def check(args):
    # Exercise this export's signs, off-diagonals and genuinely sub-binary64 input.
    k=F(1)+F(1,10**40)
    p=args.output/'input.dat-s'
    p.write_text('\n'.join(['"Global Coxeter SDP; sub-binary64 input precision fixture','1','2','2 -1','3',
        '0 1 1 2 '+decimal(-k),'1 1 1 1 1','1 1 2 2 1','1 2 1 1 1'])+'\n')
    fixture_gate=solver_gate.verify_native(p,None,fixture=True)
    info=native(p,args.output,args.precision,150,min(args.seconds,60))
    raw,phase=conic.native_vector(args.output/'native.result')
    if phase!='pdOPT' or len(raw)!=1 or abs(raw[0]-k)>F(1,10**44):
        raise ValueError('Multiprecision fixture did not recover the sub-binary64 shift')
    if abs(raw[0]-1)<F(1,2*10**40):raise ValueError('Native coefficient lost its 1e-40 shift')
    for text in ['1/3','-7/11','1/10000000000000000000000000000000000000001']:
        for exponent in (-28,0,28):
            q=rational(text,exponent);check_decimal(decimal(q),q)
    try:check_decimal('0',F(1,10**40))
    except ValueError:pass
    else:raise ValueError('Dropped tiny coefficient was accepted')
    conic.write_json(args.output/'result.json',dict(status='SDPA_GMP_EXACT_EXPORT_SUB_BINARY64_INPUT_FIXTURE_PASSED',
        expected_optimum=str(k),returned_optimum=str(raw[0]),error=str(raw[0]-k),
        native=info,input_shift='1e-40',required_absolute_accuracy='1e-44',
        fixture_gate=fixture_gate,
        exact_decimal_error_guard_and_nonzero_preservation_passed=True,
        research_problem_not_solved_by_fixture=True,source_bindings=bind([Path(__file__)])))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['check','prepare','verify','solve'])
    for name in ('model','replay','units','balanced','native-input','preflight','prepared','verified'):
        parser.add_argument('--'+name,type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--iterations',type=int,default=200)
    parser.add_argument('--seconds',type=int,default=780)
    args=parser.parse_args()
    if not os.environ.get('COXETER_CHECK_OWNER'):raise RuntimeError('Use the serial scientific watchdog')
    if args.precision not in (128,256,512) or not 1<=args.seconds<=800:
        raise ValueError('Unsupported precision or batch duration')
    args.output=args.output.resolve()
    if not args.output.is_relative_to(ROOT/'results') or args.output.exists():raise ValueError('Use a fresh results folder')
    args.output.mkdir();started=time.monotonic()
    globals()[args.action](args)
    result=read(args.output/'result.json')
    result['worker_elapsed_seconds']=time.monotonic()-started
    (args.output/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(dict(status=result['status'],output=str(args.output),elapsed_seconds=result['worker_elapsed_seconds'])),flush=True)


if __name__=='__main__':main()
