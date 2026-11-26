#!/usr/bin/env python3
"""Solve the explicit dual and reconstruct/check certificates over real fields."""
import argparse
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from research.global_formulation import conic, conic_backend as backend, balanced_dual, field_rounding as rounding, solver_gate
from scripts.solve_balanced_certificate_20261126 import solve, verified_balanced_model
from scripts.compare_conic_solvers_20261126 import safe_json


def checked_model(folder,replay):
    verified=json.loads(replay.read_text())
    if (verified.get('status')!='BALANCED_DUAL_PARAMETRIZATION_ALL_COEFFICIENTS_VERIFIED'
            or verified.get('solver_ready') is not True
            or verified['prepared_result_sha256']!=solver_gate.file_digest(folder/'result.json')):
        raise ValueError('The explicit dual has no matching completed separate replay')
    for path,row in verified['source_bindings'].items():
        if solver_gate.file_digest(path)!=row['sha256']:
            raise ValueError('A dual parametrization input or verifier source changed')
    for name,h in verified['encoding_file_sha256'].items():
        if solver_gate.file_digest(folder/name)!=h:raise ValueError('A dual model artifact changed')
    return verified


def exact_dual_identities(forms,blocks,point,field,rational,balanced,coordinates):
    """Check every reduced adjoint identity using exact trace multiplicities."""
    import numpy as np
    with gzip.open(rational/'substitutions.exact.json.gz','rt') as stream:subs=json.load(stream)
    free=np.load(rational/'free.npz',allow_pickle=False)['original_dual_free_coordinates']
    inverse={int(col):i for i,col in enumerate(free)}
    aa=[field.Q(0) for _ in free];bb=[field.Q(0) for _ in free]
    weights=[1]*len(forms)
    for block in blocks:
        n=block['size'];pos=block['start']
        for i in range(n):
            for j in range(i,n):weights[pos]=1 if i==j else 2;pos+=1
    for col,row in enumerate(forms):
        value=rounding.evaluate(row,point,field)
        if value==field.zero:continue
        a,b=field.scale(weights[col],value)
        if col in inverse:
            k=inverse[col];aa[k]+=a;bb[k]+=b
        else:
            for old,c in subs[str(col)].items():
                if int(old)<0:raise ValueError('Unexpected affine term in homogeneous certificate map')
                k=inverse[int(old)];c=field.Q(c);aa[k]+=c*a
                if b:bb[k]+=c*b
    units=json.loads((balanced/'units.exact.json').read_text())
    row={int(k):field.Q(v) for k,v in units['source_normalization'].items()}
    largest=max(abs(c) for col,c in row.items() if col>=0)
    for col,c in row.items():
        if col>=0:aa[inverse[col]]-=c/largest
    for k,(a,b) in enumerate(zip(aa,bb,strict=True)):
        if a or b:raise ValueError('Nonzero exact reduced adjoint identity '+str(k)+': '+str((a,b)))
    return dict(exact_reduced_adjoint_equations_checked=len(free),
        normalized_dual_pairing=-1,free_equality_multiplier=1,
        off_diagonal_trace_factor_two=True,
        source_reduced_certificate_problem_infeasible=True,
        original_necessary_data_primal_feasible_claimed=False)


def numeric_matrix(values,block):
    import numpy as np
    n=block['size'];M=np.zeros((n,n));pos=block['start']
    for i in range(n):
        for j in range(i,n):
            M[i,j]=M[j,i]=values[pos]/(2**.5 if i!=j else 1.);pos+=1
    return M


def proposals(forms,blocks,data,approximate,field,tolerance,progress):
    import numpy as np
    value=data['b']-data['A']@approximate
    bounds=np.abs(data['b'])+abs(data['A'])@np.maximum(1.,np.abs(approximate))
    equations=[];kernels=[];scalar=0;terms=0;limited=False
    for count,block in enumerate(blocks):
        if block['size']==1:
            pos=block['start']
            if abs(value[pos])<=tolerance*bounds[pos] and forms[pos]:
                if terms+len(forms[pos])<=8000000:
                    equations.append({col:(c,field.Q(0)) for col,c in forms[pos].items()});scalar+=1
                    terms+=len(forms[pos])
                else:limited=True
            continue
        matrix=numeric_matrix(value,block)
        vectors=rounding.kernel_vectors(matrix,field.d,tolerance)
        if vectors:
            proposed=rounding.add_kernel_equations(forms,block,vectors,field)
            added=sum(map(len,proposed))
            if terms+added<=8000000:
                kernels.append(dict(block=block['name'],vectors=vectors))
                equations.extend(proposed);terms+=added
            else:limited=True
        if count%10000==0:progress(stage='kernel_proposals',field=field.d,blocks=count,kernel_blocks=len(kernels))
    return equations,kernels,scalar,limited


def save_point(path,point,field):
    with gzip.open(path,'xt') as stream:
        json.dump(dict(schema='exact-real-quadratic-point-v1',radicand=field.d,
            embedding='sqrt(d)>0' if field.d else 'Q',values=[field.text(v) for v in point]),stream,separators=(',',':'))


def fixture(output):
    import numpy as np
    import scipy.sparse as sp
    f=rounding.Field(2)
    if f.sign(f.element(3,-2))!=1 or f.sign(f.element(1,-1))!=-1:
        raise ValueError('Incorrect exact quadratic sign')
    root=f.element(0,1)
    if f.mul(root,root)!=f.element(2) or f.div(f.one,root)!=f.element(0,'1/2'):
        raise ValueError('Incorrect exact quadratic arithmetic')
    rounding.psd([[f.one,root],[root,f.element(2)]],f)
    try:rounding.psd([[f.one,root],[root,f.one]],f)
    except ValueError:pass
    else:raise ValueError('An indefinite exact field matrix was accepted')
    forms=[{-1:f.Q(1)},{0:f.Q(1)},{1:f.Q(1)}]
    block=dict(name='irrational_rank_one',start=0,size=2)
    equations=rounding.add_kernel_equations(forms,block,[[(0,-1),(1,0)]],f)
    point,info=rounding.project(equations,[2**.5,2.],f)
    if point!=[root,f.element(2)] or not rounding.check_point(forms,[block],point,f)['valid']:
        raise ValueError('Exact algebraic kernel projection failed')
    detected=rounding.kernel_vectors([[1.,2**.5],[2**.5,2.]],2)
    if not detected or not any(b for v in detected for a,b in v):
        raise ValueError('LLL did not detect the known quadratic kernel')
    data=dict(A=sp.csc_matrix([[-1.],[1.]]),b=np.array([0.,1.]),c=np.zeros(1))
    layout=dict(zero_cones=0,scalar_cones=2,PSD_sizes=[])
    folder=output/'zero_equality_solver';folder.mkdir()
    _,native=solve(data,layout,folder,15,100,1e-11)
    if native['native_status']!='Solved':raise ValueError('Zero-equality dual solver fixture failed')
    return dict(status='EXACT_REAL_FIELD_KERNEL_PROJECTION_AND_DUAL_SOLVER_FIXTURES_VERIFIED',
        positive_embedding_sign_checked=True,indefinite_matrix_rejected=True,
        quadratic_LLL_kernel_detected=True,nonrational_point_reconstructed=True,
        all_proposed_kernel_equations_checked=True,zero_equality_solver_checked=True)


def main():
    import numpy as np
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['check','solve','attempt','verify'])
    parser.add_argument('--model',type=Path)
    parser.add_argument('--verified',type=Path)
    parser.add_argument('--balanced',type=Path)
    parser.add_argument('--rational',type=Path)
    parser.add_argument('--native-input',type=Path)
    parser.add_argument('--preflight',type=Path)
    parser.add_argument('--solver-output',type=Path)
    parser.add_argument('--candidate',type=Path)
    parser.add_argument('--seconds',type=int,default=650)
    parser.add_argument('--fields',default='0,2,5')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();started=time.monotonic()
    if not 1<=args.seconds<=750:raise ValueError('Keep the solver within the monitored batch')
    paths=[Path(__file__),ROOT/'research/global_formulation/field_rounding.py',
        ROOT/'research/global_formulation/balanced_dual.py',ROOT/'research/global_formulation/conic_backend.py',
        ROOT/'scripts/solve_balanced_certificate_20261126.py',ROOT/'research/global_formulation/conic.py',
        ROOT/'research/global_formulation/solver_gate.py']
    if args.action!='check':
        paths += [args.model/'result.json',args.model/'forms.exact.jsonl.gz',args.model/'coordinates.json',
            args.model/'blocks.json.gz',args.verified,args.balanced/'result.json',args.rational/'result.json']
        if args.action=='attempt':paths += [args.solver_output/'result.json',args.solver_output/'solution.npz']
        if args.action=='verify':paths += [args.candidate/'result.json',args.candidate/'field_point.exact.json.gz']
    def bindings():
        return {str(p):dict(bytes=p.stat().st_size,sha256=solver_gate.file_digest(p)) for p in paths}
    before=bindings();args.output.mkdir();conic.write_json(args.output/'source_bindings.json',before)
    def progress(**event):
        event['elapsed_seconds']=time.monotonic()-started
        with (args.output/'progress.jsonl').open('a') as stream:stream.write(json.dumps(event)+'\n')
    if args.action=='check':report=fixture(args.output)
    else:
        checked=checked_model(args.model,args.verified)
        verified_balanced_model(args.balanced,args.native_input,args.preflight)
        data,layout=backend.load_encoding(args.model)
        if args.action=='solve':
            solution,native=solve(data,layout,args.output,args.seconds,500,1e-11)
            np.savez_compressed(args.output/'solution.npz',**solution)
            report=dict(status='EXPLICIT_BALANCED_DUAL_CLARABEL_RETURNED_REQUIRES_FIELD_CHECK',
                native=native,diagnostics=backend.diagnostics(data,layout,solution),
                solution_sha256=solver_gate.file_digest(args.output/'solution.npz'),
                exact_certificate_obtained=False)
        else:
            meta=json.loads((args.model/'coordinates.json').read_text())
            forms=balanced_dual.read_forms(args.model/'forms.exact.jsonl.gz',meta['forms'])
            with gzip.open(args.model/'blocks.json.gz','rt') as stream:blocks=json.load(stream)
            if args.action=='verify':
                candidate=json.loads((args.candidate/'result.json').read_text())
                path=args.candidate/'field_point.exact.json.gz'
                if solver_gate.file_digest(path)!=candidate['field_point_sha256']:
                    raise ValueError('The field certificate changed')
                with gzip.open(path,'rt') as stream:raw=json.load(stream)
                field=rounding.Field(raw['radicand']);point=[field.element(a,b) for a,b in raw['values']]
                if len(point)!=len(meta['variables']):raise ValueError('Incomplete exact field point')
                test=rounding.check_point(forms,blocks,point,field,progress)
                if not test['valid']:raise ValueError('The field certificate fails a cone: '+str(test))
                identities=exact_dual_identities(forms,blocks,point,field,args.rational,args.balanced,meta)
                report=dict(status='EXACT_FIELD_CERTIFICATE_FOR_DUAL_OF_REDUCED_FEASIBILITY_VERIFIED',
                    field='Q' if not field.d else 'Q(sqrt('+str(field.d)+'))',
                    field_point_sha256=solver_gate.file_digest(path),cone_check=test,
                    exact_linear_identity_check=identities,exact_certificate_obtained=True)
            else:
                native=json.loads((args.solver_output/'result.json').read_text())
                approximate=np.load(args.solver_output/'solution.npz',allow_pickle=False)['x']
                if len(approximate)!=len(meta['variables']) or not np.isfinite(approximate).all():
                    raise ValueError('No complete finite numerical dual candidate')
                report=dict(status='FIELD_RECONSTRUCTION_ATTEMPTS_NO_EXACT_CERTIFICATE',
                    native_status=native['native']['native_status'],attempts=[],exact_certificate_obtained=False)
                accepted=None
                for denominator in [2**32,2**44,2**50]:
                    field=rounding.Field(0);point=rounding.rational_point(approximate,field,denominator)
                    test=rounding.check_point(forms,blocks,point,field,progress)
                    report['attempts'].append(dict(method='common_grid_candidate',common_denominator=denominator,check=test))
                    progress(stage='candidate_checked',method='common_grid',denominator=denominator,valid=test['valid'])
                    if test['valid']:accepted=(field,point,test);break
                if accepted is None:
                    for radicand in [int(v) for v in args.fields.split(',')]:
                        field=rounding.Field(radicand)
                        for tolerance in [1e-10,1e-12]:
                            equations,kernels,scalar,limited=proposals(forms,blocks,data,approximate,field,tolerance,progress)
                            attempt=dict(method='LLL_kernel_and_exact_affine_projection',radicand=radicand,
                                kernel_tolerance=tolerance,kernel_blocks=len(kernels),proposed_scalar_equalities=scalar,
                                candidate_proposal_term_limit_reached=limited,
                                nonrational_kernel_vectors=sum(any(b for a,b in v) for k in kernels for v in k['vectors']))
                            try:
                                point,info=rounding.project(equations,approximate,field,denominator=2**48)
                                test=rounding.check_point(forms,blocks,point,field,progress)
                                attempt.update(projection=info,check=test)
                            except ValueError as error:attempt['candidate_failure']=str(error);test={'valid':False}
                            report['attempts'].append(attempt);progress(stage='field_attempt_complete',**attempt)
                            if test['valid']:
                                accepted=(field,point,test)
                                conic.write_json(args.output/'candidate_kernel_proposals.json',kernels)
                                break
                        if accepted is not None:break
                if accepted is not None:
                    field,point,test=accepted
                    save_point(args.output/'field_point.exact.json.gz',point,field)
                    report.update(status='EXACT_FIELD_DUAL_CERTIFICATE_CANDIDATE_REQUIRES_SEPARATE_REPLAY',
                        field='Q' if not field.d else 'Q(sqrt('+str(field.d)+'))',
                        field_point_sha256=solver_gate.file_digest(args.output/'field_point.exact.json.gz'),
                        all_cones_checked=test,exact_certificate_obtained=False,separate_replay_still_required=True)
        report.update(model=str(args.model),separate_dual_parametrization_replay=str(args.verified),
            interpretation='An exact normalized dual certificate proves the reduced certificate search infeasible; it alone does not establish original primal feasibility.')
    if before!=bindings():raise ValueError('A field source, input or proof changed during the run')
    report.update(source_bindings=before,elapsed_seconds=time.monotonic()-started,
        accepted_bounds_changed=False,geometric_constraints_changed=False,original_feasibility_status='UNDECIDED')
    conic.write_json(args.output/'result.json',safe_json(report))
    print(json.dumps(safe_json({k:v for k,v in report.items() if k not in {'source_bindings','attempts'}})),flush=True)


if __name__=='__main__':main()
