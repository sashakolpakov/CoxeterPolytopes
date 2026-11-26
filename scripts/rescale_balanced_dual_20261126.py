#!/usr/bin/env python3
"""Exact dyadic cone units for the separately replayed explicit dual SDP."""
import argparse
from contextlib import redirect_stdout,redirect_stderr
import json
import math
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from research.global_formulation import conic,conic_backend as backend,solver_gate
from scripts.compare_conic_solvers_20261126 import safe_json
from scripts.solve_balanced_certificate_20261126 import verified_balanced_model


def checked_dual(folder,replay):
    result=json.loads(replay.read_text())
    if (result.get('status')!='BALANCED_DUAL_PARAMETRIZATION_ALL_COEFFICIENTS_VERIFIED'
            or result['prepared_result_sha256']!=solver_gate.file_digest(folder/'result.json')):
        raise ValueError('Missing separate exact dual replay')
    for path,row in result['source_bindings'].items():
        if solver_gate.file_digest(path)!=row['sha256']:raise ValueError('Changed dual proof source or input')
    for name,h in result['encoding_file_sha256'].items():
        if solver_gate.file_digest(folder/name)!=h:raise ValueError('Changed dual encoding')
    return result


def rescale(data,layout):
    import numpy as np
    if layout['zero_cones']:raise ValueError('This cone rescaling handles the explicit dual with no equalities')
    A=data['A'].tocsr();envelope=np.maximum(np.asarray(abs(A).max(axis=1).toarray()).ravel(),np.abs(data['b']))
    exponents=np.zeros(A.shape[0],dtype=np.int64)
    diagonals=[[0]*n for n in layout['PSD_sizes']]
    def scalar_exponent(value):return 0 if not value else -math.ceil(math.log2(value))
    def matrix_row_envelopes(values):
        all_rows=list(values[:layout['scalar_cones']]);base=layout['scalar_cones']
        for n in layout['PSD_sizes']:
            i,j=np.triu_indices(n);M=np.zeros((n,n));M[i,j]=values[base:base+len(i)];M[j,i]=M[i,j]
            all_rows.extend(np.max(M,axis=1));base+=len(i)
        return np.asarray(all_rows)
    transformed=A.copy();rhs=data['b'].copy();rounds=0
    for turn in range(12):
        current=np.maximum(np.asarray(abs(transformed).max(axis=1).toarray()).ravel(),np.abs(rhs))
        physical=matrix_row_envelopes(current);positive=physical[physical>0]
        if np.min(positive)>=.5 and np.max(positive)<=2:break
        increment=np.zeros(A.shape[0],dtype=np.int64)
        for r in range(layout['scalar_cones']):increment[r]=scalar_exponent(current[r])
        base=layout['scalar_cones']
        for block,n in enumerate(layout['PSD_sizes']):
            i,j=np.triu_indices(n);M=np.zeros((n,n));M[i,j]=current[base:base+len(i)];M[j,i]=M[i,j]
            maximum=np.max(M,axis=1)
            d=np.array([0 if not v else -round(.5*math.log2(v)) for v in maximum],dtype=np.int64)
            diagonals[block]=[a+int(b) for a,b in zip(diagonals[block],d,strict=True)]
            increment[base:base+len(i)]=d[i]+d[j];base+=len(i)
        exponents+=increment;rounds+=1
        transformed.data=np.ldexp(A.data,np.repeat(exponents,np.diff(A.indptr)))
        rhs=np.ldexp(data['b'],exponents)
    if base!=len(exponents) or np.max(np.abs(exponents))>900:raise ValueError('Unsafe cone unit exponent')
    units=np.ldexp(np.ones(len(exponents)),exponents)
    entry_exponents=np.repeat(exponents,np.diff(A.indptr))
    if (not np.isfinite(transformed.data).all() or not np.isfinite(rhs).all()
            or not np.array_equal(np.ldexp(transformed.data,-entry_exponents),A.data)
            or not np.array_equal(np.ldexp(rhs,-exponents),data['b'])):
        raise ValueError('The numerical dyadic scaling has no exact byte-value inverse')
    # Scalar units are positive. Each PSD entry uses d_i+d_j, exactly the
    # exponent of the positive diagonal congruence diag(2^d_i) M diag(2^d_i).
    if not np.all(units>0):raise ValueError('A cone unit is not positive')
    after=np.maximum(np.asarray(abs(transformed).max(axis=1).toarray()).ravel(),np.abs(rhs))
    physical_after=matrix_row_envelopes(after)
    if np.min(physical_after[physical_after>0])<.5 or np.max(physical_after)>2:
        raise ValueError('Cone row units have not reached the numerical readiness range')
    proof=dict(scalar_positive_units_verified=True,PSD_positive_diagonal_congruences_verified=True,
        dyadic_numeric_inverse_verified_for_every_entry=True,coefficients_removed_by_tolerance=0,
        scalar_cone_exponents=exponents[:layout['scalar_cones']].tolist(),PSD_diagonal_exponents=diagonals,
        smallest_input_positive_row_envelope=float(np.min(envelope[envelope>0])),
        smallest_output_positive_row_envelope=float(np.min(after[after>0])),
        maximum_output_coefficient=float(np.max(np.abs(transformed.data))),
        minimum_row_exponent=int(np.min(exponents)),maximum_row_exponent=int(np.max(exponents)))
    proof.update(congruence_rounds=rounds,
        minimum_positive_scalar_or_PSD_matrix_row_envelope=float(np.min(physical_after[physical_after>0])),
        maximum_scalar_or_PSD_matrix_row_envelope=float(np.max(physical_after)),
        PSD_row_maxima_used_instead_of_diagonals_only=True)
    return dict(A=transformed.tocsc(),b=rhs,c=data['c'].copy()),exponents,proof


def native_solve(data,layout,output,seconds):
    import clarabel
    import numpy as np
    import scipy.sparse as sp
    permutation=backend.upper_permutation(layout)
    cones=[clarabel.NonnegativeConeT(layout['scalar_cones'])]
    cones += [clarabel.PSDTriangleConeT(n) for n in layout['PSD_sizes']]
    settings=clarabel.DefaultSettings();settings.max_threads=1;settings.direct_solve_method='faer'
    settings.max_iter=500;settings.time_limit=seconds
    settings.tol_feas=settings.tol_gap_abs=settings.tol_gap_rel=1e-11
    settings.tol_infeas_abs=settings.tol_infeas_rel=1e-12
    settings.chordal_decomposition_enable=False
    def progress(info):
        record={k:getattr(info,k,None) for k in ['iterations','res_primal','res_dual','gap_abs','cost_dual','mu','solve_time']}
        with (output/'progress.jsonl').open('a') as stream:stream.write(json.dumps(safe_json(record))+'\n')
        return False
    with (output/'solver.stdout').open('x') as stream,redirect_stdout(stream),redirect_stderr(stream):
        solver=clarabel.DefaultSolver(sp.csc_matrix((len(data['c']),len(data['c']))),data['c'],
            data['A'][permutation].tocsc(),data['b'][permutation],cones,settings)
        solver.set_termination_callback(progress);result=solver.solve()
    y=np.empty(len(permutation));s=np.empty(len(permutation));y[permutation]=result.z;s[permutation]=result.s
    return dict(x=np.asarray(result.x),y=y,s=s),dict(native_status=str(result.status),
        reported_iterations=result.iterations,reported_primal_residual=result.r_prim,
        reported_dual_residual=result.r_dual,arithmetic='IEEE binary64',significand_bits=53,
        settings=str(settings),chordal_decomposition=False)


def main():
    import numpy as np
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['check','export','solve'])
    parser.add_argument('--model',type=Path);parser.add_argument('--verified',type=Path)
    parser.add_argument('--balanced',type=Path);parser.add_argument('--native-input',type=Path)
    parser.add_argument('--preflight',type=Path);parser.add_argument('--conditioned',type=Path)
    parser.add_argument('--seconds',type=int,default=750)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();started=time.monotonic()
    if not 1<=args.seconds<=800:raise ValueError('Keep room inside the 900-second monitored batch')
    paths=[Path(__file__),ROOT/'research/global_formulation/conic_backend.py',ROOT/'research/global_formulation/conic.py',
        ROOT/'research/global_formulation/solver_gate.py',ROOT/'scripts/solve_balanced_certificate_20261126.py']
    if args.action!='check':paths += [args.model/'result.json',args.verified,args.balanced/'result.json']
    if args.action=='solve':paths += [args.conditioned/'result.json',args.conditioned/'units.npz']
    def bindings():return {str(p):dict(bytes=p.stat().st_size,sha256=solver_gate.file_digest(p)) for p in paths}
    before=bindings();args.output.mkdir();conic.write_json(args.output/'source_bindings.json',before)
    if args.action=='check':
        import scipy.sparse as sp
        data=dict(A=sp.csc_matrix([[-1e-14,0.],[-1e-20,0.],[0.,-2**.5*1e-10],[0.,-1.]]),b=np.zeros(4),c=np.zeros(2))
        layout=dict(zero_cones=0,scalar_cones=1,PSD_sizes=[2])
        out,e,p=rescale(data,layout)
        if p['smallest_output_positive_row_envelope']<.1:raise ValueError('Known weak cone was not balanced')
        np.testing.assert_array_equal(np.ldexp(out['b'],-e),data['b'])
        # The same native configuration must solve a small PSD feasibility problem.
        folder=args.output/'native';folder.mkdir()
        small=dict(A=sp.csc_matrix([[-1.],[0.],[0.],[0.]]),b=np.array([0.,1.,0.,1.]),c=np.zeros(1))
        _,native=native_solve(small,layout,folder,15)
        if native['native_status']!='Solved':raise ValueError('The full-PSD native configuration failed')
        report=dict(status='DYADIC_CONE_CONGRUENCE_AND_NATIVE_FIXTURES_VERIFIED')
    else:
        verified_balanced_model(args.balanced,args.native_input,args.preflight)
        checked=checked_dual(args.model,args.verified);source,layout=backend.load_encoding(args.model)
        if args.action=='export':
            data,exponents,proof=rescale(source,layout)
            hashes=backend.save_encoding(args.output,data,layout)
            np.savez_compressed(args.output/'units.npz',cone_row_exponents=exponents)
            hashes['units.npz']=solver_gate.file_digest(args.output/'units.npz')
            conic.write_json(args.output/'congruence.exact.json',proof)
            hashes['congruence.exact.json']=solver_gate.file_digest(args.output/'congruence.exact.json')
            report=dict(status='DYADIC_EXPLICIT_DUAL_CONE_UNITS_VERIFIED_SOLVER_READY',solver_ready=True,
                encoding_file_sha256=hashes,source_model_result_sha256=solver_gate.file_digest(args.model/'result.json'),
                original_model=str(args.model),separate_dual_replay=str(args.verified),
                summary={k:v for k,v in proof.items() if k not in ['scalar_cone_exponents','PSD_diagonal_exponents']})
        else:
            receipt=json.loads((args.conditioned/'result.json').read_text())
            if receipt['status']!='DYADIC_EXPLICIT_DUAL_CONE_UNITS_VERIFIED_SOLVER_READY':raise ValueError('No cone unit proof')
            for path,row in receipt['source_bindings'].items():
                if solver_gate.file_digest(path)!=row['sha256']:raise ValueError('Changed cone-unit source')
            for name,h in receipt['encoding_file_sha256'].items():
                if solver_gate.file_digest(args.conditioned/name)!=h:raise ValueError('Changed scaled input')
            data,saved_layout=backend.load_encoding(args.conditioned)
            expected,exponents,_=rescale(source,layout)
            if (saved_layout!=layout or (data['A']!=expected['A']).nnz or not np.array_equal(data['b'],expected['b'])):
                raise ValueError('A cone-scaled coefficient differs')
            solution,native=native_solve(data,layout,args.output,args.seconds)
            np.savez_compressed(args.output/'scaled_solution.npz',**solution)
            scaled_diagnostics=backend.diagnostics(data,layout,solution)
            # Preserve solution.npz in original G units for all existing field
            # reconstruction paths; the native scaled arrays remain separately.
            solution['y']=np.ldexp(solution['y'],exponents)
            solution['s']=np.ldexp(solution['s'],-exponents)
            np.savez_compressed(args.output/'solution.npz',**solution)
            report=dict(status='CONE_BALANCED_EXPLICIT_DUAL_CLARABEL_RETURNED_REQUIRES_EXACT_CHECK',
                native=native,conditioned=str(args.conditioned),
                scaled_diagnostics=scaled_diagnostics,original_unit_diagnostics=backend.diagnostics(source,layout,solution),
                solution_sha256=solver_gate.file_digest(args.output/'solution.npz'),
                saved_solution_units='Original G units; native scaled arrays in scaled_solution.npz',exact_certificate_obtained=False)
    if before!=bindings():raise ValueError('A cone-scaling input or source changed')
    report.update(source_bindings=before,elapsed_seconds=time.monotonic()-started,
        accepted_bounds_changed=False,geometric_constraints_changed=False,original_feasibility_status='UNDECIDED')
    conic.write_json(args.output/'result.json',safe_json(report))
    print(json.dumps(safe_json({k:v for k,v in report.items() if k not in ['source_bindings','encoding_file_sha256']})),flush=True)


if __name__=='__main__':main()
