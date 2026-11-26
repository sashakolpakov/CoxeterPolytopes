#!/usr/bin/env python3
"""Run Clarabel on the exactly reduced, positively balanced certificate SDP."""
import argparse
from contextlib import redirect_stdout, redirect_stderr
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from research.global_formulation import conic, conic_backend as backend, solver_gate
from scripts.compare_conic_solvers_20261126 import safe_json


def verified_balanced_model(folder,native_input,preflight):
    gate=solver_gate.verify_native(native_input,preflight)
    receipt=json.loads((folder/'result.json').read_text())
    if (receipt.get('status')!='RATIONAL_CERTIFICATE_POSITIVE_UNITS_VERIFIED_SOLVER_READY'
            or receipt.get('solver_ready') is not True
            or receipt['source_native_sha256']!=gate['native_input_sha256']):
        raise ValueError('Missing or mismatched balanced-model proof')
    for path,row in receipt['source_bindings'].items():
        if solver_gate.file_digest(path)!=row['sha256']:
            raise ValueError('A balanced export input or source changed')
    replay=json.loads(Path(receipt['separate_rational_replay']).read_text())
    if replay.get('status')!='RATIONAL_CERTIFICATE_REDUCTION_AND_ALL_EQUATIONS_VERIFIED':
        raise ValueError('No completed separate rational replay')
    for path,row in replay['source_bindings'].items():
        if solver_gate.file_digest(path)!=row['sha256']:
            raise ValueError('An exact reduction input or verifier source changed')
    for name,h in receipt['encoding_file_sha256'].items():
        if solver_gate.file_digest(folder/name)!=h:
            raise ValueError('A balanced numerical artifact changed')
    data,layout=backend.load_encoding(folder)
    if layout['zero_cones']!=1 or layout['mode']!='rationally_reduced_normalized_certificate_feasibility':
        raise ValueError('Unexpected normalized certificate layout')
    return data,layout,gate,receipt


def solve(data,layout,output,seconds,iterations,tolerance):
    import clarabel
    import clarabel.clarabel as implementation
    import numpy as np
    import scipy.sparse as sp
    permutation=backend.upper_permutation(layout)
    cones=[clarabel.ZeroConeT(layout['zero_cones'])]
    if layout['scalar_cones']:
        cones.append(clarabel.NonnegativeConeT(layout['scalar_cones']))
    cones.extend(clarabel.PSDTriangleConeT(n) for n in layout['PSD_sizes'])
    settings=clarabel.DefaultSettings()
    settings.max_threads=1
    settings.direct_solve_method='faer'
    settings.time_limit=seconds
    settings.max_iter=iterations
    settings.tol_gap_abs=settings.tol_gap_rel=settings.tol_feas=tolerance
    settings.tol_infeas_abs=settings.tol_infeas_rel=tolerance/10
    settings.verbose=True
    events=output/'progress.jsonl'
    def progress(info):
        row={k:getattr(info,k,None) for k in ['iterations','res_primal','res_dual',
            'gap_abs','gap_rel','cost_primal','cost_dual','mu','ktratio','solve_time']}
        with events.open('a') as stream:
            stream.write(json.dumps(safe_json(row))+'\n')
        return False
    with (output/'solver.stdout').open('x') as stream,redirect_stdout(stream),redirect_stderr(stream):
        worker=clarabel.DefaultSolver(sp.csc_matrix((len(data['c']),len(data['c']))),data['c'],
            data['A'][permutation].tocsc(),data['b'][permutation],cones,settings)
        worker.set_termination_callback(progress)
        result=worker.solve()
    dual=np.empty(len(permutation));slack=np.empty(len(permutation))
    dual[permutation]=result.z
    slack[permutation]=result.s
    return dict(x=np.asarray(result.x),y=dual,s=slack),dict(
        backend='Clarabel',linear_solver='faer',version=clarabel.__version__,
        arithmetic='IEEE binary64',significand_bits=53,settings=str(settings),
        native_binary_sha256=solver_gate.file_digest(implementation.__file__),
        native_status=str(result.status),reported_iterations=result.iterations,
        reported_primal_residual=result.r_prim,reported_dual_residual=result.r_dual,
        reported_primal_objective=result.obj_val,reported_dual_objective=result.obj_val_dual,
        warm_start_used=False)


def main():
    import numpy as np
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['check','solve'])
    parser.add_argument('--model',type=Path)
    parser.add_argument('--native-input',type=Path)
    parser.add_argument('--preflight',type=Path)
    parser.add_argument('--seconds',type=int,default=650)
    parser.add_argument('--iterations',type=int,default=500)
    parser.add_argument('--tolerance',type=float,default=1e-11)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if not 1<=args.seconds<=750 or not 1e-13<=args.tolerance<=1e-8:
        raise ValueError('Unsupported batch duration or binary64 tolerance')
    started=time.monotonic()
    paths=[Path(__file__),ROOT/'research/global_formulation/conic_backend.py',
        ROOT/'research/global_formulation/solver_gate.py',ROOT/'research/global_formulation/conic.py',
        ROOT/'scripts/compare_conic_solvers_20261126.py']
    if args.action=='solve':paths += [args.model/'result.json',args.model/'units.exact.json']
    def bindings():
        return {str(p):dict(bytes=p.stat().st_size,sha256=solver_gate.file_digest(p)) for p in paths}
    before=bindings()
    args.output.mkdir()
    conic.write_json(args.output/'source_bindings.json',before)
    if args.action=='check':
        import scipy.sparse as sp
        layout=dict(zero_cones=1,scalar_cones=1,PSD_sizes=[])
        checks=[]
        for name,rhs in [('feasible',1.),('infeasible',-1.)]:
            folder=args.output/name;folder.mkdir()
            data=dict(A=sp.csc_matrix([[1.],[-1.]]),b=np.array([rhs,0.]),c=np.zeros(1))
            solution,native=solve(data,layout,folder,15,100,args.tolerance)
            expected='Solved' if rhs>0 else 'PrimalInfeasible'
            if native['native_status']!=expected:
                raise ValueError('Incorrect normalized certificate fixture status')
            if rhs<0:
                np.testing.assert_allclose(solution['y']/solution['y'][0],[1.,1.],atol=1e-9)
                np.testing.assert_allclose(data['A'].T@solution['y'],[0.],atol=1e-9)
            checks.append(dict(name=name,native_status=expected))
        report=dict(status='BALANCED_CERTIFICATE_CLARABEL_FIXTURES_VERIFIED',checks=checks)
    else:
        data,layout,gate,receipt=verified_balanced_model(args.model,args.native_input,args.preflight)
        conic.write_json(args.output/'launch.json',dict(model=str(args.model),
            shape=list(data['A'].shape),nonzeros=data['A'].nnz,verified_preflight=gate,
            requested_tolerance=args.tolerance,requested_infeasibility_tolerance=args.tolerance/10,
            time_limit_seconds=args.seconds,iteration_limit=args.iterations,
            arithmetic='IEEE binary64; tighter tolerances do not add arithmetic precision'))
        solution,native=solve(data,layout,args.output,args.seconds,args.iterations,args.tolerance)
        np.savez_compressed(args.output/'solution.npz',**solution)
        report=dict(status='BALANCED_CERTIFICATE_CLARABEL_RETURNED_REQUIRES_EXACT_CHECK',
            native=native,model=str(args.model),model_result_sha256=solver_gate.file_digest(args.model/'result.json'),
            solution_sha256=solver_gate.file_digest(args.output/'solution.npz'),
            diagnostics=backend.diagnostics(data,layout,solution),verified_preflight=gate,
            interpretation='Primal is the reduced certificate search; its infeasibility ray is a dual-certificate candidate.',
            exact_certificate_obtained=False)
    if bindings()!=before:raise ValueError('A solver input or source changed during execution')
    report.update(elapsed_seconds=time.monotonic()-started,source_bindings=before,
        accepted_bounds_changed=False,geometric_constraints_changed=False,original_feasibility_status='UNDECIDED')
    conic.write_json(args.output/'result.json',safe_json(report))
    print(json.dumps(safe_json({k:v for k,v in report.items() if k!='source_bindings'})),flush=True)


if __name__=='__main__':main()
