#!/usr/bin/env python3
"""Validate and compare independent solvers on the unchanged gated SDP."""
import argparse
from contextlib import redirect_stderr, redirect_stdout
import importlib.metadata
import json
import math
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from research.global_formulation import conic, conic_backend as backend, slater_margin, solver_gate


def bindings(paths):
    return {str(p): {'bytes': p.stat().st_size, 'sha256': solver_gate.file_digest(p)} for p in paths}


def sources():
    return [Path(__file__), ROOT/'research/global_formulation/conic_backend.py',
            ROOT/'research/global_formulation/solver_gate.py', ROOT/'research/global_formulation/conic.py',
            ROOT/'research/global_formulation/feasibility.py',ROOT/'research/global_formulation/slater_margin.py']


def safe_json(value):
    if isinstance(value, dict):
        return {str(k):safe_json(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)):
        return [safe_json(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


def checked_encoding(folder, gate):
    encoded=json.loads((folder/'result.json').read_text())
    if encoded['source_native_sha256']!=gate['native_input_sha256']:
        raise ValueError('The comparison cache describes a different frozen input')
    for p,h in encoded['encoding_file_sha256'].items():
        if solver_gate.file_digest(folder/p)!=h:
            raise ValueError('An encoding cache changed: '+p)
    return backend.load_encoding(folder)


def solve(data, layout, output, name, seconds, iterations, warm=None):
    import numpy as np
    import scipy.sparse as sp
    log = output/'solver.stdout'
    if name == 'scs':
        import scs
        import scs._scs_indirect as implementation
        cone = {'l':layout['scalar_cones'], 's':layout['PSD_sizes']}
        if layout.get('zero_cones',0):
            cone['z'] = layout['zero_cones']
        settings = {'use_indirect':True, 'verbose':True, 'eps_abs':1e-8,
                    'eps_rel':1e-8, 'eps_infeas':1e-9, 'normalize':True,
                    'max_iters':iterations, 'time_limit_secs':seconds,
                    'log_csv_filename':str(output/'iterations.csv')}
        with log.open('x') as stream, redirect_stdout(stream), redirect_stderr(stream):
            solver = scs.SCS(data, cone, **settings)
            solution = solver.solve(warm_start=warm is not None, **(warm or {}))
        info = solution['info']
        metadata = {'backend':'SCS', 'version':scs.__version__,
                    'native_binary_sha256':solver_gate.file_digest(implementation.__file__),
                    'settings':settings, 'native_status':info['status'], 'native_info':info}
        return {k:np.asarray(solution[k]) for k in ['x','y','s']}, metadata
    import clarabel
    import clarabel.clarabel as implementation
    permutation = backend.upper_permutation(layout)
    cones = []
    if layout.get('zero_cones',0):
        cones.append(clarabel.ZeroConeT(layout['zero_cones']))
    if layout['scalar_cones']:
        cones.append(clarabel.NonnegativeConeT(layout['scalar_cones']))
    cones.extend(clarabel.PSDTriangleConeT(n) for n in layout['PSD_sizes'])
    settings = clarabel.DefaultSettings()
    settings.max_iter = iterations
    settings.time_limit = seconds
    settings.max_threads = 1
    settings.direct_solve_method = 'qdldl'
    settings.tol_gap_abs = settings.tol_gap_rel = settings.tol_feas = 1e-8
    settings.tol_infeas_abs = settings.tol_infeas_rel = 1e-9
    settings.verbose = True
    with log.open('x') as stream, redirect_stdout(stream), redirect_stderr(stream):
        solver = clarabel.DefaultSolver(sp.csc_matrix((len(data['c']),len(data['c']))),
            data['c'],data['A'][permutation].tocsc(),data['b'][permutation],cones,settings)
        result = solver.solve()
    dual = np.empty(len(permutation))
    slack = np.empty(len(permutation))
    dual[permutation] = result.z
    slack[permutation] = result.s
    metadata = {'backend':'Clarabel', 'version':clarabel.__version__,
                'native_binary_sha256':solver_gate.file_digest(implementation.__file__),
                'settings':str(settings), 'native_status':str(result.status),
                'reported_iterations':result.iterations,
                'reported_primal_residual':result.r_prim,
                'reported_dual_residual':result.r_dual,
                'reported_primal_objective':result.obj_val,
                'reported_dual_objective':result.obj_val_dual,
                'upper_to_lower_svec_permutation_applied_to_both_cones':True,
                'warm_start_used':False}
    return {'x':np.asarray(result.x),'y':dual,'s':slack}, metadata


def self_check(output):
    import numpy as np
    checked = []
    # Nonuniform 3x3 data distinguish lower-column and upper-column layouts.
    P = conic.model(['x','y'], [
        {'name':'asymmetric_positions', 'entries':[
            [conic.form(4),conic.form(x=1),conic.form(y=1)],
            [conic.form(x=1),conic.form(5),conic.form(x=1,y=-2)],
            [conic.form(y=1),conic.form(x=1,y=-2),conic.form(6)]]},
        {'name':'x_lower','entries':[[conic.form(-1,x=1)]]},
        {'name':'x_upper','entries':[[conic.form(1,x=-1)]]},
        {'name':'y_lower','entries':[[conic.form(conic.F(-1,4),y=1)]]},
        {'name':'y_upper','entries':[[conic.form(conic.F(1,4),y=-1)]]}])
    fixture = output/'fixture.dat-s'
    slater_margin.write_input(fixture,slater_margin.augmented(P))
    solver_gate.verify_native(fixture,None,fixture=True)
    data,layout = backend.read_native(fixture)
    point = np.array([1.,.25,0.])
    affine = data['b']-data['A']@point
    values = np.array([[4.,1.,.25],[1.,5.,.5],[.25,.5,6.]])
    i,j = np.triu_indices(3)
    expected = values[i,j].copy()
    expected[i!=j] *= np.sqrt(2.)
    np.testing.assert_allclose(affine[layout['scalar_cones']:],expected,atol=1e-14,rtol=0)
    Z = np.array([[2.,.125,.5],[.125,3.,.75],[.5,.75,4.]])
    packed_Z = Z[i,j].copy()
    packed_Z[i!=j] *= np.sqrt(2.)
    np.testing.assert_allclose(expected@packed_Z,np.trace(values@Z),atol=1e-14,rtol=0)
    perm = backend.upper_permutation(layout)
    upper_i,upper_j = np.tril_indices(3)
    upper = values[upper_i,upper_j].copy()
    upper[upper_i!=upper_j] *= np.sqrt(2.)
    np.testing.assert_allclose(affine[perm][layout['scalar_cones']:],upper,atol=1e-14,rtol=0)
    original,original_layout = backend.original_feasibility(data,layout)
    np.testing.assert_allclose(original['b']-original['A']@point[:-1],
                              np.delete(affine,layout['scalar_cones']-1),atol=1e-14,rtol=0)
    for name in ['scs','clarabel']:
        folder = output/name
        folder.mkdir()
        candidate,meta = solve(data,layout,folder,name,15,20000 if name=='scs' else 100)
        if meta['native_status'].lower() != 'solved':
            raise ValueError('Known singular margin fixture was not solved: '+str(meta))
        np.testing.assert_allclose(candidate['x'],point,atol=2e-6,rtol=0)
        diag = backend.diagnostics(data,layout,candidate)
        if diag['primal_affine_residual_max']>2e-6 or diag['dual_adjoint_residual_max']>2e-6:
            raise ValueError('Incorrect primal or dual sign/layout')
        checked.append({'backend':name,'native_status':meta['native_status'],'diagnostics':diag})
    # An inconsistent interval has an elementary exact Farkas vector.
    infeasible = conic.model(['x'],[
        {'name':'lower','entries':[[conic.form(-2,x=1)]]},
        {'name':'upper','entries':[[conic.form(1,x=-1)]]}])
    fixture2 = output/'infeasible.dat-s'
    slater_margin.write_input(fixture2,slater_margin.augmented(infeasible))
    data2,layout2 = backend.read_native(fixture2)
    direct,direct_layout = backend.original_feasibility(data2,layout2)
    certificate,certificate_layout = backend.certificate_feasibility(data2,layout2)
    for name in ['scs','clarabel']:
        folder = output/(name+'_infeasible')
        folder.mkdir()
        _,meta = solve(direct,direct_layout,folder,name,15,20000 if name=='scs' else 100)
        if 'infeasible' not in meta['native_status'].lower():
            raise ValueError('The inconsistent interval was not detected')
        folder2 = output/(name+'_certificate')
        folder2.mkdir()
        candidate,meta2 = solve(certificate,certificate_layout,folder2,name,15,20000 if name=='scs' else 100)
        np.testing.assert_allclose(candidate['x'],[1.,1.],rtol=0,atol=2e-6)
        checked.append({'backend':name,'infeasible_fixture_status':meta['native_status'],
                        'normalized_certificate_fixture_status':meta2['native_status'],
                        'normalized_certificate':candidate['x'].tolist()})
    return {'status':'INDEPENDENT_CONIC_ENCODINGS_AND_SOLVERS_VERIFIED',
            'offdiagonal_trace_multiplicity_checked':True,
            'nonuniform_3_by_3_triangle_order_checked':True,
            'known_singular_margin_optimum_checked':True,
            'original_margin_slice_checked':True,
            'normalized_Farkas_alternative_checked':True,'checks':checked}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['check','encode','solve','diagnose'])
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--native-input',type=Path)
    parser.add_argument('--preflight',type=Path)
    parser.add_argument('--reference',type=Path)
    parser.add_argument('--reference-diagnostics',type=Path)
    parser.add_argument('--encoding',type=Path)
    parser.add_argument('--backend',choices=['scs','clarabel'])
    parser.add_argument('--mode',choices=['auxiliary','original','certificate'],default='auxiliary')
    parser.add_argument('--seconds',type=int,default=650)
    parser.add_argument('--iterations',type=int)
    args = parser.parse_args()
    if not 1<=args.seconds<=750:
        raise ValueError('Keep time for encoding and independent residual checks')
    args.output.mkdir()
    started=time.monotonic()
    source_paths=sources()
    gate=None
    if args.action!='check':
        if args.native_input is None or args.preflight is None:
            raise ValueError('All full-model operations require the original numerical preflight')
        gate=solver_gate.verify_native(args.native_input,args.preflight)
        source_paths.extend([args.native_input,args.preflight])
        source_paths.extend(p for p in [args.reference,args.reference_diagnostics] if p is not None)
        if args.encoding is not None:
            source_paths.extend([args.encoding/p for p in
                ['result.json','A.npz','vectors.npz','layout.json']])
    before=bindings(source_paths)
    conic.write_json(args.output/'source_bindings.json',before)
    if args.action=='check':
        report=self_check(args.output)
    elif args.action=='encode':
        import numpy as np
        import scipy.sparse.csgraph as graph
        from scipy import sparse
        data,layout=backend.read_native(args.native_input)
        hashes=backend.save_encoding(args.output,data,layout)
        column_norms=np.sqrt(np.asarray(data['A'].power(2).sum(axis=0)).ravel())
        report={'status':'COMPLETE_SPARSE_CONIC_ENCODING_CREATED', 'layout':layout,
            'encoding_file_sha256':hashes, 'source_native_sha256':gate['native_input_sha256'],
            'structural_column_rank':int(graph.structural_rank(data['A'])),
            'structural_rank_is_not_a_numerical_condition_estimate':True,
            'minimum_column_two_norm':float(column_norms.min()),
            'maximum_column_two_norm':float(column_norms.max()),
            'minimum_nonzero_absolute_coefficient':float(np.min(np.abs(data['A'].data))),
            'maximum_absolute_coefficient':float(np.max(np.abs(data['A'].data))),
            'uncompressed_sparse_A_bytes':sum(v.nbytes for v in [data['A'].data,data['A'].indices,data['A'].indptr]),
            'packed_PSD_Hessian_triangle_entries':sum((n*(n+1)//2)*(n*(n+1)//2+1)//2 for n in layout['PSD_sizes'])}
        if args.reference:
            point=backend.native_arrays(args.reference,layout)
            diag=backend.diagnostics(data,layout,point)
            conic.write_json(args.output/'native_crosscheck.json',diag)
            np.savez_compressed(args.output/'native_candidate.npz',**point)
            hashes['native_candidate.npz']=solver_gate.file_digest(args.output/'native_candidate.npz')
            if args.reference_diagnostics:
                previous=json.loads(args.reference_diagnostics.read_text())
                np.testing.assert_allclose(diag['primal_affine_residual_max'],previous['recomputed_primal']['maximum_absolute_affine_residual'],atol=2e-11,rtol=0)
                np.testing.assert_allclose(diag['dual_adjoint_residual_max'],previous['recomputed_dual']['ordinary_adjoint_maximum_absolute_residual'],atol=2e-12,rtol=0)
                np.testing.assert_allclose(diag['gap'],previous['recomputed_gap'],atol=2e-11,rtol=0)
                report['complete_existing_iterate_primal_dual_residual_crosscheck_passed']=True
            report['native_diagnostics']=diag
    elif args.action=='diagnose':
        import numpy as np
        if args.encoding is None or args.reference is None:
            raise ValueError('Supply the checked encoding and native checkpoint')
        data,layout=checked_encoding(args.encoding,gate)
        point=backend.native_arrays(args.reference,layout)
        auxiliary=backend.diagnostics(data,layout,point)
        original,original_layout=backend.original_feasibility(data,layout)
        x=point['x'][:-1]
        y=np.delete(point['y'],layout['scalar_cones']-1)
        original_point={'x':x,'y':y,'s':original['b']-original['A']@x}
        report={'status':'COMMON_MARGIN_AND_ORIGINAL_FEASIBILITY_DIAGNOSED',
                'auxiliary':auxiliary,'original':backend.diagnostics(original,original_layout,original_point),
                'source_checkpoint':str(args.reference),'original_feasibility_status':'UNDECIDED'}
    else:
        import numpy as np
        if args.encoding is None or args.backend is None:
            raise ValueError('Supply the checked encoding and backend')
        data,layout=checked_encoding(args.encoding,gate)
        if args.mode=='original':
            data,layout=backend.original_feasibility(data,layout)
        elif args.mode=='certificate':
            data,layout=backend.certificate_feasibility(data,layout)
        warm=None
        if args.backend=='scs' and args.mode=='auxiliary':
            if args.reference:
                warm=backend.native_arrays(args.reference,layout)
            elif (args.encoding/'native_candidate.npz').exists():
                saved=np.load(args.encoding/'native_candidate.npz',allow_pickle=False)
                warm={k:saved[k] for k in ['x','y','s']}
        options_iterations=args.iterations or (10_000_000 if args.backend=='scs' else 200)
        conic.write_json(args.output/'problem.json',{'shape':list(data['A'].shape),
            'nonzeros':data['A'].nnz,'mode':args.mode,'layout':layout,
            'warm_start_from_verified_SDPA_arrays':warm is not None})
        point,meta=solve(data,layout,args.output,args.backend,args.seconds,options_iterations,warm=warm)
        np.savez_compressed(args.output/'solution.npz',**point)
        diag=backend.diagnostics(data,layout,point)
        report={'status':'ALTERNATIVE_SOLVER_RETURNED_REQUIRES_ORIGINAL_CHECK',
            'mode':args.mode,'native':meta,'diagnostics':diag,
            'original_feasibility_status':'UNDECIDED',
            'solution_sha256':solver_gate.file_digest(args.output/'solution.npz'),
            'source_encoding':str(args.encoding)}
    if bindings(source_paths)!=before:
        raise ValueError('A pinned input or source changed during comparison')
    report.update(accepted_bounds_changed=False,geometric_constraints_changed=False,
        source_bound_geometric_review_complete=False,verified_preflight=gate,
        elapsed_seconds=time.monotonic()-started,source_bindings=before,
        installed_versions={p:importlib.metadata.version(p) for p in ['numpy','scipy','scs','clarabel']})
    conic.write_json(args.output/'result.json',safe_json(report))
    print(json.dumps(safe_json({k:v for k,v in report.items()
        if k not in ['source_bindings','checks','layout']})),flush=True)


if __name__=='__main__':
    main()
