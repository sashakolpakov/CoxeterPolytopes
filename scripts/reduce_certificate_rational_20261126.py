#!/usr/bin/env python3
"""Reduce certificate equalities over Q, then replay the reduction separately."""
import argparse
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from research.global_formulation import conic, conic_backend as backend, rational_certificate as rational, solver_gate
from scripts.compare_conic_solvers_20261126 import checked_dual_face_model


def save_exact(path, value):
    with gzip.open(path, 'xt') as stream:
        json.dump(conic.plain(value), stream, separators=(',', ':'), sort_keys=True)
        stream.write('\n')


def read_exact(path):
    with gzip.open(path, 'rt') as stream:
        value = json.load(stream)
    return {int(k): {int(c): F(v) for c, v in row.items()} for k, row in value.items()}


def fixture(output):
    # Both pivots must be substituted in reverse order; a non-leaf is rejected.
    rows = [rational.primitive({0:F(2,3), 1:F(1,7), 2:F(1,5)}),
            rational.primitive({1:F(3,11), 2:F(1,13), -1:F(-1,17)})]
    steps, remaining = rational.peel(rows, 2, cascade=True)
    if len(steps) != 2 or remaining:
        raise ValueError('The two-stage exact fixture did not eliminate both equations')
    rational.verify_peeling(rows, 2, steps, remaining)
    substitutions = rational.back_substitution(rows, steps)
    if rational.exact_remainders(rows, steps, remaining, substitutions):
        raise ValueError('An exact fixture equation remains')
    try:
        rational.verify_peeling(rows, 2, [dict(source_row=0,pivot=1)], [1])
    except ValueError:
        pass
    else:
        raise ValueError('A non-leaf pivot was accepted')
    data, layout, free = rational.encode(substitutions, {}, 3,
                                        dict(scalar_cones=3, PSD_sizes=[]))
    import numpy as np
    q = np.asarray([float(rational.substitute({i:F(1)}, substitutions).get(-1,0)
                          + rational.substitute({i:F(1)}, substitutions).get(2,0)*F(2,7))
                    for i in range(3)])
    np.testing.assert_allclose(data['b']-data['A']@np.asarray([2/7]), q, atol=1e-15)
    if free.tolist() != [2] or layout['zero_cones']:
        raise ValueError('Incorrect exact fixture free-coordinate map')
    independent_steps, independent_remaining = rational.peel(rows, 2)
    if len(independent_steps)!=1 or len(independent_remaining)!=1:
        raise ValueError('The default reduction cascaded instead of retaining the coupled core')
    # A genuinely coupled residual is retained, including its normalization.
    residual_rows=[{0:1,1:2},{0:3,1:4,-1:1}]
    s,r=rational.peel(residual_rows,2)
    if s or rational.exact_remainders(residual_rows,s,r,{})!={0:residual_rows[0],1:residual_rows[1]}:
        raise ValueError('A coupled residual was silently removed')
    # Actual matrix coordinates have rational trace coefficients, even though
    # the final numerical cone packing contains sqrt(2).
    P=conic.model(['x'],[
        dict(name='scalar',entries=[[{'x':F(2),'':F(-1)}]]),
        dict(name='matrix',entries=[
            [{'':F(4)},{'x':F(3,7)},{}],
            [{'x':F(3,7)},{'':F(5)},{}],
            [{},{},{'':F(6)}]])])
    original=dict(scalar_cones=1,PSD_sizes=[3])
    equations=rational.source_equations(P,original,np.arange(7))
    if (equations[0]!={0:F(2),2:F(6,7)}
            or equations[1]!={-1:F(1),0:F(-1),1:F(4),4:F(5),6:F(6)}):
        raise ValueError('The rational trace coefficients or normalization differ')
    packed,_,_=rational.encode({}, {}, 7, original)
    np.testing.assert_allclose(packed['b']-packed['A']@np.arange(1.,8.),
        [1.,2.,3.*np.sqrt(2.),4.*np.sqrt(2.),5.,6.*np.sqrt(2.),7.],rtol=0,atol=1e-14)
    return dict(status='RATIONAL_CERTIFICATE_ELIMINATION_FIXTURES_VERIFIED',
                primitive_integer_rows=True, triangular_back_substitution_checked=True,
                non_leaf_pivot_rejected=True, coupled_residual_retained=True,
                noncascading_default_retains_coupled_core=True,
                rational_trace_factor_two_and_numerical_svec_factor_checked=True,
                numerical_affine_encoding_checked=True)


def main():
    import numpy as np
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['check','prepare','verify','ray-diagnose'])
    parser.add_argument('--native-input',type=Path)
    parser.add_argument('--preflight',type=Path)
    parser.add_argument('--encoding',type=Path)
    parser.add_argument('--exact-model',type=Path)
    parser.add_argument('--dual-face-model',type=Path)
    parser.add_argument('--dual-face-replay',type=Path)
    parser.add_argument('--solver-output',type=Path)
    parser.add_argument('--prepared',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    started=time.monotonic()
    sources=[Path(__file__),ROOT/'research/global_formulation/rational_certificate.py',
             ROOT/'research/global_formulation/conic.py', ROOT/'research/global_formulation/conic_backend.py',
             ROOT/'scripts/compare_conic_solvers_20261126.py', ROOT/'research/global_formulation/solver_gate.py']
    if args.action!='check':
        sources += [args.exact_model,args.dual_face_replay,args.dual_face_model/'face_proof.json',
                    args.dual_face_model/'embedding.npz',args.encoding/'layout.json']
        if args.action=='verify':
            sources += [args.prepared/name for name in ['result.json','elimination.json',
                'substitutions.exact.json.gz','residual.exact.json.gz','A.npz','vectors.npz','layout.json','free.npz']]
        if args.action=='ray-diagnose':
            sources += [args.solver_output/'result.json',args.solver_output/'solution.npz',
                        ROOT/'scripts/check_reduced_certificate_result_20261126.py']
    def bindings():
        return {str(p):dict(bytes=p.stat().st_size,sha256=solver_gate.file_digest(p)) for p in sources}
    before=bindings()
    args.output.mkdir()
    conic.write_json(args.output/'source_bindings.json',before)
    if args.action=='check':
        report=fixture(args.output)
    else:
        gate=solver_gate.verify_native(args.native_input,args.preflight)
        checked_dual_face_model(args.dual_face_model,args.dual_face_replay,gate)
        face=json.loads((args.dual_face_model/'face_proof.json').read_text())
        if solver_gate.file_digest(args.exact_model)!=face['exact_model_sha256']:
            raise ValueError('A different rational source model was supplied')
        model=conic.read_model(args.exact_model)
        original_layout=json.loads((args.encoding/'layout.json').read_text())
        original_layout=dict(original_layout,scalar_cones=original_layout['scalar_cones']-1)
        maps=np.load(args.dual_face_model/'embedding.npz',allow_pickle=False)
        face_layout=json.loads((args.dual_face_model/'layout.json').read_text())
        if args.action=='ray-diagnose':
            report=diagnose_ray(model,original_layout,maps,face_layout,args.solver_output)
        else:
            equations=rational.source_equations(model,original_layout,maps['dual_coordinate_embedding'])
            del model
            source_heights=rational.heights(equations)
            rows=[rational.primitive(row) for row in equations]
            del equations
            integer_heights=rational.heights(rows)
            def progress(**event):
                event['elapsed_seconds']=time.monotonic()-started
                with (args.output/'progress.jsonl').open('a') as stream:
                    stream.write(json.dumps(event)+'\n')
            if args.action=='prepare':
                progress(stage='source_equations',equations=len(rows),source_heights=source_heights,
                         primitive_integer_heights=integer_heights)
                steps,remaining=rational.peel(rows,face_layout['scalar_cones'],progress)
            else:
                elimination=json.loads((args.prepared/'elimination.json').read_text())
                steps,remaining=elimination['steps'],elimination['remaining_source_rows']
            rational.verify_peeling(rows,face_layout['scalar_cones'],steps,remaining)
            substitutions=rational.back_substitution(rows,steps,progress=progress)
            residual=rational.exact_remainders(rows,steps,remaining,substitutions)
            data,layout,free=rational.encode(substitutions,residual,
                len(maps['dual_coordinate_embedding']),face_layout)
            if args.action=='prepare':
                save_exact(args.output/'substitutions.exact.json.gz',substitutions)
                save_exact(args.output/'residual.exact.json.gz',residual)
                conic.write_json(args.output/'elimination.json',dict(
                    schema='exact-degree-one-certificate-elimination-v1',steps=steps,
                    policy='Original degree-one scalar columns only; no cascading fill',
                    remaining_source_rows=remaining,source_model_sha256=face['exact_model_sha256'],
                    matrix_coordinates='Unscaled actual entries over Q',
                    certificate_feasibility_preserved=True))
                hashes=backend.save_encoding(args.output,data,layout)
                np.savez_compressed(args.output/'free.npz',original_dual_free_coordinates=free)
                for name in ['elimination.json','substitutions.exact.json.gz','residual.exact.json.gz','free.npz']:
                    hashes[name]=solver_gate.file_digest(args.output/name)
                status='RATIONAL_CERTIFICATE_REDUCTION_REQUIRES_SEPARATE_REPLAY'
            else:
                prepared=json.loads((args.prepared/'result.json').read_text())
                for name,h in prepared['encoding_file_sha256'].items():
                    if solver_gate.file_digest(args.prepared/name)!=h:
                        raise ValueError('A prepared rational artifact changed')
                if (read_exact(args.prepared/'substitutions.exact.json.gz')!=substitutions
                        or read_exact(args.prepared/'residual.exact.json.gz')!=residual):
                    raise ValueError('The stored rational equivalence map differs')
                saved,saved_layout=backend.load_encoding(args.prepared)
                if (saved_layout!=layout or saved['A'].shape!=data['A'].shape
                        or (saved['A']!=data['A']).nnz or not np.array_equal(saved['b'],data['b'])
                        or not np.array_equal(saved['c'],data['c'])
                        or not np.array_equal(np.load(args.prepared/'free.npz',allow_pickle=False)
                            ['original_dual_free_coordinates'],free)):
                    raise ValueError('A reduced numerical coefficient or coordinate map differs')
                hashes=prepared['encoding_file_sha256']
                status='RATIONAL_CERTIFICATE_REDUCTION_AND_ALL_EQUATIONS_VERIFIED'
            report=dict(status=status,solver_ready=args.action=='verify',
                source_native_sha256=gate['native_input_sha256'],
                exact_model_sha256=face['exact_model_sha256'],
                source_equations=len(rows),source_nonzero_equations=sum(bool(r) for r in rows),
                eliminated_equations=len(steps),remaining_equations=len(residual),
                source_certificate_variables=len(maps['dual_coordinate_embedding']),
                reduced_certificate_variables=len(free),matrix_shape=list(data['A'].shape),
                sparse_nonzeros=data['A'].nnz,source_rational_heights=source_heights,
                primitive_integer_heights=integer_heights,
                reconstruction_heights=rational.heights(substitutions.values()),
                residual_heights=rational.heights(residual.values()),
                encoding_file_sha256=hashes,certificate_feasibility_preserved=True,
                original_primal_feasibility_equivalence_claimed=False,
                all_source_equalities_verified_over_Q=True,
                separate_dual_face_replay=str(args.dual_face_replay),verified_preflight=gate)
            if args.action=='verify':
                report.update(prepared_model=str(args.prepared),
                    prepared_result_sha256=solver_gate.file_digest(args.prepared/'result.json'))
    if bindings()!=before:
        raise ValueError('An exact input or reduction source changed during execution')
    report.update(elapsed_seconds=time.monotonic()-started,source_bindings=before,
                  accepted_bounds_changed=False,geometric_constraints_changed=False,
                  original_feasibility_status='UNDECIDED')
    conic.write_json(args.output/'result.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in {'source_bindings','encoding_file_sha256'}}),flush=True)


def diagnose_ray(model,layout,maps,face_layout,solver_output):
    """Evaluate the actual binary64-derived point over Q; do not re-solve it."""
    from scripts.check_reduced_certificate_result_20261126 import projected_model
    arrays=__import__('numpy').load(solver_output/'solution.npz',allow_pickle=False)
    columns=maps['retained_ordinary_equation_columns']
    point=projected_model(model,layout,maps['dual_coordinate_embedding'],columns)
    ray=arrays['y']
    eta=F.from_float(float(ray[face_layout['zero_cones']-1]))
    values={name:-F.from_float(float(ray[k]))/eta for k,name in enumerate(point['variables'])}
    negative=[]
    cover=None
    for block in point['blocks']:
        if len(block['entries'])!=1:
            continue
        row=block['entries'][0][0]
        value=conic.value(row,values)
        if value<0:
            scale=sum(abs(c*(values[k] if k else 1)) for k,c in row.items())
            item=dict(name=block['name'],exact_value=str(value),value_float=float(value),
                absolute_term_sum_float=float(scale),
                relative_to_absolute_term_sum=float(abs(value)/scale) if scale else None,
                terms=len(row))
            if block['name']=='uniform_two_unit_triangle_cover0':
                cover=item
            negative.append(item)
    negative.sort(key=lambda item:item['value_float'])
    return dict(status='APPROXIMATE_RAY_EXACT_SCALAR_RESIDUALS_DIAGNOSED',
        native_status='PrimalInfeasible',native_status_is_not_NumericalError=True,
        exact_reconstruction_failure_is_not_disproof_of_native_classification=True,
        no_additional_precision_was_generated_by_decimal_grid_rounding=True,
        negative_scalar_count=len(negative),triangle_cover=cover,
        most_negative_scalar_rows=negative[:10],
        all_PSD_blocks_checked=False,original_point_lift_checked=False,
        exact_certificate_obtained=False)


if __name__=='__main__':
    main()
