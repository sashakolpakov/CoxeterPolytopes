#!/usr/bin/env python3
"""Prepare or separately replay exact coordinate facial reduction of the dual."""
import argparse
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from research.global_formulation import conic, conic_backend as backend, dual_faces, slater_margin, solver_gate


def equal_matrix(left,right):
    import numpy as np
    left=left.tocsc();right=right.tocsc()
    return (left.shape==right.shape and np.array_equal(left.indptr,right.indptr)
            and np.array_equal(left.indices,right.indices) and np.array_equal(left.data,right.data))


def fixture(output):
    import numpy as np
    P=conic.model(['u','v'],[
        {'name':'moment','entries':[[conic.form(1),conic.form(u=1)],
                                  [conic.form(u=1),conic.form(v=1)]]},
        {'name':'u_positive','entries':[[conic.form(u=1)]]},
        {'name':'v_positive','entries':[[conic.form(v=1)]]}])
    path=output/'fixture.dat-s'
    slater_margin.write_input(path,slater_margin.augmented(P))
    solver_gate.verify_native(path,None,fixture=True)
    data,layout=backend.read_native(path)
    data,layout=backend.original_feasibility(data,layout)
    active=np.ones(data['A'].shape[0],dtype=bool)
    try: dual_faces.exact_step(P,layout,active,{0:1})
    except ValueError: pass
    else: raise ValueError('An indefinite coefficient was accepted before its face was proved')
    steps=[]
    for _ in range(4):
        rays=dual_faces.discover(data,layout,active)
        if not rays: break
        active,step=dual_faces.exact_step(P,layout,active,rays)
        steps.append(step)
    if len(steps)!=2 or not np.array_equal(active,dual_faces.check_steps(P,layout,steps)):
        raise ValueError('The known two-round conditional facial sequence failed')
    target,target_layout,mapping=dual_faces.reduced_certificate(data,layout,active)
    if target['A'].shape!=(2,1) or target_layout['scalar_cones']!=1 or target_layout['zero_cones']!=1:
        raise ValueError('The fixture certificate cone was not reduced to the expected scalar')
    np.testing.assert_array_equal(target['A'].toarray(),[[1.],[-1.]])
    np.testing.assert_array_equal(target['b'],[-1.,0.])
    import copy
    tampered=copy.deepcopy(steps)
    tampered[0]['support_after_sha256']='0'*64
    try: dual_faces.check_steps(P,layout,tampered)
    except ValueError: pass
    else: raise ValueError('A tampered facial proof was accepted')
    conic.write_json(output/'fixture_proof.json',{'steps':steps,'layout':target_layout})
    return {'status':'EXACT_DUAL_FACIAL_REDUCTION_FIXTURE_VERIFIED',
            'conditional_rounds':2,'unproved_initial_ray_rejected':True,
            'tampered_face_hash_rejected':True,'normalized_certificate_encoding_checked':True,
            'original_feasibility_not_inferred_from_reduced_certificate_infeasibility':True}


def main():
    import numpy as np
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['check','prepare','verify'])
    parser.add_argument('--native-input',type=Path)
    parser.add_argument('--preflight',type=Path)
    parser.add_argument('--encoding',type=Path)
    parser.add_argument('--exact-model',type=Path)
    parser.add_argument('--prepared',type=Path)
    parser.add_argument('--max-rounds',type=int,default=32)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    started=time.monotonic()
    args.output.mkdir()
    source_paths=[Path(__file__),ROOT/'research/global_formulation/dual_faces.py',
                  ROOT/'research/global_formulation/conic_backend.py',ROOT/'research/global_formulation/solver_gate.py',
                  ROOT/'research/global_formulation/conic.py',ROOT/'research/global_formulation/slater_margin.py']
    source_paths.extend(p for p in [args.native_input,args.preflight,args.exact_model] if p is not None)
    if args.encoding is not None:
        source_paths.extend(args.encoding/p for p in ['result.json','A.npz','vectors.npz','layout.json'])
    if args.prepared is not None:
        source_paths.extend(args.prepared/p for p in ['result.json','face_proof.json','A.npz',
            'vectors.npz','layout.json','embedding.npz'])
    def bindings():
        return {str(p):{'bytes':p.stat().st_size,'sha256':solver_gate.file_digest(p)} for p in source_paths}
    before=bindings()
    conic.write_json(args.output/'source_bindings.json',before)
    gate=None
    if args.action=='check':
        report=fixture(args.output)
    else:
        if any(p is None for p in [args.native_input,args.preflight,args.encoding,args.exact_model]):
            raise ValueError('Require the frozen preflight, exact model and numerical encoding')
        gate=solver_gate.verify_native(args.native_input,args.preflight)
        preflight=json.loads(args.preflight.read_text())
        exact_path=args.exact_model.resolve()
        relative=str(exact_path.relative_to(ROOT))
        exact_sha=solver_gate.file_digest(exact_path)
        if preflight['bound_artifacts'].get(relative)!=exact_sha:
            raise ValueError('The exact model is not bound by the original numerical preflight')
        encoded=json.loads((args.encoding/'result.json').read_text())
        if encoded['source_native_sha256']!=gate['native_input_sha256']:
            raise ValueError('Different numerical encoding')
        for p,h in encoded['encoding_file_sha256'].items():
            if solver_gate.file_digest(args.encoding/p)!=h: raise ValueError('Changed encoding')
        data,layout=backend.load_encoding(args.encoding)
        data,layout=backend.original_feasibility(data,layout)
        with exact_path.open() as stream: problem=json.load(stream)
        if len(problem['variables'])!=data['A'].shape[1]: raise ValueError('Different source variable ordering')
        if args.action=='prepare':
            active=np.ones(data['A'].shape[0],dtype=bool)
            steps=[]
            fixed=False
            with (args.output/'progress.jsonl').open('x') as progress:
                for count in range(args.max_rounds):
                    rays=dual_faces.discover(data,layout,active)
                    if not rays:
                        fixed=True
                        break
                    active,step=dual_faces.exact_step(problem,layout,active,rays)
                    steps.append(step)
                    progress.write(json.dumps({'round':count+1,'coordinate_directions':len(rays),
                        'dual_coordinates_removed':step['dual_coordinates_removed'],
                        'elapsed_seconds':time.monotonic()-started})+'\n')
                    progress.flush()
            proof={'schema':'exact-dual-coordinate-faces-v1','source_native_sha256':gate['native_input_sha256'],
                'exact_model_sha256':exact_sha,'steps':steps,'coordinate_face_fixed_point':fixed,
                'final_dual_support_sha256':dual_faces.mask_digest(active),
                'certificate_feasibility_preserved':True,'original_primal_feasibility_equivalence_claimed':False,
                'minimal_face_or_strict_feasibility_after_reduction_claimed':False}
            conic.write_json(args.output/'face_proof.json',proof)
            reduced,reduced_layout,mapping=dual_faces.reduced_certificate(data,layout,active)
            hashes=backend.save_encoding(args.output,reduced,reduced_layout)
            np.savez_compressed(args.output/'embedding.npz',**mapping)
            hashes['embedding.npz']=solver_gate.file_digest(args.output/'embedding.npz')
            hashes['face_proof.json']=solver_gate.file_digest(args.output/'face_proof.json')
            report={'status':'EXACT_DUAL_FACE_MODEL_REQUIRES_SEPARATE_REPLAY' if fixed else 'DUAL_FACE_CHECKPOINT_NOT_FIXED_POINT',
                'solver_ready':False,'coordinate_face_fixed_point':fixed,'rounds':len(steps),
                'coordinate_directions_verified':sum(len(s['rays_0_based']) for s in steps),
                'original_dual_coordinates':data['A'].shape[0], 'retained_dual_coordinates':len(mapping['dual_coordinate_embedding']),
                'original_ordinary_equations':data['A'].shape[1],
                'retained_ordinary_equations':len(mapping['retained_ordinary_equation_columns']),
                'certificate_matrix_shape':list(reduced['A'].shape),'encoding_file_sha256':hashes,
                'source_native_sha256':gate['native_input_sha256'],'exact_model_sha256':exact_sha}
        else:
            if args.prepared is None: raise ValueError('Supply the dual face model to replay')
            recorded=json.loads((args.prepared/'result.json').read_text())
            for p,h in recorded['encoding_file_sha256'].items():
                if solver_gate.file_digest(args.prepared/p)!=h: raise ValueError('Changed prepared dual artifact')
            proof=json.loads((args.prepared/'face_proof.json').read_text())
            if proof['source_native_sha256']!=gate['native_input_sha256'] or proof['exact_model_sha256']!=exact_sha:
                raise ValueError('Proof belongs to different input')
            active=dual_faces.check_steps(problem,layout,proof['steps'])
            if dual_faces.mask_digest(active)!=proof['final_dual_support_sha256']:
                raise ValueError('Wrong final dual support')
            if not proof['coordinate_face_fixed_point'] or dual_faces.discover(data,layout,active):
                raise ValueError('The declared coordinate facial fixed point was not reached')
            expected,expected_layout,mapping=dual_faces.reduced_certificate(data,layout,active)
            saved,saved_layout=backend.load_encoding(args.prepared)
            if saved_layout!=expected_layout or not equal_matrix(saved['A'],expected['A']):
                raise ValueError('A compressed certificate coefficient differs')
            for k in ['b','c']:
                if not np.array_equal(saved[k],expected[k]): raise ValueError('A compressed vector differs')
            cached=np.load(args.prepared/'embedding.npz',allow_pickle=False)
            if any(not np.array_equal(v,cached[k]) for k,v in mapping.items()):
                raise ValueError('A dual cone embedding differs')
            report={'status':'EXACT_DUAL_FACIAL_REDUCTION_AND_ALL_REDUCED_COEFFICIENTS_VERIFIED',
                'solver_ready':True,'source_native_sha256':gate['native_input_sha256'],
                'exact_model_sha256':exact_sha,'prepared_model':str(args.prepared),
                'prepared_result_sha256':solver_gate.file_digest(args.prepared/'result.json'),
                'encoding_file_sha256':recorded['encoding_file_sha256'],
                'rounds_replayed':len(proof['steps']),'coordinate_face_fixed_point':True,
                'certificate_feasibility_preserved':True,
                'original_primal_feasibility_from_reduced_certificate_infeasibility_claimed':False}
    if bindings()!=before:
        raise ValueError('A dual face input or proof source changed during the run')
    report.update(accepted_bounds_changed=False,geometric_constraints_changed=False,
        verified_preflight=gate,elapsed_seconds=time.monotonic()-started,
        source_bindings=before)
    conic.write_json(args.output/'result.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='source_bindings'}),flush=True)


if __name__=='__main__':
    main()
