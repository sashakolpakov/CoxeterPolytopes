#!/usr/bin/env python3
"""Prepare/replay the exact dual parametrization for the reduced feasibility SDP."""
import argparse
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from research.global_formulation import conic, conic_backend as backend, balanced_dual as dual, solver_gate
from scripts.solve_balanced_certificate_20261126 import verified_balanced_model


def fixture():
    import numpy as np
    P=conic.model(['x','y'],[
        dict(name='x_unit',entries=[[{'x':F(2)}]]),
        dict(name='y_unit',entries=[[{'y':F(3)}]]),
        dict(name='matrix',entries=[[{'':F(1)},{'x':F(1)}],
                                    [{'x':F(1)},{'y':F(1)}]])])
    forms,blocks,meta=dual.build(P,dict(scalar_cones=2,PSD_sizes=[2]),np.arange(5),
        dict(scalar_cones=2,PSD_sizes=[2]),[dict(source_row=0,pivot=0),dict(source_row=1,pivot=1)],F(5))
    if forms!=[{0:F(1)},{1:F(1)},{-1:F(5)},{0:F(1,2)},{1:F(1,3)}]:
        raise ValueError('Wrong dual trace/coordinate parametrization')
    # T has rows -q01, -q11/3, q00, q01, q11.
    u0,u1=F(1),F(3)
    h=[u0,u1,F(5),u0,u1/3]
    if [h[2],-h[0]+h[3],-h[1]/3+h[4]]!=[F(5),F(0),F(0)]:
        raise ValueError('An exact normalized dual identity failed')
    data,layout=dual.encode(forms,meta)
    np.testing.assert_allclose(data['b']-data['A']@np.array([1.,3.]),
        [1.,3.,5.,np.sqrt(2.)/2,1.],rtol=0,atol=1e-15)
    return dict(status='BALANCED_DUAL_PARAMETRIZATION_FIXTURE_VERIFIED',
        exact_dual_identity_checked=True,trace_factor_two_checked=True,
        scalar_coordinate_inverse_checked=True,numerical_export_checked=True)


def main():
    import numpy as np
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['check','prepare','verify'])
    parser.add_argument('--balanced',type=Path)
    parser.add_argument('--rational',type=Path)
    parser.add_argument('--face',type=Path)
    parser.add_argument('--encoding',type=Path)
    parser.add_argument('--exact-model',type=Path)
    parser.add_argument('--native-input',type=Path)
    parser.add_argument('--preflight',type=Path)
    parser.add_argument('--prepared',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();started=time.monotonic()
    paths=[Path(__file__),ROOT/'research/global_formulation/balanced_dual.py',
        ROOT/'research/global_formulation/conic.py',ROOT/'research/global_formulation/conic_backend.py',
        ROOT/'scripts/solve_balanced_certificate_20261126.py']
    if args.action!='check':
        paths += [args.balanced/'result.json',args.balanced/'units.exact.json',
            args.rational/'elimination.json',args.rational/'result.json',args.face/'embedding.npz',
            args.face/'layout.json',args.encoding/'layout.json',args.exact_model]
        if args.action=='verify':paths += [args.prepared/name for name in
            ['result.json','forms.exact.jsonl.gz','blocks.json.gz','coordinates.json','A.npz','vectors.npz','layout.json']]
    def bindings():
        return {str(p):dict(bytes=p.stat().st_size,sha256=solver_gate.file_digest(p)) for p in paths}
    before=bindings();args.output.mkdir();conic.write_json(args.output/'source_bindings.json',before)
    if args.action=='check':report=fixture()
    else:
        _,_,gate,balanced=verified_balanced_model(args.balanced,args.native_input,args.preflight)
        if solver_gate.file_digest(args.exact_model)!=balanced['exact_model_sha256']:
            raise ValueError('Different exact primal source')
        proof=json.loads((args.rational/'elimination.json').read_text())
        if proof['source_model_sha256']!=balanced['exact_model_sha256']:
            raise ValueError('Different exact elimination source')
        # This binds the selected elimination proof to the actual separate replay.
        replay=json.loads(Path(balanced['separate_rational_replay']).read_text())
        if replay['encoding_file_sha256']['elimination.json']!=solver_gate.file_digest(args.rational/'elimination.json'):
            raise ValueError('Elimination witness differs from the replayed one')
        layout=json.loads((args.encoding/'layout.json').read_text())
        layout=dict(layout,scalar_cones=layout['scalar_cones']-1)
        face_layout=json.loads((args.face/'layout.json').read_text())
        mapping=np.load(args.face/'embedding.npz',allow_pickle=False)
        model=conic.read_model(args.exact_model)
        tau=F(balanced['exact_units']['global_free_variable_scale'])
        forms,blocks,meta=dual.build(model,layout,mapping['dual_coordinate_embedding'],face_layout,proof['steps'],tau)
        del model
        if meta['source_variable_indices']!=mapping['retained_ordinary_equation_columns'].tolist():
            raise ValueError('The dual parametrization omits a retained source variable')
        data,numeric_layout=dual.encode(forms,meta)
        if args.action=='prepare':
            dual.save_forms(args.output/'forms.exact.jsonl.gz',forms)
            with gzip.open(args.output/'blocks.json.gz','xt') as stream:json.dump(blocks,stream,separators=(',',':'))
            conic.write_json(args.output/'coordinates.json',meta)
            hashes=backend.save_encoding(args.output,data,numeric_layout)
            for name in ['forms.exact.jsonl.gz','blocks.json.gz','coordinates.json']:
                hashes[name]=solver_gate.file_digest(args.output/name)
            status='BALANCED_DUAL_MODEL_REQUIRES_SEPARATE_REPLAY'
        else:
            recorded=json.loads((args.prepared/'result.json').read_text())
            for name,h in recorded['encoding_file_sha256'].items():
                if solver_gate.file_digest(args.prepared/name)!=h:raise ValueError('Changed dual artifact')
            with gzip.open(args.prepared/'forms.exact.jsonl.gz','rt') as stream:
                count=0
                for expected,line in zip(forms,stream,strict=True):
                    found={int(k):F(v) for k,v in json.loads(line)}
                    if expected!=found:raise ValueError('A dual source coefficient differs')
                    count+=1
            with gzip.open(args.prepared/'blocks.json.gz','rt') as stream:
                if json.load(stream)!=blocks:raise ValueError('A dual source cone differs')
            if json.loads((args.prepared/'coordinates.json').read_text())!=meta:
                raise ValueError('The dual coordinate inverse differs')
            saved,saved_layout=backend.load_encoding(args.prepared)
            if (saved_layout!=numeric_layout or (saved['A']!=data['A']).nnz
                    or not np.array_equal(saved['b'],data['b']) or not np.array_equal(saved['c'],data['c'])):
                raise ValueError('A dual numerical coefficient differs')
            hashes=recorded['encoding_file_sha256'];status='BALANCED_DUAL_PARAMETRIZATION_ALL_COEFFICIENTS_VERIFIED'
        report=dict(status=status,solver_ready=args.action=='verify',
            exact_dual_linear_identities_verified=True,normalized_dual_pairing=-1,
            dual_cone_membership_still_required=True,variables=len(meta['variables']),
            scalar_cones=meta['scalar_cones'],PSD_blocks=len(meta['PSD_sizes']),
            matrix_shape=list(data['A'].shape),sparse_nonzeros=data['A'].nnz,
            encoding_file_sha256=hashes,source_native_sha256=gate['native_input_sha256'],
            exact_model_sha256=balanced['exact_model_sha256'],verified_preflight=gate,
            balanced_model=str(args.balanced),
            original_primal_feasibility_equivalence_claimed=False)
        if args.action=='verify':report.update(prepared_model=str(args.prepared),
            prepared_result_sha256=solver_gate.file_digest(args.prepared/'result.json'))
    if before!=bindings():raise ValueError('A dual proof input or source changed during execution')
    report.update(source_bindings=before,elapsed_seconds=time.monotonic()-started,
        accepted_bounds_changed=False,geometric_constraints_changed=False,original_feasibility_status='UNDECIDED')
    conic.write_json(args.output/'result.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in {'source_bindings','encoding_file_sha256'}}),flush=True)


if __name__=='__main__':main()
