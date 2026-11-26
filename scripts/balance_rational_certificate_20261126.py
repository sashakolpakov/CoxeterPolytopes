#!/usr/bin/env python3
"""Balance the remaining normalization after separately verified Q reduction.

For a homogeneous cone T x in K and c+a*x=0, use x=tau*z with
tau=abs(c)/max(abs(a)), divide the whole cone by tau, and divide the
normalization by abs(c). No cone coefficient changes and the normalization
has unit right side and maximum coefficient one. The exact map is retained.
"""
import argparse
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from research.global_formulation import conic, conic_backend as backend, solver_gate


def balance(data,layout,row,free):
    import numpy as np
    if layout['zero_cones']!=1 or np.any(data['b'][1:]!=0):
        raise ValueError('This normalization rule requires one equality and a homogeneous cone')
    constant=F(row[-1])
    largest=max(abs(F(c)) for col,c in row.items() if col>=0)
    if not constant or not largest:
        raise ValueError('No nonzero normalization to balance')
    tau=abs(constant)/largest
    A=data['A'].tocsr(copy=True)
    a,b=A.indptr[:2]
    if {int(free[k]) for k in A.indices[a:b]}!={col for col in row if col>=0}:
        raise ValueError('Normalization coordinate support differs from its exact source')
    for index in range(a,b):
        coefficient=F(row[int(free[A.indices[index]])])
        if A.data[index]!=float(coefficient):
            raise ValueError('A raw normalization coefficient differs')
        if coefficient*tau/abs(constant)!=coefficient/largest:
            raise ValueError('The exact normalization identity failed')
        A.data[index]=float(coefficient/largest)
    if data['b'][0]!=-float(constant):
        raise ValueError('Different normalization right side')
    rhs=data['b'].copy()
    rhs[0]=-float(constant/abs(constant))
    if (A[1:]!=data['A'][1:]).nnz:
        raise ValueError('A homogeneous cone coefficient changed')
    proof=dict(global_free_variable_scale=str(tau),
        free_variable_reconstruction='x_original = global_free_variable_scale * z_balanced',
        whole_cone_divisor=str(tau),normalization_divisor=str(abs(constant)),
        source_normalization={str(k):str(v) for k,v in row.items()},
        exact_positive_units_verified=True,homogeneous_cone_coefficients_unchanged=True,
        certificate_feasibility_preserved=True)
    return dict(A=A.tocsc(),b=rhs,c=data['c'].copy()),proof


def main():
    import numpy as np
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['check','export'])
    parser.add_argument('--prepared',type=Path)
    parser.add_argument('--verified',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    started=time.monotonic()
    sources=[Path(__file__),ROOT/'research/global_formulation/conic_backend.py',
             ROOT/'research/global_formulation/conic.py',ROOT/'research/global_formulation/solver_gate.py']
    if args.action=='export':
        sources += [args.verified]+[args.prepared/name for name in
            ['result.json','A.npz','vectors.npz','layout.json','free.npz',
             'substitutions.exact.json.gz','residual.exact.json.gz','elimination.json']]
    def bindings():
        return {str(p):dict(bytes=p.stat().st_size,sha256=solver_gate.file_digest(p)) for p in sources}
    before=bindings()
    args.output.mkdir()
    conic.write_json(args.output/'source_bindings.json',before)
    if args.action=='check':
        import scipy.sparse as sp
        data=dict(A=sp.csc_matrix([[2.,-1.],[-1.,0.],[0.,-1.]]),
                  b=np.array([-100.,0.,0.]),c=np.zeros(2))
        layout=dict(zero_cones=1,scalar_cones=2,PSD_sizes=[])
        balanced,proof=balance(data,layout,{-1:100,0:2,1:-1},np.arange(2))
        np.testing.assert_array_equal(balanced['A'].toarray(),[[1.,-.5],[-1.,0.],[0.,-1.]])
        np.testing.assert_array_equal(balanced['b'],[-1.,0.,0.])
        if proof['global_free_variable_scale']!='50':
            raise ValueError('The known normalization scale is wrong')
        translated=dict(data,b=np.array([-100.,1.,0.]))
        try: balance(translated,layout,{-1:100,0:2,1:-1},np.arange(2))
        except ValueError: pass
        else: raise ValueError('An unproved homogeneous-cone assumption was accepted')
        report=dict(status='RATIONAL_NORMALIZATION_BALANCE_FIXTURES_VERIFIED',
                    nonhomogeneous_cone_rejected=True,exact_reconstruction_checked=True)
    else:
        replay=json.loads(args.verified.read_text())
        if (replay.get('status')!='RATIONAL_CERTIFICATE_REDUCTION_AND_ALL_EQUATIONS_VERIFIED'
                or replay.get('solver_ready') is not True
                or replay['prepared_result_sha256']!=solver_gate.file_digest(args.prepared/'result.json')):
            raise ValueError('Missing or stale separate rational reduction replay')
        for path,row in replay['source_bindings'].items():
            if solver_gate.file_digest(path)!=row['sha256']:
                raise ValueError('A rational proof input or verifier source changed')
        for name,h in replay['encoding_file_sha256'].items():
            if solver_gate.file_digest(args.prepared/name)!=h:
                raise ValueError('A rational model artifact changed')
        with gzip.open(args.prepared/'substitutions.exact.json.gz','rt') as stream:
            expressions=json.load(stream)
        if any(F(row.get('-1','0')) for row in expressions.values()):
            raise ValueError('Exact cone reconstruction is not homogeneous')
        del expressions
        with gzip.open(args.prepared/'residual.exact.json.gz','rt') as stream:
            residual=json.load(stream)
        if len(residual)!=1:
            raise ValueError('Supply a uniform multi-equation balancing argument')
        row={int(k):F(v) for k,v in next(iter(residual.values())).items()}
        free=np.load(args.prepared/'free.npz',allow_pickle=False)['original_dual_free_coordinates']
        data,layout=backend.load_encoding(args.prepared)
        balanced,proof=balance(data,layout,row,free)
        hashes=backend.save_encoding(args.output,balanced,layout)
        conic.write_json(args.output/'units.exact.json',proof)
        hashes['units.exact.json']=solver_gate.file_digest(args.output/'units.exact.json')
        report=dict(status='RATIONAL_CERTIFICATE_POSITIVE_UNITS_VERIFIED_SOLVER_READY',
            solver_ready=True,encoding_file_sha256=hashes,
            separate_rational_replay=str(args.verified),
            source_native_sha256=replay['source_native_sha256'],
            exact_model_sha256=replay['exact_model_sha256'],
            matrix_shape=list(balanced['A'].shape),sparse_nonzeros=balanced['A'].nnz,
            exact_units=proof,normalization_right_side=float(balanced['b'][0]),
            maximum_absolute_normalization_coefficient=float(np.max(np.abs(balanced['A'].getrow(0).data))),
            certificate_feasibility_preserved=True)
    if bindings()!=before:
        raise ValueError('A normalization input or source changed during execution')
    report.update(source_bindings=before,elapsed_seconds=time.monotonic()-started,
        accepted_bounds_changed=False,geometric_constraints_changed=False,
        original_feasibility_status='UNDECIDED')
    conic.write_json(args.output/'result.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in {'source_bindings','encoding_file_sha256'}}),flush=True)


if __name__=='__main__':
    main()
