#!/usr/bin/env python3
"""Find exact coordinate recession rays exposing faces of the certificate cone.

If F_j (or -F_j) is nonzero diagonal PSD in every source block, then every
PSD certificate with zero ordinary adjoint is orthogonal to that ray. Its
positive diagonal positions, and the corresponding PSD rows, must vanish.
This is a dual facial identity, not a new geometric condition.
"""
import argparse
import csv
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from research.global_formulation import conic, conic_backend as backend, slater_margin, solver_gate


def candidates(data,layout):
    import numpy as np
    original,original_layout=backend.original_feasibility(data,layout)
    A=original['A']
    diagonal=np.zeros(A.shape[0],dtype=bool)
    diagonal[:original_layout['scalar_cones']]=True
    position=original_layout['scalar_cones']
    for n in original_layout['PSD_sizes']:
        i=np.arange(n)
        diagonal[position+n*i-i*(i-1)//2]=True
        position+=n*(n+1)//2
    answer={}
    for col in range(A.shape[1]):
        start,end=A.indptr[col:col+2]
        if start==end:
            continue
        if not np.all(diagonal[A.indices[start:end]]):
            continue
        values=A.data[start:end]
        if np.all(values<=0):
            answer[col]=1
        elif np.all(values>=0):
            answer[col]=-1
    return answer


def verify(native,rays,output):
    forced={}
    counts={col:0 for col in rays}
    with gzip.open(output/'ray_coefficients.csv.gz','wt',newline='') as stream:
        writer=csv.writer(stream)
        writer.writerow(['native_variable_1_based','ray_sign','native_block_1_based',
                         'native_diagonal_1_based','exact_positive_coefficient'])
        with native.open() as source:
            source.readline();number_variables=int(source.readline())
            number_blocks=int(source.readline())
            sizes=list(map(int,source.readline().split()))
            source.readline()
            if len(sizes)!=number_blocks or any(col>=number_variables-1 for col in rays):
                raise ValueError('Unexpected common-margin input')
            for line in source:
                fields=line.split()
                if len(fields)!=5:
                    raise ValueError('Incomplete native entry')
                k,block,i,j=map(int,fields[:4])
                if k-1 not in rays:
                    continue
                coefficient=rays[k-1]*F(fields[4])
                if i!=j or coefficient<=0:
                    raise ValueError('The proposed exact coefficient direction is not diagonal PSD')
                forced.setdefault(block,set()).add(i)
                counts[k-1]+=1
                writer.writerow([k,rays[k-1],block,i,str(coefficient)])
    if any(not c for c in counts.values()):
        raise ValueError('An alleged nonzero ray has no coefficient')
    proof={'native_input_sha256':solver_gate.file_digest(native),
           'ordinary_coordinate_rays_0_based':rays,
           'all_nonzero_coefficients_checked_exactly':True,
           'ray_nonzero_entry_counts':counts,
           'forced_dual_zero_diagonals_native_1_based':{b:sorted(v) for b,v in forced.items()},
           'native_block_sizes':sizes,
           'identity':'For each ray h, sum_b <F_b(h),Y_b>=0 by the adjoint equations. Every summand is nonnegative, so every supported dual diagonal is zero. PSD then forces its entire row and column to zero.',
           'coefficient_witness_sha256':solver_gate.file_digest(output/'ray_coefficients.csv.gz')}
    conic.write_json(output/'exact_ray_proof.json',proof)
    return proof


def check(output):
    P=conic.model(['x','y','z','u'],[
        {'name':'moment','entries':[[conic.form(1),conic.form(x=1)],
                                  [conic.form(x=1),conic.form(y=1)]]},
        {'name':'x_nonnegative','entries':[[conic.form(x=1)]]},
        {'name':'y_nonnegative','entries':[[conic.form(y=1)]]},
        {'name':'bounded_z_lower','entries':[[conic.form(z=1)]]},
        {'name':'bounded_z_upper','entries':[[conic.form(1,z=-1)]]},
        {'name':'u_nonpositive','entries':[[conic.form(u=-1)]]}])
    path=output/'fixture.dat-s'
    slater_margin.write_input(path,slater_margin.augmented(P))
    solver_gate.verify_native(path,None,fixture=True)
    data,layout=backend.read_native(path)
    rays=candidates(data,layout)
    if rays!={1:1,3:-1}:
        raise ValueError('Known recession directions or bounded/off-diagonal controls failed')
    proof=verify(path,rays,output)
    if proof['forced_dual_zero_diagonals_native_1_based'].get(1)!=[2]:
        raise ValueError('Wrong exposed PSD row')
    return {'status':'EXACT_DUAL_RECESSION_FIXTURE_VERIFIED','rays':rays,
            'bounded_coordinate_and_indefinite_coefficient_rejected':True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--native-input',type=Path)
    parser.add_argument('--preflight',type=Path)
    parser.add_argument('--encoding',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir()
    started=time.monotonic()
    if args.check:
        report=check(args.output)
    else:
        if any(p is None for p in [args.native_input,args.preflight,args.encoding]):
            raise ValueError('Require the original gated input and checked encoding')
        gate=solver_gate.verify_native(args.native_input,args.preflight)
        encoded=json.loads((args.encoding/'result.json').read_text())
        if encoded['source_native_sha256']!=gate['native_input_sha256']:
            raise ValueError('Different encoded input')
        for p,h in encoded['encoding_file_sha256'].items():
            if solver_gate.file_digest(args.encoding/p)!=h:
                raise ValueError('Changed encoding cache')
        data,layout=backend.load_encoding(args.encoding)
        rays=candidates(data,layout)
        proof=verify(args.native_input,rays,args.output)
        sizes=proof['native_block_sizes']
        forced=proof['forced_dual_zero_diagonals_native_1_based']
        report={'status':'EXACT_DIAGONAL_RECESSION_DIRECTIONS_VERIFIED' if rays else 'NO_DIAGONAL_COORDINATE_RECESSION_FOUND',
                'verified_preflight':gate,'ordinary_recession_directions':len(rays),
                'forced_PSD_dual_zero_rows':sum(len(v) for b,v in forced.items() if sizes[int(b)-1]>0),
                'forced_scalar_dual_zeros':sum(len(v) for b,v in forced.items() if sizes[int(b)-1]<0),
                'strict_certificate_cone_feasibility':'EXACTLY_EXCLUDED' if rays else 'UNDECIDED',
                'proper_dual_face_proved':bool(rays),
                'source_proof':'exact_ray_proof.json',
                'original_primal_feasibility_decided':False}
    report.update(source_sha256=solver_gate.file_digest(__file__),
                  accepted_bounds_changed=False,geometric_constraints_changed=False,
                  elapsed_seconds=time.monotonic()-started)
    conic.write_json(args.output/'result.json',report)
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    main()
