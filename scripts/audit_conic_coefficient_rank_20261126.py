#!/usr/bin/env python3
"""Certify a diagonal scalar-row minor of the full native coefficient map.

A nonzero pure-coordinate row for each column proves column independence.
This says nothing about Slater feasibility or numerical conditioning of a
weighted Newton system. Scalar coefficients are reread as exact decimal
rationals from every relevant native input row.
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
from research.global_formulation import conic, conic_backend as backend, solver_gate


def main():
    import numpy as np
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native-input',type=Path,required=True)
    parser.add_argument('--preflight',type=Path,required=True)
    parser.add_argument('--encoding',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    started=time.monotonic()
    gate=solver_gate.verify_native(args.native_input,args.preflight)
    args.output.mkdir()
    encoded=json.loads((args.encoding/'result.json').read_text())
    if encoded['source_native_sha256']!=gate['native_input_sha256']:
        raise ValueError('Different native input')
    for p,h in encoded['encoding_file_sha256'].items():
        if solver_gate.file_digest(args.encoding/p)!=h:
            raise ValueError('Changed encoding cache')
    data,layout=backend.load_encoding(args.encoding)
    A=data['A']
    scalar=A[:layout['scalar_cones']].tocsr()
    rows=np.flatnonzero(np.diff(scalar.indptr)==1)
    positions=scalar.indptr[rows]
    columns=scalar.indices[positions]
    coefficients=scalar.data[positions]
    selected=np.full(A.shape[1],-1,dtype=np.int64)
    strength=np.zeros(A.shape[1])
    for row,col,value in zip(rows,columns,coefficients,strict=True):
        if abs(value)>strength[col]:
            strength[col]=abs(value)
            selected[col]=row
    wanted={int(row):int(col) for col,row in enumerate(selected) if row>=0}
    if len(wanted)!=int(np.count_nonzero(selected>=0)):
        raise ValueError('Selected scalar rows are not distinct')
    bases={}
    total=0
    for block,size in enumerate(layout['native_block_sizes'],1):
        if size<0:
            bases[block]=total
            total-=size
    witnessed={}
    with args.native_input.open() as stream:
        for _ in range(5):
            stream.readline()
        for line in stream:
            fields=line.split()
            if len(fields)!=5:
                raise ValueError('Malformed native entry')
            k,block,i,j=map(int,fields[:4])
            if not k or block not in bases:
                continue
            row=bases[block]+i-1
            if row not in wanted:
                continue
            col=wanted[row]
            if i!=j or k!=col+1 or col in witnessed:
                raise ValueError('The proposed native row is not a unique pure-coordinate row')
            coefficient=-F(fields[4])
            if not coefficient:
                raise ValueError('Zero claimed pivot')
            witnessed[col]=(row,block,i,coefficient)
    if set(witnessed)!=set(wanted.values()):
        raise ValueError('Some proposed pivot was not found in the native input')
    with gzip.open(args.output/'scalar_minor.csv.gz','wt',newline='') as stream:
        writer=csv.writer(stream)
        writer.writerow(['column_0_based','canonical_scalar_row_0_based','native_block_1_based',
                         'native_row_1_based','exact_canonical_coefficient'])
        for col,(row,block,i,c) in sorted(witnessed.items()):
            writer.writerow([col,row,block,i,str(c)])
    full=len(witnessed)==A.shape[1]
    report={'status':'EXACT_NATIVE_SCALAR_MINOR_FULL_COLUMN_RANK_VERIFIED' if full else 'PARTIAL_EXACT_SCALAR_MINOR_VERIFIED',
        'verified_preflight':gate,'columns':A.shape[1],'independent_columns_witnessed':len(witnessed),
        'full_column_rank':full,'pure_scalar_rows_available':len(rows),
        'minimum_absolute_exact_pivot':str(min(abs(v[3]) for v in witnessed.values())),
        'maximum_absolute_exact_pivot':str(max(abs(v[3]) for v in witnessed.values())),
        'native_decimal_coefficients_checked_without_rounding':True,
        'argument':'The selected distinct scalar rows form a diagonal square minor with every diagonal entry nonzero.',
        'scope':'Rank of the coefficient map, not Slater feasibility, spectral conditioning, or geometric realizability.',
        'minor_witness':'scalar_minor.csv.gz',
        'minor_witness_sha256':solver_gate.file_digest(args.output/'scalar_minor.csv.gz'),
        'source_sha256':solver_gate.file_digest(__file__),
        'geometric_constraints_changed':False,'accepted_bounds_changed':False,
        'elapsed_seconds':time.monotonic()-started}
    conic.write_json(args.output/'result.json',report)
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    main()
