#!/usr/bin/env python3
"""Count native SDPA-GMP pair caches in the full verified explicit dual.

These are structural counts, not a solve, an infeasibility test or a prediction
of total RSS. In particular, symbolic Cholesky fill is not computed here.
"""
import argparse
import ctypes
import json
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from research.global_formulation import conic,conic_backend as backend,solver_gate
from scripts.search_balanced_dual_field_20261126 import checked_model


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model',type=Path,required=True)
    parser.add_argument('--replay',type=Path,required=True)
    parser.add_argument('--native-source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if not os.environ.get('COXETER_CHECK_OWNER'):raise RuntimeError('Use the serial scientific watchdog')
    if args.output.exists() or not args.output.resolve().is_relative_to(ROOT/'results'):
        raise ValueError('Use a fresh results folder')
    args.output.mkdir();checked_model(args.model,args.replay)
    import numpy as np
    data,layout=backend.load_encoding(args.model);A=data['A'].tocsr()
    scalar_counts=np.diff(A.indptr)[:layout['scalar_cones']]
    scalar_pairs=sum(int(k)*(int(k)+1)//2 for k in scalar_counts)
    PSD_counts=[];base=layout['scalar_cones']
    for n in layout['PSD_sizes']:
        end=base+n*(n+1)//2
        PSD_counts.append(len(np.unique(A.indices[A.indptr[base]:A.indptr[end]])))
        base=end
    if base!=A.shape[0] or layout['zero_cones']:raise ValueError('Unexpected exact-dual layout')
    PSD_pairs=sum(k*(k+1)//2 for k in PSD_counts)
    native=args.native_source.read_text()
    expected=['SDP_constraint1[l] = new int[tmp]','SDP_constraint2[l] = new int[tmp]',
        'SDP_blockIndex1[l] = new int[tmp]','SDP_blockIndex2[l] = new int[tmp]',
        'SDP_location_sparse_bMat[l] = new int[tmp]',
        'LP_constraint1[l] = new int[tmp]','LP_constraint2[l] = new int[tmp]',
        'LP_blockIndex1[l] = new int[tmp]','LP_blockIndex2[l] = new int[tmp]',
        'LP_location_sparse_bMat[l] = new int[tmp]']
    if any(q not in native for q in expected):raise ValueError('Changed native pair-cache allocation code')
    width=ctypes.sizeof(ctypes.c_int)
    report=dict(status='SDPA_GMP_FULL_MODEL_SYMBOLIC_PAIR_CACHE_COUNTS',
        variables=A.shape[1],scalar_cones=layout['scalar_cones'],PSD_blocks=len(PSD_counts),
        largest_PSD_variable_support=max(PSD_counts),largest_scalar_variable_support=int(max(scalar_counts)),
        PSD_pair_slots=PSD_pairs,scalar_pair_slots=scalar_pairs,total_pair_slots=PSD_pairs+scalar_pairs,
        native_int_bytes=width,integer_arrays_per_pair=5,
        pair_cache_integer_payload_bytes=5*width*(PSD_pairs+scalar_pairs),
        persistent_formula_enum_slots=A.shape[1]*len(PSD_counts),
        temporary_formula_integer_payload_bytes=width*A.shape[1]*len(PSD_counts),
        cache_four_derived_indices_payload_bytes=4*width*(PSD_pairs+scalar_pairs),
        possible_cache_redesign_implemented=False,Cholesky_fill_computed=False,
        total_RSS_or_actual_allocated_heap_not_inferred=True,numerical_solve_performed=False,
        original_feasibility_status='UNDECIDED',
        source_bindings={str(p):dict(sha256=solver_gate.file_digest(p)) for p in
            [Path(__file__),args.model/'result.json',args.replay,args.native_source]})
    conic.write_json(args.output/'result.json',report)
    print(json.dumps(report),flush=True)


if __name__=='__main__':main()
