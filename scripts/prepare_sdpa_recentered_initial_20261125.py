#!/usr/bin/env python3
"""Build a different numerical interior seed for the unchanged margin SDP.

Reuse saved ordinary coordinates, move the common margin safely below all
computed eigenvalues, and replace every dual block by a normalized inverse
primal slack.  The numerical checks below concern initialization only.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get('COXETER_CHECK_OWNER'):
        raise RuntimeError('Use the serial scientific watchdog')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--checkpoint', required=True)
    parser.add_argument('--source', default='results/global_edge_closure_conditioned01')
    args = parser.parse_args()
    sys.path.insert(0, str(ROOT))
    import numpy as np
    from research.global_formulation import conic

    source = (ROOT/args.source).resolve()
    native = (ROOT/args.checkpoint).resolve()/'native.result'
    model_path = source/'equilibrated.exact.json'
    out = (ROOT/args.output).resolve()
    if out.exists() or any(not p.is_relative_to(ROOT/'results') for p in (out, model_path, native)):
        raise ValueError('Use retained inputs and a fresh output')
    out.mkdir()
    started = time.monotonic()
    paths = [Path(__file__), native, model_path, source/'result.json',
             ROOT/'research/global_formulation/conic.py']
    def bindings():
        result = {}
        for path in paths:
            h = hashlib.sha256()
            with path.open('rb') as stream:
                for chunk in iter(lambda: stream.read(1024**2), b''): h.update(chunk)
            result[str(path)] = h.hexdigest()
        return result
    before = bindings()
    conic.write_json(out/'source_bindings.json', before)
    problem = conic.read_model(model_path)
    receipt = json.loads((source/'result.json').read_text())
    if conic.digest_bytes(conic.canonical(problem)) != receipt['prepared_matrix_sha256']:
        raise ValueError('Changed conditioned SDP')
    if problem['equations']:
        raise ValueError('Use the equation-free conditioned source')
    raw, phase = conic.native_vector(native)
    if len(raw) != len(problem['variables'])+1:
        raise ValueError('Wrong common-margin vector dimension')
    point = dict(zip(problem['variables'], map(float, raw[:-1]), strict=True))
    blocks = []
    smallest = float('inf')
    for block in problem['blocks']:
        A = np.asarray([[math.fsum(float(c)*(point[v] if v else 1.) for v,c in q.items())
                         for q in row] for row in block['entries']])
        if not np.all(np.isfinite(A)):
            raise ValueError('Nonfinite affine slack')
        low = float(np.linalg.eigvalsh(A)[0]) if len(A)>1 else float(A[0,0])
        smallest = min(smallest, low)
        blocks.append(A)
    rho = max(1, math.ceil(1-smallest))
    full = []
    inverses = []
    scalar = []
    inverse_traces = []
    minimum = float('inf')
    maximum_relative_centrality_error = 0.
    maximum_affine_residual = 0.
    for A in blocks:
        X = A+rho*np.eye(len(A))
        minimum = min(minimum, float(np.linalg.eigvalsh(X)[0]) if len(X)>1 else float(X[0,0]))
        maximum_affine_residual = max(maximum_affine_residual,
            float(np.max(np.abs(X-A-rho*np.eye(len(A))))))
        if len(X)==1:
            scalar.append(float(X[0,0]))
            inverse_traces.append(1/float(X[0,0]))
        else:
            np.linalg.cholesky(X)
            inverse = np.linalg.inv(X)
            inverse = (inverse+inverse.T)/2
            maximum_relative_centrality_error = max(maximum_relative_centrality_error,
                float(np.max(np.abs(X@inverse-np.eye(len(X))))))
            full.append(X)
            inverses.append(inverse)
            inverse_traces.append(float(np.trace(inverse)))
    scalar.append(float(1+rho))
    inverse_traces.append(1/float(1+rho))
    tau = 1/math.fsum(inverse_traces)
    trace_residual = abs(math.fsum(tau*t for t in inverse_traces)-1)
    if minimum <= .5 or trace_residual > 1e-12 or maximum_relative_centrality_error > 1e-6:
        raise ValueError('The proposed numerical interior seed failed validation')
    def number(c): return format(float(c), '+.17e')
    def write_matrices(stream, dual):
        stream.write('{\n')
        for X, inverse in zip(full, inverses, strict=True):
            M = tau*inverse if dual else X
            stream.write('{'+',\n'.join('{'+', '.join(number(c) for c in row)+'}' for row in M)+'},\n')
        stream.write('{'+', '.join(number(tau/c if dual else c) for c in scalar)+'}\n}\n')
    initial = out/'initial.dat'
    with initial.open('x') as stream:
        stream.write('xVec = {'+', '.join(number(v) for v in raw[:-1])+', '+number(-rho)+'}\n')
        stream.write('xMat = ')
        write_matrices(stream, False)
        stream.write('yMat = ')
        write_matrices(stream, True)
    if bindings() != before:
        raise ValueError('An initialization source changed')
    report = dict(status='NUMERICALLY_VALIDATED_RECENTERED_AUXILIARY_SEED',
        original_feasibility_claimed=False, accepted_bounds_changed=False,
        geometric_constraints_changed=False, ordinary_coordinates_reused_without_clipping=True,
        reference_native_phase=phase, old_common_margin=raw[-1], new_common_margin=-rho,
        minimum_unshifted_source_eigenvalue=smallest, minimum_shifted_eigenvalue=minimum,
        maximum_primal_affine_roundoff_residual=maximum_affine_residual,
        normalized_dual_trace_residual=trace_residual, common_complementarity=tau,
        maximum_relative_centrality_error=maximum_relative_centrality_error,
        all_nonscalar_primal_slacks_passed_numerical_Cholesky=True,
        dual_seed='normalized inverse of each new auxiliary primal slack',
        original_input_sha256=receipt['frozen_input_sha256'],
        initial_path=str(initial.relative_to(ROOT)), initial_bytes=initial.stat().st_size,
        initial_sha256=hashlib.sha256(initial.read_bytes()).hexdigest(),
        elapsed_seconds=time.monotonic()-started)
    conic.write_json(out/'result.json', report)
    print(json.dumps(conic.plain(report)), flush=True)


if __name__ == '__main__':
    main()
