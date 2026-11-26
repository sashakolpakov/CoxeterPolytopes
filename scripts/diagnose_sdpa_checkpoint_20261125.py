#!/usr/bin/env python3
"""Recompute primal/dual diagnostics from a saved common-margin SDPA state.

Floating residuals are diagnostics, never feasibility certificates.  The
separate original scalar-row check uses the exact saved decimal fractions.
Run through the serial scientific watchdog.
"""
import argparse
from fractions import Fraction as F
import gc
import hashlib
import heapq
import json
import os
from pathlib import Path
import re
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024**2), b''):
            h.update(chunk)
    return dict(sha256=h.hexdigest(), bytes=path.stat().st_size)


def largest(heap, score, index, item, limit=16):
    entry = (float(score), index, item)
    if len(heap) < limit:
        heapq.heappush(heap, entry)
    elif entry[:2] > heap[0][:2]:
        heapq.heapreplace(heap, entry)


def ordered(heap):
    return [item for _, _, item in sorted(heap, reverse=True)]


def main():
    if not os.environ.get('COXETER_CHECK_OWNER'):
        raise RuntimeError('Use the serial scientific watchdog')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--checkpoint', default='results/global_edge_closure_sdpa04')
    parser.add_argument('--source', default='results/global_edge_closure_conditioned01')
    parser.add_argument('--original', default='results/global_edge_closure_current/necessary_primal.exact.json')
    parser.add_argument('--original-candidate', default='results/global_edge_closure_terminal_check01/original_numerical_candidate.json')
    args = parser.parse_args()
    sys.path.insert(0, str(ROOT))
    import numpy as np
    from research.global_formulation import conic, feasibility

    out = (ROOT/args.output).resolve()
    source = (ROOT/args.source).resolve()
    native = (ROOT/args.checkpoint).resolve()/'native.result'
    model_path = source/'equilibrated.exact.json'
    original_path = (ROOT/args.original).resolve()
    candidate_path = (ROOT/args.original_candidate).resolve()
    paths = [Path(__file__), native, model_path, source/'result.json',
             original_path, candidate_path,
             ROOT/'research/global_formulation/conic.py',
             ROOT/'research/global_formulation/feasibility.py']
    solver_source = Path('/tmp/coxeter_sdpa_native_20261125/sdpa-7.3.21')
    paths.extend(solver_source/name for name in ('sdpa_io.cpp', 'sdpa_dataset.cpp', 'sdpa_solve.cpp'))
    if out.exists() or any(not p.is_relative_to(ROOT/'results') for p in
                           (out, native, model_path, original_path, candidate_path)):
        raise ValueError('Use retained inputs and a fresh results output')
    out.mkdir()
    started = time.monotonic()
    bindings = {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p): digest(p) for p in paths}
    conic.write_json(out/'source_bindings.json', bindings)
    print('Loading the frozen conditioned model and all three saved arrays.', flush=True)
    problem = conic.read_model(model_path)
    receipt = json.loads((source/'result.json').read_text())
    if conic.digest_bytes(conic.canonical(problem)) != receipt['prepared_matrix_sha256']:
        raise ValueError('Changed conditioned matrix')
    if bindings[str(original_path.relative_to(ROOT))]['sha256'] != receipt['frozen_input_sha256']:
        raise ValueError('Wrong original matrix')
    if problem['equations']:
        raise ValueError('The conditioned source should be equation free')
    raw, phase = conic.native_vector(native)
    variables = problem['variables']
    if len(raw) != len(variables)+1:
        raise ValueError('Wrong common-margin vector dimension')
    x = np.array([float(c) for c in raw[:-1]])
    alpha = float(raw[-1])
    index = {v: i for i, v in enumerate(variables)}
    arrays = [np.fromiter((float(s) for s in conic.NUMBER.findall(feasibility.native_array(native, field))),
                         dtype=np.float64) for field in ('xMat', 'yMat')]
    offsets = np.empty(len(problem['blocks']), dtype=np.int64)
    position = 0
    for scalar in (False, True):
        for b, block in enumerate(problem['blocks']):
            r = len(block['entries'])
            if (r == 1) == scalar:
                offsets[b] = position
                position += r*r
    if any(len(a) != position+1 for a in arrays):
        raise ValueError('Wrong cone-array size, including the final margin-upper scalar')
    lam = float(arrays[1][-1])
    upper_slack = float(arrays[0][-1])
    adjoint = np.zeros(len(x))
    adjoint_scale = np.zeros(len(x))
    constant_terms = []
    trace_terms = []
    complementarity_terms = [upper_slack*lam]
    affine_pair_terms = []
    affine_norm2 = 0.
    affine_max = 0.
    source_pair_terms = []
    maxima = dict(affine=[], trace=[], dual_condition=[], primal_condition=[],
                  negative_source=[], complementarity=[])
    cone_stats = {key: dict(minimum_eigenvalue=float('inf'), maximum_eigenvalue=float('-inf'),
                           negative_eigenvalues=0, relatively_small_eigenvalues=0,
                           worst_block=None, largest_absolute_entry=0.)
                  for key in ('saved_primal_slack', 'saved_dual_multiplier', 'unshifted_source_slack')}
    max_asymmetry = dict(xMat=0., yMat=0.)
    for b, block in enumerate(problem['blocks']):
        r = len(block['entries'])
        start = int(offsets[b])
        X, Y = (a[start:start+r*r].reshape((r, r)) for a in arrays)
        max_asymmetry['xMat'] = max(max_asymmetry['xMat'], float(np.max(np.abs(X-X.T))))
        max_asymmetry['yMat'] = max(max_asymmetry['yMat'], float(np.max(np.abs(Y-Y.T))))
        X = (X+X.T)/2
        Y = (Y+Y.T)/2
        A = np.empty((r, r))
        block_constant = 0.
        for i, row in enumerate(block['entries']):
            for j in range(i, r):
                terms = []
                weight = float(Y[i, j])*(1 if i == j else 2)
                for v, exact in row[j].items():
                    c = float(exact)
                    if v:
                        k = index[v]
                        terms.append(c*x[k])
                        adjoint[k] += c*weight
                        adjoint_scale[k] += abs(c*weight)
                    else:
                        terms.append(c)
                        block_constant += c*weight
                A[i, j] = A[j, i] = sum(terms)
        Q = A-alpha*np.eye(r)
        R = X-Q
        residual = float(np.max(np.abs(R)))
        affine_max = max(affine_max, residual)
        affine_norm2 += float(np.sum(R*R))
        affine_pair_terms.append(float(np.sum(R*Y)))
        tr = float(np.trace(Y))
        comp = float(np.sum(X*Y))
        trace_terms.append(tr)
        constant_terms.append(block_constant)
        source_pair_terms.append(float(np.sum(A*Y)))
        complementarity_terms.append(comp)
        item = dict(block=b, name=block['name'], size=r)
        largest(maxima['affine'], residual, b, dict(item, residual=residual))
        largest(maxima['trace'], tr, b, dict(item, trace=tr, constant_pairing=block_constant))
        largest(maxima['complementarity'], abs(comp), b, dict(item, complementarity=comp))
        for key, matrix in (('saved_primal_slack', X), ('saved_dual_multiplier', Y),
                            ('unshifted_source_slack', A)):
            eig = np.linalg.eigvalsh(matrix) if r > 1 else matrix[0]
            low, high = float(eig[0]), float(eig[-1])
            scale = max(abs(low), abs(high), 1e-300)
            stat = cone_stats[key]
            if low < stat['minimum_eigenvalue']:
                stat.update(minimum_eigenvalue=low, worst_block=item)
            stat['maximum_eigenvalue'] = max(stat['maximum_eigenvalue'], high)
            stat['negative_eigenvalues'] += int(np.count_nonzero(eig < 0))
            stat['relatively_small_eigenvalues'] += int(np.count_nonzero(np.abs(eig) <= 1e-10*scale))
            stat['largest_absolute_entry'] = max(stat['largest_absolute_entry'], float(np.max(np.abs(matrix))))
            if key == 'unshifted_source_slack':
                largest(maxima['negative_source'], -low, b, dict(item, minimum_eigenvalue=low,
                        minimum_eigenvalue_divided_by_spectral_scale=low/scale))
            elif r > 1:
                group = 'dual_condition' if key == 'saved_dual_multiplier' else 'primal_condition'
                largest(maxima[group], scale/max(low, 1e-300), b,
                        dict(item, minimum_eigenvalue=low, maximum_eigenvalue=high,
                             relative_minimum_eigenvalue=low/scale))
        if b and b % 50000 == 0:
            print(f'Independently evaluated {b} source cones.', flush=True)
    import math
    c0 = math.fsum(constant_terms)
    trace = math.fsum(trace_terms)
    trace_residual = trace+lam-1
    upper_residual = upper_slack-(1-alpha)
    affine_max = max(affine_max, abs(upper_residual))
    affine_norm2 += upper_residual**2
    complementarity = math.fsum(complementarity_terms)
    affine_pair = math.fsum(affine_pair_terms)+upper_residual*lam
    dot_residual = float(np.dot(x, adjoint))
    gap = -alpha+c0+lam
    worst = np.argsort(np.abs(adjoint))[-16:][::-1]
    source_id = str(model_path.relative_to(ROOT))
    report = dict(status='SAVED_ITERATE_DIAGNOSED_ORIGINAL_FEASIBILITY_UNDECIDED',
        all_residuals_recomputed_from_solution_arrays=True,
        floating_diagnostics_are_not_feasibility_or_infeasibility_certificates=True,
        original_input_sha256=receipt['frozen_input_sha256'], source=source_id, native_phase=phase,
        common_margin_exact_saved_decimal=raw[-1], variables=len(x), cones=len(problem['blocks']),
        formulation=dict(primal='minimize -alpha; P_b(x)-alpha I PSD; 1-alpha >= 0',
                         dual='maximize -sum <P_b0,Y_b>-lambda; A*(Y)=0; sum trace(Y)+lambda=1; Y PSD, lambda>=0',
                         xVec='(x,alpha)', xMat='primal auxiliary affine cone slacks',
                         yMat='dual cone multipliers, including final upper-bound multiplier'),
        recomputed_primal=dict(objective=-alpha, maximum_absolute_affine_residual=affine_max,
            affine_residual_frobenius_norm=math.sqrt(affine_norm2),
            upper_bound_slack=upper_slack, upper_bound_affine_residual=upper_residual),
        recomputed_dual=dict(objective=-(c0+lam), original_Farkas_constant=c0,
            margin_upper_multiplier=lam, source_multiplier_trace=trace,
            trace_normalization_residual=trace_residual,
            ordinary_adjoint_maximum_absolute_residual=float(np.max(np.abs(adjoint))),
            ordinary_adjoint_euclidean_residual=float(np.linalg.norm(adjoint)),
            ordinary_adjoint_maximum_residual_divided_by_one_plus_absolute_terms=
                float(np.max(np.abs(adjoint)/(1+adjoint_scale))),
            worst_coordinates=[dict(variable=variables[k], residual=float(adjoint[k]),
                                    absolute_terms=float(adjoint_scale[k]), value=float(x[k])) for k in worst]),
        recomputed_gap=gap, recomputed_complementarity=complementarity,
        complementarity_identity=dict(x_dot_ordinary_dual_residual=dot_residual,
            affine_residual_pairing=affine_pair, minus_alpha_times_trace_residual=-alpha*trace_residual,
            check=complementarity-(gap+dot_residual-alpha*trace_residual+affine_pair)),
        cone_statistics=cone_stats, cone_array_maximum_asymmetry=max_asymmetry,
        largest_free_coordinates=[dict(variable=variables[k], value=float(x[k]))
                                  for k in np.argsort(np.abs(x))[-16:][::-1]],
        block_details={key: ordered(value) for key, value in maxima.items()},
        accepted_bounds_changed=False)
    header = native.read_text().split('xVec =', 1)[0]
    report['native_cached_summary'] = {key: value for key, value in re.findall(
        r'^\s*(objValPrimal|objValDual|p\.feas\.error|d\.feas\.error|relative gap|mu)\s*=\s*(\S+)', header, re.M)}
    conic.write_json(out/'conditioned_iterate_diagnostics.json', report)
    print(json.dumps(dict(primal=report['recomputed_primal'], dual=report['recomputed_dual'],
                          cone_statistics=cone_stats, gap=gap, complementarity=complementarity)), flush=True)
    del problem, raw, arrays, X, Y, A, Q, R
    gc.collect()
    original = conic.read_model(original_path)
    point = {v: F(c) for v, c in json.loads(candidate_path.read_text()).items()}
    row_name = 'uniform_two_unit_triangle_cover0'
    matches = [b for b in original['blocks'] if b['name'] == row_name]
    if len(matches) != 1 or len(matches[0]['entries']) != 1:
        raise ValueError('Missing original scalar witness row')
    row = matches[0]['entries'][0][0]
    scalar = conic.value(row, point)
    if scalar >= 0:
        # A different initialization can violate a different source row.
        # The previous checkpoint's first failure is only a fast path.
        for block in original['blocks']:
            if len(block['entries']) != 1:
                continue
            q = block['entries'][0][0]
            proposed = conic.value(q, point)
            if proposed < 0:
                row_name, row, scalar = block['name'], q, proposed
                break
    eq_count = 0
    nonzero = []
    for k, q in enumerate(original['equations']):
        v = conic.value(q, point)
        eq_count += 1
        if v and len(nonzero) < 10:
            nonzero.append(dict(equation=k, value=v))
    report['exact_saved_original_equalities'] = dict(checked=eq_count, all_zero=not nonzero,
                                                    nonzero_examples=nonzero)
    if scalar < 0:
        report['exact_saved_original_candidate_rejection'] = dict(block=row_name, expression=row,
            exact_value=scalar, approximate_value=float(scalar), required_sign='nonnegative',
            all_original_equalities_checked_exactly=eq_count, all_original_equalities_zero=not nonzero,
            nonzero_equality_examples=nonzero,
            interpretation='The saved decimal candidate is infeasible; this does not decide whether another point exists.')
    else:
        report['exact_scalar_witness_search'] = 'No negative scalar row; original PSD verification is still required.'
    if {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p): digest(p) for p in paths} != bindings:
        raise ValueError('A diagnostic input changed')
    report['elapsed_seconds'] = time.monotonic()-started
    conic.write_json(out/'result.json', report)
    print(json.dumps(conic.plain(dict(status=report['status'],
        exact_candidate_rejection=report.get('exact_saved_original_candidate_rejection'),
        elapsed_seconds=report['elapsed_seconds']))), flush=True)


if __name__ == '__main__':
    main()
