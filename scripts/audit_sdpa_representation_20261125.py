#!/usr/bin/env python3
"""Audit exact redundancies, exposed faces and numerical units of a saved SDP."""
import argparse
from collections import Counter
from fractions import Fraction as F
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
    parser.add_argument('--source', default='results/global_edge_closure_conditioned01')
    args = parser.parse_args()
    sys.path.insert(0, str(ROOT))
    import numpy as np
    from research.global_formulation import conic
    source = (ROOT/args.source).resolve()
    model_path = source/'equilibrated.exact.json'
    units_path = source/'equilibration.exact.json'
    out = (ROOT/args.output).resolve()
    if out.exists() or any(not p.is_relative_to(ROOT/'results') for p in (out, source)):
        raise ValueError('Use retained inputs and a fresh result')
    out.mkdir()
    started = time.monotonic()
    paths = [Path(__file__), model_path, units_path, source/'result.json',
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
        raise ValueError('Changed conditioned source')
    print('Loaded the unchanged equation-free conic representation.', flush=True)
    sizes = Counter()
    zero_scalar = []
    zero_psd_blocks = []
    zero_diagonals = []
    forced_rows = {}
    scalar_rows = []
    unit_rows = {}
    false_constants = []
    positive_constant_scalar = 0
    occurrences = Counter()
    norms = Counter()
    minimum_coefficient = float('inf')
    maximum_coefficient = 0.
    for b, block in enumerate(problem['blocks']):
        M = block['entries']
        r = len(M)
        sizes[r] += 1
        if r == 1:
            q = M[0][0]
            if not q: zero_scalar.append(b)
            else: scalar_rows.append((b,q))
            if q and not any(v for v in q):
                if q[''] < 0: false_constants.append(dict(block=b, name=block['name'], value=q['']))
                else: positive_constant_scalar += 1
            if len(q) == 1 and '' not in q:
                v,c = next(iter(q.items()))
                unit_rows.setdefault(v,{})[1 if c>0 else -1] = (b,c)
        else:
            if not any(q for row in M for q in row): zero_psd_blocks.append(b)
            for i,row in enumerate(M):
                if not row[i]:
                    zero_diagonals.append((b,i))
                    for q in row:
                        if q: forced_rows.setdefault(conic.digest_bytes(conic.canonical(q)), q)
        for i,row in enumerate(M):
            for j in range(i,r):
                for v,c in row[j].items():
                    magnitude = abs(float(c))
                    minimum_coefficient = min(minimum_coefficient, magnitude)
                    maximum_coefficient = max(maximum_coefficient, magnitude)
                    if v:
                        occurrences[v] += 1
                        norms[v] += (1 if i==j else 2)*magnitude*magnitude
    seen = {}
    duplicates = 0
    opposites = 0
    opposite_examples = []
    sign_faces = []
    sign_contradictions = []
    forced_unit_variables = set()
    for b,q in scalar_rows:
        first = q[min(q)]
        normalized = conic.times(1/abs(first),q)
        key = conic.digest_bytes(conic.canonical(normalized))
        negative = conic.digest_bytes(conic.canonical(conic.times(-1,normalized)))
        if key in seen: duplicates += 1
        if negative in seen:
            opposites += 1
            if len(opposite_examples)<6:
                opposite_examples.append(dict(blocks=[seen[negative],b], normalized_row=normalized))
        seen.setdefault(key,b)
        if q.get('',F(0))>0 or not any(v for v in q): continue
        combination = []
        for v,c in q.items():
            if not v: continue
            option = unit_rows.get(v,{}).get(-1 if c>0 else 1)
            if option is None: break
            d,a = option
            combination.append(dict(block=d, weight=-c/a, variable=v))
        else:
            total = conic.plus(q,*(conic.times(item['weight'],problem['blocks'][item['block']]['entries'][0][0])
                                   for item in combination))
            expected = conic.form(q.get('',F(0)))
            if total != expected or any(item['weight']<=0 for item in combination):
                raise ValueError('Invalid exact sign-bound identity')
            item = dict(block=b, name=problem['blocks'][b]['name'],
                        constant=q.get('',F(0)), coordinate_combination=combination)
            if q.get('',F(0))<0: sign_contradictions.append(item)
            else:
                sign_faces.append(item)
                forced_unit_variables.update(v for v in q if v)
    units = json.loads(units_path.read_text())
    columns = {v:F(c) for v,c in units['column_scales'].items()}
    lengths = [math.sqrt(c) for c in norms.values() if c]
    unused = [v for v in problem['variables'] if not occurrences[v]]
    def names(indices): return [dict(block=b,name=problem['blocks'][b]['name']) for b in indices[:8]]
    report = dict(status='EXACT_REPRESENTATION_AUDIT_COMPLETE_FEASIBILITY_UNDECIDED',
        model_sha256=receipt['prepared_matrix_sha256'], original_input_sha256=receipt['frozen_input_sha256'],
        variables=len(problem['variables']), blocks=len(problem['blocks']),
        scalar_blocks=sizes[1], nonscalar_PSD_blocks=sum(c for r,c in sizes.items() if r>1),
        maximum_PSD_size=max(sizes), block_size_histogram=dict(sizes),
        identically_zero_scalar_cones=len(zero_scalar), zero_scalar_examples=names(zero_scalar),
        identically_zero_PSD_blocks=len(zero_psd_blocks),
        identically_zero_PSD_diagonals=len(zero_diagonals),
        exact_forced_PSD_kernel_equations=len(forced_rows),
        forced_kernel_examples=list(forced_rows.values())[:8],
        scalar_duplicate_positive_rays=duplicates, scalar_opposite_ray_pairs=opposites,
        opposite_ray_examples=opposite_examples,
        scalar_positive_constant_tautologies=positive_constant_scalar,
        scalar_negative_constant_contradictions=false_constants,
        exact_nonnegative_sum_zero_faces=len(sign_faces),
        exact_forced_coordinate_zeros=len(forced_unit_variables), sign_face_examples=sign_faces[:6],
        exact_sign_bound_contradictions=sign_contradictions[:6],
        exact_sign_bound_contradiction_count=len(sign_contradictions),
        unused_variables=len(unused), unused_variable_examples=unused[:12],
        column_units_equal_one=sum(c==1 for c in columns.values()), column_units_total=len(columns),
        affine_coefficient_absolute_minimum=minimum_coefficient,
        affine_coefficient_absolute_maximum=maximum_coefficient,
        numerical_operator_column_norm_quantiles=dict(zip(('min','p01','p25','median','p75','p99','max'),
            map(float,np.percentile(lengths,[0,1,25,50,75,99,100])),strict=True)),
        floating_column_norms_do_not_measure_Schur_condition_number=True,
        common_margin_has_exact_nonpositive_upper_bound=bool(zero_scalar or zero_diagonals),
        constraints_modified=False, solver_launched=False, accepted_bounds_changed=False)
    if bindings()!=before: raise ValueError('An audit source changed')
    report['elapsed_seconds']=time.monotonic()-started
    conic.write_json(out/'result.json',report)
    print(json.dumps(conic.plain({k:v for k,v in report.items() if not k.endswith('examples')})),flush=True)


if __name__=='__main__': main()
