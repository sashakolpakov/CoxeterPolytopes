#!/usr/bin/env python3
"""Check saved dual rays for exact cone-face identities, without a new SDP."""
import argparse
import hashlib
import json
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
    parser.add_argument('--checkpoint', default='results/global_edge_closure_sdpa04')
    args = parser.parse_args()
    sys.path.insert(0, str(ROOT))
    from research.global_formulation import conic, feasibility, sdpa_cli, slater_margin
    source = (ROOT/args.source).resolve()
    native = (ROOT/args.checkpoint).resolve()/'native.result'
    out = (ROOT/args.output).resolve()
    model_path = source/'equilibrated.exact.json'
    if out.exists() or any(not p.is_relative_to(ROOT/'results') for p in (out, model_path, native)):
        raise ValueError('Use retained inputs and fresh output')
    out.mkdir()
    started = time.monotonic()
    paths = [Path(__file__), native, model_path, source/'result.json', sdpa_cli.EXECUTABLE,
             *[ROOT/'research/global_formulation'/(name+'.py') for name in
               ('conic', 'feasibility', 'sdpa_cli', 'slater_margin', 'face_phase')]]
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
    sdpa_cli.self_check(ROOT, out/'output_format_and_restart_fixture')
    print('Native SDPA verified all four requested output formats and restart consistency.', flush=True)
    problem = conic.read_model(model_path)
    receipt = json.loads((source/'result.json').read_text())
    if conic.digest_bytes(conic.canonical(problem)) != receipt['prepared_matrix_sha256']:
        raise ValueError('Changed source matrix')
    # The saved auxiliary cone ordering equals the source ordering plus
    # one final scalar upper bound.  Do not duplicate the large source.
    auxiliary_shape = dict(problem, blocks=[*problem['blocks'], dict(
        name='margin_upper_bound', entries=[[{}]])])
    matrices = feasibility.native_cones(native, auxiliary_shape, field='yMat')[:-1]
    print('Searching rational PSD directions suggested by the saved dual matrices.', flush=True)
    candidate, checked = slater_margin.sparse_exposure(problem, matrices)
    report = dict(status='NO_EXACT_FACE_RECOVERED_FROM_SELECTED_SAVED_DUAL_RAYS',
        source_original_sha256=receipt['frozen_input_sha256'], accepted_bounds_changed=False,
        original_feasibility='UNDECIDED', geometric_constraints_added=False,
        selected_rays_are_a_restricted_recovery_search=True,
        output_format_and_restart_fixture='PASSED', exposure_search=checked)
    if candidate is not None:
        equations, proof = candidate
        conic.write_json(out/'verified_face.exact.json', dict(
            equations=equations, positive_zero_identity=proof,
            source_sha256=before[str(model_path)], identity_checked_exactly=True))
        report.update(status='EXACT_SOURCE_CONE_FACE_EXPOSED_FEASIBILITY_UNDECIDED',
                      forced_equations=len(equations), positive_summands=len(proof))
    if bindings() != before:
        raise ValueError('A dual-face input changed')
    report['elapsed_seconds'] = time.monotonic()-started
    conic.write_json(out/'result.json', report)
    print(json.dumps(conic.plain(report)), flush=True)


if __name__ == '__main__':
    main()
