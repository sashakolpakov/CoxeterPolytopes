#!/usr/bin/env python3
"""Check fail-closed preflight, stale-file rejection and the native restart fixture."""
import argparse
import copy
import json
import hashlib
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get('COXETER_CHECK_OWNER'):raise RuntimeError('Use the scientific watchdog')
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True)
    args=parser.parse_args();out=(ROOT/args.output).resolve()
    if out.exists() or not out.is_relative_to(ROOT/'results'):raise ValueError('Use fresh results output')
    out.mkdir();sys.path.insert(0,str(ROOT))
    from research.global_formulation import conic,feasibility,sdpa_cli,sdpa_double,solver_gate
    paths=[Path(__file__),*[ROOT/'research/global_formulation'/(n+'.py') for n in ('conic','feasibility','sdpa_cli','sdpa_double','solver_gate','checkpoint_solver')]]
    def bindings():return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    before=bindings()
    P=conic.model(['x'],[dict(name='nonnegative',entries=[[conic.form(x=1)]])])
    native=out/'fixture.dat-s';conic.write_sdpa(native,P)
    def rejects(function):
        try:function()
        except (ValueError,FileNotFoundError):return
        raise ValueError('An invalid production launch was authorized')
    rejects(lambda:sdpa_cli.call(native,out,1,1))
    rejects(lambda:sdpa_double.native_call(native,out,1))
    rejects(lambda:feasibility.solve_sdpa(None,out,1,None))
    rejects(lambda:feasibility.run(ROOT,out,native,"sdpa","necessary"))
    if (out/'solver.stdout').exists():raise ValueError('The solver launched before checking its gate')
    audit=dict(schema='conic-fixed-point-audit-v1',passed=True,counts={k:0 for k in solver_gate.AUDIT_FIELDS})
    bad=copy.deepcopy(audit);bad['counts']['constant_cones']=1
    rejects(lambda:solver_gate.require_clean_audit(bad))
    artifact=out/'bound_source.json';conic.write_json(artifact,dict(exact_model=True))
    receipt=dict(schema='solver-ready-v1',status='SOLVER_READY',presolve_fixed_point=True,
        all_exact_reduction_proofs_verified=True,all_scaling_inverse_identities_verified=True,
        fixed_point_audit=audit,scientific_contract=dict(intrinsic_feature_coordinates=26),
        frozen_input_sha256='fixture-only',native_input_sha256=solver_gate.file_digest(native),
        native_free_variables=1,native_block_sizes=[-1],
        bound_artifacts={str(artifact.relative_to(ROOT)):solver_gate.file_digest(artifact)})
    ready=out/'fixture_ready.json';conic.write_json(ready,receipt)
    solver_gate.verify_native(native,ready)
    with native.open('a') as stream:stream.write('\n')
    rejects(lambda:solver_gate.verify_native(native,ready))
    native.write_text(native.read_text()[:-1])
    artifact.write_text('{}\n')
    rejects(lambda:solver_gate.verify_native(native,ready))
    rejects(lambda:solver_gate.check_registered_input('unregistered-scientific-model'))
    oversized=out/'not_a_fixture.dat-s';oversized.write_text('"Global Coxeter SDP\n17\n1\n-1\n'+'0 '*17+'\n')
    rejects(lambda:solver_gate.verify_native(oversized,None,fixture=True))
    sdpa_cli.self_check(ROOT,out/'native_fixture')
    report=dict(status='GLOBAL_SOLVER_PREFLIGHT_REGRESSIONS_PASSED',missing_preflight_rejected_before_launch=True,all_legacy_native_backends_fail_closed=True,
        incomplete_reduction_rejected=True,changed_native_matrix_rejected=True,
        changed_bound_source_rejected=True,unregistered_constraints_rejected=True,
        large_fixture_bypass_rejected=True,native_full_precision_restart_fixture_passed=True,
        scientific_model_claimed_feasible=False,source_bindings=before)
    if bindings()!=before:raise ValueError('A gate regression source changed')
    conic.write_json(out/'result.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':main()
