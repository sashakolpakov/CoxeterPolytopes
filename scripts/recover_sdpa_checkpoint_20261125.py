#!/usr/bin/env python3
"""Check a common-margin SDPA checkpoint against the frozen original SDP."""
import argparse
from fractions import Fraction as F
import json
import os
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get("COXETER_CHECK_OWNER"): raise RuntimeError("Use the serial scientific watchdog")
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",required=True)
    parser.add_argument("--checkpoint",required=True)
    parser.add_argument("--source",required=True)
    parser.add_argument("--input",required=True)
    args=parser.parse_args();out=(ROOT/args.output).resolve();source=(ROOT/args.source).resolve()
    original_path=(ROOT/args.input).resolve();native=(ROOT/args.checkpoint).resolve()/"native.result"
    if out.exists() or any(not p.is_relative_to(ROOT/"results") for p in (out,source,original_path,native)):
        raise ValueError("Use retained results inputs and fresh output")
    sys.path.insert(0,str(ROOT))
    from research.global_formulation import conic,dual_completion,feasibility,slater_margin,solver_gate
    gate=solver_gate.verify_native(source/'input.dat-s',source/'solver_ready.json')
    checkpoint_metadata=json.loads((native.parent/'native_metadata.json').read_text())
    if checkpoint_metadata['input_sha256']!=gate['native_input_sha256']:
        raise ValueError('The checkpoint belongs to a different native input')
    started=time.monotonic();out.mkdir()
    model_path=source/"equilibrated.exact.json";subs_path=source/"equilibrated_substitutions.exact.json"
    paths=[original_path,model_path,subs_path,native,source/"result.json",Path(__file__),
        *[ROOT/"research/global_formulation"/(name+".py") for name in
            ("conic","dual_completion","feasibility","slater_margin")]]
    def binding(): return {str(p.relative_to(ROOT)):conic.digest_bytes(p.read_bytes()) for p in paths}
    before=binding();conic.write_json(out/"source_bindings.json",before)
    receipt=json.loads((source/"result.json").read_text())
    if receipt["frozen_input_sha256"]!=before[str(original_path.relative_to(ROOT))]: raise ValueError("Wrong frozen original")
    problem=conic.read_model(model_path)
    if conic.digest_bytes(conic.canonical(problem))!=receipt["prepared_matrix_sha256"]: raise ValueError("Changed prepared matrix")
    substitutions={v:{w:F(c) for w,c in q.items()} for v,q in json.loads(subs_path.read_text()).items()}
    if conic.digest_bytes(conic.canonical(substitutions))!=receipt["prepared_substitutions_sha256"]: raise ValueError("Changed reconstruction map")
    raw,phase=conic.native_vector(native)
    if len(raw)!=len(problem["variables"])+1: raise ValueError("Wrong common-margin checkpoint dimension")
    values=dict(zip(problem["variables"],raw[:-1],strict=True))
    point={v:conic.value(q,values) for v,q in substitutions.items()}
    conic.write_json(out/"original_numerical_candidate.json",point)
    original=conic.read_model(original_path)
    report=dict(status="NO_EXACT_FEASIBILITY_DECISION",native_phase=phase,native_common_margin=raw[-1],
        accepted_bounds_changed=False,source_bound_geometric_review_complete=False,
        original_numerical_residuals=feasibility.numerical_residuals(original,point))
    try:
        exact,checked=conic.reconstruct(original,problem,substitutions,raw[:-1])
        conic.write_json(out/"original_feasible_point.exact.json",exact)
        conic.write_json(out/"exact_original_point_check.json",checked)
        full=checked["check"]
        report.update(status="EXACT_ORIGINAL_SDP_FEASIBLE",original_point_sha256=full["candidate_sha256"],
            original_equalities_checked=full["equalities_checked"],original_cones_checked=len(full["blocks"]))
    except ValueError as error: report["exact_primal_recovery_failure"]=str(error)
    if report["status"]!="EXACT_ORIGINAL_SDP_FEASIBLE":
        del original,substitutions,point
        matrices=feasibility.native_cones(native,slater_margin.augmented(problem),field="yMat")[:-1]
        try: proposed,checked=dual_completion.complete(problem,problem,matrices,10)
        except ValueError as error: report["source_Farkas_recovery_failure"]=str(error)
        else:
            report["source_Farkas_completion"]=checked
            if proposed is not None:
                certificate,grams=proposed
                report["source_exact_check"]=dual_completion.verify(problem,conic.farkas(problem),certificate,grams)
                conic.write_json(out/"source_infeasibility_certificate.exact.json",certificate)
                conic.write_json(out/"source_PSD_grams.exact.json",grams)
                report.update(status="EXACT_SOURCE_SDP_INFEASIBILITY_CERTIFICATE",
                    original_infeasibility_requires_verified_reduction_chain=True)
    if binding()!=before: raise ValueError("A checkpoint recovery source changed")
    report["source_bindings"]=before;report["elapsed_seconds"]=time.monotonic()-started
    conic.write_json(out/"result.json",report)
    print(json.dumps(conic.plain({k:v for k,v in report.items() if k!="source_bindings"})),flush=True)


if __name__=="__main__": main()
