#!/usr/bin/env python3
"""Check projective cone weights as a certificate in its source SDP.

This checks the retained reduced source. An original-model claim also
requires its recorded exact reduction chain; a restricted recovery source
does not imply original infeasibility.
"""
import argparse
from fractions import Fraction as F
import json
import os
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get("COXETER_CHECK_OWNER"):raise RuntimeError("Use the serial scientific watchdog")
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",required=True)
    parser.add_argument("--source",default="results/global_shared_cusp_sdpa_equilibrated01")
    parser.add_argument("--native",default="results/global_shared_cusp_sdpa_projective01/native.result")
    args=parser.parse_args();out=(ROOT/args.output).resolve()
    if not out.is_relative_to(ROOT/"results") or out.exists():raise ValueError("Use a fresh results directory")
    out.mkdir();started=time.monotonic();sys.path.insert(0,str(ROOT))
    from research.global_formulation import conic,dual_completion,feasibility,projective_phase,sdpa_double
    model_path=ROOT/args.source/"equilibrated.exact.json";native=ROOT/args.native
    paths=[model_path,native,Path(__file__),ROOT/"research/global_formulation/conic.py",
        ROOT/"research/global_formulation/dual_completion.py",ROOT/"research/global_formulation/feasibility.py",
        ROOT/"research/global_formulation/projective_phase.py"]
    before=sdpa_double.binding(paths,ROOT);conic.write_json(out/"source_bindings.json",before)
    problem=conic.read_model(model_path)
    transformed=projective_phase.transform(problem,"shared_edge_support_m0")
    matrices=feasibility.native_cones(native,transformed,field="yMat")[:-1]
    report=dict(status="NO_EXACT_REDUCED_CERTIFICATE",source_model=str(model_path.relative_to(ROOT)),
        reduced_model_sha256=conic.digest_bytes(conic.canonical(problem)),
        original_infeasibility_requires_reduction_chain_review=True,accepted_bounds_changed=False,attempts=[])
    for digits in (6,10,14):
        proposed,checked=dual_completion.complete(problem,problem,matrices,digits)
        report["attempts"].append(checked)
        if proposed is None:continue
        point,grams=proposed;checked=dual_completion.verify(problem,conic.farkas(problem),point,grams)
        conic.write_json(out/"source_infeasibility_certificate.exact.json",point)
        conic.write_json(out/"source_PSD_grams.exact.json",grams)
        report.update(status="EXACT_SOURCE_SDP_INFEASIBILITY_CERTIFICATE",exact_check=checked)
        break
    if sdpa_double.binding(paths,ROOT)!=before:raise ValueError("A projective certificate input changed")
    report["elapsed_seconds"]=time.monotonic()-started
    conic.write_json(out/"result.json",report);print(json.dumps(conic.plain(report)),flush=True)


if __name__=="__main__":main()
