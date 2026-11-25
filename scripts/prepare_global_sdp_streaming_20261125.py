#!/usr/bin/env python3
"""Prepare one frozen SDP with bounded streaming exact cone reductions."""
import argparse
import json
import os
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get("COXETER_CHECK_OWNER"): raise RuntimeError("Use the serial scientific watchdog")
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",required=True)
    parser.add_argument("--output",required=True)
    args=parser.parse_args();source=(ROOT/args.input).resolve();out=(ROOT/args.output).resolve()
    if out.exists() or any(not p.is_relative_to(ROOT/"results") for p in (source,out)): raise ValueError("Use fresh results/ paths")
    sys.path.insert(0,str(ROOT))
    from research.global_formulation import conic,streaming_presolve
    started=time.monotonic();out.mkdir()
    paths=[source,Path(__file__),*[ROOT/"research/global_formulation"/(name+".py") for name in ("conic","feasibility","streaming_presolve")]]
    def binding(): return {str(p.relative_to(ROOT)):dict(bytes=p.stat().st_size,sha256=conic.digest_bytes(p.read_bytes())) for p in paths}
    before=binding();conic.write_json(out/"source_bindings.json",before)
    def progress(stage,data):
        item=dict(stage=stage,elapsed_seconds=time.monotonic()-started,**data)
        with (out/"progress.jsonl").open("a") as stream: stream.write(json.dumps(conic.plain(item))+"\n")
        print(json.dumps(conic.plain(item)),flush=True)
    checks=streaming_presolve.self_check()
    original=conic.read_model(source)
    original_variables=len(original["variables"])
    problem,substitutions,maps,trace=streaming_presolve.prepare(original,progress)
    del original
    conic.write_json(out/"compiled.exact.json",problem)
    conic.write_json(out/"substitutions.exact.json",substitutions)
    conic.write_json(out/"cone_maps.exact.json",maps)
    conic.write_json(out/"exact_face_reduction.json",trace)
    report=dict(status="EXACT_EQUIVALENT_FROZEN_SDP_PREPARED",orientation="necessary",
        frozen_input=str(source.relative_to(ROOT)),frozen_input_sha256=before[str(source.relative_to(ROOT))]["sha256"],
        original_variables=original_variables,compiled_variables=len(problem["variables"]),compiled_cones=len(problem["blocks"]),
        prepared_matrix_sha256=conic.digest_bytes(conic.canonical(problem)),
        prepared_substitutions_sha256=conic.digest_bytes(conic.canonical(substitutions)),
        geometric_constraints_added=False,accepted_bounds_changed=False,streaming_exact_reduction_checks=checks,
        source_bindings=before,elapsed_seconds=time.monotonic()-started)
    if binding()!=before: raise ValueError("A streaming preparation source changed")
    conic.write_json(out/"result.json",report)
    progress("completed",{k:v for k,v in report.items() if k not in ("source_bindings","elapsed_seconds")})


if __name__=="__main__": main()
