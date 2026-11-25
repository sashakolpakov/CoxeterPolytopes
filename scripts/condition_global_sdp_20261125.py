#!/usr/bin/env python3
"""Condition a frozen prepared SDP by checked, invertible changes of units."""
import argparse
from fractions import Fraction as F
import json
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get("COXETER_CHECK_OWNER"): raise RuntimeError("Use the serial scientific watchdog")
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",required=True)
    parser.add_argument("--prepared",required=True)
    parser.add_argument("--output",required=True)
    args=parser.parse_args();out=(ROOT/args.output).resolve();source=(ROOT/args.prepared).resolve()
    original=(ROOT/args.input).resolve()
    if out.exists() or any(not p.is_relative_to(ROOT/"results") for p in (out,source,original)): raise ValueError("Use results/ paths and fresh output")
    sys.path.insert(0,str(ROOT))
    from research.global_formulation import conic,conditioning,equilibration,feasibility
    paths=[original,*[source/name for name in ("compiled.exact.json","substitutions.exact.json","result.json")],Path(__file__),
        *[ROOT/"research/global_formulation"/(name+".py") for name in ("conic","conditioning","equilibration","feasibility")]]
    def binding(): return {str(p.relative_to(ROOT)):conic.digest_bytes(p.read_bytes()) for p in paths}
    before=binding();out.mkdir()
    prior=json.loads((source/"result.json").read_text())
    if prior["frozen_input_sha256"]!=before[str(original.relative_to(ROOT))] or prior["orientation"]!="necessary": raise ValueError("Wrong frozen preparation input")
    problem=conic.read_model(source/"compiled.exact.json")
    substitutions={v:{w:F(c) for w,c in q.items()} for v,q in json.loads((source/"substitutions.exact.json").read_text()).items()}
    if conic.digest_bytes(conic.canonical(problem))!=prior["prepared_matrix_sha256"]: raise ValueError("Changed prepared matrix")
    if conic.digest_bytes(conic.canonical(substitutions))!=prior["prepared_substitutions_sha256"]: raise ValueError("Changed prepared substitution")
    problem,change,_,intervals=conditioning.thin_intervals(problem)
    substitutions={v:feasibility.substitute(q,change) for v,q in substitutions.items()}
    # Reference one leaves all columns unchanged; coefficient-based
    # positive diagonal/whole-block scales balance the cone entries.
    problem,change,_,_,units=equilibration.equilibrate(problem,{v:F(1) for v in problem["variables"]})
    substitutions={v:feasibility.substitute(q,change) for v,q in substitutions.items()}
    conic.write_json(out/"equilibrated.exact.json",problem)
    conic.write_json(out/"equilibrated_substitutions.exact.json",substitutions)
    conic.write_json(out/"interval_conditioning.exact.json",intervals)
    conic.write_json(out/"equilibration.exact.json",units)
    report=dict(status="EXACT_INVERTIBLE_SDP_CONDITIONING_COMPLETE",orientation="necessary",
        frozen_input_sha256=before[str(original.relative_to(ROOT))],source_preparation=str(source.relative_to(ROOT)),
        prepared_matrix_sha256=conic.digest_bytes(conic.canonical(problem)),
        prepared_substitutions_sha256=conic.digest_bytes(conic.canonical(substitutions)),
        intervals_conditioned=len(intervals["pairs"]),compiled_variables=len(problem["variables"]),
        compiled_cones=len(problem["blocks"]),geometric_constraints_added=False,accepted_bounds_changed=False,source_bindings=before)
    if binding()!=before: raise ValueError("A conditioning source changed")
    conic.write_json(out/"result.json",report)
    print(json.dumps({k:v for k,v in report.items() if k!="source_bindings"}),flush=True)


if __name__=="__main__": main()
