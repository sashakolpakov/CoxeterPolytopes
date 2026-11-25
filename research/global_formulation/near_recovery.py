"""Candidate-only face proposals; acceptance always checks the frozen SDP."""
from fractions import Fraction as F
import json
from pathlib import Path
import time

from . import conic, feasibility


def run(root,out,source,prepared,candidate_path,cone_closure=False):
    out.mkdir();started=time.monotonic()
    paths=[source,prepared/"compiled.exact.json",prepared/"substitutions.exact.json",candidate_path,
        Path(__file__),root/"research/global_formulation/conic.py",root/"research/global_formulation/feasibility.py",
        root/"scripts/run_sdpa_near_recovery_20261125.py"]
    def bindings():return {str(p.relative_to(root)):dict(bytes=p.stat().st_size,
        sha256=conic.digest_bytes(p.read_bytes())) for p in paths}
    before=bindings();conic.write_json(out/"source_bindings.json",before)
    original=conic.read_model(source);compiled=conic.read_model(prepared/"compiled.exact.json")
    subs={v:{w:F(c) for w,c in q.items()} for v,q in json.loads((prepared/"substitutions.exact.json").read_text()).items()}
    supplied={v:F(c) for v,c in json.loads(candidate_path.read_text()).items()}
    raw={v:supplied[v] for v in compiled["variables"]}
    point={v:conic.value(q,raw) for v,q in subs.items()}
    negative=[]
    for b,block in enumerate(original["blocks"]):
        if len(block["entries"])==1:
            q=block["entries"][0][0];value=conic.value(q,point)
            if value<0:negative.append((b,q,value))
    report=dict(status="NO_EXACT_FEASIBILITY_DECISION",negative_scalar_rows=len(negative),
        original_input_unchanged=True,proposals_are_recovery_restrictions_only=True,
        conditional_cone_kernel_closure=cone_closure,attempts=[],source_bindings=before)
    print(json.dumps(dict(stage="original_candidate_loaded",negative_scalar_rows=len(negative))),flush=True)
    choices=[("all_negative_scalars",negative),
        ("negative_coordinate_rows",[(b,q,v) for b,q,v in negative if len(q)==1 and "" not in q])]
    if cone_closure:choices.reverse()
    for name,selected in choices:
        extra=[feasibility.substitute(q,subs) for _,q,_ in selected]
        try:
            if cone_closure:
                target=conic.model(compiled["variables"],compiled["blocks"],extra)
                def progress(stage,data):
                    print(json.dumps(conic.plain(dict(data,stage=stage,proposal=name))),flush=True)
                refined,updates,_,trace=feasibility.prepare(target,progress,linear_faces=False,seconds=240)
                conic.write_json(out/(name+"_conditional_cone_trace.exact.json"),trace)
            else:refined,updates=conic.eliminate(compiled,extra)
            combined={v:feasibility.substitute(q,updates) for v,q in subs.items()}
            print(json.dumps(dict(stage=name,proposed_equalities=len(extra),free_variables=len(refined["variables"]))),flush=True)
            exact,checked=conic.reconstruct(original,refined,combined,[raw[v] for v in refined["variables"]])
            conic.write_json(out/"original_feasible_point.exact.json",exact)
            conic.write_json(out/"recovery_restrictions.exact.json",extra)
            report.update(status="EXACT_ORIGINAL_SDP_FEASIBLE",exact_check=checked)
            break
        except ValueError as error:report["attempts"].append(dict(proposal=name,reason=str(error)))
    if bindings()!=before:raise ValueError("A frozen recovery source changed")
    report["elapsed_seconds"]=time.monotonic()-started
    conic.write_json(out/"result.json",report)
    print(json.dumps(conic.plain({k:v for k,v in report.items() if k not in ("source_bindings","exact_check")})),flush=True)
