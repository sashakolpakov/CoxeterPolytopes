"""Search a candidate-selected subface; only an original feasible point counts.

Zero-coordinate restrictions in this module are search heuristics, not
necessary conditions for polytopes. Their infeasibility proves nothing
about the original SDP. A successful point is checked in every original
equality and cone before reporting feasibility.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys
import time

from . import conditioning, conic, equilibration, feasibility, sdpa_double


def run(root,out,candidate_path,threshold,seconds,interval_boundaries_only=False,population_only=False):
    out.mkdir();started=time.monotonic()
    source=root/"results/global_shared_cusp_current/necessary_primal.exact.json"
    prepared=root/"results/global_shared_cusp_prepare"
    paths=[source,candidate_path,prepared/"compiled.exact.json",prepared/"substitutions.exact.json",
        Path(__file__),root/"scripts/run_sdpa_candidate_face_20261125.py",
        *[root/"research/global_formulation"/(name+".py") for name in
            ("conic","feasibility","conditioning","equilibration","sdpa_double")]]
    before=sdpa_double.binding(paths,root);conic.write_json(out/"source_bindings.json",before)
    original=conic.read_model(source);problem=conic.read_model(prepared/"compiled.exact.json")
    subs={v:{w:F(c) for w,c in q.items()} for v,q in json.loads((prepared/"substitutions.exact.json").read_text()).items()}
    candidate={v:F(c) for v,c in json.loads(candidate_path.read_text()).items()}
    extra={};selected=[]
    coordinates_to_scan={} if interval_boundaries_only else original["metadata"]["moment_coordinates"]
    for measure,coordinates in coordinates_to_scan.items():
        mass=next((candidate[name] for monomial,name in coordinates if not monomial),F(1))
        cutoff=threshold*max(abs(mass),F(1,10**12))
        for monomial,name in coordinates:
            if len(monomial)!=1 or monomial[0] in ("t","e"):continue
            if population_only and not (monomial[0].startswith("u_") or
                    (monomial[0].startswith("tail_") and monomial[0].endswith("_count"))):continue
            if candidate[name]<=cutoff and subs[name]:
                q=subs[name];extra[conic.canonical(q)]=q
                selected.append(dict(measure=measure,monomial=monomial,coordinate=name,
                    original_candidate_value=candidate[name],cutoff=cutoff))
    interval_boundaries=[]
    if interval_boundaries_only:
        names={b["name"]:b for b in problem["blocks"]}
        for block in problem["blocks"]:
            name=block["name"]
            if "_inverse_response_lower" not in name or len(block["entries"])!=1:continue
            other=names.get(name.replace("_inverse_response_lower","_inverse_response_upper"))
            if other is None:continue
            low=block["entries"][0][0];high=other["entries"][0][0]
            for v in sorted(low.keys()&high.keys()):
                if not v or low[v]*high[v]>=0:continue
                ratio=-low[v]/high[v];width=conic.plus(low,conic.times(ratio,high))
                delta=max((abs(c) for c in width.values()),default=F(0))
                magnitude=max(max(abs(c) for c in low.values()),max(abs(ratio*c) for c in high.values()))
                if not 0<delta<F(1,10**10)*magnitude:continue
                q=conic.times(1/delta,width);extra[conic.canonical(q)]=q
                interval_boundaries.append(dict(lower=name,upper=other["name"],positive_ratio=ratio,
                    width=width,positive_normalizer=delta))
                break
    conic.write_json(out/"candidate_only_restrictions.exact.json",dict(
        necessary_for_original_problem=False,relative_threshold=threshold,selected=selected,
        interval_boundaries=interval_boundaries,equations=list(extra.values())))
    report=dict(status="NO_ORIGINAL_FEASIBILITY_DECISION",candidate_only_restrictions=True,
        restricted_infeasibility_has_no_original_implication=True,selected_coordinates=len(selected),
        interval_boundaries=len(interval_boundaries),population_only=population_only,
        original_input_sha256=before[str(source.relative_to(root))]["sha256"],accepted_bounds_changed=False)
    def progress(stage,data):
        with (out/"progress.jsonl").open("a") as stream:
            stream.write(json.dumps(conic.plain(dict(stage=stage,**data)))+"\n")
    try:
        augmented=conic.model(problem["variables"],problem["blocks"],list(extra.values()))
        problem,updates,_,trace=feasibility.prepare(augmented,progress,linear_faces=False,seconds=220)
        subs={v:feasibility.substitute(q,updates) for v,q in subs.items()}
        conic.write_json(out/"conditional_reductions.exact.json",trace)
        problem,updates,_,intervals=conditioning.thin_intervals(problem)
        subs={v:feasibility.substitute(q,updates) for v,q in subs.items()}
        reference=dict(candidate)
        for pair in intervals["pairs"]:
            reference[pair["coordinate"]]=abs(conic.value(pair["width"],reference))/(2*pair["delta"])
        problem,updates,_,_,units=equilibration.equilibrate(problem,reference)
        subs={v:feasibility.substitute(q,updates) for v,q in subs.items()}
        conic.write_json(out/"equilibrated.exact.json",problem)
        conic.write_json(out/"equilibrated_substitutions.exact.json",subs)
        conic.write_json(out/"conditioning.exact.json",dict(intervals=intervals,units=units))
        report["conditional_free_variables"]=len(problem["variables"])
        report["conditional_cone_blocks"]=len(problem["blocks"])
        conic.write_sdpa(out/"input.dat-s",problem)
        command=[sys.executable,"-B",str(root/"scripts/run_sdpa_double_frozen_20261125.py"),"native",
            "--output",str(out),"--native-input",str(out/"input.dat-s"),"--iterations","80",
            "--initial-scale","0.1","--epsilon","1e-8","--search-region","10000"]
        timeout=min(seconds,max(1,int(790-(time.monotonic()-started))))
        with (out/"solver.stdout").open("x") as stream:
            try:subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,timeout=timeout,check=False)
            except subprocess.TimeoutExpired:report["native_timed_out"]=True
        if (out/"native.result").exists():
            values,phase=conic.native_vector(out/"native.result");report["native_phase"]=phase
            raw=dict(zip(problem["variables"],values,strict=True))
            point={v:conic.value(q,raw) for v,q in subs.items()}
            conic.write_json(out/"original_numerical_candidate.json",point)
            report["original_numerical_residuals"]=feasibility.numerical_residuals(original,point)
            exact,checked=conic.reconstruct(original,problem,subs,values)
            conic.write_json(out/"original_feasible_point.exact.json",exact)
            report.update(status="EXACT_ORIGINAL_SDP_FEASIBLE",exact_check=checked)
    except ValueError as error:report["recovery_or_restricted_search_failure"]=str(error)
    if sdpa_double.binding(paths,root)!=before:raise ValueError("A candidate-face input changed")
    report["elapsed_seconds"]=time.monotonic()-started
    conic.write_json(out/"result.json",report);print(json.dumps(conic.plain(report)),flush=True)
