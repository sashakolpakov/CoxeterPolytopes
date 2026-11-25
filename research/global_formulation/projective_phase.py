"""Bounded objective for homogeneous necessary-data feasibility.

For a nonnegative free coordinate g, put h=1/(1+g), y=h*x and
y_g=1-h. Maximizing h distinguishes an original feasible point (h>0)
from recession data (h=0). Each original cone is multiplied by h;
no original constraint is omitted. An h=0 answer is never a witness.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

from . import conic, feasibility, sdpa_double

H="_projective_positive_scale"


def transform(problem,pivot):
    if problem["equations"] or pivot not in problem["variables"] or H in problem["variables"]:
        raise ValueError("Use an equation-free model and a registered gauge coordinate")
    if not any(len(b["entries"])==1 and len(b["entries"][0][0])==1 and
            b["entries"][0][0].get(pivot,F(0))>0 for b in problem["blocks"]):
        raise ValueError("The gauge coordinate must have an original pure nonnegative row")
    def entry(q):
        old=q.get(pivot,F(0));constant=q.get("",F(0))
        result={v:c for v,c in q.items() if v and v!=pivot}
        if old:result[""]=old
        if constant-old:result[H]=constant-old
        # Coefficient identity: M_new = h*M_old(x), x=y/h,
        # x_pivot=(1-h)/h. This is checked for every entry.
        recovered={v:c for v,c in result.items() if v and v!=H}
        if result.get("",F(0)):recovered[pivot]=result[""]
        if result.get(H,F(0))+result.get("",F(0)):
            recovered[""]=result.get(H,F(0))+result.get("",F(0))
        if recovered!=q:raise ValueError("A homogenization coefficient identity failed")
        return result
    blocks=[dict(b,entries=[[entry(q) for q in row] for row in b["entries"]]) for b in problem["blocks"]]
    blocks.append(dict(name="projective_nonnegative_scale",entries=[[conic.form(**{H:1})]]))
    return conic.model([v for v in problem["variables"] if v!=pivot]+[H],blocks,
        metadata=dict(original_model_sha256=conic.digest_bytes(conic.canonical(problem)),gauge_coordinate=pivot,
            positive_scale_equivalence=True,zero_scale_not_an_original_point=True,
            every_affine_homogenization_identity_checked=True))


def write_input(path,problem):
    conic.write_sdpa(path,problem)
    temporary=path.with_suffix(".partial")
    with path.open() as source,temporary.open("x") as target:
        for _ in range(4):target.write(source.readline())
        source.readline()
        target.write(" ".join("-1" if v==H else "0" for v in problem["variables"])+"\n")
        shutil.copyfileobj(source,target,length=1024**2)
    temporary.replace(path)


def inverse(values,pivot):
    h=values[H]
    if h<=0:raise ValueError("A nonpositive projective scale cannot give an original point")
    result={v:c/h for v,c in values.items() if v!=H};result[pivot]=(1-h)/h
    return result


def recover(original,problem,substitutions,variables,raw,pivot):
    attempts=[]
    for digits in (6,10,16,24,36,48,60):
        grid=10**digits;values={v:F(round(c*grid),grid) for v,c in zip(variables,raw,strict=True)}
        try:
            old=inverse(values,pivot)
            point={v:conic.value(q,old) for v,q in substitutions.items()}
            checked=conic.check_point(original,point)
            return point,dict(method="exact projective inverse and all original cone checks",decimal_grid_digits=digits,
                positive_scale=values[H],check=checked)
        except ValueError as error:attempts.append(dict(digits=digits,reason=str(error)))
    raise ValueError("No exact projective recovery: "+json.dumps(attempts))


def self_check(out):
    out.mkdir()
    P=conic.model(["x"],[dict(name="nonnegative_x",entries=[[conic.form(x=1)]]),
        dict(name="PSD",entries=[[conic.form(x=1),conic.form(1)],[conic.form(1),conic.form(x=1)]]),
        dict(name="upper_x",entries=[[conic.form(3,x=-1)]])])
    Q=transform(P,"x");conic.check_point(Q,{H:F(1,3)})
    if inverse({H:F(1,3)},"x")!={"x":F(2)}:raise ValueError("Wrong projective inverse")
    write_input(out/"input.dat-s",Q);sdpa_double.native_call(out/"input.dat-s",out,80,1.,1e-9,10000.)
    raw,phase=conic.native_vector(out/"native.result")
    if abs(raw[0]-F(1,2))>F(1,10**7):raise ValueError("Wrong projective objective optimum")
    point,checked=recover(P,P,{"x":conic.form(x=1)},Q["variables"],raw,"x")
    weak=conic.model(["x"],[dict(name="nonnegative_x",entries=[[conic.form(x=1)]]),
        dict(name="weakly_infeasible_PSD",entries=[[conic.form(x=1),conic.form(1)],[conic.form(1),{}]])])
    conic.check_point(transform(weak,"x"),{H:F(0)})
    try:inverse({H:F(0)},"x")
    except ValueError:pass
    else:raise ValueError("A recession point was accepted as an original point")
    report=dict(status="EXACT_PROJECTIVE_FIXTURES_PASSED",phase=phase,
        all_affine_coefficient_identities_checked=True,known_objective_checked=True,
        zero_scale_recession_point_rejected=True,exact_original_point=point)
    conic.write_json(out/"result.json",report);print(json.dumps(conic.plain(report)),flush=True)


def run(root,out,source,seconds):
    out.mkdir();started=time.monotonic()
    original_path=root/"results/global_shared_cusp_current/necessary_primal.exact.json"
    model_path=source/"equilibrated.exact.json";subs_path=source/"equilibrated_substitutions.exact.json"
    paths=[original_path,model_path,subs_path,Path(__file__),root/"scripts/run_sdpa_projective_20261125.py",
        root/"research/global_formulation/conic.py",root/"research/global_formulation/feasibility.py",
        root/"research/global_formulation/sdpa_double.py"]
    before=sdpa_double.binding(paths,root);conic.write_json(out/"source_bindings.json",before)
    original=conic.read_model(original_path);problem=conic.read_model(model_path)
    substitutions={v:{w:F(c) for w,c in q.items()} for v,q in json.loads(subs_path.read_text()).items()}
    pivot="shared_edge_support_m0";projective=transform(problem,pivot)
    conic.write_json(out/"projective_mapping.exact.json",dict(gauge_coordinate=pivot,
        old_substitutions_sha256=conic.digest_bytes(conic.canonical(substitutions)),
        projective_model_sha256=conic.digest_bytes(conic.canonical(projective)),
        original_model_sha256=conic.digest_bytes(conic.canonical(problem)),positive_scale_required=True,
        every_affine_coefficient_identity_checked=True))
    write_input(out/"input.dat-s",projective)
    command=[sys.executable,"-B",str(root/"scripts/run_sdpa_double_frozen_20261125.py"),"native",
        "--output",str(out),"--native-input",str(out/"input.dat-s"),"--iterations","100",
        "--initial-scale","1","--epsilon","1e-9","--search-region","10000"]
    timed_out=False
    with (out/"solver.stdout").open("x") as stream:
        try:subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,timeout=seconds,check=False)
        except subprocess.TimeoutExpired:timed_out=True
    report=dict(status="NO_EXACT_FEASIBILITY_DECISION",bounded_projective_objective=True,timed_out=timed_out,
        source_model=str(source.relative_to(root)),original_input_sha256=before[str(original_path.relative_to(root))]["sha256"],
        zero_scale_not_accepted_as_original_feasibility=True,accepted_bounds_changed=False)
    if (out/"native.result").exists():
        raw,phase=conic.native_vector(out/"native.result");report["native_phase"]=phase
        if len(raw)!=len(projective["variables"]):raise ValueError("Wrong projective vector length")
        report["numerical_positive_scale"]=raw[-1]
        if raw[-1]>0:
            values=inverse(dict(zip(projective["variables"],raw,strict=True)),pivot)
            point={v:conic.value(q,values) for v,q in substitutions.items()}
            conic.write_json(out/"original_numerical_candidate.json",point)
            report["original_numerical_residuals"]=feasibility.numerical_residuals(original,point)
            try:
                exact,checked=recover(original,problem,substitutions,projective["variables"],raw,pivot)
                conic.write_json(out/"original_feasible_point.exact.json",exact)
                report.update(status="EXACT_ORIGINAL_SDP_FEASIBLE",exact_check=checked)
            except ValueError as error:report["exact_recovery_failure"]=str(error)
    if sdpa_double.binding(paths,root)!=before:raise ValueError("A projective solver input changed")
    report["elapsed_seconds"]=time.monotonic()-started
    conic.write_json(out/"result.json",report)
    print(json.dumps(conic.plain({k:v for k,v in report.items() if k!="exact_check"})),flush=True)
