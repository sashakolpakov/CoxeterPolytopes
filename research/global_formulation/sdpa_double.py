"""Standard SDPA on an unchanged frozen input, with exact original checks.

The optional SDPA Python wheel is installed outside the repository. Its
documented sparse reader and low-level SeDuMi call preserve the native free
coordinates; the cone blocks are put back in the original native order.
"""
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from . import conic, feasibility
from . import dual_completion

PACKAGE = Path.home()/".local/opt/coxeter-sdpa-python-0.2.3"


def backend():
    sys.path.insert(0,str(PACKAGE))
    import sdpap
    from sdpap.sdpacall import sdpacall
    return sdpap,sdpacall


def binding(paths,root):
    return {str(p.relative_to(root) if p.is_relative_to(root) else p):dict(
        bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths}


def certificate_input(path,problem):
    """Native primal trace equations for the normalized certificate problem.

    For F0+sum x_i Fi >=0, seek Z>=0 with <Fi,Z>=0 and
    <F0,Z>=-1. Native matrices are -Fi,-F0, and b=(0,...,0,1).
    This uses m+1 native rows, rather than eliminating the cone weights.
    """
    variables={v:i+1 for i,v in enumerate(problem["variables"])}
    last=len(variables)+1
    full=[b for b in problem["blocks"] if len(b["entries"])>1]
    scalar=[b for b in problem["blocks"] if len(b["entries"])==1]
    sizes=[len(b["entries"]) for b in full]+([-len(scalar)] if scalar else [])
    with path.open("x") as stream:
        stream.write('"Normalized Farkas trace formulation; exact original is separate\n')
        stream.write(f'{last}\n{len(sizes)}\n'+" ".join(map(str,sizes))+"\n")
        stream.write(" ".join(["0"]*len(variables)+["1"])+"\n")
        def entry(block,i,j,q):
            for v,c in q.items():
                if c:stream.write(f'{variables[v] if v else last} {block} {i} {j} {conic.decimal_text(-c)}\n')
        for b,block in enumerate(full,1):
            for i,row in enumerate(block["entries"]):
                for j in range(i,len(row)):entry(b,i+1,j+1,row[j])
        for i,block in enumerate(scalar,1):entry(len(full)+1,i,i,block["entries"][0][0])


def phase1_input(path,problem):
    """Minimize rho>=0 with every original cone shifted by rho*I.

    At rho=0 this is precisely the unshifted necessary-data problem.
    Native cone weights have trace at most one, avoiding an unbounded
    zero-objective multiplier search during the feasibility phase.
    """
    rho="_phase1_common_violation"
    if rho in problem["variables"]:raise ValueError("Reserved auxiliary coordinate")
    blocks=[dict(name=b["name"],entries=[[conic.plus(q,conic.form(**{rho:1})) if i==j else q
        for j,q in enumerate(row)] for i,row in enumerate(b["entries"])]) for b in problem["blocks"]]
    blocks.append(dict(name="phase1_nonnegative_violation",entries=[[conic.form(**{rho:1})]]))
    augmented=conic.model([*problem["variables"],rho],blocks)
    conic.write_sdpa(path,augmented)
    # The final free coordinate has objective coefficient +1. All other
    # bytes are the ordinary checked affine SDPA encoding.
    with path.open("r+b") as stream:
        for _ in range(4):stream.readline()
        objective=stream.readline()
        if not objective.endswith(b"0\n"):raise ValueError("Unexpected native objective encoding")
        stream.seek(stream.tell()-2);stream.write(b"1")
    return augmented


def native_call(input_path,out,iterations):
    import numpy as np
    sdpap,api = backend()
    A,b,c,K,J = sdpap.importsdpa(str(input_path),flipsign=True)
    if K.f or K.q or J.l or J.q or J.s:
        raise ValueError("Only the native linear/PSD product cone is intended")
    options = sdpap.param(dict(maxIteration=iterations,numThreads=1,lambdaStar=10.,
        omegaStar=2.,epsilonStar=1e-9,epsilonDash=1e-9,lowerBound=-1e12,upperBound=1e12,
        betaStar=.1,betaBar=.3,gammaStar=.9,print="display",isDimacs=False,
        sdpaResult=str(out/"backend.result"),xPrint="%+26.17e",yPrint="%+26.17e",
        sPrint="%+26.17e",infPrint="%+26.17e"))
    meta = dict(backend=api.get_backend_info(),input=str(input_path),
        input_sha256=hashlib.sha256(input_path.read_bytes()).hexdigest(),
        dimensions=dict(equalities=A.shape[0],cone_coordinates=A.shape[1],
            scalar_rows=K.l,PSD_blocks=list(K.s)),options=options,
        cold_start=True,source="https://github.com/sdpa-python/sdpa-python")
    conic.write_json(out/"backend_metadata.json",meta)
    print(json.dumps(dict(stage="standard_SDPA_start",variables=A.shape[0],
        scalar_rows=K.l,PSD_blocks=len(K.s),maximum_iterations=iterations)),flush=True)
    x,y,s,info = api.solve_sdpa(A,b,c,K,options)
    # The extension retains its C FILE stream until process teardown.
    # Flush it before inspecting the native answer in this process.
    import ctypes
    libc=ctypes.CDLL(None)
    libc.fflush.argtypes=[ctypes.c_void_p]
    if libc.fflush(None)!=0:raise OSError("The native SDPA output did not flush")
    # The native wrapper stores linear-cone entries first; original native
    # files put full PSD blocks first. x=YMat, y=XVec, s=XMat.
    def cone_text(vector):
        values = np.asarray(vector.toarray()).ravel()
        blocks=[]; position=K.l
        for size in K.s:
            M=values[position:position+size*size].reshape((size,size),order="F")
            blocks.append("{"+",\n".join("{"+", ".join(format(v,"+.17e") for v in row)+"}" for row in M)+"}")
            position += size*size
        if K.l: blocks.append("{"+", ".join(format(v,"+.17e") for v in values[:K.l])+"}")
        if position != len(values): raise ValueError("Unexpected standard SDPA cone size")
        return "{"+",\n".join(blocks)+"}"
    raw,phase=conic.native_vector(out/"backend.result")
    returned=np.asarray(y.toarray()).ravel()
    if len(raw) != len(returned) or not np.allclose([float(v) for v in raw],returned,rtol=1e-14,atol=1e-14):
        raise ValueError("The returned free vector differs from the native free coordinates")
    with (out/"native.result").open("x") as stream:
        stream.write("phase.value = "+str(phase)+"\n")
        stream.write("xVec = {"+", ".join(format(v,"+.17e") for v in returned)+"}\n")
        stream.write("xMat = "+cone_text(s)+"\nyMat = "+cone_text(x)+"\n")
    conic.write_json(out/"backend_info.json",info)
    print(json.dumps(dict(stage="standard_SDPA_completed",phase=phase,info=info)),flush=True)


def checks(root,out):
    out.mkdir()
    sdpap,api=backend()
    sources=[Path(__file__),root/"scripts/run_sdpa_double_frozen_20261125.py",
        root/"research/global_formulation/conic.py",root/"research/global_formulation/feasibility.py",
        root/"research/global_formulation/dual_completion.py"]
    sources+=sorted(PACKAGE.rglob("*.so"))+sorted(PACKAGE.rglob("METADATA"))
    before=binding(sources,root);conic.write_json(out/"source_bindings.json",before)
    toy=conic.model(["x"],[dict(name="PSD",entries=[[conic.form(2),conic.form(x=1)],
        [conic.form(x=1),conic.form(2)]]),dict(name="lower",entries=[[conic.form(-1,x=1)]]),
        dict(name="upper",entries=[[conic.form(1,x=-1)]])])
    conic.write_sdpa(out/"input.dat-s",toy)
    native_call(out/"input.dat-s",out,60)
    values,phase=conic.native_vector(out/"native.result")
    if abs(values[0]-1)>F(1,10**7): raise ValueError("Native import changed the free-vector sign")
    # A point with a unique free coordinate checks sign, a mixed scalar/PSD
    # cone checks order, and full entries check the off-diagonal convention.
    matrices=feasibility.native_cones(out/"native.result",toy,field="xMat")
    for block,M in zip(toy["blocks"],matrices,strict=True):
        expected=[[conic.value(q,{"x":values[0]}) for q in row] for row in block["entries"]]
        if max(abs(float(a-b)) for row,other in zip(M,expected,strict=True) for a,b in zip(row,other,strict=True))>1e-6:
            raise ValueError("The canonical native cone order or sign changed")
    report=dict(status="STANDARD_SDPA_NATIVE_ENCODING_CHECK_PASSED",phase=phase,
        returned_coordinate_close_to_unique_value=True,scalar_and_PSD_order_checked=True,
        off_diagonal_entries_checked=True,backend=api.get_backend_info(),source_bindings=before)
    certificate=out/"certificate_fixture";certificate.mkdir()
    P=conic.model(["x","z","w"],[
        dict(name="PSD",entries=[[{},conic.form(x=1)],[conic.form(x=1),{}]]),
        dict(name="nonnegative_z",entries=[[conic.form(z=1)]]),
        dict(name="nonnegative_w",entries=[[conic.form(w=1)]])],
        [conic.form(-1,x=1,z=-2),conic.form(2,w=1,x=-3,z=1)])
    E,subs=conic.eliminate(P)
    certificate_input(certificate/"input.dat-s",E)
    native_call(certificate/"input.dat-s",certificate,80)
    matrices=feasibility.native_cones(certificate/"native.result",E,field="yMat")
    proposed,completion=dual_completion.complete(P,E,matrices,8)
    if proposed is None:raise ValueError("Native certificate completion failed: "+str(completion))
    point,grams=proposed
    report["normalized_certificate_trace_and_affine_equality_lift"]=dual_completion.verify(P,conic.farkas(P),point,grams)
    conic.write_json(certificate/"certificate.exact.json",point)
    conic.write_json(certificate/"PSD_grams.exact.json",grams)
    phase1=out/"phase1_fixture";phase1.mkdir()
    P1=conic.model(["x"],[dict(name="PSD",entries=[[{},conic.form(x=1)],[conic.form(x=1),{}]]),
        dict(name="lower",entries=[[conic.form(-1,x=1)]]),dict(name="nonnegative_x",entries=[[conic.form(x=1)]])])
    phase1_input(phase1/"input.dat-s",P1)
    native_call(phase1/"input.dat-s",phase1,80)
    raw,phase=conic.native_vector(phase1/"native.result")
    if abs(raw[-1]-F(1,2))>F(1,10**7):raise ValueError("The phase-I objective does not equal its known optimum")
    try:conic.check_point(P1,{"x":raw[0]})
    except ValueError:pass
    else:raise ValueError("The shifted auxiliary point was accepted as an original point")
    report["phase1_known_positive_optimum_and_unshifted_rejection"]=True
    if binding(sources,root)!=before:raise ValueError("A native backend-check input changed")
    conic.write_json(out/"result.json",report)
    print(json.dumps({k:v for k,v in report.items() if k!="source_bindings"}),flush=True)


def run(root,out,input_path,prepared,input_native,iterations,seconds,orientation="necessary"):
    out.mkdir(); started=time.monotonic()
    inputs=[input_path,prepared/"compiled.exact.json",prepared/"substitutions.exact.json",prepared/"result.json",input_native,
        Path(__file__),root/"scripts/run_sdpa_double_frozen_20261125.py",
        root/"research/global_formulation/conic.py",root/"research/global_formulation/feasibility.py",
        root/"research/global_formulation/dual_completion.py",input_path.with_name("certificate_dual.exact.json")]
    inputs += sorted(PACKAGE.rglob("*.so"))+sorted(PACKAGE.rglob("METADATA"))
    before=binding(inputs,root);conic.write_json(out/"source_bindings.json",before)
    original=conic.read_model(input_path); compiled=conic.read_model(prepared/"compiled.exact.json")
    prior=json.loads((prepared/"result.json").read_text())
    if prior["orientation"]!="necessary" or prior["frozen_input_sha256"]!=hashlib.sha256(input_path.read_bytes()).hexdigest():
        raise ValueError("Wrong original frozen problem")
    if prior["prepared_matrix_sha256"]!=conic.digest_bytes(conic.canonical(compiled)):
        raise ValueError("Changed exact prepared matrix")
    substitutions={v:{w:F(c) for w,c in q.items()} for v,q in json.loads((prepared/"substitutions.exact.json").read_text()).items()}
    if prior["prepared_substitutions_sha256"]!=conic.digest_bytes(conic.canonical(substitutions)):
        raise ValueError("Changed exact substitutions")
    native_problem=compiled;native_substitutions=substitutions
    if orientation=="necessary":
        conic.write_sdpa(out/"input.dat-s",compiled)
        if (out/"input.dat-s").read_bytes()!=input_native.read_bytes():
            raise ValueError("The standard SDPA run is not the identical ordered reduced input")
    else:
        native_problem,native_substitutions=conic.eliminate(original)
        conic.write_json(out/"equality_only.exact.json",native_problem)
        conic.write_json(out/"equality_only_substitutions.exact.json",native_substitutions)
        if orientation=="phase1":phase1_input(out/"input.dat-s",native_problem)
        else:certificate_input(out/"input.dat-s",native_problem)
        input_native=out/"input.dat-s"
    command=[sys.executable,"-B",str(root/"scripts/run_sdpa_double_frozen_20261125.py"),"native",
        "--output",str(out.relative_to(root)),"--native-input",str(input_native),"--iterations",str(iterations)]
    timed_out=False;code=None
    with (out/"solver.stdout").open("x") as stream:
        try: code=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,timeout=seconds,check=False).returncode
        except subprocess.TimeoutExpired: timed_out=True
    report=dict(status="NO_EXACT_FEASIBILITY_DECISION",solver="SDPA binary64 through sdpa-python 0.2.3",
        timed_out=timed_out,returncode=code,maximum_iterations=iterations,source_bindings=before,
        original_frozen_input=str(input_path.relative_to(root)),original_input_sha256=hashlib.sha256(input_path.read_bytes()).hexdigest(),
        orientation=orientation,original_ordered_reduced_input_byte_identity_checked=orientation=="necessary",
        certificate_trace_uses_only_original_equality_elimination=orientation=="certificate",
        phase1_shifts_are_auxiliary_and_not_allowed_in_original_verification=orientation=="phase1",accepted_bounds_changed=False)
    if (out/"native.result").exists():
        values,phase=conic.native_vector(out/"native.result")
        expected=len(native_problem["variables"])+int(orientation!="necessary")
        if len(values)!=expected:raise ValueError("Wrong standard SDPA free-vector dimension")
        report["native_phase"]=phase
        if orientation=="phase1":report["phase1_violation"]=values[-1]
        geometric_values=(values if orientation=="necessary" else values[:-1] if orientation=="phase1" else
            ([v/values[-1] for v in values[:-1]] if values[-1]<0 else None))
        if geometric_values is not None:
            raw=dict(zip(native_problem["variables"],geometric_values,strict=True))
            point={v:conic.value(q,raw) for v,q in native_substitutions.items()}
            conic.write_json(out/"original_numerical_candidate.json",point)
            report["original_numerical_residuals"]=feasibility.numerical_residuals(original,point)
            try:
                exact,checked=conic.reconstruct(original,compiled,substitutions,[point[v] for v in compiled["variables"]])
                conic.write_json(out/"original_feasible_point.exact.json",exact)
                report.update(status="EXACT_ORIGINAL_SDP_FEASIBLE",exact_check=checked)
            except ValueError as error:report["exact_recovery_failure"]=str(error)
        if orientation!="necessary" and report["status"]!="EXACT_ORIGINAL_SDP_FEASIBLE":
            if orientation=="phase1":
                shapes=conic.model(native_problem["variables"],[*native_problem["blocks"],dict(name="rho",entries=[[{}]])])
                matrices=feasibility.native_cones(out/"native.result",shapes,field="yMat")[:-1]
            else:matrices=feasibility.native_cones(out/"native.result",native_problem,field="yMat")
            report["certificate_completion_attempts"]=[]
            for digits in (6,10,14):
                proposed,completion=dual_completion.complete(original,native_problem,matrices,digits)
                report["certificate_completion_attempts"].append(completion)
                if proposed is None:continue
                point,grams=proposed
                checked=dual_completion.verify(original,conic.read_model(input_path.with_name("certificate_dual.exact.json")),point,grams)
                conic.write_json(out/"original_infeasibility_certificate.exact.json",point)
                conic.write_json(out/"PSD_grams.exact.json",grams)
                report.update(status="EXACT_ORIGINAL_SDP_INFEASIBILITY_CERTIFICATE",exact_check=checked)
                break
    if binding(inputs,root)!=before:raise ValueError("A fixed standard SDPA source or input changed")
    report["elapsed_seconds"]=time.monotonic()-started
    conic.write_json(out/"result.json",report)
    print(json.dumps(conic.plain({k:v for k,v in report.items() if k not in ("source_bindings","exact_check")})),flush=True)
