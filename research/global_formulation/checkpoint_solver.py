"""Continue a frozen native SDPA input through complete restart checkpoints.

Iteration limits are checkpoint boundaries, not feasibility decisions.
The numerical model, objective and starting-state conventions are unchanged.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import re
import time

from . import conic,feasibility,sdpa_cli

TERMINAL_PHASES={"pdOPT","pdINF","pINF_dFEAS","pFEAS_dINF","pUNBD","dUNBD"}


def shape(path):
    with path.open() as stream:
        comment=stream.readline()
        variables=int(stream.readline());blocks=int(stream.readline())
        sizes=[int(s) for s in stream.readline().split()]
    if not comment.startswith('"Global Coxeter SDP') or len(sizes)!=blocks or not all(sizes):
        raise ValueError("Unexpected native input header")
    return variables,sizes


def checkpoint(path,variables,sizes):
    values,phase=conic.native_vector(path)
    if len(values)!=variables: raise ValueError("Wrong checkpoint free-vector dimension")
    expected=sum(n*n if n>0 else -n for n in sizes)
    lengths={}
    for field in ("xMat","yMat"):
        lengths[field]=len(conic.NUMBER.findall(feasibility.native_array(path,field)))
        if lengths[field]!=expected: raise ValueError("Wrong checkpoint cone dimension: "+field)
    return dict(native_phase=phase,free_variables=len(values),cone_array_entries=lengths,
        native_common_margin=str(values[-1]),complete_free_vector_and_both_cone_arrays=True,
        file_sha256=conic.digest_bytes(path.read_bytes()))


def run(root,out,input_path,initial,iterations,seconds,resume=None):
    started=time.monotonic();out.mkdir()
    paths=[input_path,initial,Path(__file__),root/"scripts/run_sdpa_checkpoint_20261125.py",
        sdpa_cli.EXECUTABLE,*[root/"research/global_formulation"/(name+".py")
        for name in ("conic","feasibility","sdpa_cli")]]
    if resume is not None: paths.extend([resume/"native_metadata.json",resume/"result.json"])
    def bindings():
        return {str(p.relative_to(root)) if p.is_relative_to(root) else str(p):
            dict(bytes=p.stat().st_size,sha256=conic.digest_bytes(p.read_bytes())) for p in paths}
    before=bindings();conic.write_json(out/"source_bindings.json",before)
    variables,sizes=shape(input_path)
    if resume is not None:
        previous=json.loads((resume/"native_metadata.json").read_text())
        if previous["input_sha256"]!=before[str(input_path.relative_to(root))]["sha256"]:
            raise ValueError("A restart requires the identical ordered SDP input")
        previous_state=json.loads((resume/"result.json").read_text())
        if not previous_state.get("checkpoint",{}).get("complete_free_vector_and_both_cone_arrays"):
            raise ValueError("The previous batch did not produce a complete restart point")
    initial_check=checkpoint(initial,variables,sizes)
    conic.write_json(out/"initial_checkpoint_check.json",initial_check)
    native=sdpa_cli.call(input_path,out,iterations,seconds,initial)
    stdout=(out/"solver.stdout").read_text()
    report=dict(status="NO_EXACT_FEASIBILITY_DECISION",native=native,
        requested_iterations=iterations,requested_native_seconds=seconds,
        identical_frozen_native_input=True,geometric_constraints_changed=False,
        auxiliary_problem="unchanged common PSD margin",accepted_bounds_changed=False,
        resume_from=str(resume.relative_to(root)) if resume is not None else None,
        native_terminal=False,source_bindings=before)
    if native["timed_out"]:
        report.update(status="TIME_LIMIT_BEFORE_CHECKPOINT_RETRY_FROM_LAST_SAVED_STATE")
    elif native["returncode"]!=0:
        report.update(status="NATIVE_SOLVER_PROCESS_ERROR",native_terminal=True)
    else:
        try: saved=checkpoint(out/"native.result",variables,sizes)
        except (ValueError,FileNotFoundError) as error:
            report.update(status="NATIVE_OUTPUT_REQUIRES_READER_DIAGNOSIS",reader_error=str(error))
        else:
            report["checkpoint"]=saved
            errors=[line for line in stdout.splitlines() if re.search(
                r"cannot move:|Strange behavior|Cholesky.*failed|\bERROR\b",line)]
            exhausted="maxIteration is reached" in stdout
            report.update(iteration_limit_reached=exhausted,numerical_error_messages=errors)
            if saved["native_phase"] in TERMINAL_PHASES:
                report.update(status="NATIVE_TERMINAL_STATUS_REQUIRES_ORIGINAL_CHECK",native_terminal=True)
            elif errors:
                report.update(status="NATIVE_SOLVER_NUMERICAL_STOP",native_terminal=True)
            elif exhausted:
                report.update(status="COMPLETE_CHECKPOINT_CONTINUE")
            else:
                report.update(status="NATIVE_NONCONVERGED_STOP",native_terminal=True)
    if bindings()!=before: raise ValueError("A frozen checkpoint source changed")
    report["elapsed_seconds"]=time.monotonic()-started
    conic.write_json(out/"result.json",report)
    print(json.dumps(conic.plain({k:v for k,v in report.items() if k!="source_bindings"})),flush=True)
