"""LP/PSD separation on a retained full SDP, with exact original dual checks."""
from pathlib import Path
import json
import time

from . import conic, separation
from .explore import linprog, exact_reconstruction


def run(root,output,input_path,seconds=650):
    if not 1 <= seconds <= 700: raise ValueError("Use at most 700 seconds inside the 900-second watchdog")
    out,source = (root/output).resolve(),(root/input_path).resolve()
    if out.exists() or not out.is_relative_to(root/"results") or not source.is_relative_to(root/"results"):
        raise ValueError("Use a retained exact input and a fresh results output")
    paths = [source,Path(__file__),root/"research/global_formulation/conic.py",
        root/"research/global_formulation/separation.py",root/"research/global_formulation/explore.py",
        root/"scripts/run_global_coxeter_formulation.py"]
    def bindings():
        return {str(p.relative_to(root)):dict(bytes=p.stat().st_size,sha256=conic.digest_bytes(p.read_bytes())) for p in paths}
    before = bindings(); out.mkdir(); started = time.monotonic()
    conic.write_json(out/"source_bindings.json",before)
    primal = conic.read_model(source); dual = conic.farkas(primal)
    print(json.dumps(dict(stage="frozen_full_SDP_loaded",variables=len(primal["variables"]),
        equalities=len(primal["equations"]),cones=len(primal["blocks"]))),flush=True)
    report = separation.explore(primal,dual,out,linprog,exact_reconstruction,seconds=seconds,
        rounds=60,retain_outer=False,iteration_seconds=120)
    report.update(input=str(source.relative_to(root)),elapsed_seconds=time.monotonic()-started,
        source_bindings=before,accepted_bounds_changed=False,
        numerical_LP_feasibility_is_not_full_SDP_feasibility=True)
    if bindings() != before: raise ValueError("A fixed search source or model changed during the run")
    conic.write_json(out/"result.json",report)
    print(json.dumps(conic.plain({k:v for k,v in report.items() if k != "source_bindings"})),flush=True)
