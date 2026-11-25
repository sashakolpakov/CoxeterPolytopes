#!/usr/bin/env python3
"""Guarded exact recovery from the normalized certificate solver's ray."""
import argparse
import os
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
if __name__=="__main__":
    if not os.environ.get("COXETER_CHECK_OWNER"):raise RuntimeError("Use the bounded serial scientific runner")
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output",required=True)
    p.add_argument("--candidate",default="results/global_shared_cusp_sdpa_certificate01/original_numerical_candidate.json")
    p.add_argument("--cone-closure",action="store_true")
    a=p.parse_args();out=(ROOT/a.output).resolve()
    if out.exists() or not out.is_relative_to(ROOT/"results"):raise ValueError("Use a fresh results directory")
    sys.path.insert(0,str(ROOT))
    from research.global_formulation.near_recovery import run
    run(ROOT,out,ROOT/"results/global_shared_cusp_current/necessary_primal.exact.json",
        ROOT/"results/global_shared_cusp_prepare",(ROOT/a.candidate).resolve(),a.cone_closure)
