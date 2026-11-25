#!/usr/bin/env python3
"""Serial candidate-only subface search; verify every original SDP row."""
import argparse
from fractions import Fraction
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get("COXETER_CHECK_OWNER"):raise RuntimeError("Use the serial scientific watchdog")
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",required=True)
    parser.add_argument("--threshold",default="1e-6")
    parser.add_argument("--interval-boundaries-only",action="store_true")
    parser.add_argument("--population-only",action="store_true")
    parser.add_argument("--seconds",type=int,default=500)
    parser.add_argument("--candidate",default="results/global_shared_cusp_sdpa_certificate01/original_numerical_candidate.json")
    args=parser.parse_args();threshold=Fraction(args.threshold)
    if not 0<=threshold<=Fraction(1,100) or not 1<=args.seconds<=550:raise ValueError("Unsupported search limits")
    out=(ROOT/args.output).resolve()
    if not out.is_relative_to(ROOT/"results") or out.exists():raise ValueError("Use a fresh results directory")
    sys.path.insert(0,str(ROOT))
    from research.global_formulation import candidate_face
    candidate_face.run(ROOT,out,(ROOT/args.candidate).resolve(),threshold,args.seconds,args.interval_boundaries_only,args.population_only)


if __name__=="__main__":main()
