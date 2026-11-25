#!/usr/bin/env python3
"""Serial standard-SDPA backend checks and solves of retained exact inputs."""
import argparse
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get("COXETER_CHECK_OWNER"):raise RuntimeError("Use the serial bounded scientific runner")
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage",choices=("check","solve","native"))
    parser.add_argument("--output",required=True)
    parser.add_argument("--input",default="results/global_shared_cusp_current/necessary_primal.exact.json")
    parser.add_argument("--prepared",default="results/global_shared_cusp_prepare")
    parser.add_argument("--native-input",default="results/global_shared_cusp_sdpa/input.dat-s")
    parser.add_argument("--iterations",type=int,default=100)
    parser.add_argument("--orientation",choices=("necessary","certificate"),default="necessary")
    parser.add_argument("--seconds",type=int,default=700)
    args=parser.parse_args()
    if not 1<=args.iterations<=500 or not 1<=args.seconds<=750:raise ValueError("Use bounded native stages")
    out=(ROOT/args.output).resolve()
    if not out.is_relative_to(ROOT/"results") or (args.stage!="native" and out.exists()):
        raise ValueError("Use a fresh retained results directory")
    sys.path.insert(0,str(ROOT))
    from research.global_formulation import sdpa_double
    if args.stage=="check":sdpa_double.checks(ROOT,out)
    elif args.stage=="native":sdpa_double.native_call((ROOT/args.native_input).resolve(),out,args.iterations)
    else:sdpa_double.run(ROOT,out,(ROOT/args.input).resolve(),(ROOT/args.prepared).resolve(),
        (ROOT/args.native_input).resolve(),args.iterations,args.seconds,args.orientation)


if __name__=="__main__":main()
