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
    parser.add_argument("--orientation",choices=("necessary","certificate","phase1"),default="necessary")
    parser.add_argument("--seconds",type=int,default=700)
    parser.add_argument("--condition-intervals",action="store_true")
    parser.add_argument("--equilibrate",action="store_true")
    parser.add_argument("--reference",default="results/global_shared_cusp_sdpa_certificate01/original_numerical_candidate.json")
    parser.add_argument("--initial-scale",type=float,default=10.)
    parser.add_argument("--epsilon",type=float,default=1e-9)
    parser.add_argument("--search-region",type=float,default=2.)
    args=parser.parse_args()
    if not 1<=args.iterations<=500 or not 1<=args.seconds<=750:raise ValueError("Use bounded native stages")
    if not 1e-12<=args.initial_scale<=1e6 or not 1e-15<=args.epsilon<=1e-4:raise ValueError("Unsupported numerical parameters")
    if not 1<=args.search_region<=1e12:raise ValueError("Unsupported numerical search region")
    out=(ROOT/args.output).resolve()
    if not out.is_relative_to(ROOT/"results") or (args.stage!="native" and out.exists()):
        raise ValueError("Use a fresh retained results directory")
    sys.path.insert(0,str(ROOT))
    from research.global_formulation import sdpa_double
    if args.stage=="check":sdpa_double.checks(ROOT,out)
    elif args.stage=="native":sdpa_double.native_call((ROOT/args.native_input).resolve(),out,args.iterations,args.initial_scale,args.epsilon,args.search_region)
    else:sdpa_double.run(ROOT,out,(ROOT/args.input).resolve(),(ROOT/args.prepared).resolve(),
        (ROOT/args.native_input).resolve(),args.iterations,args.seconds,args.orientation,
        args.condition_intervals,args.initial_scale,args.epsilon,args.equilibrate,(ROOT/args.reference).resolve(),args.search_region)


if __name__=="__main__":main()
