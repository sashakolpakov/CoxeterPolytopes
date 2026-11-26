#!/usr/bin/env python3
"""One restartable SDPA batch on the identical frozen sparse input."""
import argparse
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get("COXETER_CHECK_OWNER"): raise RuntimeError("Use the serial scientific watchdog")
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",required=True)
    parser.add_argument("--native-input",required=True)
    parser.add_argument("--initial")
    parser.add_argument("--preflight",required=True)
    parser.add_argument("--resume")
    parser.add_argument("--iterations",type=int,default=16)
    parser.add_argument("--seconds",type=int,default=800)
    args=parser.parse_args()
    if not args.resume and not args.initial:parser.error('Specify the checked initialization or a complete saved restart')
    if not 1<=args.iterations<=100 or not 1<=args.seconds<=800: raise ValueError("Leave time for checkpoint validation")
    out=(ROOT/args.output).resolve();input_path=(ROOT/args.native_input).resolve()
    resume=(ROOT/args.resume).resolve() if args.resume else None
    initial=resume/"native.result" if resume is not None else (ROOT/args.initial).resolve()
    if out.exists() or any(not p.is_relative_to(ROOT/"results") for p in (out,input_path,initial)):
        raise ValueError("Use retained results inputs and fresh output")
    sys.path.insert(0,str(ROOT))
    from research.global_formulation.checkpoint_solver import run
    run(ROOT,out,input_path,initial,args.iterations,args.seconds,resume,(ROOT/args.preflight).resolve())


if __name__=="__main__": main()
