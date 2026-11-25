#!/usr/bin/env python3
"""Serial strictly feasible margin SDP and exact source face extraction."""
import argparse
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get('COXETER_CHECK_OWNER'):raise RuntimeError('Use the serial scientific watchdog')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--output',required=True)
    parser.add_argument('--source',default='results/global_shared_cusp_sdpa_equilibrated01')
    parser.add_argument('--input',default='results/global_shared_cusp_current/necessary_primal.exact.json')
    parser.add_argument('--seconds',type=int,default=550)
    args=parser.parse_args()
    if not 1<=args.seconds<=550:raise ValueError('Leave time for exact checks and face recovery')
    out=(ROOT/args.output).resolve()
    if not out.is_relative_to(ROOT/'results') or out.exists():raise ValueError('Use a fresh results directory')
    sys.path.insert(0,str(ROOT))
    from research.global_formulation import slater_margin
    if args.check:slater_margin.self_check(ROOT,out)
    else:slater_margin.run(ROOT,out,(ROOT/args.source).resolve(),args.seconds,args.input)


if __name__=='__main__':main()
