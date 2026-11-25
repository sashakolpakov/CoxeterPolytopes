#!/usr/bin/env python3
"""Serial warm-started native SDPA with original-model exact checks."""
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
    parser.add_argument('--reference',default='results/global_shared_cusp_sdpa_certificate01/original_numerical_candidate.json')
    parser.add_argument('--resume')
    parser.add_argument('--iterations',type=int,default=60)
    parser.add_argument('--seconds',type=int,default=600)
    parser.add_argument('--primal-floor',type=float,default=1e-3)
    parser.add_argument('--complementarity',type=float,default=1e-3)
    parser.add_argument('--arithmetic',choices=('double','gmp'),default='double')
    parser.add_argument('--precision',type=int,default=128)
    args=parser.parse_args()
    if not 1<=args.iterations<=200 or not 1<=args.seconds<=750:raise ValueError('Use bounded native batches')
    if not 1e-12<=args.primal_floor<=1 or not 1e-15<=args.complementarity<=100:raise ValueError('Unsupported warm-start parameters')
    if args.precision not in (128,192,256):raise ValueError('Use a declared precision')
    out=(ROOT/args.output).resolve()
    if not out.is_relative_to(ROOT/'results') or out.exists():raise ValueError('Use a fresh results directory')
    sys.path.insert(0,str(ROOT))
    from research.global_formulation import sdpa_cli
    precision=args.precision if args.arithmetic=='gmp' else None
    if args.check:sdpa_cli.self_check(ROOT,out,precision)
    else:sdpa_cli.run(ROOT,out,(ROOT/args.source).resolve(),(ROOT/args.reference).resolve(),
        args.iterations,args.seconds,(ROOT/args.resume).resolve() if args.resume else None,args.primal_floor,args.complementarity,precision)


if __name__=='__main__':main()
