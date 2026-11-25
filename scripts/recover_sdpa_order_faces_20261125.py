#!/usr/bin/env python3
"""Serial exact scalar-face closure and original-model point recovery."""
import argparse
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get('COXETER_CHECK_OWNER'):raise RuntimeError('Use the serial scientific watchdog')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    parser.add_argument('--source',default='results/global_shared_cusp_sdpa_equilibrated01')
    parser.add_argument('--candidate-run',default='results/global_shared_cusp_sdpa_slater01')
    args=parser.parse_args();out=(ROOT/args.output).resolve()
    if not out.is_relative_to(ROOT/'results') or out.exists():raise ValueError('Use a fresh results directory')
    sys.path.insert(0,str(ROOT))
    from research.global_formulation import order_faces
    order_faces.run(ROOT,out,(ROOT/args.source).resolve(),(ROOT/args.candidate_run).resolve())


if __name__=='__main__':main()
