#!/usr/bin/env python3
"""Retained native encoding fixtures; production uses the gated checkpoint driver."""
import argparse
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get('COXETER_CHECK_OWNER'):raise RuntimeError('Use the serial scientific watchdog')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=('check',))
    parser.add_argument('--output',required=True)
    args=parser.parse_args();out=(ROOT/args.output).resolve()
    if out.exists() or not out.is_relative_to(ROOT/'results'):raise ValueError('Use fresh results output')
    sys.path.insert(0,str(ROOT))
    from research.global_formulation import sdpa_double
    sdpa_double.checks(ROOT,out)


if __name__=='__main__':main()
