#!/usr/bin/env python3
"""Retained tiny regression fixtures; production uses the gated checkpoint driver."""
import argparse
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get('COXETER_CHECK_OWNER'):raise RuntimeError('Use the serial scientific watchdog')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true',required=True)
    parser.add_argument('--output',required=True)
    args=parser.parse_args();out=(ROOT/args.output).resolve()
    if out.exists() or not out.is_relative_to(ROOT/'results'):raise ValueError('Use fresh results output')
    sys.path.insert(0,str(ROOT))
    from research.global_formulation import face_phase
    face_phase.self_check(out)


if __name__=='__main__':main()
