#!/usr/bin/env python3
"""Serial entry point for the global necessary-data LP/SDP, not diagram pricing."""
import argparse
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    if not os.environ.get("COXETER_CHECK_OWNER"):
        raise RuntimeError("Use the repository's serial bounded scientific runner")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("export", "verify-current", "audit-strengthening", "check-shared-program"))
    parser.add_argument("--threshold", type=int)
    parser.add_argument("--parity", choices=("odd", "even"))
    parser.add_argument("--output", required=True)
    parser.add_argument("--upper-dimension", type=int, choices=(57,),
        help="Use the retained n>=58 theorem to cover the complementary tail")
    parser.add_argument("--level", type=int, choices=(1, 2, 3, 4, 5, 6, 7), default=7)
    parser.add_argument("--input")
    parser.add_argument("--point")
    args = parser.parse_args()
    if args.stage in ('verify-current','audit-strengthening') and not (args.input and args.point):
        parser.error('Specify the original --input and --point explicitly')
    sys.path.insert(0,str(ROOT))
    if args.stage == "check-shared-program":
        from research.global_formulation.program_checks import run
        run(ROOT,args.output)
        return
    if args.stage == "audit-strengthening":
        from research.global_formulation.strengthening import audit
        audit(ROOT, args.output, args.input, args.point)
        return
    if args.stage == "verify-current":
        from research.global_formulation.feasibility import verify
        verify(ROOT, args.output, args.input, args.point)
        return
    if args.threshold is None or args.parity is None:
        parser.error("Model construction requires --threshold and --parity")
    from research.global_formulation.explore import run
    run(ROOT, args.output, args.threshold, args.parity, args.stage, args.upper_dimension, level=args.level)


if __name__ == "__main__":
    main()
