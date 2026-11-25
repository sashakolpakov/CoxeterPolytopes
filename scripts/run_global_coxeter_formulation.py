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
    parser.add_argument("stage", choices=("export", "explore-lp", "explore-sdpa", "explore-psd-cuts", "solve-current"))
    parser.add_argument("--threshold", required=True, type=int)
    parser.add_argument("--parity", required=True, choices=("odd", "even"))
    parser.add_argument("--output", required=True)
    parser.add_argument("--upper-dimension", type=int, choices=(57,),
        help="Use the retained n>=58 theorem to cover the complementary tail")
    parser.add_argument("--input", default="results/global_tail57_current/necessary_primal.exact.json")
    parser.add_argument("--solver", choices=("prepare", "sdpa", "clarabel", "scs", "recover"), default="sdpa")
    parser.add_argument("--orientation", choices=("necessary", "certificate"), default="necessary")
    parser.add_argument("--seconds", type=int, default=600)
    parser.add_argument("--no-facial-reduction", action="store_true")
    parser.add_argument("--resume-from")
    args = parser.parse_args()
    sys.path.insert(0, str(ROOT))
    if args.stage == "solve-current":
        from research.global_formulation.feasibility import run
        run(ROOT, args.output, args.input, args.solver, args.orientation, args.seconds,
            facial=not args.no_facial_reduction, resume_from=args.resume_from)
        return
    from research.global_formulation.explore import run
    run(ROOT, args.output, args.threshold, args.parity, args.stage, args.upper_dimension)


if __name__ == "__main__":
    main()
