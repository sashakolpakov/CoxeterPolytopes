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
    parser.add_argument("stage", choices=("export", "explore-lp", "explore-sdpa", "explore-psd-cuts", "solve-current", "verify-current", "audit-strengthening"))
    parser.add_argument("--threshold", type=int)
    parser.add_argument("--parity", choices=("odd", "even"))
    parser.add_argument("--output", required=True)
    parser.add_argument("--upper-dimension", type=int, choices=(57,),
        help="Use the retained n>=58 theorem to cover the complementary tail")
    parser.add_argument("--level", type=int, choices=(1, 2), default=1)
    parser.add_argument("--input", default="results/global_tail57_current/necessary_primal.exact.json")
    parser.add_argument("--solver", choices=("prepare", "sdpa", "recover"), default="sdpa")
    parser.add_argument("--orientation", choices=("necessary", "certificate"), default="necessary")
    parser.add_argument("--seconds", type=int, default=600)
    parser.add_argument("--no-facial-reduction", action="store_true")
    parser.add_argument("--resume-from")
    parser.add_argument("--point", default="results/global_current_exact_recovery/original_feasible_point.exact.json")
    args = parser.parse_args()
    sys.path.insert(0, str(ROOT))
    if args.stage == "audit-strengthening":
        from research.global_formulation.strengthening import audit
        audit(ROOT, args.output, args.input, args.point)
        return
    if args.stage == "verify-current":
        from research.global_formulation.feasibility import verify
        verify(ROOT, args.output, args.input, args.point)
        return
    if args.stage == "solve-current":
        from research.global_formulation.feasibility import run
        run(ROOT, args.output, args.input, args.solver, args.orientation, args.seconds,
            facial=not args.no_facial_reduction, resume_from=args.resume_from)
        return
    if args.threshold is None or args.parity is None:
        parser.error("Model construction requires --threshold and --parity")
    from research.global_formulation.explore import run
    run(ROOT, args.output, args.threshold, args.parity, args.stage, args.upper_dimension, level=args.level)


if __name__ == "__main__":
    main()
