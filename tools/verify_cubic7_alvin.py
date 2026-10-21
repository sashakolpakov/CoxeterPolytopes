#!/usr/bin/env python3
"""Command-line entry point for the independent cubic-seven AlVin replay."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from coxeter_polytopes.verify.cubic7_alvin import (  # noqa: E402
    Cubic7AlVinError,
    verify_cubic7_alvin_files,
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Exactly replay a Q(2*cos(2*pi/7)) AlVin chamber from its "
            ".coxiter graph, printed roots, and diagonal form."
        )
    )
    parser.add_argument("--coxiter", required=True, type=Path)
    parser.add_argument(
        "--roots",
        required=True,
        type=Path,
        help="generic/PARI AlVin transcript or JSON root envelope",
    )
    parser.add_argument(
        "--qf",
        required=True,
        help="actual diagonal form, e.g. '[-1,0,0],1,1,1,1'",
    )
    parser.add_argument(
        "--labelg",
        type=Path,
        help="optional explicit path to nauty labelg for exact permutation canonicalization",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="write the JSON report here instead of standard output",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    try:
        report = verify_cubic7_alvin_files(
            arguments.coxiter,
            arguments.roots,
            arguments.qf,
            labelg=arguments.labelg,
        ).to_json()
    except (Cubic7AlVinError, OSError) as exc:
        error = {
            "schema": "cubic7-alvin-independent-replay-0.1.0",
            "passes": False,
            "error": str(exc),
        }
        rendered = json.dumps(error, indent=2, sort_keys=True) + "\n"
        if arguments.output is None:
            sys.stderr.write(rendered)
        else:
            arguments.output.write_text(rendered, encoding="utf-8")
        return 1

    report["passes"] = True
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if arguments.output is None:
        sys.stdout.write(rendered)
    else:
        arguments.output.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
