#!/usr/bin/env python3
"""Freeze the exact first-shell audit for ``<-2*theta> + Lambda_24``."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from coxeter_polytopes.search.cubic7_leech import (  # noqa: E402
    first_shell_audit,
    golay_weight_distribution,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/enumerations/cubic7_leech_h24_first_shell.json"),
    )
    args = parser.parse_args()
    payload = {
        "schema": "coxeter.cubic7_leech_h24_first_shell.v1",
        "field": "Q(theta), theta^3+theta^2-2theta-1=0",
        "lattice": "<-2*theta> + Lambda_24",
        "hyperbolic_dimension": 24,
        "golay_weight_distribution": golay_weight_distribution(),
        "audit": first_shell_audit(),
        "scope": (
            "Complete mirror-norm proof, complete height-zero proof, and "
            "complete first positive control-shell enumeration only; no "
            "reflectivity, finite-volume, or compactness claim."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
