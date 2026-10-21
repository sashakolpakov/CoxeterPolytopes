#!/usr/bin/env python3
"""Enumerate and independently verify a finite rational Cartan packet box."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from coxeter_polytopes.classification.cartan_packets import (
    enumerate_q_simplex_candidates,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Exhaust connected acute Cartan cycle packets over Q for a small rank, "
            "quotient diagonal gauge and facet permutations, and replay every "
            "finite-volume simplex through the exact chamber verifier."
        )
    )
    parser.add_argument("--rank", type=int, default=3)
    parser.add_argument(
        "--edge-products",
        default="1,2,3,4",
        help=(
            "comma-separated nonzero rational length-two products; zero/label 2 "
            "is included automatically (default: 1,2,3,4 = labels 3,4,6,parallel)"
        ),
    )
    parser.add_argument("--max-assignments", type=int, default=2_000_000)
    parser.add_argument(
        "--summary-only", action="store_true", help="omit reconstructed matrices"
    )
    parser.add_argument(
        "--skip-independent-replay",
        action="store_true",
        help="skip the separate exact chamber-certificate replay",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    edge_products = tuple(
        Fraction(part.strip()) for part in args.edge_products.split(",") if part.strip()
    )
    census = enumerate_q_simplex_candidates(
        args.rank,
        edge_products,
        max_assignments=args.max_assignments,
        independent_replay=not args.skip_independent_replay,
    )
    rendered = json.dumps(
        census.to_json(include_matrices=not args.summary_only),
        indent=2,
        sort_keys=True,
    )
    if args.output:
        args.output.write_text(rendered + "\n", encoding="utf-8")
    else:
        print(rendered)


if __name__ == "__main__":
    main()
