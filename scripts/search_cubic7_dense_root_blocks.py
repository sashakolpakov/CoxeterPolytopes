#!/usr/bin/env python3
"""Freeze exact H22/H24 dense root-block Vinberg prefixes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from coxeter_polytopes.search.cubic7_e8_cubes import (
    SUPPORTED_E8_CUBE_SEEDS,
    build_h24_e8_cube_transcript,
)
from coxeter_polytopes.search.cubic7_root_blocks import (
    SUPPORTED_ROOT_BLOCK_SEEDS,
    build_h22_root_block_transcript,
)


def _write(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--family", choices=("h22", "h24", "both"), default="both")
    parser.add_argument("--height-bound", type=int, default=10)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/enumerations/cubic7_dense_root_blocks"),
    )
    args = parser.parse_args()
    if args.family in ("h22", "both"):
        payload = {
            "schema_version": "0.1.0",
            "family": "H22 <-a>+E8+E8+E6",
            "prefixes": [
                build_h22_root_block_transcript(
                    seed, args.height_bound, analyze_geometry=False
                )
                for seed in SUPPORTED_ROOT_BLOCK_SEEDS
            ],
        }
        _write(args.output_dir / f"h22_e8e8e6_h{args.height_bound}.json", payload)
    if args.family in ("h24", "both"):
        payload = {
            "schema_version": "0.1.0",
            "family": "H24 <-a>+E8^3",
            "prefixes": [
                build_h24_e8_cube_transcript(seed, args.height_bound)
                for seed in SUPPORTED_E8_CUBE_SEEDS
            ],
        }
        _write(args.output_dir / f"h24_e8cubed_h{args.height_bound}.json", payload)


if __name__ == "__main__":
    main()
