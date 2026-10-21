#!/usr/bin/env python3
"""Freeze or verify the cubic-seven Leech ``H24`` continuation audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from coxeter_polytopes.search.cubic7_leech_continuation import (  # noqa: E402
    cubic7_leech_h24_continuation_payload,
)


DEFAULT_OUTPUT = Path(
    "data/enumerations/cubic7_leech_h24_h1000_continuation.json"
)


def serialized_payload() -> str:
    return json.dumps(
        cubic7_leech_h24_continuation_payload(), indent=2, sort_keys=True
    ) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--verify",
        action="store_true",
        help="recompute and compare byte-for-byte with the frozen artifact",
    )
    arguments = parser.parse_args()
    rendered = serialized_payload()
    if arguments.verify:
        if not arguments.output.exists():
            raise SystemExit(f"missing artifact: {arguments.output}")
        if arguments.output.read_text() != rendered:
            raise SystemExit(f"artifact drift: {arguments.output}")
        print(f"PASS {arguments.output}")
        return
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(rendered)
    print(arguments.output)


if __name__ == "__main__":
    main()
