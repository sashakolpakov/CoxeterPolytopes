#!/usr/bin/env python3
"""Bounded dense/compact/inline quotient-key equivalence replay.

This is the opt-in verifier for the fixed-width quotient memo key.  It builds
three separately generated engines and completes the manageable source-orbit-0
quotient traversal.  The compact and inline builds both enable the per-state
dense-adjacency round-trip assertion.  Their parsed transcripts must equal the
dense production transcript exactly.

The verifier is intentionally not imported by the ordinary static generator
test as a replay action: compiling and running it should be scheduled only
when host memory permits.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import tempfile
from typing import Any, Mapping

from scripts.generate_h9_14_record31_a1_addition_engine import generate
from scripts.generate_h9_14_record31_a1_compact_key_engine import generate_compact
from scripts.generate_h9_14_record31_a1_inline_key_engine import generate_inline


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "h9_14_record31_a1_inline_key_equivalence.json"
)

REFERENCE_DEPTH_NODES = [1, 56, 420, 8_495, 195_250, 27_551]
REFERENCE_SEEN_COUNTS = [1, 4, 47, 1_393, 51_035, 8_396]


def _compile(source: Path, binary: Path, *, verify_compact: bool) -> None:
    command = ["c++", "-std=c++20", "-O3", "-DNDEBUG"]
    if verify_compact:
        command.append("-DRECORD31_VERIFY_COMPACT_KEY")
    command.extend(
        [
            "-I/opt/homebrew/include/nauty",
            str(source),
            "-L/opt/homebrew/lib",
            "-lnauty",
            "-o",
            str(binary),
        ]
    )
    subprocess.run(command, cwd=ROOT, check=True)


def _run(binary: Path) -> dict[str, Any]:
    completed = subprocess.run(
        [str(binary), "--quotient"],
        input="1\n0 60\n",
        text=True,
        capture_output=True,
        cwd=ROOT,
        check=True,
        timeout=600,
    )
    document = json.loads(completed.stdout)
    if not isinstance(document, dict):
        raise ArithmeticError("bounded Record-31 transcript is not an object")
    return document


def validate_transcripts(
    dense: Mapping[str, Any],
    compact: Mapping[str, Any],
    inline: Mapping[str, Any],
) -> Mapping[str, Any]:
    """Validate exact transcript identity and return the reference orbit row."""

    if dense != compact or dense != inline:
        raise ArithmeticError(
            "dense, compact, and inline bounded quotient transcripts differ"
        )
    if (
        dense.get("kind") != "h9_14_record31_a1_pivoted_addition_scout"
        or dense.get("enumeration_mode") != "nauty_partial_quotient"
        or dense.get("orbit_count") != 1
    ):
        raise ArithmeticError("bounded Record-31 quotient transcript changed")
    rows = dense.get("orbits")
    if not isinstance(rows, list) or len(rows) != 1:
        raise ArithmeticError("bounded Record-31 orbit row changed")
    row = rows[0]
    if not isinstance(row, Mapping):
        raise ArithmeticError("bounded Record-31 orbit row is not an object")
    if not (
        row.get("orbit_index") == 0
        and row.get("positive_A4_masks") == [[2, 3, 4, 5]]
        and row.get("search_node_count") == 231_773
        and row.get("partition_depth_node_counts") == REFERENCE_DEPTH_NODES
        and row.get("quotient_seen_counts") == REFERENCE_SEEN_COUNTS
        and row.get("ridge_valid_target_orbit_count") == 1
    ):
        raise ArithmeticError("bounded source-orbit-0 reference counts changed")
    return row


def verify_equivalence() -> dict[str, Any]:
    dense_source_text = generate()
    compact_source_text = generate_compact()
    inline_source_text = generate_inline()
    if len({dense_source_text, compact_source_text, inline_source_text}) != 3:
        raise ArithmeticError("one quotient-key generator made no distinct source")

    with tempfile.TemporaryDirectory(
        prefix="record31-a1-inline-key-equivalence-"
    ) as raw:
        directory = Path(raw)
        source_texts = {
            "dense": dense_source_text,
            "compact": compact_source_text,
            "inline": inline_source_text,
        }
        documents: dict[str, dict[str, Any]] = {}
        for name, source_text in source_texts.items():
            source = directory / f"{name}.cpp"
            binary = directory / name
            source.write_text(source_text, encoding="utf-8")
            _compile(source, binary, verify_compact=name != "dense")
            documents[name] = _run(binary)

    row = validate_transcripts(
        documents["dense"], documents["compact"], documents["inline"]
    )
    transcript_payload = json.dumps(
        documents["dense"], sort_keys=True, separators=(",", ":")
    ).encode()
    return {
        "schema": "h9_14_record31_a1_inline_key_equivalence_v1",
        "claim": (
            "On the complete bounded source-orbit-0 production quotient run, "
            "dense, compact-string, and zero-padded inline-24-byte quotient "
            "keys give identical exact JSON transcripts; compact and inline "
            "builds also round-trip every canonical graph adjacency bit."
        ),
        "bounded_run": {
            "source_orbit_index": 0,
            "positive_A4_mask": 60,
            "search_node_count": row["search_node_count"],
            "partition_depth_node_counts": row["partition_depth_node_counts"],
            "quotient_seen_counts": row["quotient_seen_counts"],
            "target_orbit_count": row["ridge_valid_target_orbit_count"],
        },
        "dense_compact_inline_transcripts_identical": True,
        "compact_and_inline_dense_round_trip_enabled": True,
        "inline_key_bytes": 24,
        "maximum_compact_code_bytes": 22,
        "dense_source_sha256": sha256(dense_source_text.encode()).hexdigest(),
        "compact_source_sha256": sha256(compact_source_text.encode()).hexdigest(),
        "inline_source_sha256": sha256(inline_source_text.encode()).hexdigest(),
        "common_transcript_sha256": sha256(transcript_payload).hexdigest(),
        "scope": (
            "Exact serialization-equivalence test on a complete bounded run; "
            "not a completed source-orbit-3 census and not a schedule theorem."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    payload = verify_equivalence()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not OUTPUT.is_file() or OUTPUT.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived artifact differs: {OUTPUT}")
    else:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
