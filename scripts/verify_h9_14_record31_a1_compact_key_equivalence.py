#!/usr/bin/env python3
"""Bounded exact dense/compact quotient-key equivalence replay.

The corrected Record-31 a=1 source orbit 0 is a complete, manageable quotient
run.  It traverses every search depth required by its five-addition schedule.
This verifier compiles the dense production generator and the separate compact
generator, enabling the compact key's per-state dense round-trip assertion,
and demands byte-for-byte equality of their parsed JSON transcripts.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import tempfile
from typing import Any

from scripts.generate_h9_14_record31_a1_addition_engine import generate
from scripts.generate_h9_14_record31_a1_compact_key_engine import generate_compact


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "h9_14_record31_a1_compact_key_equivalence.json"
)


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
    if (
        document.get("kind") != "h9_14_record31_a1_pivoted_addition_scout"
        or document.get("enumeration_mode") != "nauty_partial_quotient"
        or document.get("orbit_count") != 1
    ):
        raise ArithmeticError("bounded Record-31 quotient transcript changed")
    return document


def verify_equivalence() -> dict[str, Any]:
    dense_source_text = generate()
    compact_source_text = generate_compact()
    if dense_source_text == compact_source_text:
        raise ArithmeticError("compact generator did not change key serialization")

    with tempfile.TemporaryDirectory(prefix="record31-a1-key-equivalence-") as raw:
        directory = Path(raw)
        dense_source = directory / "dense.cpp"
        compact_source = directory / "compact.cpp"
        dense_binary = directory / "dense"
        compact_binary = directory / "compact"
        dense_source.write_text(dense_source_text, encoding="utf-8")
        compact_source.write_text(compact_source_text, encoding="utf-8")
        _compile(dense_source, dense_binary, verify_compact=False)
        _compile(compact_source, compact_binary, verify_compact=True)
        dense = _run(dense_binary)
        compact = _run(compact_binary)

    if dense != compact:
        raise ArithmeticError("dense and compact bounded quotient transcripts differ")
    row = dense["orbits"][0]
    if not (
        row["orbit_index"] == 0
        and row["search_node_count"] == 231_773
        and row["partition_depth_node_counts"]
        == [1, 56, 420, 8_495, 195_250, 27_551]
        and row["quotient_seen_counts"] == [1, 4, 47, 1_393, 51_035, 8_396]
        and row["ridge_valid_target_orbit_count"] == 1
    ):
        raise ArithmeticError("bounded source-orbit-0 reference counts changed")

    transcript_payload = json.dumps(
        dense, sort_keys=True, separators=(",", ":")
    ).encode()
    return {
        "schema": "h9_14_record31_a1_compact_key_equivalence_v1",
        "claim": (
            "On the complete bounded production source-orbit-0 quotient run, "
            "the dense and compact canonical keys give identical exact JSON "
            "transcripts; the compact build additionally round-trips every "
            "visited canonical graph adjacency bit."
        ),
        "bounded_run": {
            "source_orbit_index": 0,
            "positive_A4_mask": 60,
            "search_node_count": row["search_node_count"],
            "partition_depth_node_counts": row[
                "partition_depth_node_counts"
            ],
            "quotient_seen_counts": row["quotient_seen_counts"],
            "target_orbit_count": row["ridge_valid_target_orbit_count"],
        },
        "dense_and_compact_transcripts_identical": True,
        "compact_per_state_dense_round_trip_enabled": True,
        "dense_source_sha256": sha256(dense_source_text.encode()).hexdigest(),
        "compact_source_sha256": sha256(compact_source_text.encode()).hexdigest(),
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
