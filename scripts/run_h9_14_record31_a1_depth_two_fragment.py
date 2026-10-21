#!/usr/bin/env python3
"""Run and archive one terminal Record-31 ``a=1`` depth-two fragment.

The generated C++ engine is the mathematical enumerator.  This runner only
places it under the repository RSS/time watchdog, validates the requested
source/root/class headers, and writes an archive after a normal terminal exit.
No archive is written on timeout, RSS termination, malformed JSON, or a
partition-header mismatch.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess
import sys

from scripts.generate_h9_14_record31_a1_depth_two_class_engine import (
    generate_depth_two_class_engine,
)


ROOT = Path(__file__).resolve().parents[1]
WATCHDOG = ROOT / "scripts/run_with_rss_limit.py"
SOURCE_INDEX = 2
SOURCE_MASK = 6156
SOURCE_A4 = [2, 3, 11, 12]
EXPECTED_ROOT_CLASS_COUNT = 6
WATCHDOG_RE = re.compile(
    r"rss-watchdog: disposition=(\S+) scope=(\S+) "
    r"peak_gib=([0-9.]+) elapsed_seconds=([0-9.]+)"
)


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--root-class", type=int, required=True)
    parser.add_argument("--depth-two-begin", type=int, required=True)
    parser.add_argument("--depth-two-end", type=int, required=True)
    parser.add_argument("--expected-depth-two-count", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--memory-gib", type=float, default=0.6)
    parser.add_argument("--timeout", type=float, default=600.0)
    arguments = parser.parse_args()
    if not arguments.binary.is_file():
        raise SystemExit("the generated fragment binary is missing")
    if not (
        arguments.root_class >= 1
        and arguments.depth_two_begin >= 1
        and arguments.depth_two_end >= arguments.depth_two_begin
        and arguments.expected_depth_two_count >= arguments.depth_two_end
    ):
        raise SystemExit("invalid root/depth-two fragment request")
    if arguments.output.exists():
        raise SystemExit("refusing to overwrite an existing fragment archive")

    command = [
        sys.executable,
        str(WATCHDOG),
        "--memory-gib",
        str(arguments.memory_gib),
        "--timeout",
        str(arguments.timeout),
        "--",
        str(arguments.binary),
        "--matrix-quotient",
        "--depth-two-classes",
        str(arguments.root_class),
        str(arguments.depth_two_begin),
        str(arguments.depth_two_end),
    ]
    result = subprocess.run(
        command,
        input=f"1\n{SOURCE_INDEX} {SOURCE_MASK}\n",
        text=True,
        capture_output=True,
        cwd=ROOT,
    )
    if result.returncode != 0:
        sys.stderr.write(result.stderr)
        raise SystemExit(
            f"fragment did not terminate normally (exit {result.returncode}); "
            "no archive written"
        )
    try:
        document = json.loads(result.stdout)
    except json.JSONDecodeError as error:
        sys.stderr.write(result.stderr)
        raise SystemExit(f"malformed terminal fragment JSON: {error}") from error

    rows = document.get("orbits")
    if (
        document.get("kind") != "h9_14_record31_a1_pivoted_addition_scout"
        or document.get("compact_only") is not True
        or document.get("enumeration_mode")
        != "direct_binary_matrix_partial_quotient"
        or document.get("orbit_count") != 1
        or not isinstance(rows, list)
        or len(rows) != 1
    ):
        raise SystemExit("terminal fragment header changed; no archive written")
    row = rows[0]
    expected_selected = arguments.depth_two_end - arguments.depth_two_begin + 1
    if not (
        row.get("orbit_index") == SOURCE_INDEX
        and row.get("positive_A4_masks") == [SOURCE_A4]
        and row.get("root_class_range_restricted") is True
        and row.get("requested_root_class_begin") == arguments.root_class
        and row.get("requested_root_class_end") == arguments.root_class
        and row.get("root_class_discovered_count") == EXPECTED_ROOT_CLASS_COUNT
        and row.get("root_class_selected_count") == 1
        and row.get("depth_two_class_range_restricted") is True
        and row.get("requested_depth_two_class_begin")
        == arguments.depth_two_begin
        and row.get("requested_depth_two_class_end") == arguments.depth_two_end
        and row.get("depth_two_class_discovered_count")
        == arguments.expected_depth_two_count
        and row.get("depth_two_class_selected_count") == expected_selected
    ):
        raise SystemExit("terminal partition header changed; no archive written")

    matches = WATCHDOG_RE.findall(result.stderr)
    if len(matches) != 1 or matches[0][0] != "child_exit":
        sys.stderr.write(result.stderr)
        raise SystemExit("watchdog terminal certificate changed; no archive written")
    disposition, scope, peak, elapsed = matches[0]
    document["terminal_fragment_provenance"] = {
        "generated_cpp_sha256": sha256(
            generate_depth_two_class_engine().encode()
        ).hexdigest(),
        "binary_sha256": _sha256(arguments.binary),
        "watchdog_sha256": _sha256(WATCHDOG),
        "watchdog_disposition": disposition,
        "watchdog_scope": scope,
        "watchdog_peak_RSS_GiB": float(peak),
        "watchdog_elapsed_seconds": float(elapsed),
        "watchdog_memory_ceiling_GiB": arguments.memory_gib,
        "watchdog_timeout_seconds": arguments.timeout,
        "enumerator_stderr_sha256": sha256(result.stderr.encode()).hexdigest(),
    }
    payload = json.dumps(document, indent=2, sort_keys=True) + "\n"
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(payload, encoding="utf-8")
    print(f"wrote={arguments.output}")
    print(f"root_class={arguments.root_class}")
    print(
        "depth_two_range="
        f"{arguments.depth_two_begin}..{arguments.depth_two_end}"
    )
    print(f"search_nodes={row['search_node_count']}")
    print(f"literal_targets={row['ridge_valid_target_orbit_count']}")
    print(f"peak_RSS_GiB={peak}")
    print(f"elapsed_seconds={elapsed}")
    print(f"archive_sha256={sha256(payload.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
