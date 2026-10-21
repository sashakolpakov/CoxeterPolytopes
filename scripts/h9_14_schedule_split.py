"""Re-solve selected H9/14 shadows one missing-four schedule at a time."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Sequence

from coxeter_polytopes.combinatorics.facet_completion import solve_facet_selection
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "enumerations" / "h9_14_nerve_support.json"
OUTPUT = ROOT / "data" / "enumerations" / "h9_14_schedule_split.json"
TARGET_INDICES = (92, 93)
TIMEOUT_MS = 15_000
EXPECTED_SOURCE_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "2d4e6aa60912a9060c62a07a8d6997dda00fd4eec82d60e850c41fab523f2e1c"
)


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _digest(records: list[dict[str, object]]) -> str:
    critical = [
        {
            "index": record["index"],
            "branches": [
                {
                    "missing_four_count": branch["missing_four_count"],
                    "status": branch["solver"]["status"],
                    "selected_facet_complement_masks": branch["solver"][
                        "selected_facet_complement_masks"
                    ],
                }
                for branch in record["branches"]
            ],
        }
        for record in records
    ]
    return sha256(
        json.dumps(critical, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def generate() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    survivors = [
        record for record in source["records"] if record["survives_second_layer"]
    ]
    records = []
    for index in TARGET_INDICES:
        record = survivors[index]
        branches = []
        for schedule in record["dehn_sommerville_schedule"]:
            result = solve_facet_selection(
                record["graph6"],
                record["ridge_support_core_masks"],
                [schedule],
                timeout_ms=TIMEOUT_MS,
            )
            branches.append(
                {
                    "missing_four_count": schedule["missing_four_count"],
                    "solver": result.to_record(),
                }
            )
            print(
                f"index={index} a={schedule['missing_four_count']} "
                f"status={result.status} seconds={result.elapsed_seconds:.3f}",
                flush=True,
            )
        statuses = {branch["solver"]["status"] for branch in branches}
        aggregate = "unsat" if statuses == {"unsat"} else "unknown"
        records.append(
            {
                "index": index,
                "graph6": record["graph6"],
                "branches": branches,
                "aggregate_backend_status": aggregate,
                "theorem_status": "undecided",
            }
        )
    manifest = {
        "schema_version": "1.0.0",
        "object": "bounded exact schedule-split backend dispositions",
        "source": {
            "path": "data/enumerations/h9_14_nerve_support.json",
            "sha256": _sha256(SOURCE),
        },
        "timeout_ms_per_branch": TIMEOUT_MS,
        "nonclaim": (
            "Backend UNSAT is not promoted to theorem exclusion without an "
            "independent proof certificate."
        ),
        "records": records,
    }
    manifest["certificate_sha256"] = _digest(records)
    return manifest


def verify(manifest: dict[str, object]) -> None:
    if manifest.get("source", {}).get("sha256") != _sha256(SOURCE):
        raise ValueError("schedule-split source mismatch")
    if manifest.get("source", {}).get("sha256") != EXPECTED_SOURCE_SHA256:
        raise ValueError("schedule-split source is not frozen")
    records = manifest.get("records")
    if not isinstance(records, list) or [record.get("index") for record in records] != list(TARGET_INDICES):
        raise ValueError("schedule-split target coverage mismatch")
    if manifest.get("certificate_sha256") != _digest(records):
        raise ValueError("schedule-split certificate digest mismatch")
    if manifest.get("certificate_sha256") != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("schedule-split certificate is not frozen")
    for record in records:
        statuses = {branch["solver"]["status"] for branch in record["branches"]}
        expected = "unsat" if statuses == {"unsat"} else "unknown"
        if record.get("aggregate_backend_status") != expected:
            raise ValueError("schedule-split aggregate mismatch")
        if record.get("theorem_status") != "undecided":
            raise ValueError("backend result was promoted without a certificate")


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    if args.verify_only:
        manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    else:
        manifest = generate()
        OUTPUT.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    verify(manifest)
    print(json.dumps({"output": str(OUTPUT), "records": manifest["records"]}, indent=2))


if __name__ == "__main__":
    main()
