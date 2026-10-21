#!/usr/bin/env python3
"""Certify the lossless Record-31 ``a=1`` depth-two partition.

The calibration source is completed source orbit 0.  A fresh replay compares
its ordinary terminal quotient with (i) four whole-root fragments and (ii)
all 22 individual depth-two classes below the unique target-bearing root.
Only canonical target-family unions are compared: fragment node counts and
internal multiplicities are intentionally nonadditive.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import tempfile
from typing import Mapping, Sequence

from scripts.generate_h9_14_record31_a1_depth_two_class_engine import (
    generate_depth_two_class_engine,
)
from scripts.explore_h9_14_record31_a1_downstream import canonical_family_key


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (
    ROOT
    / "data/enumerations/"
    "h9_14_record31_a1_depth_two_partition_equivalence.json"
)
FULL_SOURCE0 = ROOT / "data/enumerations/h9_14_record31_a1_orbit0_corrected.json"
SOURCE_MASK = 60
EXPECTED_ROOT_CLASS_COUNT = 4
EXPECTED_TARGET_ROOT = 4
EXPECTED_TARGET_ROOT_DEPTH_TWO_CLASS_COUNT = 22


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"malformed {name}")
    return value


def _digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _mask(vertices: Sequence[object]) -> int:
    return sum(1 << int(vertex) for vertex in vertices)


def _target_keys(document: Mapping[str, object]) -> set[tuple[int, ...]]:
    rows = document.get("orbits")
    if not isinstance(rows, list) or len(rows) != 1:
        raise ValueError("partition calibration must contain one source row")
    row = _mapping(rows[0], "partition row")
    families = row.get("ridge_valid_target_minimal_nonfaces")
    if not isinstance(families, list):
        raise ValueError("partition target-family list changed")
    return {
        canonical_family_key(tuple(_mask(face) for face in family))
        for family in families
    }


def _summary(document: Mapping[str, object]) -> dict[str, object]:
    rows = document.get("orbits")
    if not isinstance(rows, list) or len(rows) != 1:
        raise ValueError("partition transcript row count changed")
    row = _mapping(rows[0], "partition transcript row")
    keys = sorted(_target_keys(document))
    return {
        "search_node_count": int(row["search_node_count"]),
        "root_class_discovered_count": int(row["root_class_discovered_count"]),
        "root_class_selected_count": int(row["root_class_selected_count"]),
        "depth_two_class_discovered_count": int(
            row["depth_two_class_discovered_count"]
        ),
        "depth_two_class_selected_count": int(
            row["depth_two_class_selected_count"]
        ),
        "ridge_valid_target_orbit_count": int(
            row["ridge_valid_target_orbit_count"]
        ),
        "canonical_target_family_keys": [list(key) for key in keys],
    }


def _run(
    binary: Path, root_class: int, depth_two_begin: int, depth_two_end: int
) -> Mapping[str, object]:
    result = subprocess.run(
        [
            str(binary),
            "--quotient",
            "--depth-two-classes",
            str(root_class),
            str(depth_two_begin),
            str(depth_two_end),
        ],
        input=f"1\n0 {SOURCE_MASK}\n",
        text=True,
        capture_output=True,
        check=True,
        cwd=ROOT,
    )
    return _mapping(json.loads(result.stdout), "fresh partition transcript")


def fresh() -> dict[str, object]:
    full = _mapping(json.loads(FULL_SOURCE0.read_text()), "source-0 transcript")
    full_keys = _target_keys(full)
    if len(full_keys) != 1:
        raise ValueError("source-0 calibration target changed")

    with tempfile.TemporaryDirectory(prefix="record31-a1-depth-two-") as directory:
        temporary = Path(directory)
        source = temporary / "enumerate.cpp"
        binary = temporary / "enumerate"
        source.write_text(generate_depth_two_class_engine(), encoding="utf-8")
        subprocess.run(
            [
                "c++",
                "-std=c++20",
                "-O3",
                "-DNDEBUG",
                "-I/opt/homebrew/include/nauty",
                str(source),
                "-L/opt/homebrew/lib",
                "-lnauty",
                "-o",
                str(binary),
            ],
            check=True,
            cwd=ROOT,
        )

        root_documents = [
            _run(binary, root_class, 1, 1_000_000)
            for root_class in range(1, EXPECTED_ROOT_CLASS_COUNT + 1)
        ]
        root_summaries = [_summary(document) for document in root_documents]
        root_union = set().union(*map(_target_keys, root_documents))
        if root_union != full_keys:
            raise ValueError("whole-root target union differs from source 0")

        target_root_summary = root_summaries[EXPECTED_TARGET_ROOT - 1]
        if (
            target_root_summary["depth_two_class_discovered_count"]
            != EXPECTED_TARGET_ROOT_DEPTH_TWO_CLASS_COUNT
        ):
            raise ValueError("target-root depth-two class count changed")
        depth_two_documents = [
            _run(binary, EXPECTED_TARGET_ROOT, class_id, class_id)
            for class_id in range(1, EXPECTED_TARGET_ROOT_DEPTH_TWO_CLASS_COUNT + 1)
        ]
        depth_two_summaries = [
            _summary(document) for document in depth_two_documents
        ]
        depth_two_union = set().union(*map(_target_keys, depth_two_documents))
        if depth_two_union != full_keys:
            raise ValueError("depth-two target union differs from source 0")

    document: dict[str, object] = {
        "kind": "h9_14_record31_a1_depth_two_partition_equivalence",
        "compact_only": True,
        "calibration_source_orbit_index": 0,
        "calibration_positive_A4_mask": SOURCE_MASK,
        "full_source_search_node_count": int(full["orbits"][0]["search_node_count"]),
        "full_source_target_family_keys": [list(key) for key in sorted(full_keys)],
        "root_fragment_summaries": root_summaries,
        "root_fragment_target_union_equals_full": True,
        "target_bearing_root_class": EXPECTED_TARGET_ROOT,
        "target_root_individual_depth_two_summaries": depth_two_summaries,
        "individual_depth_two_target_union_equals_full": True,
        "fragment_internal_counts_are_nonadditive": True,
        "full_source_sha256": sha256(FULL_SOURCE0.read_bytes()).hexdigest(),
        "generated_cpp_sha256": sha256(
            generate_depth_two_class_engine().encode()
        ).hexdigest(),
        "scope": (
            "Exact bounded equivalence calibration of the partition only; "
            "no new source exclusion and no finite-volume inference."
        ),
    }
    document["certificate_sha256"] = _digest(document)
    verify(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record31_a1_depth_two_partition_equivalence":
        raise ValueError("wrong depth-two equivalence certificate kind")
    if document.get("compact_only") is not True:
        raise ValueError("partition calibration scope changed")
    if document.get("calibration_source_orbit_index") != 0:
        raise ValueError("partition calibration source changed")
    if document.get("calibration_positive_A4_mask") != SOURCE_MASK:
        raise ValueError("partition calibration mask changed")
    if document.get("full_source_sha256") != sha256(FULL_SOURCE0.read_bytes()).hexdigest():
        raise ValueError("completed source-0 transcript changed")
    if document.get("generated_cpp_sha256") != sha256(
        generate_depth_two_class_engine().encode()
    ).hexdigest():
        raise ValueError("depth-two generated source changed")
    if document.get("root_fragment_target_union_equals_full") is not True:
        raise ValueError("root-fragment target equality is not certified")
    if document.get("individual_depth_two_target_union_equals_full") is not True:
        raise ValueError("depth-two target equality is not certified")
    if document.get("fragment_internal_counts_are_nonadditive") is not True:
        raise ValueError("partition count semantics changed")

    root_rows = document.get("root_fragment_summaries")
    depth_rows = document.get("target_root_individual_depth_two_summaries")
    if not isinstance(root_rows, list) or len(root_rows) != EXPECTED_ROOT_CLASS_COUNT:
        raise ValueError("root-fragment coverage changed")
    if (
        not isinstance(depth_rows, list)
        or len(depth_rows) != EXPECTED_TARGET_ROOT_DEPTH_TWO_CLASS_COUNT
    ):
        raise ValueError("depth-two fragment coverage changed")
    for row in root_rows:
        summary = _mapping(row, "root summary")
        if summary.get("root_class_discovered_count") != EXPECTED_ROOT_CLASS_COUNT:
            raise ValueError("root-class discovery count changed")
        if summary.get("root_class_selected_count") != 1:
            raise ValueError("root fragment selected more than one class")
    for row in depth_rows:
        summary = _mapping(row, "depth-two summary")
        if (
            summary.get("depth_two_class_discovered_count")
            != EXPECTED_TARGET_ROOT_DEPTH_TWO_CLASS_COUNT
            or summary.get("depth_two_class_selected_count") != 1
        ):
            raise ValueError("individual depth-two fragment coverage changed")

    full_keys = document.get("full_source_target_family_keys")
    root_keys = {
        tuple(key)
        for row in root_rows
        for key in _mapping(row, "root summary")["canonical_target_family_keys"]
    }
    depth_keys = {
        tuple(key)
        for row in depth_rows
        for key in _mapping(row, "depth-two summary")["canonical_target_family_keys"]
    }
    if not isinstance(full_keys, list) or root_keys != {tuple(k) for k in full_keys}:
        raise ValueError("frozen root target union changed")
    if depth_keys != {tuple(k) for k in full_keys}:
        raise ValueError("frozen depth-two target union changed")
    if document.get("certificate_sha256") != _digest(document):
        raise ValueError("depth-two certificate digest changed")
    return {
        "root_fragment_count": len(root_rows),
        "individual_depth_two_fragment_count": len(depth_rows),
        "canonical_target_union_count": len(root_keys),
        "certificate_sha256": document["certificate_sha256"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh", action="store_true")
    arguments = parser.parse_args()
    if arguments.fresh:
        document = fresh()
        ARCHIVE.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    else:
        document = _mapping(json.loads(ARCHIVE.read_text()), "archive")
    print(json.dumps(verify(document), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
