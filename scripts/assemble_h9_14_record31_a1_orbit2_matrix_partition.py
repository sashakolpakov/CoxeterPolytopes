#!/usr/bin/env python3
"""Assemble exact target unions from partitioned source-orbit-2 transcripts.

Fragment search-node counts and literal target multiplicities are deliberately
not added.  For each requested root class, the assembler first proves that the
selected depth-two ranges cover every discovered depth-two class, then unions
targets by the complete ``C2 x S9`` canonical family key.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping

from scripts.explore_h9_14_record31_a1_downstream import canonical_family_key


ROOT = Path(__file__).resolve().parents[1]
SOURCE_INDEX = 2
SOURCE_A4 = [2, 3, 11, 12]


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"malformed {name}")
    return value


def _mask(vertices: list[object]) -> int:
    return sum(1 << int(vertex) for vertex in vertices)


def _vertices(mask: int) -> list[int]:
    return [vertex for vertex in range(14) if mask >> vertex & 1]


def assemble(paths: list[Path], expected_roots: set[int]) -> dict[str, object]:
    if not paths:
        raise ValueError("at least one partition transcript is required")
    coverage: dict[int, set[int]] = {}
    discovered: dict[int, int] = {}
    representatives: dict[tuple[int, ...], tuple[int, ...]] = {}
    fragment_rows: list[dict[str, object]] = []
    raw_target_count = 0

    for path in sorted(path.resolve() for path in paths):
        document = _mapping(json.loads(path.read_text()), path.name)
        if (
            document.get("kind") != "h9_14_record31_a1_pivoted_addition_scout"
            or document.get("compact_only") is not True
            or document.get("enumeration_mode")
            != "direct_binary_matrix_partial_quotient"
            or document.get("orbit_count") != 1
        ):
            raise ValueError(f"wrong matrix-fragment header: {path.name}")
        rows = document.get("orbits")
        if not isinstance(rows, list) or len(rows) != 1:
            raise ValueError(f"wrong matrix-fragment row count: {path.name}")
        row = _mapping(rows[0], path.name)
        if (
            row.get("orbit_index") != SOURCE_INDEX
            or row.get("positive_A4_masks") != [SOURCE_A4]
            or row.get("root_class_range_restricted") is not True
            or row.get("root_class_selected_count") != 1
            or row.get("root_class_discovered_count") != 6
            or row.get("depth_two_class_range_restricted") is not True
        ):
            raise ValueError(f"matrix-fragment source/partition changed: {path.name}")
        root_begin = int(row["requested_root_class_begin"])
        root_end = int(row["requested_root_class_end"])
        if root_begin != root_end:
            raise ValueError(f"fragment selects multiple root classes: {path.name}")
        root_class = root_begin
        depth_count = int(row["depth_two_class_discovered_count"])
        if root_class in discovered and discovered[root_class] != depth_count:
            raise ValueError(f"depth-two discovery count disagrees in root {root_class}")
        discovered[root_class] = depth_count
        begin = int(row["requested_depth_two_class_begin"])
        end = int(row["requested_depth_two_class_end"])
        selected_classes = set(range(max(1, begin), min(end, depth_count) + 1))
        if int(row["depth_two_class_selected_count"]) != len(selected_classes):
            raise ValueError(f"selected depth-two count changed: {path.name}")
        overlap = coverage.setdefault(root_class, set()) & selected_classes
        if overlap:
            raise ValueError(f"overlapping depth-two partition in root {root_class}")
        coverage[root_class].update(selected_classes)

        families = row.get("ridge_valid_target_minimal_nonfaces")
        if not isinstance(families, list):
            raise ValueError(f"matrix-fragment target list changed: {path.name}")
        if int(row["ridge_valid_target_orbit_count"]) != len(families):
            raise ValueError(f"matrix-fragment target count changed: {path.name}")
        raw_target_count += len(families)
        for family in families:
            masks = tuple(_mask(face) for face in family)
            key = canonical_family_key(masks)
            literal = tuple(sorted(masks))
            if key not in representatives or literal < representatives[key]:
                representatives[key] = literal
        fragment_rows.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": sha256(path.read_bytes()).hexdigest(),
                "root_class": root_class,
                "depth_two_class_begin": begin,
                "depth_two_class_end": end,
                "depth_two_class_discovered_count": depth_count,
                "search_node_count": int(row["search_node_count"]),
                "pure_120_facet_target_count": int(
                    row["pure_120_facet_target_count"]
                ),
                "literal_target_orbit_count": len(families),
            }
        )

    if set(coverage) != expected_roots:
        raise ValueError(
            f"root coverage {sorted(coverage)} differs from requested "
            f"{sorted(expected_roots)}"
        )
    for root_class in sorted(expected_roots):
        expected = set(range(1, discovered[root_class] + 1))
        if coverage[root_class] != expected:
            missing = sorted(expected - coverage[root_class])
            raise ValueError(
                f"root {root_class} has incomplete depth-two coverage: {missing}"
            )

    ordered = sorted(representatives.items())
    target_families = [
        [_vertices(mask) for mask in representative]
        for _, representative in ordered
    ]
    return {
        "kind": "h9_14_record31_a1_pivoted_addition_scout",
        "compact_only": True,
        "enumeration_mode": "partitioned_direct_binary_matrix_partial_quotient",
        "orbit_count": 1,
        "orbits": [
            {
                "orbit_index": SOURCE_INDEX,
                "positive_A4_masks": [SOURCE_A4],
                "root_class_ids": sorted(expected_roots),
                "fragment_count": len(fragment_rows),
                "fragment_internal_counts_combined": False,
                "raw_fragment_target_orbit_count": raw_target_count,
                "ridge_valid_target_orbit_count": len(ordered),
                "ridge_valid_target_minimal_nonfaces": target_families,
                "canonical_target_family_keys": [list(key) for key, _ in ordered],
                "fragments": fragment_rows,
            }
        ],
        "scope": (
            "Exact compact-incidence canonical target union over complete "
            "partitioned root/depth-two classes; fragment internal counts are "
            "nonadditive and no finite-volume inference is made."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected-roots", required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("transcripts", type=Path, nargs="+")
    arguments = parser.parse_args()
    expected_roots = {
        int(value) for value in arguments.expected_roots.split(",") if value
    }
    document = assemble(arguments.transcripts, expected_roots)
    payload = json.dumps(document, indent=2, sort_keys=True) + "\n"
    if arguments.output is not None:
        arguments.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
