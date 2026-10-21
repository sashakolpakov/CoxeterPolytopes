#!/usr/bin/env python3
"""Assemble streamed Record-1 quotient class logs without overclaiming.

Completed seed classes are durable proof boundaries.  This assembler accepts
one or more JSONL logs, independently reconstructs the 21-class inventory,
checks every complete class terminal and every streamed target, and emits a
whole-census terminal only when all 21 classes are present and complete.

Cutoff logs and interrupted logs may be supplied: they are retained in the
diagnostic counts but can never satisfy a class boundary.  This script does
not resume within a class; the durable restart boundary is a completed class.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from typing import Any, Iterable, Sequence

from scripts.independent_h9_14_record1_a1_orbit0_quotient_audit import (
    EDGE_MASKS,
    FOUR,
    FULL,
    N,
    SOURCE_BLOCKS,
    independent_candidates,
    local_interval,
    mask,
    signature,
)


PROTOCOL = "h9_14_record1_a1_orbit0_incidence_quotient_jsonl_v1"
EXPECTED_CLASSES = 21
EXPECTED_ADDITIONS = 19
EXPECTED_FACETS = 158
EXPECTED_RIDGES = 711


def hits_all(value: int, clutter: Sequence[int]) -> bool:
    return all(value & row for row in clutter)


def minimal_transversals(clutter: Sequence[int]) -> tuple[int, ...]:
    answer = []
    for value in range(1, FULL + 1):
        if not hits_all(value, clutter):
            continue
        if all(
            not hits_all(value ^ (1 << vertex), clutter)
            for vertex in range(N)
            if value >> vertex & 1
        ):
            answer.append(value)
    return tuple(answer)


def expected_inventory(candidates: Sequence[int]) -> tuple[dict[str, int], ...]:
    groups: dict[tuple[int, ...], list[int]] = defaultdict(list)
    for value in candidates:
        groups[signature(value, SOURCE_BLOCKS)].append(value)
    rows = sorted((min(values), len(values)) for values in groups.values())
    return tuple(
        {
            "class_index": index,
            "seed_mask": seed,
            "candidate_orbit_size": size,
            "seed_stabilizer_order": 60480 // size,
        }
        for index, (seed, size) in enumerate(rows)
    )


def independent_target_audit(
    target: dict[str, Any], candidates: frozenset[int], seed: int
) -> dict[str, int]:
    if target.get("source_A4_mask") != FOUR:
        raise ArithmeticError("streamed target source A4 changed")
    additions = tuple(int(value) for value in target["additional_A5_masks"])
    if len(additions) != EXPECTED_ADDITIONS or len(set(additions)) != len(additions):
        raise ArithmeticError("streamed target addition family changed shape")
    if seed not in additions or any(value not in candidates for value in additions):
        raise ArithmeticError("streamed target escaped its seed/candidate ground")
    if target.get("scheduled_facets") != EXPECTED_FACETS:
        raise ArithmeticError("streamed target facet schedule changed")
    if target.get("active_ridges") != EXPECTED_RIDGES:
        raise ArithmeticError("streamed target ridge schedule changed")
    canonical_key = str(target.get("canonical_full_key", ""))
    if (
        len(canonical_key) != 84
        or not canonical_key.startswith("17")
        or any(character not in "0123456789abcdef" for character in canonical_key)
    ):
        raise ArithmeticError("streamed target canonical-key encoding changed")

    for size in (1, 2, 3):
        for values in combinations(range(N), size):
            subset = mask(values)
            if any(subset & edge == edge for edge in EDGE_MASKS):
                continue
            lower, upper = local_interval(subset)
            incidence = sum(value & subset == subset for value in additions)
            if not lower <= incidence <= upper:
                raise ArithmeticError("streamed target violates a local interval")

    shadow = tuple(
        value
        for values in combinations(range(N), 5)
        if (value := mask(values)) & FOUR == FOUR
        and all(value & edge != edge for edge in EDGE_MASKS)
    )
    clutter = (*EDGE_MASKS, FOUR, *shadow, *additions)
    blocker = minimal_transversals(clutter)
    if len(blocker) != EXPECTED_FACETS or any(value.bit_count() != 5 for value in blocker):
        raise ArithmeticError("streamed target is not a pure 158-facet blocker")
    ridge_degrees: Counter[int] = Counter()
    for facet in blocker:
        outside = FULL ^ facet
        for vertex in range(N):
            if outside >> vertex & 1:
                ridge_degrees[facet | (1 << vertex)] += 1
    if len(ridge_degrees) != EXPECTED_RIDGES or set(ridge_degrees.values()) != {2}:
        raise ArithmeticError("streamed target fails the exact ridge law")
    return {"facet_count": len(blocker), "active_ridge_count": len(ridge_degrees)}


def parse_log(path: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    events = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(events) < 3 or events[0].get("event") != "run_header":
        raise ArithmeticError(f"{path}: missing run header")
    if events[0].get("protocol") != PROTOCOL:
        raise ArithmeticError(f"{path}: protocol changed")
    inventories = [row for row in events if row.get("event") == "class_inventory"]
    if len(inventories) != 1:
        raise ArithmeticError(f"{path}: expected one class inventory")

    segments: list[dict[str, Any]] = []
    active: dict[str, Any] | None = None
    for event in events[2:]:
        kind = event.get("event")
        if kind == "class_started":
            if active is not None:
                raise ArithmeticError(f"{path}: nested class segment")
            active = {"started": event, "targets": [], "checkpoints": []}
        elif kind == "target":
            if active is None or event.get("class_index") != active["started"].get("class_index"):
                raise ArithmeticError(f"{path}: target outside its class segment")
            active["targets"].append(event)
        elif kind == "checkpoint":
            if active is None or event.get("class_index") != active["started"].get("class_index"):
                raise ArithmeticError(f"{path}: checkpoint outside its class segment")
            active["checkpoints"].append(event)
        elif kind == "class_terminal":
            if active is None or event.get("class_index") != active["started"].get("class_index"):
                raise ArithmeticError(f"{path}: terminal outside its class segment")
            active["terminal"] = event
            segments.append(active)
            active = None
        elif kind == "run_terminal":
            if active is not None:
                raise ArithmeticError(f"{path}: run ended inside a class")
        else:
            raise ArithmeticError(f"{path}: unknown event {kind!r}")
    # An active segment means interruption.  It is deliberately ignored as a
    # proof boundary, but retained in the caller's incomplete-segment count.
    if active is not None:
        active["interrupted"] = True
        segments.append(active)
    return inventories, segments


def assemble(paths: Sequence[Path]) -> dict[str, Any]:
    if not paths:
        raise ValueError("at least one JSONL log is required")
    candidates_tuple = independent_candidates()
    candidates = frozenset(candidates_tuple)
    expected = expected_inventory(candidates_tuple)
    expected_by_class = {row["class_index"]: row for row in expected}

    complete_segments: dict[int, list[dict[str, Any]]] = defaultdict(list)
    incomplete_segment_count = 0
    source_digests = {}
    for path in paths:
        source_digests[str(path)] = sha256(path.read_bytes()).hexdigest()
        inventories, segments = parse_log(path)
        production_rows = inventories[0].get("classes")
        reduced = tuple(
            {
                "class_index": int(row["class_index"]),
                "seed_mask": int(row["seed_mask"]),
                "candidate_orbit_size": int(row["candidate_orbit_size"]),
                "seed_stabilizer_order": int(row["seed_stabilizer_order"]),
            }
            for row in production_rows
        )
        if reduced != expected:
            raise ArithmeticError(f"{path}: class inventory differs from direct replay")
        for segment in segments:
            terminal = segment.get("terminal")
            if (
                terminal is None
                or terminal.get("terminal_status") != "complete_exhaustive"
                or terminal.get("class_proof_complete") is not True
            ):
                incomplete_segment_count += 1
                continue
            class_index = int(terminal["class_index"])
            if class_index not in expected_by_class:
                raise ArithmeticError(f"{path}: completed unknown class")
            if int(segment["started"]["seed_mask"]) != expected_by_class[class_index]["seed_mask"]:
                raise ArithmeticError(f"{path}: class seed changed")
            depth_nodes = tuple(int(value) for value in terminal["depth_nodes"])
            if len(depth_nodes) != 20 or sum(depth_nodes) != int(terminal["nodes"]):
                raise ArithmeticError(f"{path}: class node accounting changed")
            keys = []
            for target in segment["targets"]:
                independent_target_audit(
                    target,
                    candidates,
                    expected_by_class[class_index]["seed_mask"],
                )
                keys.append(str(target["canonical_full_key"]))
            if len(keys) != len(set(keys)):
                raise ArithmeticError(f"{path}: duplicate target stream within class")
            if len(keys) != int(terminal["full_target_orbit_count_within_class"]):
                raise ArithmeticError(f"{path}: target stream/terminal count mismatch")
            segment["target_keys"] = tuple(sorted(keys))
            complete_segments[class_index].append(segment)

    canonical_complete: dict[int, dict[str, Any]] = {}
    for class_index, segments in complete_segments.items():
        reference = segments[0]
        reference_payload = {
            "terminal": reference["terminal"],
            "target_keys": reference["target_keys"],
        }
        for duplicate in segments[1:]:
            payload = {
                "terminal": duplicate["terminal"],
                "target_keys": duplicate["target_keys"],
            }
            if payload != reference_payload:
                raise ArithmeticError(
                    f"independent complete logs disagree on class {class_index}"
                )
        canonical_complete[class_index] = reference

    missing = sorted(set(range(EXPECTED_CLASSES)) - set(canonical_complete))
    global_keys = sorted(
        {
            key
            for segment in canonical_complete.values()
            for key in segment["target_keys"]
        }
    )
    whole = not missing
    class_summaries = {
        str(index): {
            "seed_mask": expected_by_class[index]["seed_mask"],
            "nodes": int(segment["terminal"]["nodes"]),
            "quotient_seen_counts": segment["terminal"]["quotient_seen_counts"],
            "ridge_valid_seed_stabilizer_target_count": int(
                segment["terminal"]["ridge_valid_seed_stabilizer_target_count"]
            ),
            "full_target_orbit_count_within_class": len(segment["target_keys"]),
        }
        for index, segment in sorted(canonical_complete.items())
    }
    proof_payload = json.dumps(
        {
            "class_summaries": class_summaries,
            "global_target_keys": global_keys,
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return {
        "schema": "h9_14_record1_a1_orbit0_quotient_assembly_v1",
        "terminal_status": (
            "complete_exhaustive" if whole else "incomplete_missing_classes"
        ),
        "enumeration_complete": whole,
        "whole_census_complete": whole,
        "complete_class_count": len(canonical_complete),
        "required_class_count": EXPECTED_CLASSES,
        "missing_class_indices": missing,
        "incomplete_or_cutoff_segment_count": incomplete_segment_count,
        "global_target_orbit_count": len(global_keys),
        "global_target_canonical_keys": global_keys,
        "class_summaries": class_summaries,
        "proof_payload_sha256": sha256(proof_payload).hexdigest(),
        "source_log_sha256": source_digests,
        "scope": (
            "Whole Record-1 a=1 source-orbit-0 necessary compact incidence "
            "quotient only; no label, polytopality, or finite-volume inference."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logs", nargs="+", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--require-complete", action="store_true")
    arguments = parser.parse_args()
    payload = assemble(arguments.logs)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    if arguments.require_complete and not payload["whole_census_complete"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
