#!/usr/bin/env python3
"""Freeze the bounded Record-1 ``a=1`` target-orbit compact no-go.

The earlier exploratory command used ``--max-models 10`` and left no archived
UNSAT transcript.  Its report's phrase "exactly nine targets" is impossible:
the first frozen target has a 60,480-element labelled orbit under the fixed
source stabilizer.  This certifier proves and freezes that orbit theorem.  It
does *not* claim a complete quotient incidence census for the source.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

from scripts.explore_h9_14_open_positive_a4 import mask, vertices
from scripts.explore_h9_14_open_positive_completion import (
    audit_compact_labels,
    build_problem,
)
from scripts.explore_h9_14_record14_fixed_orbit_z3 import minimal_transversals
from scripts.independent_h9_14_record1_a1_orbit0_target_orbit_audit import (
    replay as independent_replay,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_record1_a1_orbit0_target_orbit_no_go.json"
PRODUCTION_SOURCE = ROOT / "scripts/explore_h9_14_open_positive_completion.py"
INDEPENDENT_ORBIT_SOURCE = (
    ROOT / "scripts/independent_h9_14_record1_a1_orbit0_target_orbit_audit.py"
)
INDEPENDENT_LABEL_SOURCE = (
    ROOT / "scripts/independent_h9_14_record1_a1_target0_label_audit.py"
)
LABEL_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py"
EXPECTED_CERTIFICATE_DIGEST = (
    "199a9cbfa94d32cf96371d8dec1dd1411ba6d2dbd1afc5d0d7e4fc2e93fb4eae"
)
N = 14

REPRESENTATIVE_TARGET = (
    (1, 3, 4, 7, 8),
    (2, 3, 4, 7, 8),
    (3, 4, 5, 7, 8),
    (1, 3, 4, 8, 9),
    (0, 4, 5, 9, 10),
    (1, 4, 5, 9, 10),
    (0, 5, 6, 9, 10),
    (1, 2, 7, 8, 11),
    (1, 2, 8, 9, 11),
    (1, 4, 8, 9, 11),
    (1, 4, 9, 10, 11),
    (1, 6, 9, 10, 11),
    (0, 3, 5, 7, 12),
    (0, 2, 6, 7, 12),
    (2, 3, 7, 8, 12),
    (0, 5, 6, 10, 12),
    (2, 6, 7, 11, 12),
    (2, 6, 10, 11, 12),
    (6, 9, 10, 11, 12),
)


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _family_bytes(family: Sequence[int]) -> bytes:
    return (",".join(map(str, sorted(family))) + "\n").encode("ascii")


def _move(value: int, action: Sequence[int]) -> int:
    return mask(action[vertex] for vertex in vertices(value))


def _actions() -> Iterator[tuple[int, ...]]:
    # The graph/source signatures give singleton blocks {13},{0} and free
    # blocks {1,2}, {3,4,5}, {6,...,12}.
    for leaves in permutations((1, 2)):
        for inside in permutations((3, 4, 5)):
            for outside in permutations((6, 7, 8, 9, 10, 11, 12)):
                action = list(range(N))
                action[1], action[2] = leaves
                action[3:6] = inside
                action[6:13] = outside
                yield tuple(action)


def _swap(left: int, right: int) -> tuple[int, ...]:
    action = list(range(N))
    action[left], action[right] = action[right], action[left]
    return tuple(action)


def fresh_record() -> dict[str, object]:
    problem = build_problem(1, 1, 0)
    target = tuple(sorted(mask(face) for face in REPRESENTATIVE_TARGET))
    if problem["fours"] != (mask((0, 3, 4, 5)),):
        raise AssertionError("Record-1 source representative changed")
    if (
        len(problem["shadow"]),
        problem["addition_count"],
        problem["scheduled_facets"],
        len(problem["candidates"]),
        len(problem["small"]),
        len(problem["fives"]),
    ) != (9, 19, 158, 1278, 259, 634):
        raise AssertionError("Record-1 fixed-source problem shape changed")
    if len(target) != 19 or not all(face in problem["candidates"] for face in target):
        raise AssertionError("representative target left the candidate ground set")
    degrees = tuple(
        sum(bool(face & (1 << vertex)) for face in target)
        for vertex in range(N)
    )
    if any(
        not low <= degree <= high
        for degree, low, high in zip(
            degrees, problem["lower"], problem["upper"], strict=True
        )
    ):
        raise AssertionError("representative violates a local addition interval")

    blockers = minimal_transversals((*problem["base"], *target))
    blocker_histogram = Counter(value.bit_count() for value in blockers)
    if blocker_histogram != {5: 158}:
        raise AssertionError("representative is not a literal pure target")
    ridge_histogram: Counter[int] = Counter()
    for six in combinations(range(N), 6):
        six_mask = mask(six)
        ridge_histogram[
            sum(facet & six_mask == facet for facet in blockers)
        ] += 1
    if any(degree not in (0, 2) for degree in ridge_histogram):
        raise AssertionError("representative violates a ridge equation")

    # It suffices to check invariance on adjacent-transposition generators.
    generators = (
        _swap(1, 2),
        _swap(3, 4),
        _swap(4, 5),
        *(_swap(vertex, vertex + 1) for vertex in range(6, 12)),
    )
    candidate_set = frozenset(problem["candidates"])
    edge_set = frozenset(problem["edge_rows"])
    shadow_set = frozenset(problem["shadow"])
    for action in generators:
        if frozenset(_move(row, action) for row in problem["edge_rows"]) != edge_set:
            raise AssertionError("source generator moved the missing graph")
        if tuple(_move(row, action) for row in problem["fours"]) != problem["fours"]:
            raise AssertionError("source generator moved the A4 source")
        if frozenset(_move(row, action) for row in problem["shadow"]) != shadow_set:
            raise AssertionError("source generator moved the forced shadow")
        if frozenset(_move(row, action) for row in problem["candidates"]) != candidate_set:
            raise AssertionError("source generator moved the candidate ground set")

    action_count = 0
    fixed_target_actions = 0
    stream = sha256()
    for action in _actions():
        action_count += 1
        moved = tuple(sorted(_move(face, action) for face in target))
        fixed_target_actions += int(moved == target)
        stream.update(_family_bytes(moved))
    if action_count != 60_480 or fixed_target_actions != 1:
        raise AssertionError("production orbit-stabilizer count changed")

    label = audit_compact_labels(problem["edge_rows"], problem["fours"], target)
    label_counts = [stage["passing_state_count"] for stage in label["stages"]]
    if label_counts != [108, 84, 8, 0] or not label["unrestricted_compact_label_no_go"]:
        raise AssertionError("production compact-label obstruction changed")

    independent = independent_replay()
    comparison = {
        "source_stabilizer_order": action_count,
        "representative_target_stabilizer_order": fixed_target_actions,
        "labelled_target_orbit_size": action_count // fixed_target_actions,
        "action_image_stream_sha256": stream.hexdigest(),
        "representative_additional_A5_sha256": sha256(
            _family_bytes(target)
        ).hexdigest(),
        "representative_minimal_transversal_size_histogram": {
            str(size): count for size, count in sorted(blocker_histogram.items())
        },
        "representative_ridge_degree_histogram": {
            str(degree): count for degree, count in sorted(ridge_histogram.items())
        },
        "representative_compact_label_stage_counts": label_counts,
    }
    if any(independent.get(key) != value for key, value in comparison.items()):
        raise AssertionError("production/independent orbit replay mismatch")

    document: dict[str, object] = {
        "kind": "h9_14_record1_a1_orbit0_target_orbit_compact_no_go",
        "record_index": 1,
        "graph6": problem["graph6"],
        "schedule_a": 1,
        "source_positive_A4_orbit_index": 0,
        "representative_A4s": [list(vertices(row)) for row in problem["fours"]],
        "source_stabilizer": {
            "signature_partition": [[13], [0], [1, 2], [3, 4, 5], list(range(6, 13))],
            "structure": "S2 x S3 x S7",
            "order": action_count,
            "generator_count": len(generators),
        },
        "incidence_problem": {
            "forced_A5_shadow_count": len(problem["shadow"]),
            "additional_A5_count": problem["addition_count"],
            "scheduled_facet_count": problem["scheduled_facets"],
            "candidate_count": len(problem["candidates"]),
            "base_small_transversal_count": len(problem["small"]),
            "base_order_five_transversal_count": len(problem["fives"]),
        },
        "representative_target": {
            "additional_A5s": [list(face) for face in REPRESENTATIVE_TARGET],
            "additional_A5_sha256": comparison[
                "representative_additional_A5_sha256"
            ],
            "minimal_transversal_size_histogram": comparison[
                "representative_minimal_transversal_size_histogram"
            ],
            "ridge_degree_histogram": comparison[
                "representative_ridge_degree_histogram"
            ],
            "compact_label_stage_counts": label_counts,
            "compact_label_terminal_state_count": 0,
        },
        "orbit_evidence": {
            "representative_target_stabilizer_order": fixed_target_actions,
            "labelled_target_orbit_size": action_count // fixed_target_actions,
            "action_image_stream_sha256": stream.hexdigest(),
            "production_independent_exact_equality": True,
            "independent_replay": independent,
        },
        "cutoff_correction": {
            "retracted_statement": "exactly nine literal pure targets",
            "reason": (
                "The exploratory command used --max-models 10, preserved no "
                "UNSAT transcript, and the single frozen target already has "
                "60,480 distinct source-stabilizer images."
            ),
            "minimum_literal_pure_target_count": 60_480,
            "complete_quotient_incidence_census_available": False,
        },
        "source_files": {
            str(path.relative_to(ROOT)): _sha(path)
            for path in (
                PRODUCTION_SOURCE,
                INDEPENDENT_ORBIT_SOURCE,
                INDEPENDENT_LABEL_SOURCE,
                LABEL_SOURCE,
            )
        },
        "theorem": {
            "whole_60480_member_labelled_target_orbit_compact_label_excluded": True,
            "quotient_incidence_census_complete": False,
            "whole_source_orbit_compact_no_go": False,
            "whole_record1_a1_compact_no_go": False,
            "finite_volume_noncompact_chamber_excluded": False,
        },
    }
    document["certificate_digest"] = _digest(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record1_a1_orbit0_target_orbit_compact_no_go":
        raise ValueError("wrong Record-1 target-orbit certificate kind")
    sources = document.get("source_files")
    if not isinstance(sources, Mapping):
        raise ValueError("malformed Record-1 source hash map")
    for path in (
        PRODUCTION_SOURCE,
        INDEPENDENT_ORBIT_SOURCE,
        INDEPENDENT_LABEL_SOURCE,
        LABEL_SOURCE,
    ):
        if sources.get(str(path.relative_to(ROOT))) != _sha(path):
            raise ValueError(f"Record-1 dependency changed: {path.name}")
    theorem = document.get("theorem")
    orbit = document.get("orbit_evidence")
    correction = document.get("cutoff_correction")
    if not all(isinstance(value, Mapping) for value in (theorem, orbit, correction)):
        raise ValueError("malformed Record-1 theorem evidence")
    assert isinstance(theorem, Mapping)
    assert isinstance(orbit, Mapping)
    assert isinstance(correction, Mapping)
    if (
        orbit.get("labelled_target_orbit_size") != 60_480
        or orbit.get("representative_target_stabilizer_order") != 1
        or not orbit.get("production_independent_exact_equality")
    ):
        raise ValueError("Record-1 orbit-stabilizer evidence changed")
    if (
        not theorem.get("whole_60480_member_labelled_target_orbit_compact_label_excluded")
        or theorem.get("quotient_incidence_census_complete")
        or theorem.get("whole_source_orbit_compact_no_go")
        or theorem.get("finite_volume_noncompact_chamber_excluded")
    ):
        raise ValueError("Record-1 theorem scope changed")
    if correction.get("minimum_literal_pure_target_count") != 60_480:
        raise ValueError("Record-1 cutoff correction changed")
    stored = document.get("certificate_digest")
    body = dict(document)
    body.pop("certificate_digest", None)
    digest = _digest(body)
    if stored != digest:
        raise ValueError("Record-1 certificate digest mismatch")
    if EXPECTED_CERTIFICATE_DIGEST not in ("TO_BE_FROZEN", digest):
        raise ValueError("Record-1 target-orbit certificate is not frozen")
    return {
        "record_index": 1,
        "schedule_a": 1,
        "source_orbit_index": 0,
        "labelled_target_orbit_size": 60_480,
        "target_orbit_compact_label_no_go": True,
        "quotient_incidence_census_complete": False,
        "whole_source_orbit_compact_no_go": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "certificate_digest": digest,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh", action="store_true")
    parser.add_argument("--write", type=Path)
    args = parser.parse_args()
    if args.fresh:
        document = fresh_record()
        if args.write is not None:
            args.write.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
        print(json.dumps(document, indent=2, sort_keys=True))
        return
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    print(json.dumps(verify(document), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
