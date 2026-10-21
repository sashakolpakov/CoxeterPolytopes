#!/usr/bin/env python3
"""Independent direct replay of the five frozen record-37 blockers."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "data/enumerations/h9_14_record37_certificate.json"
EXPECTED_CERTIFICATE_SHA256 = (
    "4da8617bd4dfc77889c3c3f987e7988efd87d6d366d7efd29c7403370536cc5c"
)
VERTEX_COUNT = 14
FULL_MASK = (1 << VERTEX_COUNT) - 1


def _canonical_sha256(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _minimal_transversals(hyperedges: tuple[int, ...]) -> tuple[int, ...]:
    answers = []
    for candidate in range(1, FULL_MASK + 1):
        if not all(candidate & edge for edge in hyperedges):
            continue
        bit = 1
        minimal = True
        while bit <= candidate:
            if candidate & bit and all((candidate ^ bit) & edge for edge in hyperedges):
                minimal = False
                break
            bit <<= 1
        if minimal:
            answers.append(candidate)
    return tuple(answers)


def _packed_face_table_sha256(facet_complements: tuple[int, ...]) -> str:
    table = [False] * (1 << VERTEX_COUNT)
    for complement in facet_complements:
        facet = FULL_MASK ^ complement
        subset = facet
        while True:
            table[subset] = True
            if subset == 0:
                break
            subset = (subset - 1) & facet
    packed = bytes(
        sum(bool(table[index + offset]) << offset for offset in range(8))
        for index in range(0, len(table), 8)
    )
    return sha256(packed).hexdigest()


def _ridge_data(
    facet_complements: tuple[int, ...],
) -> tuple[dict[int, int], int]:
    facets = tuple(FULL_MASK ^ complement for complement in facet_complements)
    ridge_to_facets: dict[int, list[int]] = {}
    for index, facet in enumerate(facets):
        remaining = facet
        while remaining:
            bit = remaining & -remaining
            ridge_to_facets.setdefault(facet ^ bit, []).append(index)
            remaining ^= bit
    histogram = dict(
        sorted(Counter(len(indices) for indices in ridge_to_facets.values()).items())
    )
    adjacency = [set() for _ in facets]
    # Match the usual dual graph: only a ridge incident to exactly two facets
    # contributes a dual edge.  Degree-one and degree-three ridges are already
    # the pseudomanifold obstruction being certified here.
    for indices in ridge_to_facets.values():
        if len(indices) == 2:
            first, second = indices
            adjacency[first].add(second)
            adjacency[second].add(first)
    unseen = set(range(len(facets)))
    components = 0
    while unseen:
        components += 1
        stack = [unseen.pop()]
        while stack:
            current = stack.pop()
            new = adjacency[current] & unseen
            unseen.difference_update(new)
            stack.extend(new)
    return histogram, components


def replay() -> dict[str, object]:
    manifest = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    critical = dict(manifest)
    embedded = critical.pop("certificate_sha256")
    digest = _canonical_sha256(critical)
    if embedded != digest or digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("record-37 certificate digest mismatch")
    _, edges = decode_graph6(manifest["graph6"])
    edge_masks = tuple(sum(1 << vertex for vertex in edge) for edge in edges)
    target = manifest["target_count_state_audits"]
    multiplicities = Counter(
        item["family_index"] for item in target["raw_to_unique_family_map"]
    )
    if multiplicities != {0: 3, 1: 2, 2: 2, 3: 3, 4: 2}:
        raise AssertionError("record-37 12-to-5 deduplication map changed")
    family_replays = []
    for family in target["families"]:
        missing_fives = tuple(family["missing_five_masks"])
        blocker = _minimal_transversals(edge_masks + missing_fives)
        histogram = dict(sorted(Counter(mask.bit_count() for mask in blocker).items()))
        if histogram != {5: 102}:
            raise AssertionError("a frozen record-37 blocker lost purity")
        if list(blocker) != family["facet_complement_masks"]:
            raise AssertionError("a frozen record-37 blocker list changed")
        blocker_digest = sha256(
            b"".join(mask.to_bytes(2, "little") for mask in blocker)
        ).hexdigest()
        if blocker_digest != family["facet_complement_sha256"]:
            raise AssertionError("record-37 blocker digest mismatch")
        face_digest = _packed_face_table_sha256(blocker)
        if face_digest != family["packed_face_table_sha256"]:
            raise AssertionError("record-37 face-table digest mismatch")
        ridge_histogram, component_count = _ridge_data(blocker)
        expected_ridges = {
            int(key): value for key, value in family["ridge_degree_histogram"].items()
        }
        if ridge_histogram != expected_ridges:
            raise AssertionError("record-37 ridge histogram mismatch")
        if component_count != family["dual_component_count"]:
            raise AssertionError("record-37 dual-component count mismatch")
        if ridge_histogram == {2: 459}:
            raise AssertionError("a record-37 target family unexpectedly became closed")
        family_replays.append(
            {
                "missing_five_masks": list(missing_fives),
                "blocker_sha256": blocker_digest,
                "face_table_sha256": face_digest,
                "ridge_degree_histogram": ridge_histogram,
                "dual_component_count": component_count,
            }
        )
    return {
        "certificate_sha256": digest,
        "raw_target_state_count": 12,
        "unique_family_count": len(family_replays),
        "passing_family_count": 0,
        "families": family_replays,
    }


def main() -> None:
    print(json.dumps(replay(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
