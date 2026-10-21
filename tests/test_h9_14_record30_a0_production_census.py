from __future__ import annotations

from collections import Counter, defaultdict
from hashlib import sha256
import json
from pathlib import Path

import pytest

from scripts.independent_h9_14_record30_a0_target_labels import replay


ROOT = Path(__file__).resolve().parents[1]
CENSUS = (
    ROOT
    / "data/enumerations/h9_14_record30_a0_connected_ridge_census.json"
)
LABELS = (
    ROOT / "data/enumerations/h9_14_record30_a0_all_target_labels.json"
)
SOURCE = (
    ROOT / "scripts/enumerate_h9_14_record30_a0_ridge_components.cpp"
)
EXPECTED_CENSUS_SHA256 = (
    "3de6293068012c6213fcae03fb824d773cab6ead16c7820242482d5752cfbf0c"
)
EXPECTED_LABEL_FILE_SHA256 = (
    "3bebba2d96461e2d46151f42c2258e1818267e9b6c7c13beeb1d30e1875e805b"
)
EXPECTED_SOURCE_SHA256 = (
    "e6d7c4a6b42ce3eb0c4b59fc45000c18d8299e614c6ba7e276c90e249c4af430"
)
EXPECTED_F_VECTOR = (14, 88, 328, 806, 1364, 1596, 1248, 594, 132)
EXPECTED_STAGE_COUNTS = (
    (420, 126, 18, 0),
    (420, 96, 6, 0),
    (420, 96, 6, 0),
    (420, 96, 30, 0),
    (420, 126, 18, 0),
    (420, 96, 6, 0),
    (420, 96, 6, 0),
    (420, 96, 6, 0),
)
N = 14
FULL = (1 << N) - 1
EDGE_MASKS = ((1 << 0) | (1 << 11), (1 << 1) | (1 << 12), (1 << 2) | (1 << 13))


def _mask(vertices: list[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def _hits_all(value: int, clutter: tuple[int, ...]) -> bool:
    return all(value & member for member in clutter)


def _minimal_transversals(clutter: tuple[int, ...]) -> tuple[int, ...]:
    result = []
    for value in range(1, 1 << N):
        if not _hits_all(value, clutter):
            continue
        if all(
            not _hits_all(value ^ (1 << vertex), clutter)
            for vertex in range(N)
            if value >> vertex & 1
        ):
            result.append(value)
    return tuple(result)


def _face_table(facet_complements: tuple[int, ...]) -> set[int]:
    result = {0}
    for complement in facet_complements:
        facet = FULL ^ complement
        subset = facet
        while True:
            result.add(subset)
            if subset == 0:
                break
            subset = (subset - 1) & facet
    return result


def _minimal_nonfaces(faces: set[int]) -> tuple[int, ...]:
    result = []
    for value in range(1, 1 << N):
        if value in faces:
            continue
        if all(
            value ^ (1 << vertex) in faces
            for vertex in range(N)
            if value >> vertex & 1
        ):
            result.append(value)
    return tuple(result)


@pytest.fixture(scope="module")
def census() -> dict[str, object]:
    return json.loads(CENSUS.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def labels() -> dict[str, object]:
    return json.loads(LABELS.read_text(encoding="utf-8"))


def test_frozen_production_transcript_and_scope(census) -> None:
    assert sha256(CENSUS.read_bytes()).hexdigest() == EXPECTED_CENSUS_SHA256
    assert sha256(SOURCE.read_bytes()).hexdigest() == EXPECTED_SOURCE_SHA256
    assert census["kind"] == "h9_14_record30_a0_connected_ridge_census"
    assert census["canonical_seed"] == [0, 1, 2, 3, 4]
    assert census["candidate_complement_count"] == 326
    assert census["ridge_group_count"] == 833
    assert census["search_nodes"] == 28_039_413
    assert census["target_paths"] == 19
    assert census["target_orbit_count"] == 8
    assert len(census["target_orbit_representatives"]) == 8
    assert "compact a=0 incidence census only" in census["scope"].lower()
    assert "no finite-volume" in census["scope"].lower()


def test_all_eight_targets_pass_literal_schedule_and_ridge_audit(census) -> None:
    for rows in census["target_orbit_representatives"]:
        a5s = tuple(_mask(row) for row in rows)
        assert len(a5s) == len(set(a5s)) == 8
        clutter = (*EDGE_MASKS, *a5s)
        blocker = _minimal_transversals(clutter)
        assert Counter(value.bit_count() for value in blocker) == {5: 132}

        faces = _face_table(blocker)
        assert tuple(
            sum(value.bit_count() == size for value in faces)
            for size in range(1, 10)
        ) == EXPECTED_F_VECTOR
        assert tuple(sorted(_minimal_nonfaces(faces))) == tuple(sorted(clutter))

        ridges: dict[int, list[int]] = defaultdict(list)
        for index, complement in enumerate(blocker):
            outside = FULL ^ complement
            for vertex in range(N):
                if outside >> vertex & 1:
                    ridges[complement | (1 << vertex)].append(index)
        assert len(ridges) == 594
        assert Counter(map(len, ridges.values())) == {2: 594}

        adjacency = [set() for _ in blocker]
        for left, right in ridges.values():
            adjacency[left].add(right)
            adjacency[right].add(left)
        assert Counter(map(len, adjacency)) == {9: 132}
        reached = {0}
        stack = [0]
        while stack:
            current = stack.pop()
            for other in adjacency[current]:
                if other not in reached:
                    reached.add(other)
                    stack.append(other)
        assert len(reached) == 132


def test_frozen_all_target_label_replay_is_exact(census, labels) -> None:
    assert sha256(LABELS.read_bytes()).hexdigest() == EXPECTED_LABEL_FILE_SHA256
    assert replay(census) == labels
    assert labels["summary_sha256"] == (
        "3b6162151fcf6da220e063099da955ab2c0590f88682d643fbb710b6835ce9be"
    )
    assert labels["target_family_sha256"] == (
        "9f8d03f5dbf2d6500a2b052a0fd45c82feb28641eeb697d4329063d4a5f691ac"
    )
    assert labels["target_orbit_count"] == 8
    assert labels["terminal_sparse_state_count_histogram"] == {"0": 8}
    assert labels["unrestricted_compact_label_no_go_target_count"] == 8
    assert labels["sparse_label_survivor_count"] == 0
    assert labels[
        "all_nonempty_input_targets_have_unrestricted_compact_label_no_go"
    ] is True
    assert labels["finite_volume_noncompact_chamber_excluded"] is False
    assert labels["new_hyperbolic_coxeter_polytope"] is False
    assert tuple(
        tuple(stage["passing_state_count"] for stage in audit["stages"])
        for audit in labels["target_audits"]
    ) == EXPECTED_STAGE_COUNTS
