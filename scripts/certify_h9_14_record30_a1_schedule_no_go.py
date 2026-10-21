#!/usr/bin/env python3
"""Build and verify the compact Record-30 ``a=1`` no-go certificate.

The certificate joins three logically separate computations:

* the nauty-quotiented discovery census;
* a structurally independent colored-incidence exhaustion with no nauty; and
* the standalone complete finite Lanner/actual-vertex label gate.

The two incidence programs are compared as exact labelled target sets after
expansion under the fixed-A4 stabilizer.  SHA-256 values are transcript
checksums only; the equality assertion is performed on the literal Python
sets before any digest is recorded.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import permutations, product
import json
import os
from pathlib import Path
import subprocess
import tempfile
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record30_a1_schedule_no_go.json"
PRODUCTION_SOURCE = ROOT / "scripts/explore_h9_14_record30_a1_nauty.cpp"
INDEPENDENT_SOURCE = (
    ROOT / "scripts/independent_h9_14_record30_a1_exhaustive.cpp"
)
LABEL_SOURCE = ROOT / "scripts/independent_h9_14_record30_a1_labels.py"

EDGES = ((0, 11), (1, 12), (2, 13))
ISOLATES = tuple(range(3, 11))
EXPECTED_STATE_COUNTS = (
    (1, 11, 256, 9_229, 309_175),
    (1, 9, 153, 4_123, 159_275, 1_399_930),
    (0, 0, 0, 0, 0, 0, 0),
)
EXPECTED_TARGET_ORBIT_SIZES = ((48,), (80, 40), ())
EXPECTED_LABEL_TRANSCRIPTS = ((108, 84, 0), (108, 84, 36, 108, 0), (108, 84, 0))


def _json_bytes(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def _digest(value: object) -> str:
    return sha256(_json_bytes(value)).hexdigest()


def _mask(vertices: Iterable[int]) -> int:
    answer = 0
    for vertex in vertices:
        answer |= 1 << vertex
    return answer


def _transform(mask: int, action: tuple[int, ...]) -> int:
    answer = 0
    for vertex in range(14):
        if mask & (1 << vertex):
            answer |= 1 << action[vertex]
    return answer


def _active_actions() -> tuple[tuple[int, ...], ...]:
    answer: list[tuple[int, ...]] = []
    for component_order in permutations(range(3)):
        for flips in product(range(2), repeat=3):
            action = list(range(14))
            for source, target in enumerate(component_order):
                flip = flips[source]
                action[EDGES[source][0]] = EDGES[target][flip]
                action[EDGES[source][1]] = EDGES[target][1 - flip]
            answer.append(tuple(action))
    assert len(answer) == 48
    return tuple(answer)


ACTIVE_ACTIONS = _active_actions()


def _fixed_four_stabilizer(four_vertices: Iterable[int]) -> tuple[tuple[int, ...], ...]:
    four = _mask(four_vertices)
    inside = tuple(vertex for vertex in ISOLATES if four & (1 << vertex))
    outside = tuple(vertex for vertex in ISOLATES if not four & (1 << vertex))
    answer: list[tuple[int, ...]] = []
    for active in ACTIVE_ACTIONS:
        if _transform(four & sum(1 << vertex for edge in EDGES for vertex in edge), active) != (
            four & sum(1 << vertex for edge in EDGES for vertex in edge)
        ):
            continue
        for inside_image in permutations(inside):
            for outside_image in permutations(outside):
                action = list(active)
                for source, target in zip(inside, inside_image, strict=True):
                    action[source] = target
                for source, target in zip(outside, outside_image, strict=True):
                    action[source] = target
                transformed_four = _transform(four, tuple(action))
                if transformed_four != four:
                    raise AssertionError("constructed action did not fix A4")
                answer.append(tuple(action))
    expected = {70: 27_648, 336: 5_760}
    full_orbit_size = 1_935_360 // len(answer)
    assert len(answer) == expected[full_orbit_size]
    return tuple(answer)


def _target_orbits(
    row: dict[str, object],
) -> tuple[frozenset[tuple[int, ...]], ...]:
    representatives = row["target_addition_representatives"]
    assert isinstance(representatives, list)
    fixed_four = row["fixed_A4"]
    assert isinstance(fixed_four, list)
    stabilizer = _fixed_four_stabilizer(fixed_four)
    answer: set[tuple[int, ...]] = set()
    orbits: list[frozenset[tuple[int, ...]]] = []
    for representative in representatives:
        assert isinstance(representative, list)
        masks = tuple(_mask(face) for face in representative)
        orbit = {
            tuple(sorted(_transform(face, action) for face in masks))
            for action in stabilizer
        }
        if answer.intersection(orbit):
            raise AssertionError("target orbits overlap literally")
        answer.update(orbit)
        orbits.append(frozenset(orbit))
    assert len(answer) == row["literal_target_count"]
    return tuple(orbits)


def _families(row: dict[str, object]) -> frozenset[tuple[int, ...]]:
    return frozenset().union(*_target_orbits(row))


def _literal_sha256(orbit_index: int, families: frozenset[tuple[int, ...]]) -> str:
    transcript = [[orbit_index, *family] for family in sorted(families)]
    return _digest(transcript)


COMMON_COUNTERS = (
    "fixed_A4",
    "corrected_forced_nonminimal_A5_shadow_count",
    "additional_minimal_A5_count",
    "scheduled_missing_independent_five_count",
    "scheduled_facet_count",
    "explicit_stabilizer_order",
    "base_small_transversal_count",
    "initial_order_five_cover_count",
    "candidate_count",
    "raw_unordered_addition_count",
    "strict_fixed_link_exclusion",
    "search_node_count",
    "canonical_state_counts",
    "isomorphic_prune_count",
    "facet_lower_bound_prune_count",
    "uncoverable_prune_count",
    "scheduled_incidence_orbit_count",
    "scheduled_blocker_count_histogram",
    "pure_blocker_orbit_count",
    "face_table_reject_count",
    "ridge_reject_count",
    "necessary_link_reject_count",
    "necessary_incidence_target_orbit_count",
    "literal_target_count",
    "target_orbit_size_histogram",
    "target_orbit_sizes",
)


def _compare_incidence(
    production: list[dict[str, object]], independent: list[dict[str, object]]
) -> dict[str, object]:
    if len(production) != 3 or len(independent) != 3:
        raise AssertionError("expected exactly three fixed-A4 branches")
    literal_hashes: list[str] = []
    all_literal: list[list[int]] = []
    for orbit_index, (left, right) in enumerate(zip(production, independent, strict=True)):
        if left["positive_A4_orbit_index"] != orbit_index:
            raise AssertionError("production orbit order changed")
        if right["positive_A4_orbit_index"] != orbit_index:
            raise AssertionError("independent orbit order changed")
        for key in COMMON_COUNTERS:
            if left[key] != right[key]:
                raise AssertionError(f"orbit {orbit_index} counter {key} disagrees")
        left_families = _families(left)
        right_families = _families(right)
        # This is the proof-relevant comparison.  Digests below merely make
        # the already-established literal equality compact to audit.
        if left_families != right_families:
            raise AssertionError(f"orbit {orbit_index} literal target sets disagree")
        literal_hashes.append(_literal_sha256(orbit_index, left_families))
        all_literal.extend(
            [orbit_index, *family] for family in sorted(left_families)
        )
    return {
        "common_counter_names": list(COMMON_COUNTERS),
        "production_and_independent_common_counters_identical": True,
        "production_and_independent_literal_target_sets_identical": True,
        "fixed_A4_orbit_count": 3,
        "necessary_incidence_target_orbit_count": 3,
        "literal_target_count": 168,
        "literal_target_set_sha256_by_A4_orbit": literal_hashes,
        "literal_target_set_sha256_all_orbits": _digest(all_literal),
        "hash_role": "Reproducibility checksum after exact literal-set equality; not a substitute for equality.",
    }


def _compile(
    source: Path, executable: Path, *, nauty_include: Path | None = None,
    nauty_library: Path | None = None,
) -> None:
    compiler = os.environ.get("CXX", "clang++")
    command = [
        compiler,
        "-std=c++20",
        "-O3",
        "-Wall",
        "-Wextra",
        "-Werror",
    ]
    if nauty_include is not None:
        command.extend([f"-I{nauty_include}"])
    command.extend([str(source), "-o", str(executable)])
    if nauty_library is not None:
        command.append(str(nauty_library))
    subprocess.run(command, cwd=ROOT, check=True)


def _nauty_paths() -> tuple[Path, Path]:
    configured = os.environ.get("NAUTY_PREFIX")
    prefixes = ([Path(configured)] if configured else []) + [
        Path("/opt/homebrew/opt/nauty"),
        Path("/usr/local/opt/nauty"),
    ]
    for prefix in prefixes:
        header_directory = prefix / "include/nauty"
        library = prefix / "lib/libnauty.a"
        if header_directory.joinpath("nauty.h").exists() and library.exists():
            return header_directory, library
    raise RuntimeError("nauty headers/static library not found; set NAUTY_PREFIX")


def _run_json(command: list[str]) -> dict[str, object]:
    completed = subprocess.run(
        command, cwd=ROOT, check=True, text=True, capture_output=True
    )
    return json.loads(completed.stdout)


def build() -> dict[str, object]:
    include, library = _nauty_paths()
    with tempfile.TemporaryDirectory(prefix="record30-a1-certificate-") as directory:
        temporary = Path(directory)
        production_executable = temporary / "production"
        independent_executable = temporary / "independent"
        _compile(
            PRODUCTION_SOURCE,
            production_executable,
            nauty_include=include,
            nauty_library=library,
        )
        _compile(INDEPENDENT_SOURCE, independent_executable)
        production = [
            _run_json([str(production_executable), str(index)])
            for index in range(3)
        ]
        independent = [
            _run_json([str(independent_executable), str(index)])
            for index in range(3)
        ]
    cross_audit = _compare_incidence(production, independent)
    labels = _run_json([os.environ.get("PYTHON", "python3"), str(LABEL_SOURCE)])

    incidence_target_orbits = [
        _target_orbits(row) for row in production
    ]
    label_orbits: set[tuple[int, int]] = set()
    for audit in labels["label_audits"]:
        family = tuple(sorted(_mask(face) for face in audit["additional_A5s"]))
        memberships = [
            (orbit_index, target_index)
            for orbit_index, branch in enumerate(incidence_target_orbits)
            for target_index, orbit in enumerate(branch)
            if family in orbit
        ]
        if len(memberships) != 1:
            raise AssertionError("label target does not identify one incidence orbit")
        label_orbits.add(memberships[0])
    expected_label_orbits = {
        (orbit_index, target_index)
        for orbit_index, branch in enumerate(incidence_target_orbits)
        for target_index in range(len(branch))
    }
    if label_orbits != expected_label_orbits:
        raise AssertionError("label replay did not cover every incidence target orbit")

    source_hashes = {
        str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest()
        for path in (PRODUCTION_SOURCE, INDEPENDENT_SOURCE, LABEL_SOURCE)
    }
    certificate: dict[str, object] = {
        "kind": "h9_14_record30_a1_compact_schedule_no_go_certificate",
        "record_index": 30,
        "fixed_positive_A4_count": 1,
        "compact_only": True,
        "schedule": {
            "f_vector": [14, 88, 328, 805, 1358, 1582, 1232, 585, 130],
            "missing_independent_five_count": 14,
            "facet_count": 130,
            "corrected_shadow_counts": [10, 9, 8],
            "additional_minimal_A5_counts": [4, 5, 6],
        },
        "production_nauty_incidence_replay": production,
        "production_nauty_transcript_sha256": _digest(production),
        "independent_explicit_incidence_replay": independent,
        "independent_explicit_incidence_transcript_sha256": _digest(independent),
        "incidence_cross_audit": cross_audit,
        "independent_complete_label_replay": labels,
        "independent_complete_label_replay_sha256": _digest(labels),
        "source_sha256": source_hashes,
        "theorem_scope": {
            "record30_a1_compact_schedule_excluded": True,
            "finite_volume_noncompact_chamber_excluded": False,
            "whole_record30_shadow_excluded": False,
            "new_hyperbolic_coxeter_polytope": False,
        },
        "proof_chain": [
            "The three fixed-A4 orbits exhaust the a=1 action quotient.",
            "The two-active orbit is impossible by the strict edge-link g3 theorem.",
            "Nauty and a no-nauty colored-incidence traversal independently exhaust the other two branches and have exactly the same 168 labelled targets.",
            "The complete order-four/order-five Lanner gate leaves zero finite label states on all three target orbits.",
        ],
        "scope": "Compact Record-30 a=1 schedule only. No assertion about ideal vertices, noncompact finite volume, or existence of a new hyperbolic Coxeter polytope.",
    }
    certificate["certificate_sha256"] = _digest(certificate)
    verify(certificate)
    return certificate


def verify(certificate: dict[str, object]) -> dict[str, object]:
    stored_digest = certificate["certificate_sha256"]
    without_digest = dict(certificate)
    del without_digest["certificate_sha256"]
    if stored_digest != _digest(without_digest):
        raise AssertionError("certificate digest changed")
    production = certificate["production_nauty_incidence_replay"]
    independent = certificate["independent_explicit_incidence_replay"]
    assert isinstance(production, list) and isinstance(independent, list)
    comparison = _compare_incidence(production, independent)
    if comparison != certificate["incidence_cross_audit"]:
        raise AssertionError("archived incidence cross-audit changed")
    for index, rows in enumerate((production, independent)):
        del index
        state_counts = tuple(tuple(row["canonical_state_counts"]) for row in rows)
        if state_counts != EXPECTED_STATE_COUNTS:
            raise AssertionError("canonical state transcript changed")
        orbit_sizes = tuple(tuple(row["target_orbit_sizes"]) for row in rows)
        if orbit_sizes != EXPECTED_TARGET_ORBIT_SIZES:
            raise AssertionError("target orbit sizes changed")
    labels = certificate["independent_complete_label_replay"]
    transcripts = tuple(
        tuple(row["passing_state_count"] for row in audit["stages"])
        for audit in labels["label_audits"]
    )
    if transcripts != EXPECTED_LABEL_TRANSCRIPTS:
        raise AssertionError("finite-label transcript changed")
    if labels["terminal_finite_label_state_count_across_targets"] != 0:
        raise AssertionError("a target survived the complete finite-label gate")
    scope = certificate["theorem_scope"]
    if scope != {
        "record30_a1_compact_schedule_excluded": True,
        "finite_volume_noncompact_chamber_excluded": False,
        "whole_record30_shadow_excluded": False,
        "new_hyperbolic_coxeter_polytope": False,
    }:
        raise AssertionError("theorem scope changed")
    return {
        "record_index": 30,
        "fixed_positive_A4_count": 1,
        "necessary_incidence_target_orbit_count": 3,
        "literal_target_count": 168,
        "terminal_finite_label_state_count": 0,
        "record30_a1_compact_schedule_excluded": True,
        "finite_volume_noncompact_chamber_excluded": False,
        "new_hyperbolic_coxeter_polytope": False,
        "certificate_sha256": stored_digest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--recompute", action="store_true")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    if arguments.recompute:
        certificate = build()
        if arguments.write:
            ARCHIVE.write_text(
                json.dumps(certificate, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
    else:
        certificate = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    print(json.dumps(verify(certificate), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
