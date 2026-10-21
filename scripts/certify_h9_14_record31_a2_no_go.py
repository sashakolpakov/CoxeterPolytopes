#!/usr/bin/env python3
"""Freeze and verify the compact-only Record 31 a=2 schedule no-go."""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import permutations
import json
import os
from pathlib import Path
import subprocess
import tempfile
from typing import Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record31_a2_compact_no_go.json"
CHECKPOINT = ROOT / "data/enumerations/h9_14_record31_a2_working_checkpoint.json"
POSITIVE_ARCHIVE = ROOT / "data/enumerations/h9_14_record31_positive_a5.json"
CPP_SOURCE = ROOT / "scripts/explore_h9_14_record31_a2_additions.cpp"
PRODUCTION_DOWNSTREAM = ROOT / "scripts/explore_h9_14_record31_a2_downstream.py"
INDEPENDENT_LINK = ROOT / "scripts/independent_h9_14_record31_a2_link_audit.py"
INDEPENDENT_DOWNSTREAM = ROOT / "scripts/independent_h9_14_record31_a2_downstream_audit.py"
CERTIFIER = Path(__file__).resolve()

GRAPH6 = "M?????????O?O?o??"
EXPECTED_F_VECTOR = [14, 87, 319, 769, 1275, 1463, 1126, 531, 118]
EXPECTED_FAILURES = [
    4, 7, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 31, 32, 36, 38, 39, 41, 43, 45, 46, 48,
]
EXPECTED_RESIDUALS = [index for index in range(49) if index not in EXPECTED_FAILURES]
EXPECTED_TARGET_DISTRIBUTION = {0: 1, 2: 1, 3: 1, 5: 1, 34: 6, 40: 2}
ACTIVE = (0, 1, 11, 12, 13)
ISOLATES = tuple(range(2, 11))
ACTIVE_FLIP = {0: 1, 1: 0, 11: 12, 12: 11, 13: 13}


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"Record 31 a=2 {name} is malformed")
    return value


def _source(path: Path) -> dict[str, str]:
    return {
        "path": str(path.relative_to(ROOT)),
        "sha256": sha256(path.read_bytes()).hexdigest(),
    }


def _digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _mask(vertices: Iterable[int]) -> int:
    return sum(1 << int(vertex) for vertex in vertices)


def canonical_family_key(raw_family: Sequence[Sequence[int]]) -> tuple[int, ...]:
    """Complete direct invariant for C2 x S9 with A4/A5 row colors."""

    family = tuple(_mask(face) for face in raw_family)
    groups = {
        order: tuple(face for face in family if face.bit_count() == order)
        for order in (4, 5)
    }
    if len(groups[4]) != 2 or not 3 <= len(groups[5]) <= 6:
        raise ValueError("Record 31 a=2 target family shape changed")
    active_mask = _mask(ACTIVE)
    best: tuple[int, ...] | None = None
    for flip in (False, True):
        moved: dict[int, tuple[int, ...]] = {}
        for order, rows in groups.items():
            values = []
            for face in rows:
                if flip:
                    face = (face & ~active_mask) | _mask(
                        ACTIVE_FLIP[vertex]
                        for vertex in ACTIVE
                        if face >> vertex & 1
                    )
                values.append(face)
            moved[order] = tuple(values)
        for fours in permutations(moved[4]):
            for fives in permutations(moved[5]):
                rows = fours + fives
                active_rows = tuple(
                    sum(((face >> vertex) & 1) << index
                        for index, vertex in enumerate(ACTIVE))
                    for face in rows
                )
                isolate_columns = tuple(
                    sorted(
                        sum(((face >> vertex) & 1) << index
                            for index, face in enumerate(rows))
                        for vertex in ISOLATES
                    )
                )
                candidate = active_rows + isolate_columns
                if best is None or candidate < best:
                    best = candidate
    if best is None:
        raise AssertionError("Record 31 a=2 canonicalization failed")
    return best


def _target_rows(census: Mapping[str, object]) -> list[Sequence[Sequence[int]]]:
    rows = census.get("orbits")
    if not isinstance(rows, list) or len(rows) != 49:
        raise ValueError("Record 31 a=2 incidence orbit coverage changed")
    return [
        family
        for row in rows
        for family in _mapping(row, "incidence orbit").get(
            "ridge_valid_target_minimal_nonfaces", []
        )
    ]


def _target_distribution(census: Mapping[str, object]) -> dict[int, int]:
    rows = census.get("orbits")
    if not isinstance(rows, list):
        raise ValueError("Record 31 a=2 incidence rows changed")
    return {
        int(row["orbit_index"]): int(row["ridge_valid_target_orbit_count"])
        for row in rows
        if int(row["ridge_valid_target_orbit_count"])
    }


def _failure_indices(census: Mapping[str, object]) -> list[int]:
    rows = census.get("orbits")
    if not isinstance(rows, list):
        raise ValueError("Record 31 a=2 incidence rows changed")
    return [
        int(row["orbit_index"])
        for row in rows
        if bool(row["fixed_local_failure"])
    ]


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record31_a2_compact_schedule_no_go":
        raise ValueError("wrong Record 31 a=2 certificate kind")
    if document.get("graph6") != GRAPH6:
        raise ValueError("Record 31 graph source changed")
    if document.get("schedule") != {
        "a": 2,
        "missing_five_count": 21,
        "facet_count": 118,
        "f_vector": EXPECTED_F_VECTOR,
    }:
        raise ValueError("Record 31 a=2 schedule changed")
    expected_scope = {
        "compact_only": True,
        "entire_a2_schedule_excluded": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "new_hyperbolic_coxeter_polytope": False,
    }
    if document.get("scope") != expected_scope:
        raise ValueError("Record 31 a=2 scope changed or was over-promoted")

    production = _mapping(document.get("production"), "production census")
    independent = _mapping(document.get("independent_incidence"), "independent census")
    if production.get("enumeration_mode") != "nauty_partial_quotient":
        raise ValueError("Record 31 a=2 production quotient changed")
    if independent.get("enumeration_mode") != "direct_binary_matrix_partial_quotient":
        raise ValueError("Record 31 a=2 independent quotient changed")
    if production.get("orbit_count") != 49 or independent.get("orbit_count") != 49:
        raise ValueError("Record 31 a=2 positive-A4 coverage changed")
    if _failure_indices(production) != EXPECTED_FAILURES:
        raise ValueError("production strict face-link partition changed")
    if _failure_indices(independent) != EXPECTED_FAILURES:
        raise ValueError("independent-incidence strict face-link partition changed")
    for census in (production, independent):
        rows = census["orbits"]
        for expected_index, raw in enumerate(rows):
            row = _mapping(raw, "incidence orbit")
            if int(row["orbit_index"]) != expected_index:
                raise ValueError("Record 31 a=2 incidence row order changed")
            witnesses = row.get("fixed_local_failure_witnesses")
            if not isinstance(witnesses, list):
                raise ValueError("Record 31 a=2 local witness list changed")
            if bool(row["fixed_local_failure"]) != bool(witnesses):
                raise ValueError("Record 31 a=2 local witness flag changed")
            for witness in witnesses:
                witness = _mapping(witness, "local witness")
                interval = witness.get("additional_A5_incidence_interval")
                if not isinstance(interval, list) or len(interval) != 2 or interval[0] <= interval[1]:
                    raise ValueError("Record 31 a=2 local witness is not contradictory")
                if witness.get("g_prefix", [])[-1] != 0 or interval != [0, -1]:
                    raise ValueError("Record 31 a=2 frozen strict-link witness changed")
    if _target_distribution(production) != EXPECTED_TARGET_DISTRIBUTION:
        raise ValueError("production Record 31 a=2 target distribution changed")
    if _target_distribution(independent) != EXPECTED_TARGET_DISTRIBUTION:
        raise ValueError("independent Record 31 a=2 target distribution changed")
    production_keys = {canonical_family_key(family) for family in _target_rows(production)}
    independent_keys = {canonical_family_key(family) for family in _target_rows(independent)}
    if len(production_keys) != 12 or production_keys != independent_keys:
        raise ValueError("Record 31 a=2 independent target-orbit equality failed")

    link = _mapping(document.get("independent_link_audit"), "independent link audit")
    if (
        link.get("kind") != "independent_h9_14_record31_a2_strict_face_link_audit"
        or link.get("strict_face_link_excluded_orbit_indices") != EXPECTED_FAILURES
        or link.get("residual_orbit_indices") != EXPECTED_RESIDUALS
    ):
        raise ValueError("independent Record 31 a=2 link audit changed")
    link_rows = link.get("orbits")
    if not isinstance(link_rows, list) or len(link_rows) != 49:
        raise ValueError("independent Record 31 a=2 link row coverage changed")
    for index in range(49):
        production_witnesses = production["orbits"][index]["fixed_local_failure_witnesses"]
        independent_witnesses = link_rows[index]["witnesses"]
        # The two traversals intentionally choose different, symmetry-related
        # subsets.  Freeze the complete multiset of numerical contradictions,
        # not an accidental first-witness ordering.
        def numerical_signature(witness: Mapping[str, object]) -> str:
            payload = {
                key: value for key, value in witness.items() if key != "subset"
            }
            return json.dumps(payload, sort_keys=True, separators=(",", ":"))

        if sorted(map(numerical_signature, production_witnesses)) != sorted(
            map(numerical_signature, independent_witnesses)
        ):
            raise ValueError("production/independent strict-link numerical witnesses differ")
        for raw_witness in independent_witnesses:
            witness = _mapping(raw_witness, "independent local witness")
            subset = witness.get("subset")
            interval = witness.get("additional_A5_incidence_interval")
            if (
                not isinstance(subset, list)
                or not 1 <= len(subset) <= 3
                or len(set(subset)) != len(subset)
                or not isinstance(interval, list)
                or interval != [0, -1]
                or witness.get("g_prefix", [])[-1] != 0
            ):
                raise ValueError("independent Record 31 a=2 local witness changed")

    checkpoint_targets = _mapping(document.get("checkpoint"), "checkpoint").get(
        "ridge_valid_target_orbits"
    )
    if not isinstance(checkpoint_targets, list):
        raise ValueError("Record 31 a=2 checkpoint targets changed")
    checkpoint_keys = {canonical_family_key(family) for family in checkpoint_targets}
    if checkpoint_keys != production_keys:
        raise ValueError("Record 31 a=2 checkpoint lost a target orbit")

    downstream = _mapping(document.get("downstream"), "downstream audit")
    if (
        downstream.get("kind") != "h9_14_record31_a2_downstream_compact_audit"
        or downstream.get("input_target_orbit_count") != 12
        or downstream.get("connected_mod2_homology_sphere_orbit_count") != 12
        or downstream.get("unrestricted_compact_label_no_go_orbit_count") != 12
        or downstream.get("compact_label_survivor_orbit_count") != 0
    ):
        raise ValueError("Record 31 a=2 production downstream result changed")
    downstream_rows = downstream.get("target_orbits")
    if not isinstance(downstream_rows, list) or len(downstream_rows) != 12:
        raise ValueError("Record 31 a=2 production downstream rows changed")
    if any(
        row["sparse_compact_label_gate"]["terminal_sparse_state_count"] != 0
        for row in downstream_rows
    ):
        raise ValueError("a production Record 31 a=2 compact label survived")

    hostile = _mapping(document.get("independent_downstream"), "independent downstream")
    if (
        hostile.get("kind") != "independent_h9_14_record31_a2_downstream_audit"
        or hostile.get("classical_lanner_placement_counts") != {"4": 108, "5": 420}
        or hostile.get("literal_target_orbit_count") != 12
        or hostile.get("compact_label_survivor_count") != 0
        or hostile.get("all_targets_have_unrestricted_compact_label_no_go") is not True
    ):
        raise ValueError("independent Record 31 a=2 downstream result changed")
    hostile_rows = hostile.get("ridge_survivor_label_audits")
    if not isinstance(hostile_rows, list) or len(hostile_rows) != 12:
        raise ValueError("independent Record 31 a=2 label rows changed")
    if any(row["terminal_sparse_state_count"] != 0 for row in hostile_rows):
        raise ValueError("an independent Record 31 a=2 compact label survived")

    provenance = document.get("source_provenance")
    expected_provenance = [
        _source(path)
        for path in (
            CPP_SOURCE,
            PRODUCTION_DOWNSTREAM,
            INDEPENDENT_LINK,
            INDEPENDENT_DOWNSTREAM,
            POSITIVE_ARCHIVE,
            CHECKPOINT,
            CERTIFIER,
        )
    ]
    if provenance != expected_provenance:
        raise ValueError("Record 31 a=2 source provenance changed")
    if document.get("certificate_sha256") != _digest(document):
        raise ValueError("Record 31 a=2 certificate digest changed")
    return {
        "strict_face_link_excluded_orbit_count": 31,
        "enumerated_residual_orbit_count": 18,
        "ridge_target_orbit_count": 12,
        "compact_label_survivor_count": 0,
        "certificate_sha256": document["certificate_sha256"],
    }


def assemble(
    production_path: Path,
    independent_path: Path,
    link_path: Path,
    downstream_path: Path,
    independent_downstream_path: Path,
) -> dict[str, object]:
    document: dict[str, object] = {
        "kind": "h9_14_record31_a2_compact_schedule_no_go",
        "claim": "Every completion of H9/14 graph-shadow Record 31 in the compact Dehn--Sommerville schedule a=2 fails an exact necessary compact Coxeter condition.",
        "graph6": GRAPH6,
        "schedule": {
            "a": 2,
            "missing_five_count": 21,
            "facet_count": 118,
            "f_vector": EXPECTED_F_VECTOR,
        },
        "scope": {
            "compact_only": True,
            "entire_a2_schedule_excluded": True,
            "whole_shadow_compact_coxeter_theorem": False,
            "finite_volume_noncompact_chamber_excluded": False,
            "new_hyperbolic_coxeter_polytope": False,
        },
        "production": json.loads(production_path.read_text(encoding="utf-8")),
        "independent_incidence": json.loads(independent_path.read_text(encoding="utf-8")),
        "independent_link_audit": json.loads(link_path.read_text(encoding="utf-8")),
        "checkpoint": json.loads(CHECKPOINT.read_text(encoding="utf-8")),
        "downstream": json.loads(downstream_path.read_text(encoding="utf-8")),
        "independent_downstream": json.loads(
            independent_downstream_path.read_text(encoding="utf-8")
        ),
        "source_provenance": [
            _source(path)
            for path in (
                CPP_SOURCE,
                PRODUCTION_DOWNSTREAM,
                INDEPENDENT_LINK,
                INDEPENDENT_DOWNSTREAM,
                POSITIVE_ARCHIVE,
                CHECKPOINT,
                CERTIFIER,
            )
        ],
        "nonclaims": [
            "No noncompact finite-volume chamber is excluded.",
            "No parabolic or Vinberg endpoint census is performed.",
            "The topology audit is necessary and does not itself prove convex polytopality.",
            "No new hyperbolic Coxeter polytope is constructed.",
            "The compact schedules a=0,1 of Record 31 remain open.",
        ],
    }
    document["certificate_sha256"] = _digest(document)
    verify(document)
    return document


def _input_projection() -> str:
    document = json.loads(POSITIVE_ARCHIVE.read_text(encoding="utf-8"))
    branch = next(row for row in document["production"]["branches"] if row["a"] == 2)
    lines = [str(len(branch["orbits"]))]
    for index, row in enumerate(branch["orbits"]):
        masks = [_mask(face) for face in row["representative"]]
        lines.append(f"{index} {masks[0]} {masks[1]}")
    return "\n".join(lines) + "\n"


def fresh() -> dict[str, object]:
    with tempfile.TemporaryDirectory(prefix="record31-a2-") as directory:
        directory_path = Path(directory)
        binary = directory_path / "enumerate"
        subprocess.run(
            [
                "c++", "-std=c++20", "-O3", "-DNDEBUG",
                "-I/opt/homebrew/include/nauty", str(CPP_SOURCE),
                "-L/opt/homebrew/lib", "-lnauty", "-o", str(binary),
            ],
            check=True,
            cwd=ROOT,
        )
        projection = _input_projection()
        paths = {
            "production": directory_path / "production.json",
            "independent": directory_path / "independent.json",
            "link": directory_path / "link.json",
            "downstream": directory_path / "downstream.json",
            "independent_downstream": directory_path / "independent_downstream.json",
        }
        for mode, key in (("--quotient", "production"), ("--matrix-quotient", "independent")):
            result = subprocess.run(
                [str(binary), mode], input=projection, text=True,
                capture_output=True, check=True, cwd=ROOT,
            )
            paths[key].write_text(result.stdout, encoding="utf-8")
        subprocess.run(
            ["python3", str(INDEPENDENT_LINK), str(POSITIVE_ARCHIVE)],
            check=True, cwd=ROOT, stdout=paths["link"].open("w", encoding="utf-8"),
        )
        # The checked-in checkpoint is exactly the canonical target handoff;
        # verify() compares it with both fresh incidence enumerators.
        subprocess.run(
            ["python3", str(PRODUCTION_DOWNSTREAM), str(CHECKPOINT)],
            check=True, cwd=ROOT, stdout=paths["downstream"].open("w", encoding="utf-8"),
            env={**os.environ, "PYTHONPATH": f"{ROOT}:{ROOT / 'src'}"},
        )
        subprocess.run(
            ["python3", str(INDEPENDENT_DOWNSTREAM), str(CHECKPOINT)],
            check=True, cwd=ROOT,
            stdout=paths["independent_downstream"].open("w", encoding="utf-8"),
            env={**os.environ, "PYTHONPATH": f"{ROOT}:{ROOT / 'src'}"},
        )
        return assemble(
            paths["production"], paths["independent"], paths["link"],
            paths["downstream"], paths["independent_downstream"],
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh", action="store_true")
    parser.add_argument("--assemble", action="store_true")
    parser.add_argument("--production", type=Path)
    parser.add_argument("--independent", type=Path)
    parser.add_argument("--link", type=Path)
    parser.add_argument("--downstream", type=Path)
    parser.add_argument("--independent-downstream", type=Path)
    args = parser.parse_args()
    if args.fresh:
        document = fresh()
        ARCHIVE.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    elif args.assemble:
        required = (
            args.production, args.independent, args.link, args.downstream,
            args.independent_downstream,
        )
        if any(path is None for path in required):
            parser.error("--assemble requires all five input paths")
        document = assemble(*required)  # type: ignore[arg-type]
        ARCHIVE.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    else:
        document = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    print(json.dumps(verify(document), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
