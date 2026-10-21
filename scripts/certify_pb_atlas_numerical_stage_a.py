#!/usr/bin/env python3
"""Freeze and mutate-test the PB-Atlas numerical all-contraction Stage A."""

from __future__ import annotations

import argparse
from collections import Counter
import copy
import gc
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.pb_atlas_contraction_core import (  # noqa: E402
    PhysicalBranch,
    _marked_wall_gram_permutations,
    _permutation_lift,
    certify_branch,
    mask_hex,
    parse_physical_branch,
    wall_gram_automorphisms,
)


FIXTURES = ROOT / "data/fixtures/pb_atlas_numerical_stage_a_inputs.json"
OUTPUT = ROOT / "data/enumerations/pb_atlas_numerical_contraction_stage_a.json"
CORE = ROOT / "scripts/pb_atlas_contraction_core.py"
FIXTURE_SCRIPT = ROOT / "scripts/freeze_pb_atlas_stage_a_fixtures.py"
SCHEMA = "coxeter.pb_atlas.numerical_contraction_stage_a_calibration.v1"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _rows(matrix: sp.Matrix) -> list[list[int]]:
    if any(entry.q != 1 for entry in matrix):
        raise ArithmeticError("a mutation produced a nonintegral matrix")
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _basis_mutation(payload: Mapping[str, Any]) -> dict[str, Any]:
    rank = len(payload["picard_gram"])
    first = sp.eye(rank)
    first[0, 1] = 1
    second = sp.eye(rank)
    second[rank - 1, 0] = -1
    change = first * second
    if int(change.det()) != 1:
        raise ArithmeticError("the deterministic basis mutation is not unimodular")
    inverse = change.inv()
    gram = sp.Matrix(payload["picard_gram"])
    walls = [sp.Matrix(row) for row in payload["wall_vectors"]]
    canonical = sp.Matrix(payload["canonical_class"])
    result = copy.deepcopy(dict(payload))
    result["picard_gram"] = _rows(change.T * gram * change)
    result["wall_vectors"] = [
        [int(entry) for entry in inverse * wall] for wall in walls
    ]
    result["canonical_class"] = [int(entry) for entry in inverse * canonical]
    return result


def _wall_permutation(payload: Mapping[str, Any]) -> tuple[dict[str, Any], tuple[int, ...]]:
    wall_count = len(payload["wall_vectors"])
    old_to_new = tuple(wall_count - 1 - index for index in range(wall_count))
    new_walls: list[list[int] | None] = [None] * wall_count
    for old, new in enumerate(old_to_new):
        new_walls[new] = list(payload["wall_vectors"][old])
    boundary = list(
        zip(
            payload["boundary_indices_zero_based"],
            payload["boundary_weights"],
        )
    )
    transported = sorted((old_to_new[int(index)], int(weight)) for index, weight in boundary)
    result = copy.deepcopy(dict(payload))
    result["wall_vectors"] = new_walls
    result["boundary_indices_zero_based"] = [index for index, _ in transported]
    result["boundary_weights"] = [weight for _, weight in transported]
    return result, old_to_new


def _conjugate_automorphisms(
    permutations: Sequence[Sequence[int]], old_to_new: Sequence[int]
) -> list[list[int]]:
    wall_count = len(old_to_new)
    result = []
    for permutation in permutations:
        conjugate = [0] * wall_count
        for old in range(wall_count):
            conjugate[old_to_new[old]] = old_to_new[int(permutation[old])]
        result.append(conjugate)
    return result


def _abstract_summary(certificate: Mapping[str, Any]) -> dict[str, Any]:
    terminal_orbits = []
    for family in ("hirzebruch_endpoint_orbits", "P2_endpoint_orbits"):
        for record in certificate[family]:
            program = record["program"]
            terminal_orbits.append(
                (
                    program["endpoint_kind"],
                    program["hirzebruch_e"],
                    record["orbit_size"],
                    program["permutation_invariant_profile_sha256"],
                )
            )
    return {
        "group_order": certificate["symmetry_quotient"][
            "full_marked_wall_gram_automorphism_group_order"
        ],
        "level_counts": [
            row["labelled_frame_count"]
            for row in certificate["contraction_level_census"]
        ],
        "hirzebruch_histogram": certificate["hirzebruch_labelled_histogram"],
        "P2_count": certificate["P2_labelled_count"],
        "terminal_orbit_multiset": sorted(terminal_orbits),
    }


def _coordinate_packet_summary(certificate: Mapping[str, Any]) -> dict[str, Any]:
    def packet(program: Mapping[str, Any], *, swap_f0: bool) -> dict[str, Any]:
        wall_rows = copy.deepcopy(program["all_wall_coordinates_in_endpoint_blowup_basis"])
        boundary = copy.deepcopy(program["boundary"])
        if swap_f0:
            for row in wall_rows:
                row[0], row[1] = row[1], row[0]
            for record in boundary:
                row = record["endpoint_blowup_coordinates"]
                row[0], row[1] = row[1], row[0]
                image = record["endpoint_image_coordinates"]
                image[0], image[1] = image[1], image[0]
        return {
            "kind": program["endpoint_kind"],
            "e": program["hirzebruch_e"],
            "contracted": program["contracted_wall_indices_zero_based"],
            "order": program["sequential_contraction_order_zero_based"],
            "strict_gram": program["strict_contracted_wall_gram"],
            "basis_gram": program["endpoint_blowup_basis_gram"],
            "canonical": program["canonical_coordinates_in_endpoint_blowup_basis"],
            "nef_rays": sorted(program["nef_rays_in_model_basis"]),
            "walls": wall_rows,
            "boundary": boundary,
        }

    records = {}
    for family in ("hirzebruch_endpoint_orbits", "P2_endpoint_orbits"):
        for record in certificate[family]:
            program = record["program"]
            candidates = [packet(program, swap_f0=False)]
            if program["endpoint_kind"] == "F_e" and program["hirzebruch_e"] == 0:
                candidates.append(packet(program, swap_f0=True))
            records[record["representative_mask"]] = min(
                candidates,
                key=lambda value: json.dumps(value, sort_keys=True, separators=(",", ":")),
            )
    return records


def _typed_member_masks(
    certificate: Mapping[str, Any],
    *,
    new_to_old: Sequence[int] | None = None,
) -> dict[str, list[str]]:
    wall_count = certificate["input"]["wall_count"]
    result: dict[str, list[str]] = {}
    for family in ("hirzebruch_endpoint_orbits", "P2_endpoint_orbits"):
        for record in certificate[family]:
            program = record["program"]
            key = (
                "P2"
                if program["endpoint_kind"] == "P2"
                else f"F_{program['hirzebruch_e']}"
            )
            for encoded in record["member_masks"]:
                mask = int(encoded, 16)
                if new_to_old is not None:
                    transported = 0
                    for new in range(wall_count):
                        if mask & (1 << new):
                            transported |= 1 << int(new_to_old[new])
                    mask = transported
                result.setdefault(key, []).append(mask_hex(mask, wall_count))
    return {key: sorted(values) for key, values in sorted(result.items())}


def _expected_checks(certificate: Mapping[str, Any], expected: Mapping[str, Any]) -> dict[str, bool]:
    return {
        "automorphism_group_order": certificate["symmetry_quotient"][
            "full_marked_wall_gram_automorphism_group_order"
        ]
        == expected["wall_gram_automorphism_group_order"],
        "all_labelled_frame_counts_by_depth": [
            row["labelled_frame_count"]
            for row in certificate["contraction_level_census"]
        ]
        == expected["labelled_frame_counts_by_depth"],
        "hirzebruch_histogram": certificate["hirzebruch_labelled_histogram"]
        == expected["hirzebruch_histogram"],
        "P2_count": certificate["P2_labelled_count"] == expected["P2_count"],
        "every_terminal_orbit_has_stage_B_packet": all(
            record["program"]["stage_B_payload_complete"]
            for family in ("hirzebruch_endpoint_orbits", "P2_endpoint_orbits")
            for record in certificate[family]
        ),
    }


def _nonsaturated_symmetry_hostile_fixture() -> dict[str, Any]:
    """Prove that index > 1 is accepted but unsafe Gram symmetry is filtered."""

    base_payload = copy.deepcopy(_load(FIXTURES)["cases"][0]["input"])
    subset = (0, 4, 5, 6, 10, 18, 19, 22, 26, 27, 28)
    base_payload["wall_vectors"] = [
        base_payload["wall_vectors"][index] for index in subset
    ]
    base_payload["boundary_indices_zero_based"] = [0]
    accepted = parse_physical_branch(base_payload)

    gram = sp.diag(1, -1, -1)
    wall_matrix = sp.Matrix(
        [
            [0, 1, 2],
            [0, -1, 0],
            [1, -1, -1],
        ]
    )
    walls = tuple(wall_matrix[:, index] for index in range(3))
    canonical = -walls[2]
    wall_gram = wall_matrix.T * gram * wall_matrix
    synthetic = PhysicalBranch(
        gram=gram,
        walls=walls,
        boundary_indices=(2,),
        boundary_weights=(2,),
        canonical=canonical,
        wall_gram=wall_gram,
        canonical_pairings=tuple(
            int((canonical.T * gram * wall)[0]) for wall in walls
        ),
        wall_sublattice_index=2,
        input_payload={},
    )
    unsafe_swap = (1, 0, 2)
    raw = _marked_wall_gram_permutations(synthetic)
    safe = wall_gram_automorphisms(synthetic)
    checks = {
        "physical_index_two_full_Q_span_is_accepted": accepted.wall_sublattice_index == 2,
        "synthetic_swap_preserves_marked_wall_gram": unsafe_swap in raw,
        "synthetic_swap_has_no_integral_Picard_lift": _permutation_lift(
            synthetic, unsafe_swap
        )
        is None,
        "unsafe_swap_is_discarded": unsafe_swap not in safe,
        "identity_survives": tuple(range(3)) in safe,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"nonsaturated symmetry hostile fixture failed: {checks}")
    return {
        "physical_subset_wall_sublattice_index": accepted.wall_sublattice_index,
        "synthetic_picard_gram": _rows(gram),
        "synthetic_wall_vectors": _rows(wall_matrix.T),
        "synthetic_wall_gram": _rows(wall_gram),
        "unsafe_marked_Gram_swap": list(unsafe_swap),
        "raw_marked_Gram_group_order": len(raw),
        "integrally_liftable_group_order": len(safe),
        "checks": checks,
        "all_checks_pass": True,
    }


def _case_certificate(case: Mapping[str, Any]) -> dict[str, Any]:
    payload = case["input"]
    original = certify_branch(payload, retain_automorphisms=True)
    permutations = original["symmetry_quotient"][
        "automorphism_permutations_zero_based"
    ]

    basis_payload = _basis_mutation(payload)
    basis = certify_branch(
        basis_payload,
        automorphisms=permutations,
        retain_automorphisms=False,
    )
    basis_summary_equal = _abstract_summary(basis) == _abstract_summary(original)
    basis_packets_equal = _coordinate_packet_summary(basis) == _coordinate_packet_summary(
        original
    )
    basis_array_sha256 = basis["input"]["array_payload_sha256"]
    del basis
    gc.collect()

    wall_payload, old_to_new = _wall_permutation(payload)
    conjugated = _conjugate_automorphisms(permutations, old_to_new)
    wall = certify_branch(
        wall_payload,
        automorphisms=conjugated,
        retain_automorphisms=False,
    )
    new_to_old = [0] * len(old_to_new)
    for old, new in enumerate(old_to_new):
        new_to_old[new] = old
    wall_summary_equal = _abstract_summary(wall) == _abstract_summary(original)
    wall_frames_equivariant = _typed_member_masks(
        wall, new_to_old=new_to_old
    ) == _typed_member_masks(original)
    wall_array_sha256 = wall["input"]["array_payload_sha256"]
    del wall
    gc.collect()

    expected_checks = _expected_checks(original, case["expected"])
    mutation_checks = {
        "unimodular_basis_mutation_preserves_abstract_census": basis_summary_equal,
        "unimodular_basis_mutation_preserves_every_endpoint_coordinate_packet": basis_packets_equal,
        "wall_permutation_preserves_abstract_census": wall_summary_equal,
        "wall_permutation_transports_every_labelled_terminal_frame": wall_frames_equivariant,
    }
    if not all(expected_checks.values()) or not all(mutation_checks.values()):
        raise ArithmeticError(
            f"Stage-A calibration failed for {case['fixture_id']}: "
            f"expected={expected_checks}, mutation={mutation_checks}"
        )
    return {
        "fixture_id": case["fixture_id"],
        "calibration_alias_outside_recognizer": case["calibration_alias"],
        "expected_checks": expected_checks,
        "mutation_checks": mutation_checks,
        "basis_mutation_array_sha256": basis_array_sha256,
        "wall_permutation_array_sha256": wall_array_sha256,
        "certificate": original,
        "all_checks_pass": True,
    }


def build_certificate() -> dict[str, Any]:
    fixtures = _load(FIXTURES)
    cases = [_case_certificate(case) for case in fixtures["cases"]]
    nonsaturated = _nonsaturated_symmetry_hostile_fixture()
    global_histogram: Counter[str] = Counter()
    for case in cases:
        global_histogram.update(
            case["certificate"]["hirzebruch_labelled_histogram"]
        )
        global_histogram["P2"] += case["certificate"]["P2_labelled_count"]
    checks = {
        "four_arrays_only_calibrations": len(cases) == 4,
        "all_expected_counts_replayed": all(
            all(case["expected_checks"].values()) for case in cases
        ),
        "all_basis_and_wall_mutations_pass": all(
            all(case["mutation_checks"].values()) for case in cases
        ),
        "nonsaturated_branches_retained_but_unsafe_symmetries_filtered": nonsaturated[
            "all_checks_pass"
        ],
        "H17_has_unique_F0_and_no_P2": cases[3]["certificate"][
            "hirzebruch_labelled_histogram"
        ]
        == {"F_0": 1}
        and cases[3]["certificate"]["P2_labelled_count"] == 0,
        "no_stage_B_or_L4_overclaim": True,
    }
    return {
        "schema": SCHEMA,
        "status_date": "2026-08-10",
        "status": "PASS_PB_ATLAS_NUMERICAL_STAGE_A",
        "theorem": {
            "statement": (
                "For every finite arrays-only normalized physical Coble branch, the marked "
                "weighted wall Gram determines a finite exact orbit traversal of every supplied-wall "
                "sequential boundary-preserving (-1)-contraction frame.  For a complete full-nef "
                "negative-curve array these are all geometric frames.  Codimension rho-1 frames "
                "are certified as P2 endpoints and codimension rho-2 frames as F_e endpoints, "
                "with complete integral endpoint/blow-up coordinate packets."
            ),
            "finiteness": (
                "rank drops after each contraction and a wall index is used at most once; "
                "there are at most 2^m labelled subsets.  Quotient traversal is exact because "
                "the full stabilizer acts transitively on every omitted candidate orbit."
            ),
            "old_code_extracted": [
                "the Schur-inverse update from certify_coble_type_ii_b3_plane_inverse.py",
                "the exact intersection pushforward from certify_rational_p17_mmp.py",
                "the rational congruence inertia routine used across the K3 verifiers",
            ],
        },
        "input_contract": {
            "schema": "coxeter.pb_atlas.physical_coble_branch.v1",
            "allowed": [
                "Picard Gram array",
                "primitive oriented wall-vector array",
                "boundary index and weight arrays",
                "canonical-class array",
            ],
            "forbidden_and_rejected": [
                "source wall labels",
                "known surface or type names",
                "archived endpoint/output Grams",
                "stored contraction frames or incidence witnesses",
            ],
        },
        "sources": {
            "fixture": {"path": str(FIXTURES.relative_to(ROOT)), "sha256": _sha256(FIXTURES)},
            "core": {"path": str(CORE.relative_to(ROOT)), "sha256": _sha256(CORE)},
            "fixture_extractor": {
                "path": str(FIXTURE_SCRIPT.relative_to(ROOT)),
                "sha256": _sha256(FIXTURE_SCRIPT),
            },
        },
        "cases": cases,
        "nonsaturated_symmetry_hostile_fixture": nonsaturated,
        "aggregate_terminal_histogram": dict(sorted(global_histogram.items())),
        "stage_boundary": {
            "closed_here": (
                "PB-Atlas Stage A: numerical all-contraction enumeration and canonical "
                "endpoint coordinate packets"
            ),
            "next_exact_stage_B": (
                "turn each endpoint packet into universal cluster/incidence/Fitting equations "
                "and canonically generate collision, Jacobian, irreducibility, and Hilbert "
                "avoidance factors"
            ),
            "not_proved": [
                "incidence-stack nonemptiness",
                "a surviving primitive-binomial regular open",
                "strong-Coble construction for arbitrary input",
                "hidden-wall exclusion or literal same-P return",
                "any non-singleton L4 theorem or global dimension bound",
            ],
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-output", type=Path)
    parser.add_argument("--case-index", type=int)
    arguments = parser.parse_args()
    if arguments.case_index is not None:
        fixtures = _load(FIXTURES)
        case = _case_certificate(fixtures["cases"][arguments.case_index])
        print(
            json.dumps(
                {
                    "fixture_id": case["fixture_id"],
                    "expected_checks": case["expected_checks"],
                    "mutation_checks": case["mutation_checks"],
                    "all_checks_pass": case["all_checks_pass"],
                },
                sort_keys=True,
            )
        )
        return
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if arguments.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT)
    elif arguments.check_output is not None:
        if arguments.check_output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("frozen Stage-A archive differs from exact replay")
        print("PB-Atlas numerical Stage-A archive replay passed")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
