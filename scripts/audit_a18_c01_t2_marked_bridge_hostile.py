#!/usr/bin/env python3
"""Hostile exact marked bridge from the Gram-only T2 compiler to C01.

This audit is deliberately separate from the Gram-only lattice precompiler.
It starts from that frozen hostile archive, then reads the archived (18,4,0)
leading-root model, the C01 factory marking, and the existing strong-promotion
certificate.  It accepts the datum bridge only after constructing an explicit
integral unimodular isometry which transports the form, involution, ordered
24-root marking, and all fourteen source-wall orbits exactly.

The result closes exact marked recovery R for this one datum and permits the
already-certified C01 cone-completeness C and Vinberg V conclusions to be
imported.  It does not prove forward universality U for a variable H9/14
class, and it does not add a new polytope.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form


ROOT = Path(__file__).resolve().parents[1]
HOSTILE = (
    ROOT
    / "data"
    / "enumerations"
    / "a18_c01_t2_lattice_precompiler_hostile_audit.json"
)
FACTORY = ROOT / "data" / "enumerations" / "alexeev_k3_h9_14_factory.json"
CONSTRUCTION = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_remaining_high_fixed_section_construction.json"
)
PROMOTION = (
    ROOT / "data" / "enumerations" / "a18_c01_c07_k3_nef_promotions.json"
)
PROMOTION_SCRIPT = ROOT / "scripts" / "certify_a18_c01_c07_k3_nef_promotions.py"
OUTPUT = (
    ROOT / "data" / "enumerations" / "a18_c01_t2_marked_bridge_hostile_audit.json"
)

SCHEMA = "coxeter.a18_c01.t2_marked_bridge.hostile_audit.v1"
STATUS_DATE = "2026-08-08"
CASE_ID = "A18_4_0_C01"
SELECTED_SHEET_MASK = 831
EXPECTED_MATCHING_SHEET_MASKS = [192, 440, 583, 831]
EXPECTED_ACTIVE_ORBIT_INDICES = [1, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 16, 17, 20]
EXPECTED_ORDERED_SOURCE_WALLS = [
    1,
    9,
    3,
    8,
    4,
    7,
    5,
    6,
    17,
    11,
    16,
    12,
    15,
    13,
    14,
    18,
    24,
    19,
    25,
    32,
    22,
    31,
    23,
    30,
]
EXPECTED_MARKED_ISOMETRY_SHA256 = (
    "c1c7e49891a9467985e4df99d44ac41dff8bd21afa39c166dfc168a322ba1e14"
)
EXPECTED_SOURCE_SHA256 = {
    "hostile_ic0_archive": (
        "deaca02f71afe2e7fa6e9418ba628dd0772c5a970a0bb09adc390bb40dc34b08"
    ),
    "factory": "b1ee6fd769527ae6dab01a2bb58099e7c7d1776d0417b6a8d710316697e8a775",
    "construction": (
        "c339b56a6d0b68c0939dab85ba6390c7429a30e56eeafec822ac0d5373ec0e91"
    ),
    "strong_promotion_archive": (
        "2ab7acc379e48ac55944e92b8d614c0992d7c482af1300c670f2bf8e899ba1b1"
    ),
    "strong_promotion_script": (
        "cea0682c9acd10748fa29cab9b35ca65b0a1963704e0e2022230a89a84fa159f"
    ),
}

Matrix = sp.Matrix


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _matrix(rows: Sequence[Sequence[Any]]) -> Matrix:
    return Matrix([[sp.Rational(str(value)) for value in row] for row in rows])


def _matrix_payload(matrix: Matrix) -> list[list[str]]:
    return [[str(value) for value in matrix.row(index)] for index in range(matrix.rows)]


def _matrix_digest(matrix: Matrix) -> str:
    raw = json.dumps(_matrix_payload(matrix), separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def _integer_matrix(matrix: Matrix) -> Matrix:
    if not all(value.q == 1 for value in matrix):
        raise ArithmeticError("expected an integral matrix")
    return matrix.applyfunc(int)


def _same_column_lattice(left: Matrix, right: Matrix) -> bool:
    if left.shape != right.shape:
        return False
    return hermite_normal_form(_integer_matrix(left)) == hermite_normal_form(
        _integer_matrix(right)
    )


def _basis_transition(target_basis: Matrix, source_basis: Matrix) -> Matrix:
    """Return C with target_basis*C=source_basis, using independent rows."""

    if target_basis.rows != source_basis.rows or target_basis.cols != source_basis.cols:
        raise ArithmeticError("basis dimensions differ")
    pivot_rows = list(target_basis.T.rref()[1])
    if len(pivot_rows) != target_basis.cols:
        raise ArithmeticError("target columns are dependent")
    transition = target_basis[pivot_rows, :].inv() * source_basis[pivot_rows, :]
    if target_basis * transition != source_basis:
        raise ArithmeticError("basis transition did not reconstruct")
    return transition


def _ordered_source_ids(
    active_records: Sequence[dict[str, Any]], mask: int
) -> tuple[list[int], list[list[int]]]:
    ordered: list[int] = []
    blocks: list[list[int]] = []
    pair_index = 0
    for record in active_records:
        block = [int(value) for value in record["source_wall_indices"]]
        if len(block) == 2:
            if (mask >> pair_index) & 1:
                block.reverse()
            pair_index += 1
        elif len(block) != 1:
            raise ArithmeticError("a C01 active orbit is not one- or two-sheeted")
        blocks.append(block)
        ordered.extend(block)
    if pair_index != 10:
        raise ArithmeticError("the C01 active marking no longer has ten pairs")
    return ordered, blocks


def audit() -> dict[str, Any]:
    paths = {
        "hostile_ic0_archive": HOSTILE,
        "factory": FACTORY,
        "construction": CONSTRUCTION,
        "strong_promotion_archive": PROMOTION,
        "strong_promotion_script": PROMOTION_SCRIPT,
    }
    actual_hashes = {label: _sha256(path) for label, path in paths.items()}
    if actual_hashes != EXPECTED_SOURCE_SHA256:
        raise ArithmeticError("a pinned marked-bridge source changed")

    hostile = json.loads(HOSTILE.read_text(encoding="utf-8"))
    factory = json.loads(FACTORY.read_text(encoding="utf-8"))
    construction = json.loads(CONSTRUCTION.read_text(encoding="utf-8"))
    promotion = json.loads(PROMOTION.read_text(encoding="utf-8"))
    if not (
        hostile["all_checks_pass"]
        and construction["all_checks_pass"]
        and promotion["all_checks_pass"]
        and promotion["source_file_sha256"]["factory"] == actual_hashes["factory"]
        and promotion["source_file_sha256"]["construction"]
        == actual_hashes["construction"]
    ):
        raise ArithmeticError("an upstream hostile/factory/promotion gate is not green")

    candidate = hostile["unique_candidate"]
    hostile_rank18 = candidate["rank18"]
    hostile_fold = candidate["fold_back"]
    hostile_gram = _matrix(hostile_rank18["gram"])
    hostile_action = _matrix(hostile_rank18["involution"])
    hostile_roots = _matrix(hostile_fold["source_roots"])
    hostile_folded = _matrix(hostile_fold["folded_roots"])
    hostile_orbits = hostile_fold["source_root_orbits_one_based"]

    inventory = next(
        row
        for row in construction["inventory"]["examples"]
        if row["key"] == [18, 4, 0]
    )
    archive_gram = _matrix(inventory["ambient"]["gram"])
    wall_rows = _matrix(inventory["walls"]["coordinates_in_leading_root_basis"])
    factory_case = next(
        row for row in factory["constructions"] if row["id"] == CASE_ID
    )
    section = factory_case["section"]
    archive_row_action = _matrix(
        section["case"]["representative_nonidentity_primal_actions"][0]
    )
    archive_column_action = archive_row_action.T
    active_records = [
        row
        for row in section["restricted_mirrors"]["records"]
        if row["hyperbolic_status"] == "active_facet"
    ]
    promotion_case = next(row for row in promotion["cases"] if row["id"] == CASE_ID)

    wall_lookup = {
        tuple(wall_rows.row(index)): index + 1 for index in range(wall_rows.rows)
    }
    if len(wall_lookup) != 48:
        raise ArithmeticError("the archived source chamber has duplicate walls")
    action_permutation = [
        wall_lookup[tuple(wall_rows.row(index) * archive_row_action)]
        for index in range(48)
    ]
    full_chamber_action_exact = (
        sorted(action_permutation) == list(range(1, 49))
        and all(action_permutation[action_permutation[index] - 1] == index + 1 for index in range(48))
        and archive_row_action**2 == sp.eye(18)
        and archive_row_action * archive_gram * archive_row_action.T == archive_gram
    )

    if [int(row["orbit_index"]) for row in active_records] != EXPECTED_ACTIVE_ORBIT_INDICES:
        raise ArithmeticError("the ordered list of fourteen active C01 orbits changed")
    if len(active_records) != len(hostile_orbits):
        raise ArithmeticError("hostile and factory facet markings have different lengths")

    matching_masks: list[int] = []
    hostile_root_gram = hostile_roots.T * hostile_gram * hostile_roots
    for mask in range(1 << 10):
        source_ids, _ = _ordered_source_ids(active_records, mask)
        archive_roots = wall_rows[[index - 1 for index in source_ids], :].T
        if archive_roots.T * archive_gram * archive_roots == hostile_root_gram:
            matching_masks.append(mask)
    if matching_masks != EXPECTED_MATCHING_SHEET_MASKS:
        raise ArithmeticError("the exhaustive sheet-order isometry search changed")

    ordered_source_ids, ordered_source_blocks = _ordered_source_ids(
        active_records, SELECTED_SHEET_MASK
    )
    if ordered_source_ids != EXPECTED_ORDERED_SOURCE_WALLS:
        raise ArithmeticError("the selected ordered source-wall marking changed")
    archive_roots = wall_rows[[index - 1 for index in ordered_source_ids], :].T

    pivot_columns = list(hostile_roots.rref()[1])
    if pivot_columns != list(range(18)):
        raise ArithmeticError("the frozen root-coordinate pivot set changed")
    marked_isometry = (
        archive_roots[:, pivot_columns] * hostile_roots[:, pivot_columns].inv()
    )
    marked_isometry_inverse = marked_isometry.inv()
    if _matrix_digest(marked_isometry) != EXPECTED_MARKED_ISOMETRY_SHA256:
        raise ArithmeticError("the frozen mask-831 marked isometry changed")

    archive_folded_columns: list[Matrix] = []
    orbit_rows: list[dict[str, Any]] = []
    flat_position = 0
    orbit_transport_exact = True
    for facet_index, (record, hostile_orbit, ordered_block) in enumerate(
        zip(active_records, hostile_orbits, ordered_source_blocks), start=1
    ):
        raw_block = [int(value) for value in record["source_wall_indices"]]
        positions = list(range(flat_position + 1, flat_position + len(raw_block) + 1))
        flat_position += len(raw_block)
        if positions != [int(value) for value in hostile_orbit]:
            orbit_transport_exact = False
        source_rows = [wall_rows.row(index - 1) for index in ordered_block]
        folded_row = sum(source_rows, sp.zeros(1, 18))
        archived_sum = _matrix([record["ambient_orbit_sum_primal_coordinates"]])
        action_images = [
            wall_lookup[tuple(wall_rows.row(index - 1) * archive_row_action)]
            for index in ordered_block
        ]
        expected_images = ordered_block if len(ordered_block) == 1 else list(reversed(ordered_block))
        this_orbit_exact = (
            folded_row == archived_sum
            and action_images == expected_images
            and int(record["presentation_scale_from_primitive"]) == 1
            and int(record["square"]) == (-2 if len(ordered_block) == 1 else -4)
        )
        orbit_transport_exact &= this_orbit_exact
        archive_folded_columns.append(folded_row.T)
        orbit_rows.append(
            {
                "hostile_facet_one_based": facet_index,
                "factory_orbit_index": int(record["orbit_index"]),
                "hostile_source_root_positions_one_based": positions,
                "factory_source_wall_indices_unordered": raw_block,
                "factory_source_wall_indices_ordered": ordered_block,
                "archive_action_images_in_order": action_images,
                "presentation_scale_from_primitive": int(
                    record["presentation_scale_from_primitive"]
                ),
                "folded_square": int(record["square"]),
                "orbit_and_folded_sum_transport_exact": this_orbit_exact,
            }
        )
    archive_folded = Matrix.hstack(*archive_folded_columns)

    hostile_fixed_basis = _matrix(hostile_rank18["fixed_basis_columns"])
    hostile_anti_basis = _matrix(hostile_rank18["anti_basis_columns"])
    archive_fixed_basis = _matrix(
        section["saturated_fixed_lattice"]["basis_primal_rows"]
    ).T
    archive_anti_basis = _matrix(
        promotion_case["ambient_NS"]["coinvariant_lattice"]["basis_in_N"]
    ).T
    transported_fixed_basis = marked_isometry * hostile_fixed_basis
    transported_anti_basis = marked_isometry * hostile_anti_basis
    fixed_transition = _basis_transition(archive_fixed_basis, transported_fixed_basis)
    anti_transition = _basis_transition(archive_anti_basis, transported_anti_basis)
    eigensublattice_transport_exact = (
        _same_column_lattice(transported_fixed_basis, archive_fixed_basis)
        and _same_column_lattice(transported_anti_basis, archive_anti_basis)
        and abs(fixed_transition.det()) == 1
        and abs(anti_transition.det()) == 1
    )

    exact_form_transport = marked_isometry.T * archive_gram * marked_isometry == hostile_gram
    exact_action_transport = (
        marked_isometry * hostile_action
        == archive_column_action * marked_isometry
    )
    exact_ordered_root_transport = marked_isometry * hostile_roots == archive_roots
    exact_folded_transport = (
        marked_isometry * hostile_folded == archive_folded
        and archive_folded.T * archive_gram * archive_folded
        == _matrix(section["restricted_mirrors"]["facet_gram"])
    )
    selected_root_norms = list((archive_roots.T * archive_gram * archive_roots).diagonal())
    full_wall_norms = list((wall_rows * archive_gram * wall_rows.T).diagonal())

    cone_claim = promotion_case["nef_and_quotient_cone"]
    period_claim = promotion_case["period_and_torelli"]
    exact_cone = factory_case["exact_cone"]
    action_digest_matches_promotion = (
        _matrix_digest(archive_row_action)
        == promotion_case["ambient_NS"]["action_sha256"]
    )
    existing_C_importable = (
        action_digest_matches_promotion
        and period_claim["picard_lattice_is_exactly_N"] is True
        and period_claim["source_48_wall_chamber_preserved"] is True
        and period_claim["invariant_kahler_chamber_exists"] is True
        and cone_claim["all_restricted_long_walls_redundant"] is True
        and cone_claim["invariant_nef_section_equality"]
        == "factory cone = P_r intersect Fix(g) = Nef(X) intersect Fix(g)"
    )
    existing_V_importable = (
        all(exact_cone["checks"].values())
        and exact_cone["decision"]["finite_volume"] is True
        and exact_cone["decision"]["compact"] is False
        and cone_claim["finite_volume_noncompact"] is True
        and cone_claim["archived_dimension_facets"] == [9, 14]
    )

    checks = {
        "all_pinned_sources_exact": actual_hashes == EXPECTED_SOURCE_SHA256,
        "upstream_hostile_ic0_green": hostile["all_checks_pass"] is True,
        "upstream_strong_promotion_green": promotion["all_checks_pass"] is True,
        "full_48_wall_chamber_action_exact": full_chamber_action_exact,
        "ten_pair_sheet_search_exhaustive": matching_masks
        == EXPECTED_MATCHING_SHEET_MASKS,
        "selected_mask_831_order_exact": ordered_source_ids
        == EXPECTED_ORDERED_SOURCE_WALLS,
        "marked_isometry_integral": all(value.q == 1 for value in marked_isometry),
        "marked_isometry_unimodular": marked_isometry.det() == 1,
        "marked_isometry_inverse_integral": all(
            value.q == 1 for value in marked_isometry_inverse
        ),
        "rank18_form_congruence_exact": exact_form_transport,
        "involution_action_conjugacy_exact": exact_action_transport,
        "ordered_24_root_transport_exact": exact_ordered_root_transport,
        "fourteen_orbit_source_records_exact": orbit_transport_exact,
        "fourteen_folded_sums_and_gram_exact": exact_folded_transport,
        "selected_24_source_roots_all_minus_two": selected_root_norms
        == [-2] * 24,
        "full_source_wall_norm_histogram_exact": Counter(map(int, full_wall_norms))
        == {-4: 24, -2: 24},
        "fixed_and_anti_primitive_lattices_transport_exact": (
            eigensublattice_transport_exact
        ),
        "existing_C01_cone_completeness_importable": existing_C_importable,
        "existing_C01_terminal_V_importable": existing_V_importable,
    }
    if not all(checks.values()):
        failed = [label for label, passed in checks.items() if not passed]
        raise ArithmeticError(f"marked bridge failed: {failed}")

    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "case_id": CASE_ID,
        "source_file_sha256": actual_hashes,
        "scope": {
            "claim": "exact datum-specific marked IC0/R bridge to archived C01",
            "new_polytope_claim": False,
            "class_wide_IC2_or_IC4_claim": False,
            "variable_H9_14_forward_universality_U_claim": False,
        },
        "sheet_order_search": {
            "two_sheet_orbit_count": 10,
            "searched_mask_count": 1024,
            "exact_matching_masks": matching_masks,
            "selected_mask": SELECTED_SHEET_MASK,
            "selected_mask_bits_low_to_high": [
                (SELECTED_SHEET_MASK >> index) & 1 for index in range(10)
            ],
            "selected_ordered_factory_source_wall_indices_one_based": (
                ordered_source_ids
            ),
        },
        "marked_isometry": {
            "coordinate_convention": (
                "columns; x_archive=P*x_hostile; archived action is transposed "
                "from the factory row convention"
            ),
            "matrix_P": _matrix_payload(marked_isometry),
            "matrix_P_sha256": _matrix_digest(marked_isometry),
            "determinant": int(marked_isometry.det()),
            "inverse_integral": True,
            "form_equality": "P^T*G_archive*P=G_hostile",
            "action_equality": "P*T_hostile=A_archive^T*P",
            "ordered_root_equality": "P*R_hostile=R_archive_ordered",
            "all_three_equalities_exact": (
                exact_form_transport
                and exact_action_transport
                and exact_ordered_root_transport
            ),
        },
        "source_chamber": {
            "source_wall_count": 48,
            "action_permutation_one_based": action_permutation,
            "action_preserves_full_48_wall_set": full_chamber_action_exact,
            "source_wall_norm_histogram": {"-4": 24, "-2": 24},
            "active_source_root_count": 24,
            "active_source_roots_all_square_minus_two": True,
        },
        "active_marked_orbits": orbit_rows,
        "eigensublattice_transport": {
            "fixed_rank": 10,
            "anti_rank": 8,
            "fixed_basis_change": _matrix_payload(fixed_transition),
            "fixed_basis_change_determinant": int(fixed_transition.det()),
            "anti_basis_change": _matrix_payload(anti_transition),
            "anti_basis_change_determinant": int(anti_transition.det()),
            "primitive_column_lattices_equal": eigensublattice_transport_exact,
        },
        "imported_existing_C01_conclusions": {
            "strong_promotion_action_sha256": promotion_case["ambient_NS"][
                "action_sha256"
            ],
            "action_digest_matches_exact_bridge_target": (
                action_digest_matches_promotion
            ),
            "C_complete_invariant_nef_section": existing_C_importable,
            "C_statement": cone_claim["invariant_nef_section_equality"],
            "V_finite_volume_noncompact": existing_V_importable,
            "V_finite_and_ideal_vertices": cone_claim[
                "finite_and_ideal_vertices"
            ],
            "surface_realization_statement": period_claim[
                "global_torelli_conclusion"
            ],
        },
        "obligation_ledger": {
            "IC0": "closed for the single pinned 14-facet Gram by the separate hostile precompiler",
            "R": "closed for this labelled C01 datum by the explicit integral marked isometry P",
            "C": "imported for this datum from the pinned existing C01 strong promotion through P",
            "V": "imported for this datum from the pinned exact C01 factory cone through P",
            "U": "open for any nontrivial variable H9/14 class",
            "class_wide_surface_softening_theorem": False,
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    payload = audit()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if (
            not arguments.output.is_file()
            or arguments.output.read_text(encoding="utf-8") != rendered
        ):
            raise SystemExit("hostile C01 marked-bridge archive differs from exact replay")
        print("archive_check=PASS")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print("selected_sheet_mask=831")
    print("marked_isometry_det=1")
    print("exact_form_action_ordered_roots=true")
    print("full_48_wall_action_and_14_orbits=true")
    print("datum_R_C_V=closed")
    print("variable_class_U=open")
    print(f"payload_sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
