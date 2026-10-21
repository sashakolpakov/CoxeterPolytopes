"""Certify the exact Gamma16 blowup obstruction to complement bounds.

The theorem in the accompanying report is symbolic: a Gamma-stable negative
exceptional lattice receives no nonzero equivariant projection from the full
finite-covolume trace representation.  The computation here certifies the
independent sharpness example.  Tumarkin's reflection-maximal Gamma16 K3
lattice is enlarged by an arbitrary diagonal exceptional block.  Every wall
reflection extends by the identity, fixes the blowup canonical vector, and
leaves the original wall Gram unchanged, while the orthogonal-complement rank
is the freely chosen number of blowups.  The replay also constructs the
divisor action of tGamma_pi o Z_i o Gamma_pi - sum(E_j x E_j) and checks that
it is exactly the claimed block reflection, including the exceptional sign.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.verify_k3_high_rank_nef_chambers import build_gram, exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATH = ROOT / "data" / "known_examples" / "k3_high_rank_nef_chambers.json"
MAXIMALITY_PATH = (
    ROOT
    / "data"
    / "enumerations"
    / "wang_h16_maximal_reflection_no_correction.json"
)
CANONICAL_CALIBRATION_PATH = (
    ROOT
    / "data"
    / "enumerations"
    / "canonical_orthogonal_k3_cover_calibration.json"
)
SIX_CARRIER_PATH = (
    ROOT
    / "data"
    / "enumerations"
    / "six_carrier_integral_packet_rank_calibration.json"
)
COBLE_BLOWDOWN_PATH = (
    ROOT
    / "data"
    / "enumerations"
    / "coble_polytope_first_blowdown_incidence_bridge.json"
)
ARCHIVE_PATH = (
    ROOT
    / "data"
    / "enumerations"
    / "equivariant_minimalization_trace_descent.json"
)

SAMPLE_BLOWUP_COUNTS = (0, 1, 2, 5, 11)

MINIMAL_ZERO_DEFECT_IDS = (
    "wehler_h2_3",
    "bugaenko_h8_11",
    "coble_type_ii_h10_30",
    "coble_type_i_h11_26",
    "coble_type_vi_char3_h14_60",
    "coble_type_mii_char3_h17_48",
)

# The rank of a full total-transform exceptional frame for the displayed
# rational model.  These rows are defect controls only: the archived weak/full
# nef carrier does not certify the canonical-preserving operator hypotheses of
# the descent theorem.
RATIONAL_DEFECT_SPECS = {
    "toric_h6_9": {"minimal_model": "P2", "minimal_rank": 1, "blowups": 6},
    "coble_type_ii_h10_30": {
        "minimal_model": "P2",
        "minimal_rank": 1,
        "blowups": 10,
    },
    "coble_type_i_h11_26": {
        "minimal_model": "P2",
        "minimal_rank": 1,
        "blowups": 11,
    },
    "coble_type_vi_char3_h14_60": {
        "minimal_model": "P2",
        "minimal_rank": 1,
        "blowups": 14,
    },
    "coble_type_mii_char3_h17_48": {
        "minimal_model": "P1xP1",
        "minimal_rank": 2,
        "blowups": 16,
    },
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_rank17_case() -> dict[str, Any]:
    payload = json.loads(SOURCE_PATH.read_text())
    return next(case for case in payload["cases"] if int(case["rho"]) == 17)


def _diag_block(left: sp.Matrix, right: sp.Matrix) -> sp.Matrix:
    return sp.diag(left, right)


def _reflection(root: sp.Matrix, gram: sp.Matrix) -> sp.Matrix:
    norm = (root.T * gram * root)[0]
    return sp.eye(gram.rows) - sp.Rational(2, norm) * root * root.T * gram


def _sample_row(
    blowups: int,
    basis_gram: sp.Matrix,
    wall_coordinates: sp.Matrix,
) -> dict[str, Any]:
    exceptional = -sp.eye(blowups)
    enlarged_gram = _diag_block(basis_gram, exceptional)
    canonical = sp.Matrix([0] * basis_gram.rows + [1] * blowups)

    # With covariant correspondence conventions, Gamma_pi acts by pushforward
    # and its transpose by pullback.  Thus
    #   tGamma_pi o Z_i o Gamma_pi
    # is diag(s_i, 0) on divisors.  The cycle -E_j x E_j acts by
    # -E_j * (E_j^T G), which is +1 on E_j because E_j^2=-1.
    exceptional_cycle_correction = sp.zeros(enlarged_gram.rows)
    for offset in range(blowups):
        exceptional_class = sp.zeros(enlarged_gram.rows, 1)
        exceptional_class[basis_gram.rows + offset, 0] = 1
        exceptional_cycle_correction -= (
            exceptional_class * exceptional_class.T * enlarged_gram
        )

    extended_reflections: list[sp.Matrix] = []
    reflection_checks: list[bool] = []
    correspondence_formula_checks: list[bool] = []
    for column in range(wall_coordinates.cols):
        root = wall_coordinates[:, column]
        reflection = _reflection(root, basis_gram)
        extended = _diag_block(reflection, sp.eye(blowups))
        pull_reflect_push = _diag_block(reflection, sp.zeros(blowups))
        correspondence_action = pull_reflect_push + exceptional_cycle_correction
        extended_reflections.append(extended)
        correspondence_formula_checks.append(correspondence_action == extended)
        reflection_checks.append(
            all(entry.q == 1 for entry in reflection)
            and reflection * reflection == sp.eye(reflection.rows)
            and reflection.T * basis_gram * reflection == basis_gram
            and extended.T * enlarged_gram * extended == enlarged_gram
            and extended * canonical == canonical
        )

    trace_inclusion = sp.Matrix.vstack(
        sp.eye(basis_gram.rows), sp.zeros(blowups, basis_gram.rows)
    )
    complement_inclusion = sp.Matrix.vstack(
        sp.zeros(basis_gram.rows, blowups), sp.eye(blowups)
    )
    trace_pairing = trace_inclusion.T * enlarged_gram * trace_inclusion
    complement_pairing = (
        complement_inclusion.T * enlarged_gram * complement_inclusion
    )
    orthogonality = trace_inclusion.T * enlarged_gram * complement_inclusion

    checks = {
        "enlarged_signature_is_1_16_plus_blowups": exact_inertia(enlarged_gram)
        == (1, 16 + blowups, 0),
        "trace_block_unchanged": trace_pairing == basis_gram,
        "exceptional_block_is_minus_identity": complement_pairing == exceptional,
        "trace_and_exceptional_blocks_orthogonal": orthogonality
        == sp.zeros(basis_gram.rows, blowups),
        "canonical_square_is_minus_blowups": (canonical.T * enlarged_gram * canonical)[
            0
        ]
        == -blowups,
        "all_19_reflections_integral_isometric_involutions_fixing_K": all(
            reflection_checks
        ),
        "all_19_Gamma16_correspondence_formulas_equal_block_reflections": all(
            correspondence_formula_checks
        ),
    }
    return {
        "blowups": blowups,
        "picard_rank": basis_gram.rows + blowups,
        "trace_rank": basis_gram.rows,
        "orthogonal_complement_rank": blowups,
        "canonical_square": -blowups,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def _defect_calibration() -> dict[str, Any]:
    canonical_payload = json.loads(CANONICAL_CALIBRATION_PATH.read_text())
    canonical_by_id = {row["id"]: row for row in canonical_payload["rows"]}
    minimal_rows: list[dict[str, Any]] = []
    for row_id in MINIMAL_ZERO_DEFECT_IDS:
        source = canonical_by_id[row_id]
        if "K3" not in source["minimal_carrier"]:
            raise ArithmeticError(f"{row_id} is no longer a K3 carrier")
        if source["canonical_orthogonality"] != "automatic because K_X=0":
            raise ArithmeticError(f"{row_id} lost canonical orthogonality")
        minimal_rows.append(
            {
                "id": row_id,
                "carrier": source["minimal_carrier"],
                "trace_rank": int(source["trace_rank"]),
                "characteristic": int(source["characteristic"]),
                "exceptional_frame_rank": 0,
                "D_pi_rank": 0,
                "D_pi_zero": True,
                "reason": "pi is the identity on the minimal K3 carrier, so C_pi has zero columns",
                "canonical_operator_scope": True,
                "complex_rank_20_scope": int(source["characteristic"]) == 0,
                "control_strength": "vacuous identity-contraction zero control",
            }
        )

    six_payload = json.loads(SIX_CARRIER_PATH.read_text())
    six_by_id = {row["id"]: row for row in six_payload["rows"]}
    blowdown_payload = json.loads(COBLE_BLOWDOWN_PATH.read_text())
    vi_blowdowns = blowdown_payload["type_VI_plane_compiler"][
        "boundary_preserving_P2_blowdowns"
    ]["blowdowns"]
    mii_blowdown = blowdown_payload["type_MII_raw_compiler"][
        "P_compiled_blowdown"
    ]
    blowdown_source_checks = {
        "Type_VI_has_five_determinant_one_P2_frames_of_rank_14": len(vi_blowdowns)
        == 5
        and all(
            len(row["active_wall_indices_zero_based"]) == 14
            and int(row["exceptional_basis_determinant_absolute"]) == 1
            for row in vi_blowdowns
        ),
        "Type_MII_has_P1xP1_frame_of_rank_16": mii_blowdown["endpoint"]
        == "P1xP1"
        and len(mii_blowdown["exceptional_wall_indices_zero_based"]) == 16
        and mii_blowdown["exceptionals_pairwise_disjoint"] is True,
    }
    if not all(blowdown_source_checks.values()):
        raise ArithmeticError("the frozen Type-VI/MII blowdown frames changed")

    rational_rows: list[dict[str, Any]] = []
    for row_id, spec in RATIONAL_DEFECT_SPECS.items():
        source = six_by_id[row_id]
        rho = int(source["carrier"]["rho"])
        trace_rank = int(source["trace"]["rank_d_times_n_plus_1"])
        blowups = int(spec["blowups"])
        if trace_rank != rho:
            raise ArithmeticError(f"{row_id} is no longer a full-NS trace row")
        if rho - int(spec["minimal_rank"]) != blowups:
            raise ArithmeticError(f"{row_id} blowup rank no longer matches rho")
        if spec["minimal_model"] == "P2":
            gram = sp.diag(1, *([-1] * blowups))
        elif spec["minimal_model"] == "P1xP1":
            gram = _diag_block(sp.Matrix([[0, 1], [1, 0]]), -sp.eye(blowups))
        else:
            raise ArithmeticError(f"unknown minimal model for {row_id}")
        trace_inclusion = sp.eye(rho)
        exceptional_frame = sp.Matrix.vstack(
            sp.zeros(int(spec["minimal_rank"]), blowups), sp.eye(blowups)
        )
        defect = trace_inclusion.T * gram * exceptional_frame
        defect_rank = int(defect.rank())
        if defect_rank != blowups:
            raise ArithmeticError(f"{row_id} defect matrix lost full column rank")
        rational_rows.append(
            {
                "id": row_id,
                "carrier": source["carrier"]["surface"],
                "minimal_model": spec["minimal_model"],
                "trace_rank": trace_rank,
                "exceptional_frame_rank": blowups,
                "D_pi_shape": [defect.rows, defect.cols],
                "D_pi_matrix": [
                    [int(defect[row, column]) for column in range(defect.cols)]
                    for row in range(defect.rows)
                ],
                "D_pi_rank": defect_rank,
                "D_pi_zero": defect == sp.zeros(defect.rows, defect.cols),
                "reason": "T_Q=Num(Y)_Q, so J is invertible; G_Y is nondegenerate and a full exceptional-frame matrix C_pi has independent columns",
                "canonical_operator_scope": False,
                "calibration_kind": "exact standard-basis matrix replay from the frozen full-NS and blowdown model",
                "scope_warning": "defect-only negative control; canonical-preserving simple-reflection extensions are not certified on this rational model",
            }
        )

    checks = {
        "six_minimal_K3_rows_have_zero_column_defect": all(
            row["exceptional_frame_rank"] == 0
            and row["D_pi_zero"]
            and row["canonical_operator_scope"]
            for row in minimal_rows
        ),
        "five_rational_full_NS_rows_have_full_rank_nonzero_defect": all(
            row["D_pi_rank"] == row["exceptional_frame_rank"] > 0
            and not row["D_pi_zero"]
            for row in rational_rows
        ),
        "rational_rows_are_not_promoted_to_canonical_operator_controls": all(
            not row["canonical_operator_scope"] for row in rational_rows
        ),
        "Type_VI_MII_blowdown_sources_pass": all(blowdown_source_checks.values()),
    }
    return {
        "formula": "D_pi = J^t G_Y C_pi",
        "minimal_K3_and_K3_cover_zero_controls": minimal_rows,
        "rational_full_NS_nonzero_controls": rational_rows,
        "blowdown_source_checks": blowdown_source_checks,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
def build_certificate() -> dict[str, Any]:
    case = _load_rank17_case()
    wall_gram = build_gram(case)
    basis_indices = tuple(int(index) - 1 for index in case["basis_indices"])
    basis_gram = wall_gram.extract(basis_indices, basis_indices)
    wall_coordinates = basis_gram.inv() * wall_gram.extract(
        basis_indices, range(wall_gram.cols)
    )

    if wall_gram.rank() != 17 or exact_inertia(wall_gram) != (1, 16, 2):
        raise ArithmeticError("the frozen Gamma16 wall Gram changed")
    if basis_gram.det() != 2:
        raise ArithmeticError("the frozen Gamma16 wall basis changed")
    if not all(entry.q == 1 for entry in wall_coordinates):
        raise ArithmeticError("a Gamma16 wall is no longer integral")
    if wall_coordinates.T * basis_gram * wall_coordinates != wall_gram:
        raise ArithmeticError("wall coordinates no longer reconstruct the Gram")

    rows = [
        _sample_row(blowups, basis_gram, wall_coordinates)
        for blowups in SAMPLE_BLOWUP_COUNTS
    ]
    defect_calibration = _defect_calibration()
    maximality = json.loads(MAXIMALITY_PATH.read_text())
    maximality_checks = {
        "source_archive_passes": maximality["all_checks_pass"] is True,
        "source_is_reflection_maximal": maximality["maximality"]["conclusion"]
        == "Gamma16 is maximal among discrete/arithmetic reflection overgroups",
        "source_dimension": int(maximality["polytope"]["dimension"]) == 16,
        "source_facets": int(maximality["polytope"]["facets"]) == 19,
        "source_trace_lattice_rank": basis_gram.rows == 17,
    }

    checks = {
        "Gamma16_wall_rank_17": wall_gram.rank() == 17,
        "Gamma16_wall_inertia_1_16_2": exact_inertia(wall_gram) == (1, 16, 2),
        "Gamma16_basis_determinant_2": basis_gram.det() == 2,
        "all_wall_coordinates_integral": all(
            entry.q == 1 for entry in wall_coordinates
        ),
        "all_sample_blowup_rows_pass": all(row["all_checks_pass"] for row in rows),
        "defect_calibration_passes": defect_calibration["all_checks_pass"],
        "maximality_source_passes": all(maximality_checks.values()),
    }

    return {
        "schema_version": 1,
        "artifact_id": "equivariant_minimalization_trace_descent",
        "sources": {
            "Gamma16_wall_gram": str(SOURCE_PATH.relative_to(ROOT)),
            "Gamma16_wall_gram_sha256": _sha256(SOURCE_PATH),
            "Gamma16_maximality": str(MAXIMALITY_PATH.relative_to(ROOT)),
            "Gamma16_maximality_sha256": _sha256(MAXIMALITY_PATH),
            "canonical_K3_cover_calibration": str(
                CANONICAL_CALIBRATION_PATH.relative_to(ROOT)
            ),
            "canonical_K3_cover_calibration_sha256": _sha256(
                CANONICAL_CALIBRATION_PATH
            ),
            "six_carrier_calibration": str(SIX_CARRIER_PATH.relative_to(ROOT)),
            "six_carrier_calibration_sha256": _sha256(SIX_CARRIER_PATH),
            "Coble_VI_MII_blowdown_calibration": str(
                COBLE_BLOWDOWN_PATH.relative_to(ROOT)
            ),
            "Coble_VI_MII_blowdown_calibration_sha256": _sha256(
                COBLE_BLOWDOWN_PATH
            ),
        },
        "theorem_boundary": {
            "symbolic_family": "Bl_r(X), r >= 0, for an exact-Picard Gamma16 K3 X",
            "trace_lattice": "pi^*(U + E8(-1) + E7(-1))",
            "exceptional_complement": "<-1>^r",
            "descent_invariant": "a finite numerical operator orbit of one geometric exceptional lattice",
            "extended_simple_reflections": "s_alpha direct_sum identity_r",
            "correspondence_formula": "tGamma_pi o Z_i o Gamma_pi - sum_j(E_j x E_j)",
            "canonical_vector": "E_1 + ... + E_r",
            "conclusion": "reflection maximality does not bound nonminimal orthogonal-complement rank in the intrinsic correspondence category",
            "not_claimed": "the blown-up full nef cone has only the original 19 walls",
        },
        "maximality_checks": maximality_checks,
        "D_pi_calibration": defect_calibration,
        "sample_rows": rows,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-archive", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if arguments.write_archive:
        ARCHIVE_PATH.write_text(rendered)
    if arguments.check_archive:
        if not ARCHIVE_PATH.exists() or ARCHIVE_PATH.read_text() != rendered:
            raise SystemExit("recomputed certificate differs from archive")
    print(rendered, end="")


if __name__ == "__main__":
    main()
