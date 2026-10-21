#!/usr/bin/env python3
"""Certify the rational Mori-boundary/K3-cover dimension dichotomy."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
CANONICAL = ROOT / "data/enumerations/canonical_orthogonal_k3_cover_calibration.json"
SIX = ROOT / "data/enumerations/six_carrier_integral_packet_rank_calibration.json"
TORIC = ROOT / "data/enumerations/toric_full_nef_coxeter_factory.json"
COBLE_II = ROOT / "data/enumerations/coble_one_boundary_full_ns_h10.json"
COBLE_I = ROOT / "data/enumerations/coble_type_i_two_boundary_full_ns_h11.json"
COBLE_VI = (
    ROOT / "data/enumerations/coble_type_vi_char3_five_boundary_full_ns_h14.json"
)
COBLE_MII = ROOT / "data/enumerations/coble_type_mii_full_ns_h17.json"
ARCHIVE = ROOT / "data/enumerations/rational_mori_boundary_k3_dichotomy.json"

PINS = {
    CANONICAL: "70b393759109e07a5ecbb9e8cc7cc87844ecb63ff266d8fd2898dcb331ba4082",
    SIX: "e36f558b0d40e817c718ccb4760687240dab89a3e2f8a12980f5cd28f816eabd",
    TORIC: "7bdd361d5d349b10110fa908d65e7508fc17cc653a097afe91be803da064a389",
    COBLE_II: "1e90fcd5f9b1ba8538937b91d2f26ea401f6952e8484a597624939edaf2a873e",
    COBLE_I: "bcfa2de1e96bf4f41282345a64f8ee6f3a75e715287bab4962d7082584b4d985",
    COBLE_VI: "8247ad85f7970d3bf1f08fe15517b70376409c17cc2784b8f14b129aaecc2ede",
    COBLE_MII: "6ed6fba936c362e02793375fa3adfea0279d312c853a8bb4bab30c23b4b8c79c",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_intersection(self_square: int) -> int:
    """K.C for a smooth rational curve C by adjunction."""
    return -2 - self_square


def adjusted_intersection(self_square: int, boundary_degree: int) -> int:
    """Twice (K+B/2).C, kept integral."""
    return 2 * canonical_intersection(self_square) + boundary_degree


def picard_ceiling(branch: str, characteristic: int) -> int:
    if branch == "weak_anticanonical":
        return 10
    if branch != "index_two_boundary":
        raise ValueError(branch)
    if characteristic == 2:
        raise ValueError("the smooth double-cover theorem excludes characteristic two")
    return 20 if characteristic == 0 else 22


def dimension_ceiling(branch: str, characteristic: int) -> int:
    return picard_ceiling(branch, characteristic) - 1


def _square_histogram(gram: sp.Matrix) -> dict[str, int]:
    squares = [int(gram[index, index]) for index in range(gram.rows)]
    return {str(square): squares.count(square) for square in sorted(set(squares))}


def _weak_alphabet_replay(
    gram_rows: list[list[int]],
    picard_rank: int,
    source_full_nef: bool,
) -> dict[str, Any]:
    gram = sp.Matrix(gram_rows)
    checks = {
        "wall_gram_is_square_symmetric": gram.rows == gram.cols
        and gram == gram.T,
        "every_mori_generator_has_square_minus_one_or_minus_two": all(
            gram[index, index] in (-1, -2) for index in range(gram.rows)
        ),
        "physical_walls_span_N1": int(gram.rank()) == picard_rank,
        "source_certifies_complete_mori_generation": source_full_nef,
    }
    return {
        "curve_count": gram.rows,
        "picard_rank": picard_rank,
        "wall_gram_rank": int(gram.rank()),
        "square_histogram": _square_histogram(gram),
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def _index_two_alphabet_replay(
    gram_rows: list[list[int]],
    boundary_components: int,
    picard_rank: int,
    source_full_nef: bool,
    boundary_scale: int = 1,
) -> dict[str, Any]:
    archived_gram = sp.Matrix(gram_rows)
    scales = [boundary_scale] * boundary_components + [1] * (
        archived_gram.rows - boundary_components
    )
    change = sp.diag(*scales)
    gram = change * archived_gram * change
    squares = [int(gram[index, index]) for index in range(gram.rows)]
    boundary_degrees = [
        int(sum(gram[index, column] for index in range(boundary_components)))
        for column in range(boundary_components, gram.cols)
    ]
    nonboundary_pairs = list(
        zip(squares[boundary_components:], boundary_degrees, strict=True)
    )
    checks = {
        "physical_gram_is_integral_square_symmetric": gram.rows == gram.cols
        and gram == gram.T
        and all(entry.q == 1 for entry in gram),
        "boundaries_are_pairwise_disjoint_square_minus_four": all(
            gram[index, index] == -4 for index in range(boundary_components)
        )
        and all(
            gram[left, right] == 0
            for left in range(boundary_components)
            for right in range(boundary_components)
            if left != right
        ),
        "other_generators_have_exact_root_or_crossing_alphabet": all(
            pair in ((-2, 0), (-1, 2)) for pair in nonboundary_pairs
        ),
        "physical_walls_span_N1": int(gram.rank()) == picard_rank,
        "source_certifies_complete_mori_generation": source_full_nef,
    }
    return {
        "curve_count": gram.rows,
        "boundary_components": boundary_components,
        "archived_boundary_scale_to_physical_curve": boundary_scale,
        "picard_rank": picard_rank,
        "physical_wall_gram_rank": int(gram.rank()),
        "physical_square_histogram": _square_histogram(gram),
        "other_generator_pair_histogram": {
            "(-2,0)": nonboundary_pairs.count((-2, 0)),
            "(-1,2)": nonboundary_pairs.count((-1, 2)),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def build_certificate() -> dict[str, Any]:
    for path, expected in PINS.items():
        actual = sha256(path)
        if actual != expected:
            raise RuntimeError(f"source changed: {path}: {actual}")

    canonical = json.loads(CANONICAL.read_text(encoding="utf-8"))
    six = json.loads(SIX.read_text(encoding="utf-8"))
    toric_source = json.loads(TORIC.read_text(encoding="utf-8"))
    coble_ii_source = json.loads(COBLE_II.read_text(encoding="utf-8"))
    coble_i_source = json.loads(COBLE_I.read_text(encoding="utf-8"))
    coble_vi_source = json.loads(COBLE_VI.read_text(encoding="utf-8"))
    coble_mii_source = json.loads(COBLE_MII.read_text(encoding="utf-8"))
    canonical_rows = {row["id"]: row for row in canonical["rows"]}
    six_rows = {row["id"]: row for row in six["rows"]}

    coble_ii_case = coble_ii_source["cases"]["II"]
    alphabet_sources = {
        "coble_type_ii_h10_30": {
            "path": COBLE_II,
            "gram": coble_ii_case["negative_wall_gram"],
            "boundary_scale": 2,
            "source_full_nef": coble_ii_source["all_checks_pass"] is True
            and coble_ii_case["all_checks_pass"] is True
            and "full_nef" in coble_ii_case["classification"]
            and coble_ii_case["construction"]["wall_order"].startswith("b=-K"),
            "source_note": (
                "the archive stores b=-K as the first nef-wall normal; "
                "doubling it recovers the physical boundary B=-2K"
            ),
        },
        "coble_type_i_h11_26": {
            "path": COBLE_I,
            "gram": coble_i_source["facets"]["negative_wall_gram"],
            "boundary_scale": 1,
            "source_full_nef": coble_i_source["all_checks_pass"] is True
            and "full_nef" in coble_i_source["classification"],
            "source_note": "the archive uses the physical divisor Gram",
        },
        "coble_type_vi_char3_h14_60": {
            "path": COBLE_VI,
            "gram": coble_vi_source["facets"]["negative_wall_gram"],
            "boundary_scale": 1,
            "source_full_nef": coble_vi_source["all_checks_pass"] is True
            and "full_nef" in coble_vi_source["classification"],
            "source_note": "the archive uses the physical divisor Gram",
        },
        "coble_type_mii_char3_h17_48": {
            "path": COBLE_MII,
            "gram": coble_mii_source["facets"]["physical_wall_gram"],
            "boundary_scale": 1,
            "source_full_nef": coble_mii_source["all_checks_pass"] is True
            and coble_mii_source["curve_exhaustion"]["full_nef_conclusion"]
            is True,
            "source_note": "the archive uses the physical divisor Gram",
        },
    }

    coble_ids = (
        "coble_type_ii_h10_30",
        "coble_type_i_h11_26",
        "coble_type_vi_char3_h14_60",
        "coble_type_mii_char3_h17_48",
    )
    coble_rows = []
    for row_id in coble_ids:
        row = canonical_rows[row_id]
        six_row = six_rows[row_id]
        boundary_components = int(row["boundary_components"])
        picard_rank = int(six_row["carrier"]["rho"])
        source = alphabet_sources[row_id]
        alphabet = _index_two_alphabet_replay(
            source["gram"],
            boundary_components,
            picard_rank,
            bool(source["source_full_nef"]),
            int(source["boundary_scale"]),
        )
        bound = picard_ceiling("index_two_boundary", int(row["characteristic"]))
        row_checks = {
            "alphabet_replay_passes": alphabet["all_checks_pass"],
            "rational_surface_rank_formula": picard_rank
            == 10 + boundary_components,
            "canonical_calibration_rank_matches_surface": int(
                row["carrier_rank_or_embedded_rank"]
            )
            == picard_rank,
            "full_nef_dimension_matches_picard_rank": int(row["n"])
            == picard_rank - 1,
            "facet_count_matches_physical_alphabet": int(row["facets"])
            == alphabet["curve_count"],
            "K3_ceiling_holds": picard_rank <= bound,
        }
        coble_rows.append(
            {
                "id": row_id,
                "characteristic": row["characteristic"],
                "dimension": row["n"],
                "boundary_components": boundary_components,
                "embedded_rank": row["carrier_rank_or_embedded_rank"],
                "rank_formula": 10 + boundary_components,
                "carrier_bound": bound,
                "alphabet_source": str(source["path"].relative_to(ROOT)),
                "alphabet_source_note": source["source_note"],
                "alphabet_replay": alphabet,
                "checks": row_checks,
                "passes": all(row_checks.values()),
            }
        )

    toric = six_rows["toric_h6_9"]
    toric_rho = int(toric["carrier"]["rho"])
    toric_alphabet = _weak_alphabet_replay(
        toric_source["h6_endpoint"]["wall_gram"],
        toric_rho,
        toric_source["all_checks_pass"] is True
        and toric_source["h6_endpoint"]["checks"][
            "fan_is_smooth_complete_with_expected_cycle"
        ]
        is True,
    )
    toric_checks = {
        "alphabet_replay_passes": toric_alphabet["all_checks_pass"],
        "full_nef_dimension_matches_picard_rank": int(toric["n"])
        == toric_rho - 1,
        "weak_anticanonical_ceiling_holds": toric_rho
        <= picard_ceiling("weak_anticanonical", 0),
    }
    toric_row = {
        "id": "toric_h6_9",
        "dimension": toric["n"],
        "picard_rank": toric_rho,
        "weak_anticanonical_ceiling": picard_ceiling("weak_anticanonical", 0),
        "alphabet_source": str(TORIC.relative_to(ROOT)),
        "alphabet_replay": toric_alphabet,
        "checks": toric_checks,
        "passes": all(toric_checks.values()),
    }

    identities = {
        "minus_one_anticanonical_degree": -canonical_intersection(-1),
        "minus_two_anticanonical_degree": -canonical_intersection(-2),
        "boundary_adjusted_intersection": adjusted_intersection(-4, -4),
        "root_adjusted_intersection": adjusted_intersection(-2, 0),
        "crossing_adjusted_intersection": adjusted_intersection(-1, 2),
    }
    checks = {
        "weak_alphabet_makes_minus_K_nonnegative": identities[
            "minus_one_anticanonical_degree"
        ] == 1
        and identities["minus_two_anticanonical_degree"] == 0,
        "index_two_alphabet_annihilates_K_plus_half_B": all(
            identities[key] == 0
            for key in (
                "boundary_adjusted_intersection",
                "root_adjusted_intersection",
                "crossing_adjusted_intersection",
            )
        ),
        "weak_picard_and_dimension_ceilings": picard_ceiling(
            "weak_anticanonical", 0
        ) == 10
        and dimension_ceiling("weak_anticanonical", 0) == 9,
        "complex_cover_ceilings": picard_ceiling("index_two_boundary", 0) == 20
        and dimension_ceiling("index_two_boundary", 0) == 19,
        "odd_characteristic_cover_ceilings": picard_ceiling(
            "index_two_boundary", 3
        ) == 22
        and dimension_ceiling("index_two_boundary", 3) == 21,
        "four_coble_alphabet_and_ceiling_controls_pass": all(
            row["passes"] for row in coble_rows
        ),
        "toric_alphabet_and_ceiling_control_passes": toric_row["passes"],
        "global_maximality_bridge_not_claimed": True,
    }

    return {
        "schema": "coxeter.rational_mori_boundary_k3_dichotomy.v1",
        "status_date": "2026-08-10",
        "scope": (
            "Carrier-relative AG dimension theorem. It does not assert that every "
            "reflection-maximal chamber has such a carrier."
        ),
        "source_pins": {str(path.relative_to(ROOT)): digest for path, digest in PINS.items()},
        "symbolic_identities": identities,
        "theorem": {
            "weak_anticanonical": {
                "curve_squares": [-1, -2],
                "conclusion": "-K is nef; rho<=10; full-nef dimension<=9",
            },
            "index_two_boundary": {
                "boundary": "pairwise disjoint square-minus-four curves",
                "other_curve_types": [[-2, 0], [-1, 2]],
                "forced_identity": "B=-2K",
                "conclusion": (
                    "in characteristic not two, the canonical double cover is "
                    "smooth, connected, and K3"
                ),
                "complex_picard_ceiling": 20,
                "odd_characteristic_picard_ceiling": 22,
            },
            "trace_pushback": "[k:Q]*(n+1) <= rho(S)",
            "mori_span_reason": (
                "curve classes span N1 by definition, so any complete finite "
                "generating set for the closed Mori cone spans N1"
            ),
            "numerical_to_linear_reason": (
                "2K+B is an integral numerically trivial divisor and Pic(S) is "
                "torsion-free with Pic^0(S)=0 for rational S"
            ),
        },
        "controls": {"toric": toric_row, "coble": coble_rows},
        "open_gate": (
            "Construct from every raw maximal Coxeter datum either a carrier with "
            "one chosen geometric exceptional lattice having finite numerical "
            "operator orbit, or a complete rational Mori-alphabet carrier."
        ),
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def certify(payload: dict[str, Any]) -> None:
    if not payload["all_checks_pass"]:
        failed = [key for key, value in payload["checks"].items() if not value]
        raise ArithmeticError(f"failed checks: {failed}")
    if "every raw maximal" in payload["scope"]:
        raise ArithmeticError("carrier-relative theorem was promoted to a global theorem")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    certify(payload)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if ARCHIVE.read_text(encoding="utf-8") != rendered:
            raise SystemExit("archive mismatch")
        print("PASS: rational Mori-boundary/K3 dichotomy archive matches")
    elif args.write:
        ARCHIVE.write_text(rendered, encoding="utf-8")
        print(ARCHIVE)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
