from __future__ import annotations

from functools import lru_cache
import hashlib
import json

from scripts.certify_enriques_type_ii_downstairs_anti_gram_census import (
    BRANDHORST_SHIMADA,
    BRANDHORST_SHIMADA_SHA256,
    MARTIN_PDF,
    MARTIN_PDF_SHA256,
    OUTPUT,
    analyze,
)


ARCHIVE_SHA256 = "04d4ac3b122d3ef0cf0b9feeb90e03f06e0692aa6a982df9d52a89df090e72e0"


@lru_cache(maxsize=1)
def _payload() -> dict:
    return analyze()


def test_archive_replays_from_all_pinned_primary_sources() -> None:
    assert hashlib.sha256(BRANDHORST_SHIMADA.read_bytes()).hexdigest() == (
        BRANDHORST_SHIMADA_SHA256
    )
    assert hashlib.sha256(MARTIN_PDF.read_bytes()).hexdigest() == MARTIN_PDF_SHA256
    payload = _payload()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert hashlib.sha256(OUTPUT.read_bytes()).hexdigest() == ARCHIVE_SHA256
    assert payload["all_checks_pass"] is True
    assert payload["scope"]["geometric_status"].startswith(
        "a wallwise forced realization of the known Type-II"
    )


def test_complete_sheet_switch_and_exact_inertia_census() -> None:
    payload = _payload()
    alphabet = payload["downstairs_alphabet"]
    assert len(alphabet["m1_edges_zero_based"]) == 15
    assert alphabet["m1_values"] == [-2, 2]
    assert alphabet["raw_assignment_count"] == 2**15 == 32768
    assert alphabet["wall_span_rank_over_Q"] == 10
    assert alphabet["wall_matrix_rank_mod_2"] == 9

    quotient = payload["sheet_switch_quotient"]
    assert quotient["connected_support_cycle_rank"] == 4
    assert quotient["predicted_orbit_count_from_cycle_signs"] == 16
    assert quotient["orbit_count"] == 16
    assert quotient["orbit_size_histogram"] == {"2048": 16}
    assert quotient["all_canonical_keys_sha256"] == (
        "999d0031c19e4e7626b18258199fbe535c5d8cac38a0e0718bd355cbed19f4f0"
    )

    census = payload["exact_inertia_census"]
    assert census["raw_assignment_histogram"] == {
        "p0_n9_z3": 4096,
        "p1_n10_z1": 12288,
        "p1_n9_z2": 2048,
        "p2_n10_z0": 14336,
    }
    assert census["sheet_switch_orbit_histogram"] == {
        "p0_n9_z3": 2,
        "p1_n10_z1": 6,
        "p1_n9_z2": 1,
        "p2_n10_z0": 7,
    }
    assert sum(
        item["raw_assignment_multiplicity"] for item in census["orbit_records"]
    ) == 32768


def test_full_half_sum_saturation_forces_the_D9_minus_two_orbit() -> None:
    theorem = _payload()["projection_theorem"]
    assert theorem["surviving_orbit_count"] == 1
    assert theorem["survivor_canonical_key_sha256"] == (
        "0b7a0c40c9cd879907ba18053a577c679171ab176131831fb396b366124ef835"
    )
    survivor = theorem["survivor_generated_anti_lattice"]
    assert survivor["rank"] == 9
    assert survivor["absolute_determinant"] == 2048
    assert survivor["smith_diagonal"] == [2] * 8 + [8]
    assert survivor["minimum_test"]["signed_vector_count"] == 0
    assert survivor["minimum_test"]["rootless_through_norm_two"] is True
    classification = survivor["scaled_lattice_classification"]
    assert classification["isometry_class"] == "D9(-2)"
    assert classification["signed_norm_four_vector_count"] == 144
    assert classification["basis_determinant_absolute"] == 1
    assert classification["dynkin_arm_lengths"] == [1, 1, 6]


def test_other_signature_orbit_is_cubic_and_has_explicit_glue_roots() -> None:
    theorem = _payload()["projection_theorem"]
    assert theorem["eliminated_orbit_canonical_key_sha256"] == (
        "38dafd994b4d29e6aa57b55febdaa876cc90e0730c48e29731c4b53e1626d832"
    )
    killed = theorem["eliminated_generated_anti_lattice"]
    assert killed["absolute_determinant"] == 512
    assert killed["smith_diagonal"] == [2] * 9
    assert killed["scaled_lattice_classification"]["isometry_class"] == "<-2>^9"
    assert killed["scaled_lattice_classification"]["basis_determinant_absolute"] == 1
    assert killed["minimum_test"]["signed_vector_count"] == 18
    assert len(killed["minus_two_roots_modulo_sign"]) == 9

    explicit = theorem["eliminated_explicit_ternary_witness_search"]
    assert explicit["witness_count_modulo_sign"] == 36
    assert explicit["minimum_support_size"] == 9
    assert explicit["minimum_support_witness_count_modulo_sign"] == 36
    assert explicit["all_witness_coefficients_sha256"] == (
        "6c7a35b3bafc46c5569742bcbe83c915681e7ef641175e0112a350654e357415"
    )
    witness = explicit["first_three_minimum_support_witnesses"][0]
    assert all(value % 2 == 0 for value in witness["downstairs_sum"])
    assert witness["four_times_square"] == -8
    assert witness["square"] == -2


def test_primary_source_calibration_and_sheet_invariance_are_scoped() -> None:
    payload = _payload()
    calibration = payload["imported_geometric_calibration"]
    assert "D9(-2)" in calibration["brandhorst_shimada"]["imported_statement"]
    assert "12B" in calibration["brandhorst_shimada"]["imported_statement"]
    assert calibration["martin"]["role"].endswith("no lift Gram imported")
    invariance = payload["sheet_switch_quotient"]["rootlessness_invariance"]
    assert invariance == {
        "coefficient_transport": "z_i -> epsilon_i*z_i",
        "norm_reason": "for D'=EDE and z'=zE, z'D'z'^t=zDz^t",
        "parity_reason": "epsilon_i=1 mod 2, hence zW=0 mod 2 is preserved",
    }
    assert "new hyperbolic Coxeter polytope" in payload["scope"]["not_claimed"][0]
