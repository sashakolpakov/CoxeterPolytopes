from __future__ import annotations

import json
from pathlib import Path

from scripts.audit_rational_weak_log_finite_iff_compiler import (
    ARCHIVE,
    EXPECTED_OWNER_SHA256,
    build,
    cusp_certificate,
    delta246_descent,
    delta246_integral_physical_assignments,
    fujino_multiplier,
    hirzebruch_endpoint_pairs,
    owner_hashes,
    sha256_json,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = (
    ROOT
    / "results/rational_weak_log_finite_iff_compiler_independent_hostile_audit.md"
)


def frozen() -> dict[str, object]:
    return json.loads(ARCHIVE.read_text())


def test_live_audit_matches_frozen_archive() -> None:
    assert build() == frozen()
    payload = frozen().copy()
    digest = payload.pop("content_sha256")
    assert digest == sha256_json(payload)


def test_audited_owner_bytes_are_pinned_after_endpoint_repair() -> None:
    assert owner_hashes() == EXPECTED_OWNER_SHA256
    assert EXPECTED_OWNER_SHA256["report"].startswith("48cf57f67dea")


def test_hirzebruch_endpoint_omission_is_closed_exactly() -> None:
    assert hirzebruch_endpoint_pairs() == [
        [2, 0],
        [2, 2],
        [3, 1],
        [4, 0],
    ]
    row = frozen()["mmp_audit"]  # type: ignore[index]
    assert row["finite_hirzebruch_endpoints"] == [0, 2, 3, 4]
    assert row["automorphism_quotient_needed_for_finite_presentation"] is False


def test_affine_cusp_relation_is_independent_and_positive() -> None:
    row = cusp_certificate()
    assert row["coefficients"] == [2, 4, 6, 5, 4, 3, 2, 1, 3]
    assert row["target"] == row["reconstruction"]
    assert row["checks"] == {
        "target_is_nonzero_null": True,
        "root_coordinate_rank": 9,
        "root_gram_rank": 8,
        "all_roots_have_square_minus_two": True,
        "off_diagonal_intersections_nonnegative": True,
        "all_roots_incident_to_target": True,
        "all_relation_coefficients_positive": True,
        "positive_relation_reconstructs_target": True,
    }


def test_fujino_gate_uses_L_equals_2D_and_multiplier_96() -> None:
    assert fujino_multiplier(2, 0) == 96
    row = frozen()["nef_audit"]["big_branch"]  # type: ignore[index]
    assert row["L_D"] == "2D=-2K-B is nef Cartier on an actual return"
    assert row["fixed_system"] == "|96L_D|=|192D|"
    assert row["global_generation_is_necessary"] is True
    assert row["global_generation_is_sufficient_for_D_nef"] is True


def test_246_rational_and_physical_host_obstructions_are_exhaustive() -> None:
    assert delta246_descent()["rational_form_is_anisotropic"] is True
    assert delta246_integral_physical_assignments() == []
    row = frozen()["host_audit"]  # type: ignore[index]
    assert row["real_signature_alone_sufficient"] is False
    assert row["delta246_integral_physical_assignments"] == []


def test_constructibility_claim_does_not_claim_implementation() -> None:
    row = frozen()["constructibility_audit"]  # type: ignore[index]
    assert row["basepoint_free_locus_constructible"] is True
    assert row["algebraically_closed_nonemptiness_decidable_by_elimination"] is True
    assert row["generic_degree_192_elimination_implemented"] is False


def test_empty_and_property4_scopes_are_not_conflated() -> None:
    row = frozen()["return_and_scope_audit"]  # type: ignore[index]
    assert row["empty_means_no_carrier_in_this_weak_log_species"] is True
    assert row["empty_means_no_carrier_in_any_possible_species"] is False
    assert row["universal_property4_nonemptiness_proved"] is False
    text = REPORT.read_text()
    for required in (
        "PASS after one necessary endpoint-finiteness",
        "initially omitted a real finiteness issue",
        "|96L_D|=|192D|",
        "generic degree-`192D`",
        "`EMPTY` has the narrow, correct meaning",
    ):
        assert required in text
