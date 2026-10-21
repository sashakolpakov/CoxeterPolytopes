"""Focused regressions for the independent projected-K3 finiteness audit."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import re
import subprocess
import sys

from scripts.audit_projected_k3_maximal_abstract_chamber_finiteness import (
    AUDIT_ARCHIVE,
    AUDIT_REPORT,
    AUDIT_REPORT_SHA256,
    DELTA_246_GRAM,
    OWNER_ARCHIVE,
    OWNER_REPORT,
    OWNER_SCRIPT,
    OWNER_SHA256,
    OWNER_TEST,
    PRIMARY_NIKULIN_SHA256,
    build_certificate,
    cartan_matrix,
    digest,
    encoded_certificate,
    finite_index_replay,
    normalized_squared_gram,
    rank_budget,
    scale_matrix,
    sha256_file,
)


ROOT = Path(__file__).resolve().parents[1]
CERTIFIER = (
    ROOT / "scripts/audit_projected_k3_maximal_abstract_chamber_finiteness.py"
)
FROZEN_OWNER = {
    OWNER_REPORT: OWNER_SHA256["report"],
    OWNER_SCRIPT: OWNER_SHA256["script"],
    OWNER_ARCHIVE: OWNER_SHA256["archive"],
    OWNER_TEST: OWNER_SHA256["test"],
}
FROZEN_AUDIT = {
    AUDIT_REPORT: "edc663c0e7f35880123149d81d832812d5f28837c837df16eaa0f1daf725916a",
    CERTIFIER: "c8be6882b15bdfca48929c80176fb863d0ee9493f7990515d2ee0b3ba53c579a",
    AUDIT_ARCHIVE: "4f7ee617fcdb8e034bfef01c91935a476062e13221d27557dbd413268b0faa39",
}


def test_frozen_owner_and_independent_packages() -> None:
    for path, expected in FROZEN_OWNER.items():
        assert sha256_file(path) == expected
    for path, expected in FROZEN_AUDIT.items():
        assert sha256_file(path) == expected
    assert sha256_file(AUDIT_REPORT) == AUDIT_REPORT_SHA256


def test_primary_source_pin_and_exact_locator() -> None:
    payload = build_certificate()
    source = payload["primary_source"]
    assert source["locator"] == "Theorem 5.2.1, pp. 119-122"
    assert source["scope_matches_owner"] is True
    assert source["torsion_free_projective_modules_allowed"] is True
    assert source["pdf_sha256"] == PRIMARY_NIKULIN_SHA256

    # This local primary artifact is present in the research workspace but is
    # not required in a clean checkout; when present, verify its bytes.
    primary = Path("/private/tmp/nikulin_1981_eng.pdf")
    if primary.exists():
        assert sha256_file(primary) == PRIMARY_NIKULIN_SHA256


def test_complex_k3_rank_budget_is_exact() -> None:
    assert rank_budget() == [
        {
            "field_degree": 1,
            "maximum_native_rank": 20,
            "maximum_hyperbolic_dimension": 19,
        },
        {
            "field_degree": 2,
            "maximum_native_rank": 10,
            "maximum_hyperbolic_dimension": 9,
        },
        {
            "field_degree": 3,
            "maximum_native_rank": 6,
            "maximum_hyperbolic_dimension": 5,
        },
        {
            "field_degree": 4,
            "maximum_native_rank": 5,
            "maximum_hyperbolic_dimension": 4,
        },
        {
            "field_degree": 5,
            "maximum_native_rank": 4,
            "maximum_hyperbolic_dimension": 3,
        },
        {
            "field_degree": 6,
            "maximum_native_rank": 3,
            "maximum_hyperbolic_dimension": 2,
        },
    ]


def test_finite_index_and_positive_cone_bridge() -> None:
    replay = finite_index_replay(256)
    assert replay == {
        "maximum_tested_positive_index_W": 256,
        "factorizations_tested": 2932,
        "maximum_resulting_full_index_R": 512,
    }
    bridge = build_certificate()["reflective_lattice_bridge"]
    assert bridge["chain"] == "W(P) <= R(L) <= O+(L) <= O(L)"
    assert bridge["positive_cone_index"] == "[O(L):O+(L)] <= 2"
    assert bridge["nikulin_index"] == "[O(L):R(L)] < infinity"


def test_similarity_replay_preserves_reflections_and_normalized_gram() -> None:
    expected_cartan = (
        (Fraction(2), Fraction(-1), Fraction(0)),
        (Fraction(-2), Fraction(2), Fraction(-1)),
        (Fraction(0), Fraction(-3), Fraction(2)),
    )
    expected_normalized = (
        (Fraction(4), Fraction(2), Fraction(0)),
        (Fraction(2), Fraction(4), Fraction(3)),
        (Fraction(0), Fraction(3), Fraction(4)),
    )
    assert cartan_matrix(DELTA_246_GRAM) == expected_cartan
    assert normalized_squared_gram(DELTA_246_GRAM) == expected_normalized
    for scale in range(1, 257):
        gram = scale_matrix(DELTA_246_GRAM, scale)
        assert cartan_matrix(gram) == expected_cartan
        assert normalized_squared_gram(gram) == expected_normalized


def test_physical_star_bound_is_local_and_unused_globally() -> None:
    payload = build_certificate()
    physical = payload["physical_star"]
    assert physical == {
        "strict_real_upper_bound": 14,
        "nonnegative_integral_upper_bound": 13,
        "used_in_global_chamber_finiteness": False,
    }
    for numerator in range(201):
        a = Fraction(numerator, 100)
        assert 6 * a + 2 <= 14


def test_maximality_and_complete_chamber_orbit_are_both_required() -> None:
    section = build_certificate()["maximality_and_chambers"]
    assert section["R_is_discrete_reflection_overgroup"] is True
    assert section["R_is_arithmetic_when_arithmetic_only_maximality_is_used"] is True
    assert section["maximality_forces_W_equals_R"] is True
    assert section["complete_mirror_arrangement_locally_finite"] is True
    assert section["generic_compact_path_crosses_finitely_many_mirrors"] is True
    assert section["successive_crossing_reflections_give_chamber_transitivity"] is True
    assert section[
        "unmarked_normalized_chambers_per_oriented_similarity_host_at_most"
    ] == 1


def test_scale_glue_and_nonmaximal_controls_have_opposite_scope() -> None:
    controls = build_certificate()["hostile_controls"]
    scaling = controls["projected_scale_glue_family"]
    assert scaling["ambient_hosts_infinite"] is True
    assert scaling["native_similarity_host_fixed"] is True
    assert scaling["abstract_chamber_fixed"] is True
    redoubling = controls["redoubling_family"]
    assert redoubling["abstract_chambers_infinite"] is True
    assert redoubling["reflection_maximal"] is False


def test_exact_scope_and_nonclaims_remain_explicit() -> None:
    payload = build_certificate()
    scope = payload["exact_scope_gate"]
    assert scope["full_native_trace_module"] is True
    assert scope["admissible_native_form"] is True
    assert scope["exact_complete_facet_operator_return"] is True
    assert scope["merely_semialgebraic_NS_slice_suffices"] is False
    assert scope["universal_projected_K3_existence_proved"] is False
    conclusion = payload["conclusion"]
    assert conclusion["finite_unmarked_normalized_maximal_chambers"] is True
    assert conclusion["finite_ambient_NS_hosts"] is False
    assert conclusion["full_ABSW_without_exhaustive_carrier_stack"] is False
    assert conclusion["AG_only_no_Nikulin_import"] is False


def test_connected_local_star_does_not_replace_global_input() -> None:
    connected = build_certificate()["hostile_controls"]["connected_frame"]
    assert connected["bounded_spanning_local_star_may_be_disconnected"] is True
    assert connected["right_angled_pentagon_control"] is True
    assert connected["local_star_replaces_Nikulin_fixed_data_theorem"] is False


def test_reports_have_no_known_stripped_tex_tokens() -> None:
    combined = OWNER_REPORT.read_text(encoding="utf-8") + AUDIT_REPORT.read_text(
        encoding="utf-8"
    )
    assert re.search(r"(?<!\\)qquad", combined) is None
    assert "This is not yet an AG-only derivation of ABSW" in combined
    assert "universal complex-K3 return" in combined
    assert "[O(L):O^+(L)]" in combined


def test_archive_is_exact_replay() -> None:
    archived = json.loads(AUDIT_ARCHIVE.read_text(encoding="utf-8"))
    assert archived["payload"] == build_certificate()
    assert archived["payload_sha256"] == digest(archived["payload"])
    assert AUDIT_ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()


def test_cli_check() -> None:
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check-output"],
        cwd=ROOT,
        check=True,
    )
