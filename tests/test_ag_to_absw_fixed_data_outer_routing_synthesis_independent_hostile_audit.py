import json
import subprocess
import sys

from scripts.audit_ag_to_absw_fixed_data_outer_routing_synthesis import (
    DOCUMENT_PINS,
    EVIDENCE_PINS,
    OUTPUT,
    OWNER_PACKAGE_PINS,
    ROOT,
    build_certificate,
    build_payload,
    digest,
    encoded_certificate,
    local_link_audit,
    sha256_file,
    STRIPPED_TEX_RE,
)


AUDIT_REPORT = (
    ROOT
    / "results/ag_to_absw_fixed_data_outer_routing_synthesis_independent_hostile_audit.md"
)
AUDIT_REPORT_SHA256 = (
    "60f0b50fa78ba2d301104991aaf03985e4747c7e8f89b374e217608063e8d5c9"
)


def test_all_fifteen_frozen_inputs_are_pinned() -> None:
    assert len(DOCUMENT_PINS) == 6
    assert len(OWNER_PACKAGE_PINS) == 3
    assert len(EVIDENCE_PINS) == 6
    assert len({*DOCUMENT_PINS, *OWNER_PACKAGE_PINS, *EVIDENCE_PINS}) == 15


def test_each_synthesis_document_is_individually_audited() -> None:
    documents = build_payload()["document_audits"]
    assert set(documents) == set(DOCUMENT_PINS)
    for audit in documents.values():
        assert audit["all_checks_pass"] is True
        assert all(audit["checks"].values())
        assert audit["checks"]["all_six_frozen_evidence_reports_linked"] is True
        assert audit["all_local_links_resolve"] is True
        assert audit["broken_local_links"] == []
        assert audit["links_escaping_repository"] == []
        assert audit["stripped_tex_hits"] == []
        assert audit["defective_x3_patterns"] == []
        assert audit["defective_field_scope_patterns"] == []


def test_finite_pair_scope_is_conditional_and_classification_importing() -> None:
    scope = build_payload()["scope"]
    assert "reflection maximality" in scope["finite_pair_theorem"]
    assert scope["finite_pair_step_imports_nikulin_1982_classification"] is True
    assert scope["finite_pair_step_is_ag_only"] is False
    assert scope["global_absw_proved"] is False
    assert (
        scope[
            "field_finiteness_requires_exact_admissible_full_trace_order_module_k3_return"
        ]
        is True
    )


def test_x3_quantifiers_are_not_reversed() -> None:
    scope = build_payload()["scope"]
    assert scope["x3_is_one_counterexample_to_all_finite_volume_universality"] is True
    assert scope["x3_obstructs_every_finite_volume_input"] is False
    assert scope["x3_proves_a_maximal_only_no_go"] is False
    assert scope["x3_absolute_reflection_maximality_proved"] is False
    assert scope["x3_mixed_divisor_complex_k3_positive"] is True


def test_h21_and_projected_host_scopes_are_literal() -> None:
    scope = build_payload()["scope"]
    assert scope["h21_absolute_reflection_maximality"] is True
    assert scope["h21_current_complex_stack_empty"] is True
    assert scope["h21_every_ag_carrier_empty"] is False
    assert scope["h21_characteristic_two_mixed_curve_cremona_positive"] is True
    assert scope["projected_ambient_host_finiteness_is_true_in_general"] is False
    assert (
        scope["projected_ambient_host_finiteness_needed_for_abstract_count"]
        is False
    )
    assert scope["delta_246_projected_return_proved_for_one_packet"] is True
    assert scope["every_exact_projected_complex_k3_return_has_n_at_most_19"] is True


def test_independent_arithmetic_replay_hits_both_endpoint_controls() -> None:
    replay = build_payload()["hostile_arithmetic_replay"]
    assert replay["x3_forced_square_class"] == 3
    assert replay["x3_three_is_rational_square"] is False
    assert replay["h21_complex_k3_rank_gate_passes"] is False
    assert replay["fano_diagonal_point_determinant"] == -2
    assert replay["fano_configuration_possible_over_characteristic_zero"] is False
    assert replay["fano_determinant_mod_two"] == 0
    assert replay["native_reflection_coefficient_similarity_invariant"] is True


def test_owner_archive_is_cross_checked_without_importing_owner_producer() -> None:
    checks = build_payload()["owner_archive_cross_checks"]
    assert all(checks.values())


def test_independent_report_links_and_tex_transport_are_clean() -> None:
    assert sha256_file(AUDIT_REPORT) == AUDIT_REPORT_SHA256
    text = AUDIT_REPORT.read_text(encoding="utf-8")
    links = local_link_audit(AUDIT_REPORT, text)
    assert links["all_local_links_resolve"] is True
    assert links["broken_local_links"] == []
    assert links["links_escaping_repository"] == []
    assert list(STRIPPED_TEX_RE.finditer(text)) == []


def test_archive_is_an_exact_stable_replay_and_does_not_self_pin() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    live = build_certificate()
    assert archived == live
    assert archived["payload_sha256"] == digest(archived["payload"])
    assert OUTPUT.read_text(encoding="utf-8") == encoded_certificate()
    assert "audit_archive" not in archived["payload"]["input_sha256"]
    assert str(OUTPUT.relative_to(ROOT)) not in archived["payload"]["input_sha256"]


def test_cli_check_output() -> None:
    subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.audit_ag_to_absw_fixed_data_outer_routing_synthesis",
            "--check-output",
        ],
        cwd=ROOT,
        check=True,
    )
