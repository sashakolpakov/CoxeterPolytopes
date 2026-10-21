"""Independent hostile tests for the compact maximal (2,5,6) K3 no-go."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.audit_k3_compact_maximal_256_phys_empty_independent_hostile import (
    AUDIT_ARCHIVE,
    AUDIT_REPORT,
    Kphi,
    OWNER_HASHES,
    PHI,
    build_certificate,
    compact_reflection_area_maximality,
    derive_codifferent_shell,
    digest,
    encoded_certificate,
    owner_hash_audit,
    reconstruct_raw_packet,
    singerman_maximality_crosscheck,
    valuation_obstruction,
)


AUDITOR = (
    ROOT / "scripts/audit_k3_compact_maximal_256_phys_empty_independent_hostile.py"
)


def test_owner_bytes_are_frozen_and_auditor_is_independent() -> None:
    assert owner_hash_audit() == OWNER_HASHES
    source = AUDITOR.read_text(encoding="utf-8")
    assert "certify_k3_compact_maximal_256_phys_empty import" not in source
    assert "from scripts.certify_k3" not in source


def test_independent_phi_arithmetic_and_raw_reconstruction() -> None:
    one = Kphi(Fraction(1), Fraction(0))
    assert PHI * PHI == PHI + one
    assert PHI.trace() == 1
    assert PHI.norm() == -1
    assert PHI.sign() == 1
    assert PHI.sign(conjugate=True) == -1

    raw = reconstruct_raw_packet()
    assert raw["basis"] == "1,phi with phi^2=phi+1"
    assert raw["determinant"] == ["0", "6"]
    assert raw["determinant_text"] == "6*phi"
    assert raw["LDL_identity_signs"] == [-1, -1, 1]
    assert raw["LDL_conjugate_signs"] == [-1, -1, -1]
    assert raw["cartan_entries_in_Ok"] is True
    assert raw["compact_angle_sum"] == "13/15"
    assert raw["trace_rank"] == 6
    assert raw["cycle_order_equals_maximal_order"] is True


def test_direct_area_proof_exhausts_every_proper_reflection_overgroup() -> None:
    proof = compact_reflection_area_maximality()
    assert proof["P_area_over_pi"] == "2/15"
    assert proof["proper_overgroup_chamber_area_upper_bound_over_pi"] == "1/15"
    assert proof["compact_coxeter_polygon_s_ge_4_area_floor_over_pi"] == "1/6"
    assert proof["overgroup_chamber_must_be_triangle"] is True
    assert proof["compact_triangle_area_floor_over_pi"] == "1/42"
    assert proof["derived_index_upper_bound"] == 5
    assert proof["solutions_by_index"] == {
        "2": [[2, 3, 10]],
        "3": [],
        "4": [],
        "5": [],
    }
    assert proof["index_two_means_adjacent_chamber_double"] is True
    assert [case["result"] for case in proof["doubling_cases"]] == [
        "triangle (3,3,5)",
        "triangle angles (pi/10,pi/10,2pi/3)",
        "quadrilateral",
    ]
    assert not any(
        case["equals_target_2_5_6"] for case in proof["doubling_cases"]
    )
    assert proof["proper_reflection_overgroup_exists"] is False
    assert proof["absolute_full_reflection_maximality"] is True
    assert proof["classification_dependency"] is False


def test_published_triangle_classification_is_only_a_crosscheck() -> None:
    crosscheck = singerman_maximality_crosscheck()
    assert crosscheck["target"] == [2, 5, 6]
    assert crosscheck["target_in_complete_nonmaximal_patterns"] is False
    assert crosscheck["absolute_full_reflection_maximality"] is True
    assert crosscheck["published_dependency_not_reproved"] is True


def test_codifferent_is_derived_and_trace_shell_has_no_cutoff() -> None:
    shell = derive_codifferent_shell()
    assert shell["trace_gram"] == [[2, 1], [1, 3]]
    assert shell["inverse_trace_gram"] == [
        [[3, 5], [-1, 5]],
        [[-1, 5], [2, 5]],
    ]
    assert shell["trace_minus_two_parameter"] == "A=-2"
    assert shell["negativity_equivalence"] == "|B+1|<sqrt(5)"
    assert shell["integer_boundary_proof"] == "2^2<5<3^2"
    assert shell["total_negative_B_values"] == [-3, -2, -1, 0, 1]
    assert [row["c=B+1"] for row in shell["rows"]] == [-2, -1, 0, 1, 2]
    assert [row["trace"] for row in shell["rows"]] == ["-2"] * 5
    assert [row["norm"] for row in shell["rows"]] == [
        "1/5",
        "4/5",
        "1",
        "4/5",
        "1/5",
    ]
    assert shell["excluded_c_signs"] == {"-3": [-1, 1], "3": [1, -1]}
    assert shell["cardinality"] == 5


def test_inert_prime_parity_kills_all_twenty_five_pairs() -> None:
    obstruction = valuation_obstruction(derive_codifferent_shell())
    assert obstruction["minimal_polynomial_mod_3_values"] == [2, 2, 1]
    assert obstruction["three_inert"] is True
    assert obstruction["energy_shell_all_p_units"] is True
    assert obstruction["m6_norm_ratio_square_class"] == "[3]"
    assert obstruction["form_similarity_scale_cancels"] is True
    assert obstruction["line_representative_changes_contribute_only_squares"] is True
    assert obstruction["ordered_pair_count"] == 25
    assert obstruction["square_count"] == 0
    assert all(row["v_p(u_left/u_right)"] == 0 for row in obstruction["pairs"])
    assert all(
        row["v_p(u_left/(3*u_right))"] == -1
        for row in obstruction["pairs"]
    )
    assert not any(row["square_possible"] for row in obstruction["pairs"])


def test_all_declared_freedoms_and_scope_are_closed_without_overclaim() -> None:
    payload = build_certificate()
    assert payload["maximality"]["classification_dependency"] is False
    assert payload["published_classification_crosscheck"][
        "published_dependency_not_reproved"
    ] is True

    freedoms = payload["freedom_audit"]
    assert freedoms["common_form_scale_lambda"] == "cancels in u1/u3"
    assert freedoms["facet_representative_scaling"] == (
        "changes h1/h3 by a square"
    )
    assert freedoms["coefficient_order"] == (
        "cycle order O0=Ok, so no order lane remains"
    )
    assert freedoms["all_freedoms_close"] is True

    embedding = payload["embedding_scope"]
    assert embedding["trace_rank"] == 6
    assert embedding["K3_complement_rank"] == 16
    assert embedding["discriminant_length_at_most"] == 6
    assert embedding["stable_embedding_after_PHYS"] is True
    assert embedding["failure_is_before_EMB"] is True

    scope = payload["claim_scope"]
    assert scope["full_trace"] is True
    assert scope["full_Picard_NS_equals_trace_lattice"] is True
    assert scope["same_vector_physical_roots"] is True
    assert scope["projected_root_larger_Picard_excluded"] is False
    assert scope["alternative_AG_carrier_excluded"] is False
    assert scope["global_Property4_disjunction_refuted"] is False
    assert scope["universal_full_Picard_K3_nonemptiness_refuted"] is True
    assert payload["verdict"] == "PASS"


def test_audit_report_states_source_caution_proof_and_scope() -> None:
    text = AUDIT_REPORT.read_text(encoding="utf-8")
    for needle in (
        "PASS, with one source-presentation caution",
        "Classification-free absolute reflection maximality",
        "\\operatorname{area}(P)",
        "N& (a,b,c)",
        "(2,3,10)",
        "None is \\((2,5,6)\\)",
        "The codifferent shell is complete",
        "Every scaling and representative freedom cancels or squares",
        "v_{\\mathfrak p}(3)=1",
        "full-trace, full-Picard, same-vector",
        "does **not** exclude a larger Picard lattice",
        "did not edit any owner byte",
    ):
        assert needle in text
    assert text.count("\\[") == text.count("\\]")
    assert "../scripts/audit_k3_compact" in text
    assert "../data/enumerations/k3_compact" in text
    assert "../tests/test_k3_compact" in text


def test_archive_is_an_exact_replay() -> None:
    archived = json.loads(AUDIT_ARCHIVE.read_text(encoding="utf-8"))
    assert archived["payload"] == build_certificate()
    assert archived["payload_sha256"] == digest(archived["payload"])
    assert AUDIT_ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()


def test_cli_check() -> None:
    subprocess.run(
        [sys.executable, str(AUDITOR), "--check-output"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
