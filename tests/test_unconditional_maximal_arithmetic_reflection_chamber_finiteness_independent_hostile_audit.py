"""Regressions for the independent unconditional chamber-finiteness audit."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
import re
import subprocess
import sys

from scripts.audit_unconditional_maximal_arithmetic_reflection_chamber_finiteness import (
    AUDIT_ARCHIVE,
    AUDIT_REPORT,
    AUDIT_REPORT_SHA256,
    CONJUGATE_HEIGHT,
    FIELD_DEGREE_MAX,
    OWNER_ARCHIVE,
    OWNER_REPORT,
    OWNER_SCRIPT,
    OWNER_SHA256,
    OWNER_TEST,
    RANK_MAX,
    SCOPE_NO_GO_REPORT,
    SCOPE_NO_GO_SHA256,
    SOURCE_PINS,
    build_certificate,
    coefficient_bounds,
    cyclic_integrality_countercontrol,
    decompose_closed_walk,
    digest,
    directed_edges,
    encoded_certificate,
    sha256_file,
    simple_cycle_slot_bound,
)


ROOT = Path(__file__).resolve().parents[1]
CERTIFIER = (
    ROOT
    / "scripts/audit_unconditional_maximal_arithmetic_reflection_chamber_finiteness.py"
)
OWNER_CERTIFIER = (
    ROOT
    / "scripts/certify_unconditional_maximal_arithmetic_reflection_chamber_finiteness.py"
)
FROZEN_OWNER = {
    OWNER_REPORT: OWNER_SHA256["report"],
    OWNER_SCRIPT: OWNER_SHA256["script"],
    OWNER_ARCHIVE: OWNER_SHA256["archive"],
    OWNER_TEST: OWNER_SHA256["test"],
}
FROZEN_AUDIT = {
    AUDIT_REPORT: "1f1592ef207580b765d56bdca9461429456ef9fa46eb180b0ae65d69dc68ef3c",
    CERTIFIER: "e6fc5f7543093b6167e6d46612b64c7247007cc16b06f081e14260a283193400",
    AUDIT_ARCHIVE: "87ec814ae05d29b04ac068296c4f995144264fed353fa02f607eb22841cff734",
    SCOPE_NO_GO_REPORT: SCOPE_NO_GO_SHA256,
}


def test_frozen_owner_audit_and_scope_bytes() -> None:
    for path, expected in FROZEN_OWNER.items():
        assert sha256_file(path) == expected
    for path, expected in FROZEN_AUDIT.items():
        assert sha256_file(path) == expected
    assert sha256_file(AUDIT_REPORT) == AUDIT_REPORT_SHA256


def test_optional_primary_source_pins_match_when_present() -> None:
    for record in SOURCE_PINS.values():
        path = Path(record["path"])
        if path.exists():
            assert sha256_file(path) == record["sha256"]


def test_terminal_scope_is_classical_unconditional_not_ag_only() -> None:
    payload = build_certificate()
    terminal = payload["terminal_scope"]
    assert terminal["maximum_dimension"] == 29
    assert terminal["maximum_ground_field_degree"] == 25
    assert terminal["finite_ground_field_set"] is True
    assert terminal["finite_native_pair_set"] is True
    assert terminal[
        "finite_unmarked_normalized_chambers_for_reflection_maximal_inputs"
    ] is True
    assert terminal["AG_carrier_existence_hypothesis"] is False
    assert payload["verdict"] == "PASS_UNCONDITIONAL_CLASSICAL_NOT_AG_ONLY"


def test_noncompact_refinement_is_scoped_as_a_literature_import() -> None:
    refinement = build_certificate()["noncompact_refinement"]
    assert refinement["ground_field"] == "Q"
    assert refinement["dimensions"] == "n<=19 or n=21"
    assert refinement["cusp_to_k_isotropy"] is True
    assert refinement["nondistinguished_definiteness_forces_degree_one"] is True
    assert refinement["Esselmann_integral_reflective_lattice_theorem_imported"] is True
    assert refinement["owner_machine_certificate_claims_to_replay_this_import"] is False


def test_cyclic_integrality_does_not_require_entrywise_integrality() -> None:
    control = cyclic_integrality_countercontrol()
    assert control["symmetrizer"] == ["36", "9", "1"]
    assert control["nonintegral_individual_entries"] == ["1/2", "1/3", "1/6"]
    assert control["integral_two_cycles"] == ["1", "1", "1"]
    assert control["forward_three_cycle"] == "1"
    assert control["reverse_three_cycle"] == "1"
    assert control["cycle_square"] == control["product_of_two_cycles"] == "1"


def test_closed_walks_factor_into_simple_cycles_with_same_edges() -> None:
    walk = (0, 1, 2, 1, 3, 0)
    pieces = decompose_closed_walk(walk)
    assert pieces == [(1, 2, 1), (0, 1, 3, 0)]
    source = Counter(directed_edges(walk))
    target: Counter[tuple[int, int]] = Counter()
    for piece in pieces:
        assert piece[0] == piece[-1]
        assert len(set(piece[:-1])) == len(piece) - 1
        target.update(directed_edges(piece))
    assert target == source

    replay = build_certificate()["finite_pair_factorization"][
        "closed_walk_decomposition_replay"
    ]
    assert replay == {
        "vertices": 4,
        "maximum_edges": 7,
        "closed_walks_checked": 3276,
        "simple_cycle_factors_checked": 8664,
        "maximum_factors_in_one_walk": 3,
    }


def test_simple_cycle_slot_count_is_independently_fixed() -> None:
    assert simple_cycle_slot_bound(2) == 1
    assert simple_cycle_slot_bound(3) == 5
    assert simple_cycle_slot_bound(4) == 20
    assert RANK_MAX == 30
    assert simple_cycle_slot_bound(RANK_MAX) == 24893912605687593731774060679329


def test_vieta_box_is_exact_finite_and_matches_owner_archive() -> None:
    assert FIELD_DEGREE_MAX == 25
    assert CONJUGATE_HEIGHT == 62**30
    bounds = coefficient_bounds(FIELD_DEGREE_MAX, CONJUGATE_HEIGHT)
    assert len(bounds) == 25
    assert bounds[0] == 25 * CONJUGATE_HEIGHT
    assert bounds[-1] == CONJUGATE_HEIGHT**25
    assert all(isinstance(value, int) and value > 0 for value in bounds)

    box = build_certificate()["finite_algebraic_integer_box"]
    assert box["uniform_Vieta_bounds"] == [str(value) for value in bounds]
    assert (
        box["uniform_Vieta_bounds_sha256"]
        == "d8b4b1eea61f46b202b2a4fce6813e57866bbfbf79406e0ced26e45b45f5c3a1"
    )
    assert box["coefficient_box_union_upper_bound_decimal_digits"] == 17602
    assert (
        box["coefficient_box_union_upper_bound_sha256"]
        == "1b9cf329c1ca0200fcb965ca5422726958f4c47e68bb33b67bc5bd933fbd1f27"
    )
    assert box["finite"] is True


def test_field_generation_bridges_and_noncircular_order_are_explicit() -> None:
    factorization = build_certificate()["finite_pair_factorization"]
    assert factorization["uses_ABSW_finiteness_conclusion"] is False
    assert factorization["bounded_frame_taken_from_complete_RL_chamber"] is True
    assert factorization["individual_Cartan_entry_integrality_required"] is False
    assert factorization["cyclic_Cartan_product_integrality_required"] is True
    assert factorization["full_rank_frame_ground_field_inheritance"] is True
    assert factorization["finite_index_minimal_Vinberg_field_invariance"] is True
    assert factorization["field_generation"] == "k=Q(simple cyclic frame products)"
    assert factorization["fixed_data_invoked_only_after_finite_native_pairs"] is True


def test_fixed_data_maximality_and_complete_chamber_orbit_all_remain() -> None:
    closure = build_certificate()["fixed_data_and_chambers"]
    assert closure[
        "Nikulin_5_2_1_fixed_field_fixed_rank_similarity_classes"
    ] is True
    assert closure["projective_modules_in_scope"] is True
    assert closure["similarity_preserves_reflections_and_normalized_Gram"] is True
    assert closure["reflection_maximality_forces_W_equals_R"] is True
    assert closure[
        "arithmetic_only_and_absolute_reflection_maximality_agree"
    ] is True
    assert closure["complete_mirror_chambers_one_group_orbit"] is True


def test_weak_surface_corollary_is_not_misreported_as_leverage() -> None:
    corollary = build_certificate()["weak_rational_surface_corollary"]
    assert corollary["maximum_blowups"] == 2997
    assert corollary["effective_or_nef_return_claimed"] is False
    assert corollary["genuine_reflection_automorphisms_claimed"] is False
    assert corollary["independent_AG_finiteness_leverage"] is False

    ag_scope = build_certificate()["AG_scope"]
    assert ag_scope["unconditional_weak_numerical_carrier_exists"] is True
    assert ag_scope["weak_carrier_nonemptiness_restricts_inputs"] is False
    assert ag_scope["universal_strong_intrinsic_carrier_proved"] is False
    assert ag_scope["positive_AG_only_ABSW_proof_claimed"] is False
    assert ag_scope["only_presently_closed_positive_global_route"] == (
        "classical arithmetic"
    )


def test_replay_boundary_does_not_claim_to_machine_prove_literature() -> None:
    boundary = build_certificate()["replay_boundary"]
    assert boundary["owner_and_scope_report_bytes_verified"] is True
    assert boundary["owner_archive_payload_digest_verified"] is True
    assert boundary["finite_field_box_recomputed_independently"] is True
    assert boundary["temporary_source_hashes_checked_when_present"] is True
    assert boundary["temporary_source_files_required"] is False
    assert boundary["published_literature_machine_proved"] is False


def test_reports_preserve_exact_scope_and_tex_transport() -> None:
    combined = (
        OWNER_REPORT.read_text(encoding="utf-8")
        + AUDIT_REPORT.read_text(encoding="utf-8")
        + SCOPE_NO_GO_REPORT.read_text(encoding="utf-8")
    )
    assert re.search(r"(?<!\\)qquad", combined) is None
    for marker in (
        "PASS as an unconditional classical theorem",
        "Only cyclic Cartan products are required to be integral",
        "The selected cycles generate the original ground field",
        "No honest positive **AG-only** theorem",
        "A prime wall cannot itself be its regular reflection action",
    ):
        assert marker in combined


def test_archive_is_exact_replay() -> None:
    archived = json.loads(AUDIT_ARCHIVE.read_text(encoding="utf-8"))
    assert archived["payload"] == build_certificate()
    assert archived["payload_sha256"] == digest(archived["payload"])
    assert AUDIT_ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()


def test_owner_and_independent_cli_checks() -> None:
    subprocess.run(
        [sys.executable, str(OWNER_CERTIFIER), "--check-output"],
        cwd=ROOT,
        check=True,
    )
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check-output"],
        cwd=ROOT,
        check=True,
    )
