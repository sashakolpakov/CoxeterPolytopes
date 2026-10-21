"""Independent hostile tests for the projected (2,4,6) K3 eigenslice."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.audit_compact_246_projected_k3_eigenslice_carrier_independent_hostile import (
    AUDIT_ARCHIVE,
    AUDIT_REPORT,
    OWNER_HASHES,
    ambient_packet,
    build_certificate,
    correspondence_packet,
    digest,
    encoded_certificate,
    hidden_root_window,
    k3_and_chamber_packet,
    owner_hash_audit,
    visible_root_reduction,
)


AUDITOR = (
    ROOT
    / "scripts/audit_compact_246_projected_k3_eigenslice_carrier_"
    "independent_hostile.py"
)


def test_owner_bytes_are_terminal_and_auditor_is_independent() -> None:
    assert owner_hash_audit() == OWNER_HASHES
    source = AUDITOR.read_text(encoding="utf-8")
    assert "certify_compact_246_projected_k3_eigenslice_carrier import" not in source
    assert "import sympy" not in source


def test_ambient_lattice_snf_primitive_trace_and_projector() -> None:
    packet = ambient_packet()
    assert packet["ambient_even_integral"] is True
    assert packet["raw_LDL_pivots"] == [-12, -3, 1]
    assert packet["ambient_signature"] == [1, 5]
    assert packet["ambient_determinant"] == -393336
    assert packet["determinantal_divisors"] == [1, 1, 1, 1, 18, 393336]
    assert packet["smith_invariants"] == [1, 1, 1, 1, 18, 21852]
    assert packet["discriminant_length"] == 2
    assert packet["trace_primitive_direct_summand"] is True
    assert packet["projector_idempotent"] is True
    assert packet["projector_self_adjoint"] is True
    assert packet["kernel_negative_definite"] is True
    assert packet["projected_delta_vectors"] == [
        ["1/15", 0, 0, 0, 0, 0],
        [0, "1/15", 0, 0, 0, 0],
        [0, 0, "1/15", 0, 0, 0],
    ]


def test_visible_root_sign_ellipsoid_and_coordinate_bound_are_exhaustive() -> None:
    reduction = visible_root_reduction()
    sign = reduction["visible_root_implies_negative_Gy"]
    assert "negative definite" in sign["y_positive"]
    assert "boundary" in sign["y_null"]
    assert sign["closed_compact_P_points_are_positive"] is True
    assert sign["conclusion"].endswith("G(y,y)<0")
    assert reduction["a_matrix_leading_minors"] == [18, 396, 10926]
    assert reduction["dual_coordinate_square_bounds"] == [
        "1105/607",
        "840/607",
        "660/607",
    ]
    assert {tuple(row["a"]) for row in reduction["ellipsoid_lattice_points"]} == {
        (-1, 0, 0),
        (1, 0, 0),
        (0, -1, 0),
        (0, 1, 0),
        (0, 0, -1),
        (0, 0, 1),
        (0, 0, 0),
    }
    assert reduction["completed_square_identity_checked_on_729_pairs"] is True
    assert "15*G(x,x)=-2" in reduction["a_zero_root_impossible"]
    assert "multiple of 15 equal to 2" in reduction[
        "zero_projection_root_impossible"
    ]


def test_qp_circumradius_window_is_exact_and_finite() -> None:
    window = hidden_root_window()
    assert window["interior_point"] == ["11/12", 2, "5/2"]
    assert window["integral_multiple_12p"] == [11, 24, 30]
    assert window["p_square"] == "65/12"
    assert window["wall_slacks"] == [1, 1, 1]
    assert [row["cosh_twice_distance"] for row in window["vertices"]] == [
        "177/65",
        "79/65",
        "17/13",
    ]
    assert window["max_cosh_twice_radius"] == "177/65"
    assert window["Qp_leading_minors"] == ["804/65", "612/13", 36]
    assert "d*cosh" in window["hyperplane_distance_identity"]
    assert window["universal_Qp_bound"] == "2124/65"
    assert window["coordinate_square_bounds"] == [
        "31329/4225",
        "111864/4225",
        "36108/845",
    ]
    assert window["integer_coordinate_bounds"] == [2, 5, 6]
    assert window["box_size"] == 715
    assert window["root_window_count"] == 22


def test_residue_filter_is_exhaustive_and_roots_above_facets_are_unique() -> None:
    window = hidden_root_window()
    observed = {
        (tuple(row["a"]), tuple(row["y"]), tuple(row["x"]))
        for row in window["residue_survivors"]
    }
    expected = {
        ((-1, 0, 0), (-1, 0, 0), (0, 0, 0)),
        ((1, 0, 0), (1, 0, 0), (0, 0, 0)),
        ((0, -1, 0), (0, -1, 0), (0, 0, 0)),
        ((0, 1, 0), (0, 1, 0), (0, 0, 0)),
        ((0, 0, -1), (0, 0, -1), (0, 0, 0)),
        ((0, 0, 1), (0, 0, 1), (0, 0, 0)),
    }
    assert observed == expected
    assert window["residue_exhaustive_because_bounds_below_15"] is True
    assert "rho=+/-delta_i" in window["unique_roots_above_facets"]
    assert "generic facet point" in window["actual_ambient_wall_argument"]
    assert window["all_visible_root_mirrors_on_closed_P_are_facets"] is True


def test_raw_actions_are_only_rational_and_correspondences_are_exact() -> None:
    ambient = ambient_packet()
    packet = correspondence_packet(ambient["_ambient"], ambient["_projector"])
    assert packet["raw_reflections_integral_on_T"] is True
    assert packet["extensions_are_involutive_isometries"] is True
    assert packet["extensions_self_adjoint"] is True
    assert packet["extensions_commute_with_projector"] is True
    assert packet["extensions_preserve_full_NS_integrally"] is False
    assert all(
        row["ambient_delta_reflection_preserves_T"] is False
        and row["projected_action_equals_raw"] is False
        for row in packet["physical_root_reflection_comparison"]
    )
    assert "D_i^vee x f(D_i)" in packet["divisor_product_formula"]
    assert packet[
        "projector_and_extensions_are_rational_algebraic_correspondences"
    ] is True


def test_stable_k3_gate_walls_irreducibility_and_eigen_nef_equality() -> None:
    ambient = ambient_packet()
    window = hidden_root_window()
    k3 = k3_and_chamber_packet(ambient, window)
    assert k3["complement_signature"] == [2, 14]
    assert k3["complement_rank"] == 16
    assert k3["discriminant_length"] == 2
    assert k3["stable_embedding_inequality"] == "16>=2+2"
    assert k3["primitive_K3_embedding"] is True
    assert k3["generic_period_exact_NS"] is True
    assert k3["no_zero_projection_roots"] is True
    assert k3["p_regular_for_ambient_root_arrangement"] is True
    assert k3["facet_roots_effective_simple"] is True
    assert "primitivity" in k3["facet_irreducibility_argument"]
    assert k3["facet_curves_smooth_rational"] is True
    assert "Delta(2,4,6)" in k3["exact_eigen_nef_equality"]
    assert "n<=19" in k3["carrier_relative_dimension_bound"]


def test_scope_is_positive_but_not_a_universal_fourth_branch() -> None:
    scope = build_certificate()["scope"]
    assert scope["one_projected_K3_carrier_proved"] is True
    assert scope["outside_original_k_Q_literal_full_nef_Property2"] is True
    assert scope["Picard_rank"] == 6
    assert scope["eigenslice_rank"] == 3
    assert scope["raw_reflections_extend_rationally"] is True
    assert scope["raw_reflections_extend_integrally"] is False
    assert scope["ambient_root_reflections_return_raw_Coxeter_action"] is False
    assert scope["universal_fourth_branch_nonempty"] is False
    assert scope["finite_host_control_for_fourth_branch"] is False
    assert scope["arbitrary_non_K3_eigenslice_dimension_bound"] is False
    assert scope["K3_carrier_relative_n_bound"] == 19


def test_report_contains_every_hostile_gate_and_scope_boundary() -> None:
    text = AUDIT_REPORT.read_text(encoding="utf-8")
    for needle in (
        "PASS at the stated projected-root scope",
        "must satisfy",
        "G(y,y)<0",
        "all \\(27\\) triples",
        "all \\(715\\) points",
        "Exactly \\(22\\) remain",
        "unique ambient roots",
        "actual walls",
        "stable primitive-embedding",
        "effective simple roots",
        "irreducible, not merely effective",
        "two inclusions in the chamber equality",
        "divisor-product",
        "n\\le19",
        "universal nonemptiness or finite-host control",
        "No owner byte",
    ):
        assert needle in text
    assert text.count("\\[") == text.count("\\]")


def test_archive_replay_and_cli() -> None:
    archived = json.loads(AUDIT_ARCHIVE.read_text(encoding="utf-8"))
    assert archived["payload"] == build_certificate()
    assert archived["payload_sha256"] == digest(archived["payload"])
    assert AUDIT_ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()
    subprocess.run(
        [sys.executable, str(AUDITOR), "--check-output"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
