"""Hostile regression tests for the rational isotropic period gate."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_rational_isotropic_weak_log_period_gate import (
    CANONICAL,
    FIBRE,
    GRAM,
    build_certificate,
    digest,
    dot,
    encoded_certificate,
    exceptional,
    lattice_record,
    nontorsion_recurrence,
    normalize_orders,
    riemann_roch_euler,
    torsion_recurrence,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_isotropic_weak_log_period_gate.md"
CERTIFIER = (
    ROOT / "scripts/certify_rational_isotropic_weak_log_period_gate.py"
)
ARCHIVE = (
    ROOT
    / "data/enumerations/rational_isotropic_weak_log_period_gate.json"
)
AUDIT = (
    ROOT / "results/rational_isotropic_weak_log_period_gate_hostile_audit.md"
)


def independent_dot(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return left[0] * right[0] - sum(
        a * b for a, b in zip(left[1:], right[1:])
    )


def test_marked_blowup_lattice_is_exactly_I_1_9():
    assert len(GRAM) == 10
    assert all(len(row) == 10 for row in GRAM)
    assert [GRAM[index][index] for index in range(10)] == [1] + [-1] * 9
    assert all(
        GRAM[row][column] == 0
        for row in range(10)
        for column in range(10)
        if row != column
    )
    assert FIBRE == tuple(-entry for entry in CANONICAL)
    assert independent_dot(FIBRE, FIBRE) == 0
    assert independent_dot(CANONICAL, CANONICAL) == 0
    assert independent_dot(CANONICAL, FIBRE) == 0
    assert dot(FIBRE, FIBRE) == independent_dot(FIBRE, FIBRE)
    assert lattice_record()["determinant"] == -1
    assert lattice_record()["F_primitive"] is True


def test_exceptional_pairings_and_invalid_indices():
    for index in range(1, 10):
        e_i = exceptional(index)
        assert dot(e_i, e_i) == -1
        assert dot(FIBRE, e_i) == 1
    with pytest.raises(ValueError):
        exceptional(0)
    with pytest.raises(ValueError):
        exceptional(10)


def test_riemann_roch_polynomial_is_identically_one_on_large_range():
    for r in range(0, 2001):
        # Independent expansion: F^2=K.F=0, so the correction term vanishes.
        assert r * r * independent_dot(FIBRE, FIBRE) - r * independent_dot(
            FIBRE, CANONICAL
        ) == 0
        assert riemann_roch_euler(r) == 1


def test_torsion_profiles_have_exact_floor_formula_and_first_move():
    for order in range(1, 41):
        rows = torsion_recurrence(order, 400)
        for row in rows:
            r = row["r"]
            expected = r // order + 1
            assert row["h0"] == expected
            assert row["h1"] == expected - 1
            assert row["h2"] == 0
            assert row["chi"] == 1
            if r:
                assert row["normal_power_trivial"] is (r % order == 0)
                assert row["restriction_map_rank"] == int(r % order == 0)
        assert all(row["h0"] == 1 for row in rows[1:order])
        assert rows[order]["h0"] == 2


def test_non_torsion_profile_never_moves():
    rows = nontorsion_recurrence(1000)
    assert all(row["h0"] == 1 for row in rows)
    assert all(row["h1"] == row["h2"] == 0 for row in rows)
    assert all(row["chi"] == 1 for row in rows)
    assert all(row["normal_power_trivial"] is False for row in rows[1:])
    assert all(row["restriction_map_rank"] == 0 for row in rows[1:])


def test_archive_is_a_byte_exact_replay_with_digest():
    expected = encoded_certificate()
    assert ARCHIVE.read_text(encoding="utf-8") == expected
    envelope = json.loads(expected)
    assert envelope["payload"] == build_certificate()
    assert envelope["payload_sha256"] == digest(envelope["payload"])


def test_archive_scope_forbids_the_dangerous_upgrades():
    payload = build_certificate(horizon=8, orders=(1, 4))
    scope = payload["scope"]
    assert scope["smooth_genus_one_representative_assumed"] is True
    assert scope["same_marked_NS_K_F_and_nefness_proved"] is True
    for forbidden_claim in (
        "every_isotropic_weak_log_class_covered",
        "same_complete_nef_cone_claimed",
        "prescribed_Coxeter_polytope_realized",
        "complete_nef_chamber_returned",
        "all_physical_facets_returned",
        "finite_hyperbolic_volume_claimed",
        "reflection_saturation_claimed",
        "reflection_maximality_claimed",
        "raw_polytope_exhaustivity_claimed",
        "computer_replay_replaces_geometric_proof",
    ):
        assert scope[forbidden_claim] is False


def test_report_contains_the_exact_theorem_proof_and_two_examples():
    text = REPORT.read_text(encoding="utf-8")
    normalized = " ".join(text.split())
    required = (
        r"N:=\mathcal O_F(F)\in\operatorname {Pic}^0(F)",
        r"\left\lfloor\frac r m\right\rfloor+1",
        r"\left\lfloor\frac r m\right\rfloor",
        r"0\longrightarrow\mathcal O_S((r-1)F)",
        r"H^2(S,\mathcal O_S(rF))",
        r"the section \(u^q\) restricts nontrivially",
        "base-point-free pencil",
        r"[F]=3H-E_1-\cdots-E_9=[-K_S]",
        r"\mathcal O_C(3)(-p_1-\cdots-p_9)",
        "base locus of the resulting cubic pencil",
        "translate of the Abel--Jacobi isomorphism",
        r"same-\((\operatorname {NS},K,F)\) obstruction",
        "does not assume that the pencil morphism already has connected fibres",
        "does it silently use a Stein factorization",
        r"\(F\sim-K_S\) linearly, not merely numerically",
    )
    for phrase in required:
        assert phrase in normalized


def test_report_states_every_material_nonclaim_explicitly():
    normalized = " ".join(REPORT.read_text(encoding="utf-8").split())
    exclusions = (
        "every isotropic weak-log class has a smooth genus-one representative",
        "realize a prescribed Coxeter polytope",
        "a complete nef chamber or all effective facet curves",
        "finite hyperbolic volume",
        "reflection group is saturated or maximal",
        "exhaustive for raw arithmetic polytopes",
        "not a claim that the two complete nef cones are marked-isomorphic",
        "finite replay does not replace the cohomological or moduli arguments",
    )
    for phrase in exclusions:
        assert phrase in normalized


def test_package_imports_no_reflective_classification():
    combined = REPORT.read_text(encoding="utf-8") + CERTIFIER.read_text(
        encoding="utf-8"
    )
    for forbidden in ("Nikulin", "Vinberg", "ABSW"):
        assert forbidden not in combined


def test_manual_hostile_audit_closes_the_exact_proof_attacks_only():
    normalized = " ".join(AUDIT.read_text(encoding="utf-8").split())
    required = (
        "PASS AFTER REPAIR",
        "manual internal proof audit",
        "not represented as an independent external referee report",
        "does not covertly use either one",
        "silently assume connected fibres or a trivial Stein factorization",
        "equality of divisor classes",
        "linear-equivalence statement",
        "opposite answers to semiampleness",
        "does **not** assert that their complete nef cones are marked-isomorphic",
        "finite hyperbolic volume",
        "reflection saturation or maximality",
        "exhaustivity for raw arithmetic polytopes",
    )
    for phrase in required:
        assert phrase in normalized
    for placeholder in (
        "REPORT_SHA256",
        "CERTIFIER_SHA256",
        "ARCHIVE_SHA256",
        "TEST_SHA256",
    ):
        assert placeholder not in normalized


def test_invalid_replay_parameters_are_rejected():
    with pytest.raises(ValueError):
        riemann_roch_euler(-1)
    with pytest.raises(ValueError):
        torsion_recurrence(0, 4)
    with pytest.raises(ValueError):
        torsion_recurrence(2, -1)
    with pytest.raises(ValueError):
        nontorsion_recurrence(-1)
    with pytest.raises(ValueError):
        normalize_orders(())
    with pytest.raises(ValueError):
        normalize_orders((1, -2))
    with pytest.raises(ValueError):
        build_certificate(horizon=0)
