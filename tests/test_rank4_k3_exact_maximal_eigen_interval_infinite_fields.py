"""Regression tests for the exact maximal rank-four K3 eigen-intervals."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_rank4_k3_exact_maximal_eigen_interval_infinite_fields import (
    ambient_facet_record,
    build_certificate,
    companion,
    digest,
    endpoint_record,
    endpoint_reflections,
    enumerate_sector_roots,
    f,
    field_form,
    field_record,
    finite_root_audit,
    fm_multiply,
    fm_transpose,
    int_multiply,
    int_transpose,
    lift_field_matrix,
    maximality_trace_candidates,
    original_gram,
    original_tau,
    q_pair,
    reflection_record,
    square_equation_solutions,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = (
    ROOT
    / "results/rank4_k3_exact_maximal_eigen_interval_infinite_fields.md"
)
CERTIFIER = (
    ROOT
    / "scripts/certify_rank4_k3_exact_maximal_eigen_interval_infinite_fields.py"
)
ARCHIVE = (
    ROOT
    / "data/enumerations/rank4_k3_exact_maximal_eigen_interval_infinite_fields.json"
)


def test_unimodular_change_gives_the_diagonal_lattice():
    # Columns send (Y,X,U,Z) to (x1,x2,y1,y2).
    change = (
        (0, 1, 0, 0),
        (1, 0, 0, 0),
        (0, 0, 0, 1),
        (0, 0, 1, -1),
    )
    for a in range(2, 301):
        d = a - 1
        pullback = int_multiply(
            int_multiply(int_transpose(change), original_gram(a)),
            change,
        )
        assert pullback == (
            (2, 0, 0, 0),
            (0, -2, 0, 0),
            (0, 0, -2, 0),
            (0, 0, 0, -4 * d),
        )


def test_endpoint_roots_and_uniform_polarization_are_exact():
    for a in range(2, 501):
        row = endpoint_record(a)
        assert row["endpoint_root_q_norms"] == [1, 1]
        assert row["exact_open_component"] == "0<s<1"
        assert row["polarization_square"] == 14
        assert row["polarized_correspondence_height"] == 32 * a - 10


def test_bounded_root_shell_has_no_hidden_cut_in_the_interval():
    for a in range(2, 31):
        row = finite_root_audit(a, max_height=80)
        assert row["p0_strict_for_positive_height_roots"] is True
        assert row["p1_weak_for_positive_height_roots"] is True
        assert row["p1_equalities"] == [[1, -1, 1, 0]]
        assert row["square14_polarization_strict"] is True


def test_endpoint_field_reflections_are_integral_isometries_and_centralizers():
    for a in range(2, 301):
        row = reflection_record(a)
        assert row["field_reflection_checks"] == [True, True]
        assert row["lifts_preserve_original_gram"] is True
        assert row["lifts_commute_with_tau"] is True
        assert row["translation_half_trace"] == ["-1", "2"]
        assert row["cyclic_Cartan_product"] == ["0", "4"]

        form = field_form(a)
        for reflection in endpoint_reflections(a):
            assert fm_multiply(reflection, reflection, a) == (
                (f(1), f(0)),
                (f(0), f(1)),
            )
            assert fm_multiply(
                fm_multiply(fm_transpose(reflection), form, a),
                reflection,
                a,
            ) == form
            lift = lift_field_matrix(reflection, a)
            assert int_multiply(
                int_multiply(int_transpose(lift), original_gram(a)),
                lift,
            ) == original_gram(a)
            assert int_multiply(lift, original_tau(a)) == int_multiply(
                original_tau(a), lift
            )


def test_translation_has_the_claimed_field_matrix():
    for a in range(2, 101):
        s0, s1 = endpoint_reflections(a)
        translation = fm_multiply(s1, s0, a)
        assert translation == (
            (f(-1, 2), f(-2, 2)),
            (f(0, 2), f(-1, 2)),
        )


def test_maximality_trace_window_has_only_two_candidate_half_traces():
    for a in range(2, 501):
        assert maximality_trace_candidates(a, bound=40) == [[-1, 1], [0, 1]]


def test_both_remaining_maximality_square_equations_are_empty():
    for a in range(2, 151):
        result = square_equation_solutions(a, bound=100)
        assert result == {"P_equals_t": [], "P_equals_t_minus_one": []}


def test_ambient_root_family_supplies_unbounded_distinct_facets():
    for a in range(2, 41):
        d = a - 1
        seen = set()
        for k in range(1, 21):
            row = ambient_facet_record(a, k)
            assert row["root_q_norm"] == 1
            assert row["support_root_pairing"] == 0
            assert row["support_q_norm"] == -1 - row["height"] ** 2
            seen.add(tuple(row["root"]))
            support = tuple(row["support_point"])
            height = row["height"]
            for root in enumerate_sector_roots(d, min(30, height + 5)):
                assert q_pair(support, root, d) <= 0
        assert len(seen) == 20


def test_ambient_missing_long_reflections_are_integral():
    for d in range(1, 101):
        swap = (0, 1, -1, 0)
        z_sign = (0, 0, 0, 1)
        assert q_pair(swap, swap, d) == 2
        assert q_pair(z_sign, z_sign, d) == 2 * d
        # The norm-two reflection swaps X and U.
        vector = (7, 2, -3, 5)
        pairing = q_pair(vector, swap, d)
        image = tuple(
            entry - pairing * root_entry
            for entry, root_entry in zip(vector, swap)
        )
        assert image == (7, -3, 2, 5)
        # The second reflection changes only Z.
        pairing = q_pair(vector, z_sign, d)
        coefficient = pairing // d
        assert pairing % d == 0
        image = tuple(
            entry - coefficient * root_entry
            for entry, root_entry in zip(vector, z_sign)
        )
        assert image == (7, 2, -3, -5)


def test_Pell_field_image_is_not_parameter_counting():
    fields = [field_record(a) for a in range(2, 513)]
    distinct = {row["squarefree_field_parameter"] for row in fields}
    assert len(distinct) > 200
    assert field_record(2)["squarefree_field_parameter"] == 3
    assert field_record(7)["squarefree_field_parameter"] == 3
    for row in fields:
        assert (
            row["a"] ** 2
            - row["squarefree_field_parameter"] * row["Pell_y"] ** 2
            == 1
        )


def test_archive_replays_exactly_and_keeps_scope_literal():
    envelope = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = build_certificate(max_a=64)
    assert envelope["payload"] == payload
    assert envelope["payload_sha256"] == digest(payload)
    conclusions = payload["conclusions"]
    assert conclusions["finite_volume_eigen_interval"] is True
    assert conclusions["exact_prescribed_eigen_nef_return"] is True
    assert conclusions["complete_restricted_root_saturation"] is True
    assert conclusions["reflection_maximal_in_declared_Z_t_host"] is True
    assert conclusions["fixed_ample_polarization_square"] == 14
    assert conclusions["bounded_polarized_correspondence_height"] is False
    assert conclusions["finite_field_set"] is False
    assert conclusions["ambient_full_root_chamber_finite_volume"] is False
    scope = payload["scope"]
    assert scope["n_equals_one_only"] is True
    assert scope["n_at_least_two_claimed"] is False
    assert scope["maximal_across_overorders_claimed"] is False
    assert scope["maximal_across_commensurable_overlattices_claimed"] is False
    assert scope["maximal_in_full_O_L_claimed"] is False


def test_report_states_the_positive_theorem_and_every_nonclaim():
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    required = [
        "hyperbolic-dimension-one counterfamily",
        r"P_a=\mathbf P\bigl(\operatorname {Nef}(X_a)\cap E_t^+\bigr)",
        "complete restricted arrangement has one",
        "maximal among reflection groups preserving the",
        r"\cosh^2\ell_a=t",
        "minimal Vinberg field is exactly",
        "primitive ample with square \\(14\\)",
        r"h.\tau_a(h)=32a-10",
        "complete ambient K3 root chamber has infinite volume",
        "does not settle the required \\(n\\ge2\\) field theorem",
        "No maximality across arbitrary overorders or commensurable overlattices",
    ]
    for phrase in required:
        assert phrase in text


def test_no_reflective_classification_is_imported():
    combined = REPORT.read_text(encoding="utf-8") + CERTIFIER.read_text(
        encoding="utf-8"
    )
    for forbidden in (
        "Nikulin--Vinberg",
        "Vinberg classification",
        "ABSW theorem",
    ):
        assert forbidden not in combined


def test_invalid_inputs_are_rejected():
    with pytest.raises(ValueError):
        reflection_record(1)
    with pytest.raises(ValueError):
        ambient_facet_record(1, 2)
    with pytest.raises(ValueError):
        ambient_facet_record(2, 0)
    with pytest.raises(ValueError):
        build_certificate(max_a=1)
