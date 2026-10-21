"""Regression tests for the projected-root K3 carrier of Delta(2,4,6)."""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import sympy as sp

from scripts.certify_compact_246_projected_k3_eigenslice_carrier import (
    G,
    H,
    OUTPUT,
    SCALE,
    ambient_gram,
    build_certificate,
    orthogonal_projector,
    rational_reflection_extension,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/compact_246_projected_k3_eigenslice_carrier.md"
EXPECTED_CERTIFICATE_SHA256 = (
    "57bae0a4f93ffe188d9dbd5c950249a22ee16c39e9767484009ba694ad007c3a"
)


def test_certificate_passes_and_is_archived() -> None:
    payload = build_certificate()
    assert payload["all_checks_pass"]
    assert payload["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_rank_six_lattice_and_primitive_trace_block() -> None:
    gram = ambient_gram()
    payload = build_certificate()["ambient_lattice"]
    assert SCALE == 15
    assert gram == (15 * G).row_join(G).col_join(G.row_join(H))
    assert payload["signature"] == [1, 5]
    assert payload["determinant"] == -393336
    assert payload["smith_diagonal"] == [1, 1, 1, 1, 18, 21852]
    assert payload["discriminant_length"] == 2
    assert payload["trace_embedding_primitive"] is True
    assert all(int(gram[i, i]) % 2 == 0 for i in range(6))


def test_projector_and_projected_physical_roots() -> None:
    gram = ambient_gram()
    projector = orthogonal_projector()
    assert projector * projector == projector
    assert projector.T * gram == gram * projector
    assert projector.rank() == 3
    assert gram[3:, 3:] == -2 * sp.eye(3)
    for i in range(3):
        delta = sp.eye(6)[:, 3 + i]
        expected = sp.zeros(6, 1)
        expected[i] = sp.Rational(1, 15)
        assert projector * delta == expected


def test_raw_reflections_have_rational_orthogonal_extensions() -> None:
    gram = ambient_gram()
    projector = orthogonal_projector()
    for i in range(3):
        extension = rational_reflection_extension(i)
        assert extension * extension == sp.eye(6)
        assert extension.T * gram * extension == gram
        assert extension.T * gram == gram * extension
        assert extension * projector == projector * extension
        assert any(value.q == 15 for value in extension)


def test_visible_root_a_ellipsoid_is_exhaustive() -> None:
    payload = build_certificate()["visible_root_reduction"]
    assert payload["ellipsoid_leading_principal_minors"] == [18, 396, 10926]
    expected = {
        (-1, 0, 0),
        (1, 0, 0),
        (0, -1, 0),
        (0, 1, 0),
        (0, 0, -1),
        (0, 0, 1),
        (0, 0, 0),
    }
    assert {tuple(row["a"]) for row in payload["ellipsoid_lattice_points"]} == expected

    matrix = 30 * sp.eye(3) + G
    # The archived dual bounds prove every eligible coordinate is in this box;
    # this loop independently filters all 27 points.
    survivors = set()
    for entries in itertools.product(range(-1, 2), repeat=3):
        vector = sp.Matrix(entries)
        if (vector.T * matrix * vector)[0] < 30:
            survivors.add(entries)
    assert survivors == expected


def test_finite_hidden_root_window_and_residue_filter() -> None:
    window = build_certificate()["finite_hidden_root_window"]
    assert window["integer_coordinate_bounds"] == [2, 5, 6]
    assert window["enumerated_candidate_count"] == 22
    assert {tuple(row["y"]) for row in window["mod_15_residue_survivors"]} == {
        (-1, 0, 0),
        (1, 0, 0),
        (0, -1, 0),
        (0, 1, 0),
        (0, 0, -1),
        (0, 0, 1),
    }
    assert window["conclusion"].endswith("three original facet lines")


def test_hidden_root_window_rebuilt_independently() -> None:
    p = G.inv() * sp.ones(3, 1)
    p_square = (p.T * G * p)[0]
    qp = sp.Rational(2, 1) / p_square * (G * p) * (G * p).T - G
    max_cosh_2r = sp.Rational(177, 65)
    candidates = []
    for entries in itertools.product(range(-2, 3), range(-5, 6), range(-6, 7)):
        y = sp.Matrix(entries)
        square = (y.T * G * y)[0]
        if -square not in (2, 6, 12):
            continue
        if (y.T * qp * y)[0] <= max_cosh_2r * (-square):
            candidates.append((entries, int(square)))
    assert len(candidates) == 22

    survivors = set()
    for index, positive_norm in enumerate((12, 6, 2)):
        for sign in (-1, 1):
            residue = [0, 0, 0]
            residue[index] = sign
            for entries, square in candidates:
                if square != -positive_norm:
                    continue
                if all((entries[j] - residue[j]) % 15 == 0 for j in range(3)):
                    survivors.add(entries)
    assert survivors == {
        (-1, 0, 0),
        (1, 0, 0),
        (0, -1, 0),
        (0, 1, 0),
        (0, 0, -1),
        (0, 0, 1),
    }


def test_k3_gate_and_scope_are_exact() -> None:
    payload = build_certificate()
    k3 = payload["K3_realization"]
    assert k3["K3_complement_signature"] == [2, 14]
    assert k3["discriminant_length"] == 2
    assert k3["stable_embedding_inequality"] == "16>=2+2"
    assert k3["primitive_embedding"] is True
    assert "Delta(2,4,6)" in k3["exact_return"]
    assert "n<=19" in k3["carrier_relative_dimension_bound"]
    scope = payload["scope"]
    assert "proper rational projector/eigen-nef" in scope["proved"]
    assert "literal full-nef" in scope["outside_current_contract"]
    assert any("universal existence" in item for item in scope["not_proved"])
    assert any("integral extension" in item for item in scope["not_proved"])


def test_report_states_positive_result_and_frontier() -> None:
    text = REPORT.read_text(encoding="utf-8")
    for token in (
        "positive fourth-category calibration",
        "proper rational projector/eigen-nef",
        "ambient square-\\((-2)\\) mirror meeting the **closed** restricted",
        "a^t(30I_3+G)a<30",
        "Q_p(y)\\le d\\cosh(2R)",
        "y=\\pm e_1,\\quad\\pm e_2",
        "preserve \\(S\\) integrally",
        "is an actual wall",
        "are then effective simple roots",
        "hence}\\qquad n\\le19",
        "does not obstruct algebraic geometry itself",
    ):
        assert token in text
    assert "or a replacement dimension bound" not in text
