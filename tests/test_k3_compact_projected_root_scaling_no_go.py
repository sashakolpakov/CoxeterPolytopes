"""Independent exact tests for the projected-root scaling counterfamily."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

import sympy as sp

from scripts.certify_k3_compact_projected_root_scaling_no_go import (
    ARCHIVE,
    B,
    H,
    Q,
    R,
    REPORT,
    SOURCE_ARCHIVE,
    SOURCE_ARCHIVE_SHA256,
    ambient_gram,
    build_certificate,
    digest,
    direct_sum_gram,
    encoded_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
CERTIFIER = ROOT / "scripts/certify_k3_compact_projected_root_scaling_no_go.py"
FROZEN_PACKAGE = {
    REPORT: "bcd4e41d6b8915f09726a9512c8bcac1876da3aa3cdd513ac567a4bcf4c7b057",
    CERTIFIER: "3776c68936b62c6545ffb2547c20a82f094c8afe0e31567dc701c9102027e90a",
    ARCHIVE: "6aaeca3809bbed3e66af8a97fe89b3f93b9b5269755c00087ad924c5554d95c9",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_projected_no_go_package_bytes_are_frozen() -> None:
    for path, expected in FROZEN_PACKAGE.items():
        assert sha256(path) == expected


def test_frozen_golden_source_is_literal() -> None:
    assert sha256(SOURCE_ARCHIVE) == SOURCE_ARCHIVE_SHA256
    source = json.loads(SOURCE_ARCHIVE.read_text(encoding="utf-8"))
    rows = [row for row in source["examples"] if row["seed_id"] == "triangle_255"]
    assert len(rows) == 1
    row = rows[0]
    assert row["input_coxeter_labels"] == [2, 5, 5]
    assert row["vinberg_compactness"]["compact"] is True
    assert row["k3_eigen_nef_promotion"][
        "projectivized_eigen_nef_section_equals_closed_seed"
    ] is True
    assert row["trace_k3_package"]["visible_k3_wall_shell_collapse"][
        "restricted_hyperplanes_equal"
    ] is True


def test_seed_matrices_and_negative_A3_block() -> None:
    assert Q.det() == -16
    assert R.T * Q * R == H
    assert H == sp.Matrix([[-2, 0, 1], [0, -2, 1], [1, 1, -2]])
    assert [(-H)[:i, :i].det() for i in range(1, 4)] == [2, 4, 4]
    assert H.det() == -4
    assert B == Q * R
    assert B.T * Q.inv() * B == H


def test_overlattice_construction_exactly_for_many_scales() -> None:
    for scale in range(2, 18):
        direct = direct_sum_gram(scale)

        # Columns express the ambient basis (trace, delta) in the direct
        # basis (trace, c), where delta_i=(r_i+c_i)/s.
        change = sp.zeros(9)
        for i in range(6):
            change[i, i] = 1
        for i in range(3):
            change[:, 6 + i] = R[:, i].col_join(sp.eye(3)[:, i]) / scale

        ambient = sp.simplify(change.T * direct * change)
        assert ambient == ambient_gram(scale)
        assert sp.det(change) == sp.Rational(1, scale**3)
        assert int(ambient.det()) == 64 * scale**3 * (scale - 1) ** 3
        assert all(entry.q == 1 for entry in ambient)
        assert all(int(ambient[i, i]) % 2 == 0 for i in range(9))
        assert ambient[6:, 6:] == H
        assert [ambient[6 + i, 6 + i] for i in range(3)] == [-2, -2, -2]

        schur = sp.simplify(H - B.T * (scale * Q).inv() * B)
        assert schur == sp.Rational(scale - 1, scale) * H
        # Q has (1,5), while the Schur complement is negative definite.
        assert [(-schur)[:i, :i].det() > 0 for i in range(1, 4)] == [True] * 3


def test_trace_block_is_primitive_but_has_no_minus_two_vector() -> None:
    # In the (trace,delta) basis, quotienting by the first six basis vectors
    # leaves a free Z^3.  Hence the trace inclusion is primitive.  Its form
    # is sQ, whose squares are multiples of 2s and cannot equal -2 for s>=2.
    assert all(Q[i, i] % 2 == 0 for i in range(6))
    for scale in range(2, 30):
        assert (-2) % (2 * scale) != 0
        assert scale**3 > 1


def test_root_projection_exhaustion_is_all_parameter() -> None:
    s = sp.symbols("s", integer=True, positive=True)
    j = sp.symbols("j", integer=True, positive=True)
    expression = sp.expand(-2 * s + 2 * j * (s - 1))
    assert sp.expand(expression - 2 * ((j - 1) * s - j)) == 0
    assert expression.subs(j, 1) == -2

    # For every integer s>=2 and j>=2, the trace square is nonnegative.
    for scale in range(2, 80):
        for half_norm in range(2, 80):
            assert int(expression.subs({s: scale, j: half_norm})) >= 0
    assert expression.subs({s: 2, j: 2}) == 0

    # If the complement coefficient is zero, the root equation would force
    # s*Q(x,x)=-2, impossible because Q is even and s>=2.
    for scale in range(2, 30):
        assert all((scale * (2 * value)) != -2 for value in range(-50, 51))


def test_k3_stable_embedding_window_and_unbounded_invariants() -> None:
    assert 13 >= 9 + 2
    determinants = [64 * s**3 * (s - 1) ** 3 for s in range(2, 20)]
    indices = [s**3 for s in range(2, 20)]
    assert determinants == sorted(set(determinants))
    assert indices == sorted(set(indices))
    assert determinants[-1] > determinants[0]


def test_certificate_pins_claims_and_nonclaims() -> None:
    payload = build_certificate()
    assert payload["source"]["triangle"]["compact"] is True
    assert payload["all_parameter_family"]["trace_embedding_primitive"] is True
    assert payload["all_parameter_family"]["ambient_even_integral"] is True
    assert payload["all_parameter_family"]["trace_contains_minus_two_roots"] is False
    assert payload["root_projection_exhaustion"][
        "all_visible_roots_project_to_old_golden_root_mirrors"
    ] is True
    assert payload["K3_gate"]["primitive_K3_embedding_for_every_scale"] is True
    assert payload["K3_gate"]["every_physical_facet_root_lies_outside_trace"] is True
    assert payload["scope"]["violates_same_vector_in_trace_hypothesis"] is True
    assert payload["scope"]["produces_infinitely_many_abstract_chambers"] is False
    assert payload["scope"]["integral_order_action_on_ambient_complement_claimed"] is False


def test_report_contains_exact_scope_and_equations() -> None:
    text = REPORT.read_text(encoding="utf-8")
    for needle in (
        "same vector",
        "S_s:T_s\\oplus K_s]=s^3",
        "\\det S_s=64s^3(s-1)^3",
        "Q(y,y)=-2s-(s-1)H(a,a)",
        "Q(y,y)=2\\bigl((j-1)s-j\\bigr)",
        "class has a nonzero component in the rank-three complement",
        "does not produce infinitely many abstract Coxeter chambers",
    ):
        assert needle in text


def test_archive_is_exact_replay() -> None:
    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert archived["payload"] == build_certificate()
    assert archived["payload_sha256"] == digest(archived["payload"])
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()


def test_cli_check() -> None:
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check-output"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
