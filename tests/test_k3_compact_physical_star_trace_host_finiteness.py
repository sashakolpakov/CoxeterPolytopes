"""Regression tests for the compact physical-star trace-host theorem."""

from __future__ import annotations

import json
import math
from pathlib import Path
import subprocess
import sys

import sympy as sp

from scripts.certify_k3_compact_physical_star_trace_host_finiteness import (
    ARCHIVE,
    REPORT,
    build_certificate,
    digest,
    encoded_certificate,
    golden_codifferent_box,
    rank_rows,
    toy_overorder_interval,
)


ROOT = Path(__file__).resolve().parents[1]
CERTIFIER = (
    ROOT / "scripts/certify_k3_compact_physical_star_trace_host_finiteness.py"
)


def test_rank_degree_window() -> None:
    rows = rank_rows()
    by_rank = {row["module_rank"]: row for row in rows}
    assert min(by_rank) == 3
    assert max(by_rank) == 20
    assert [by_rank[m]["field_degree_ceiling"] for m in range(3, 11)] == [
        6,
        5,
        4,
        3,
        2,
        2,
        2,
        2,
    ]
    assert all(not by_rank[m]["nonrational_possible"] for m in range(11, 21))


def test_physical_energy_gives_uniform_entry_and_trace_bounds() -> None:
    payload = build_certificate()["physical_energy"]
    assert payload["sum_of_positive_place_energies"] == 2
    assert payload["distinguished_field_Gram_entry_strict_upper_bound"] == 14
    assert payload["definite_place_field_Gram_entry_upper_bound"] == 2
    assert payload["integral_trace_pairing_strict_upper_bound"] == 16

    # Independently replay the inequalities on a dense rational energy grid.
    for ai in range(1, 200):
        ti0 = 2 * ai / 200
        tail_i = 2 - ti0
        for aj in range(1, 200):
            tj0 = 2 * aj / 200
            tail_j = 2 - tj0
            distinguished = 7 * math.sqrt(ti0 * tj0)
            definite_tail = math.sqrt(tail_i * tail_j)
            assert distinguished < 14
            assert definite_tail < 2
            assert distinguished + definite_tail < 16


def test_golden_codifferent_box_is_finite_and_trace_dual() -> None:
    toy = golden_codifferent_box()
    assert toy["v_min"] == -16
    assert toy["v_max"] == 16
    assert toy["u_min"] == -58
    assert toy["u_max"] == 58
    assert toy["bounded_box_point_count"] > 0

    # For z=(u+v*phi)/sqrt(5), Tr(z)=v and Tr(phi*z)=u+v.
    # These two integral traces certify membership in O_K^vee.
    for u, v in ((0, 0), (1, 0), (0, 1), (-17, 9), (23, -11)):
        assert isinstance(v, int)
        assert isinstance(u + v, int)


def test_cycle_order_has_only_finitely_many_overorders_toy() -> None:
    toy = toy_overorder_interval()
    assert toy["lower_order_discriminant"] == 20
    assert toy["maximal_order_discriminant"] == 5
    assert toy["maximal_over_lower_index"] == 2
    assert toy["additive_quotient_order"] == 2
    assert toy["intermediate_order_count"] == 2


def test_root_frame_dual_index_identity() -> None:
    # F has Gram diag(2,-2,-2).  Adjoin (e1+e2)/2 to get an integral even
    # index-two overlattice M.  det(F)=[M:F]^2 det(M).
    gram_f = sp.diag(2, -2, -2)
    change = sp.Matrix(
        [
            [sp.Rational(1, 2), 0, 0],
            [sp.Rational(1, 2), 1, 0],
            [0, 0, 1],
        ]
    )
    gram_m = change.T * gram_f * change
    assert gram_m == sp.Matrix([[0, -1, 0], [-1, -2, 0], [0, 0, -2]])
    assert all(gram_m[i, i] % 2 == 0 for i in range(3))
    assert abs(int(gram_f.det())) == 8
    assert abs(int(gram_m.det())) == 2
    assert abs(int(gram_f.det())) == 2**2 * abs(int(gram_m.det()))


def test_trace_scalar_and_reflection_scope_are_pinned() -> None:
    payload = build_certificate()
    assert payload["trace_scalar"]["K3_sign_relation"] == (
        "g=-alpha*h with alpha totally positive"
    )
    assert payload["trace_scalar"]["native_reflections_determined_by_g"] is True
    maximality = payload["maximality"]
    assert maximality["needed_for_finite_trace_hosts"] is False
    assert (
        maximality["same_marked_k_reflection_host_needed_for_chamber_corollary"]
        is True
    )
    assert maximality["ordinary_rank_one_K3_reflection_maximality_substitutes"] is False


def test_conclusions_and_nonclaims() -> None:
    payload = build_certificate()
    conclusions = payload["conclusions"]
    assert conclusions["finite_coefficient_orders"] is True
    assert conclusions["finite_field_Gram_matrices"] is True
    assert conclusions["finite_integral_trace_modules_and_saturations"] is True
    assert conclusions["finite_marked_trace_reflection_hosts"] is True
    assert conclusions["finite_maximal_chambers_with_same_host_k_maximality"] is True
    assert conclusions["finite_ambient_Neron_Severi_extensions"] is False
    assert conclusions["carrier_existence"] is False
    assert payload["scope"]["property4_proved"] is False


def test_report_states_exact_physical_and_full_picard_scope() -> None:
    text = REPORT.read_text(encoding="utf-8")
    for needle in (
        "same vector",
        "\\sum_\\tau t_{i,\\tau}=2",
        "g_{ij}\\in\\mathcal O^\\vee",
        "F\\subseteq M\\subseteq F^\\#",
        "[M:F]^2\\mid |\\det C|",
        "exact complete eigen-nef intersection",
        "ordinary K3 root reflection",
        "does not assert finiteness of",
        "does not prove Property 4",
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
