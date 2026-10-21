"""Regression tests for the finite-volume K3 simplex field theorem."""

from __future__ import annotations

from fractions import Fraction
import json

from scripts.certify_k3_finite_volume_simplex_field_bound import (
    ARCHIVE,
    REPORT,
    build_certificate,
    digest,
    encoded_certificate,
    rank_rows,
    sha256_file,
)


def test_rank_rows_and_uniform_exponent_are_exact() -> None:
    rows = rank_rows()
    by_rank = {row["module_rank"]: row for row in rows}
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
    assert [by_rank[m]["binary_discriminant_exponent"] for m in range(3, 11)] == [
        240,
        180,
        120,
        66,
        24,
        26,
        28,
        30,
    ]
    assert all(by_rank[m]["rational_only"] for m in range(11, 21))
    assert max(row["binary_discriminant_exponent"] for row in rows) == 240


def test_symmetric_cycle_square_identity() -> None:
    gram = [
        [Fraction(2), Fraction(-1, 2), Fraction(-1, 3), Fraction(-1, 5)],
        [Fraction(-1, 2), Fraction(3), Fraction(-1, 7), Fraction(-1, 11)],
        [Fraction(-1, 3), Fraction(-1, 7), Fraction(5), Fraction(-1, 13)],
        [Fraction(-1, 5), Fraction(-1, 11), Fraction(-1, 13), Fraction(7)],
    ]

    def cartan(i: int, j: int) -> Fraction:
        return 2 * gram[i][j] / gram[i][i]

    for cycle in [(0, 1), (0, 1, 2), (0, 1, 2, 3)]:
        directed = Fraction(1)
        pair_product = Fraction(1)
        for pos, i in enumerate(cycle):
            j = cycle[(pos + 1) % len(cycle)]
            directed *= cartan(i, j)
            pair_product *= cartan(i, j) * cartan(j, i)
        assert directed * directed == pair_product


def test_archive_replays_and_scope_is_exact() -> None:
    envelope = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = build_certificate()
    assert envelope["payload"] == payload
    assert envelope["payload_sha256"] == digest(payload)
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()
    assert payload["report_sha256"] == sha256_file(REPORT)
    assert payload["uniform_binary_discriminant_exponent"] == 240
    assert int(payload["uniform_field_discriminant_bound"]) == 2**240
    assert payload["uniform_field_discriminant_bound_expression"] == "256^30"
    assert payload["conclusions"] == {
        "finite_field_set": True,
        "compact_simplices_covered": True,
        "ideal_finite_volume_simplices_covered": True,
        "imposed_mixed_height_needed": False,
        "angle_classification_used": False,
        "reflection_maximality_used_for_bound": False,
        "carrier_existence_proved": False,
        "essential_surjectivity_proved": False,
        "general_higher_facet_polytope_closed": False,
    }


def test_report_contains_the_cycle_proof_and_nonclaims() -> None:
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    for phrase in [
        r"m[k:\mathbf Q]\le20",
        r"c(C)^2=\prod_{r=1}^{\ell}p_{i_ri_{r+1}}",
        r"|\sigma(c(C))|\le2^\ell\le2^m",
        r"(m+5)d(d-1)",
        r"2^{240}=256^{30}",
        "Every closed directed walk decomposes into simple directed cycles",
        "does not construct the K3 carrier",
        "prove essential surjectivity",
        "general higher-facet polytopes",
    ]:
        assert phrase in text
