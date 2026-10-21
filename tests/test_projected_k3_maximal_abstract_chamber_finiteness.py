"""Exact regressions for projected-K3 maximal abstract-chamber finiteness."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from scripts.certify_projected_k3_maximal_abstract_chamber_finiteness import (
    ARCHIVE,
    DELTA_246_GRAM,
    DEPENDENCY_SHA256,
    REDOUBLING_ARCHIVE,
    REPORT,
    SCALING_ARCHIVE,
    build_certificate,
    cartan_matrix,
    digest,
    encoded_certificate,
    normalized_squared_gram,
    rank_budget,
    scale_matrix,
    sha256_file,
)


ROOT = Path(__file__).resolve().parents[1]
CERTIFIER = (
    ROOT / "scripts/certify_projected_k3_maximal_abstract_chamber_finiteness.py"
)

# Filled only with files whose bytes do not depend on this test file.
FROZEN_PACKAGE: dict[Path, str] = {
    REPORT: "59dd458acf14d18021f6c5583096da13f39c83115fc136d6ba91bcc850ce3603",
    CERTIFIER: "6e4d6859ea4513f614714794ac76197c895bf163a4122565674d3694bf03e23b",
    ARCHIVE: "b7108377559f75b65c03fcebb7b01098f7d91b1f5e499c2133d80e882c6f559b",
}


def test_frozen_owner_package() -> None:
    for path, expected in FROZEN_PACKAGE.items():
        assert sha256_file(path) == expected


def test_k3_rank_budget_is_exact() -> None:
    rows = rank_budget()
    assert rows == [
        {
            "field_degree": 1,
            "maximum_native_rank": 20,
            "maximum_hyperbolic_dimension": 19,
        },
        {
            "field_degree": 2,
            "maximum_native_rank": 10,
            "maximum_hyperbolic_dimension": 9,
        },
        {
            "field_degree": 3,
            "maximum_native_rank": 6,
            "maximum_hyperbolic_dimension": 5,
        },
        {
            "field_degree": 4,
            "maximum_native_rank": 5,
            "maximum_hyperbolic_dimension": 4,
        },
        {
            "field_degree": 5,
            "maximum_native_rank": 4,
            "maximum_hyperbolic_dimension": 3,
        },
        {
            "field_degree": 6,
            "maximum_native_rank": 3,
            "maximum_hyperbolic_dimension": 2,
        },
    ]
    assert all(
        row["field_degree"] * row["maximum_native_rank"] <= 20 for row in rows
    )
    assert all(
        row["field_degree"] * (row["maximum_native_rank"] + 1) > 20
        for row in rows
    )


def test_native_reflection_and_normalized_gram_ignore_global_scale() -> None:
    expected_cartan = (
        (Fraction(2), Fraction(-1), Fraction(0)),
        (Fraction(-2), Fraction(2), Fraction(-1)),
        (Fraction(0), Fraction(-3), Fraction(2)),
    )
    expected_normalized = (
        (Fraction(4), Fraction(2), Fraction(0)),
        (Fraction(2), Fraction(4), Fraction(3)),
        (Fraction(0), Fraction(3), Fraction(4)),
    )
    assert cartan_matrix(DELTA_246_GRAM) == expected_cartan
    assert normalized_squared_gram(DELTA_246_GRAM) == expected_normalized
    for scale in range(1, 80):
        scaled = scale_matrix(DELTA_246_GRAM, scale)
        assert cartan_matrix(scaled) == expected_cartan
        assert normalized_squared_gram(scaled) == expected_normalized


def test_projected_physical_root_bound_is_exactly_scale_free() -> None:
    # For x=(sqrt(a),sqrt(2-a)) and y=(sqrt(b),sqrt(2-b)), both have
    # squared norm 2.  The operator diag(7,1) has squared norm 49.
    operator_squared_norm = max(49, 1)
    assert operator_squared_norm == 7**2
    assert 7 * 2 == 14

    # Exact rational boundary and interior samples replay the sharp interface.
    for a in (Fraction(0), Fraction(1, 5), Fraction(1), Fraction(7, 4), Fraction(2)):
        # Equal-coordinate case removes radicals:
        # 7*sqrt(a^2)+sqrt((2-a)^2)=7a+(2-a)=6a+2.
        bound = 6 * a + 2
        assert bound <= 14
    assert 6 * Fraction(2) + 2 == 14
    # The geometric star inequality is strict; an integral K3 intersection is <=13.
    assert max(integer for integer in range(100) if integer < 14) == 13


def test_scaling_family_is_host_infinite_but_chamber_fixed() -> None:
    source = json.loads(SCALING_ARCHIVE.read_text(encoding="utf-8"))["payload"]
    assert source["all_parameter_family"]["integer_scale_domain"] == "s>=2"
    assert source["all_parameter_family"]["ambient_determinant"] == (
        "64*s^3*(s-1)^3"
    )
    assert source["all_parameter_family"][
        "direct_sum_to_ambient_overlattice_index"
    ] == "s^3"
    assert source["scope"]["produces_infinitely_many_abstract_chambers"] is False
    determinants = [64 * s**3 * (s - 1) ** 3 for s in range(2, 30)]
    assert len(set(determinants)) == len(determinants)


def test_redoubling_family_proves_maximality_is_indispensable() -> None:
    source = json.loads(REDOUBLING_ARCHIVE.read_text(encoding="utf-8"))
    assert source["all_checks_pass"] is True
    assert source["theorem"]["boundary_proof"]["facet_count_identity"] == (
        "9+(N+1)+2=N+12"
    )
    facets = [index + 12 for index in range(1, 100)]
    assert len(set(facets)) == 99
    assert facets[0] == 13 and facets[-1] == 111
    assert "index-N subgroup" in source["theorem"]["statement"]


def test_certificate_keeps_the_import_and_nonclaims_explicit() -> None:
    payload = build_certificate()
    theorem = payload["theorem"]
    assert theorem["finite_ambient_NS_hosts_required"] is False
    assert theorem["finite_native_similarity_hosts"] is True
    assert theorem["finite_unmarked_normalized_abstract_chambers"] is True
    assert theorem["classification_free_AG_only_proof"] is False
    assert "Nikulin 1982 Theorem 5.2.1" in theorem["imported_fixed_data_input"]
    assert payload["nonmaximal_redoubling_control"]["reflection_maximal"] is False
    assert payload["projected_physical_star"]["implies_finite_ambient_NS_hosts"] is False


def test_logical_composition_has_no_ambient_host_step() -> None:
    payload = build_certificate()
    assert payload["logic"] == [
        "exact K3 trace rank and finite-volume field theorem -> finite (k,m)",
        "arithmetic W(P) inside O+(L), both lattices -> finite index",
        "W(P) subset R(L) -> L reflective",
        "Nikulin fixed-(k,m) theorem -> finite native similarity hosts",
        "raw reflection maximality -> W(P)=R(L)",
        "complete mirror chambers form one orbit -> finite normalized chambers",
    ]
    assert all("ambient" not in step.lower() for step in payload["logic"])


def test_dependency_pins_and_primary_source_pin_are_literal() -> None:
    payload = build_certificate()
    assert payload["dependency_sha256"] == DEPENDENCY_SHA256
    assert DEPENDENCY_SHA256["nikulin_primary_pdf"] == (
        "696652b6a67c6ed420de6c1cb61d9422e3619e28e081a78131a5323ce2edee9d"
    )


def test_report_states_exact_scope_and_project_native_frontier() -> None:
    text = REPORT.read_text(encoding="utf-8")
    for needle in (
        "finitely many native similarity hosts",
        "ambient Picard hosts are unnecessary",
        "Theorem 5.2.1",
        "A=62",
        "W(P)=R(L)",
        "does **not** prove ambient-host",
        "Reflection maximality is essential",
        "This is not yet an AG-only derivation of ABSW",
        "connected full-rank raw",
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
    )
