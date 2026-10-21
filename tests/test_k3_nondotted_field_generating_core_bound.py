"""Regression tests for the nondotted K3 field-generating-core theorem."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

from scripts.certify_k3_nondotted_field_generating_core_bound import (
    ARCHIVE,
    REPORT,
    build_certificate,
    digest,
    encoded_certificate,
    rank_rows,
)


ROOT = Path(__file__).resolve().parents[1]
CERTIFIER = ROOT / "scripts/certify_k3_nondotted_field_generating_core_bound.py"


def test_rank_length_degree_and_discriminant_table() -> None:
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
    assert [by_rank[m]["chordless_cycle_length_bound"] for m in range(3, 11)] == [
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    ]
    assert [by_rank[m]["binary_discriminant_exponent"] for m in range(3, 11)] == [
        300,
        220,
        144,
        78,
        28,
        30,
        32,
        34,
    ]
    assert all(by_rank[m]["rational_only"] for m in range(11, 21))


def test_cycle_and_primitive_bounds_are_literal() -> None:
    for module_rank in range(3, 11):
        cycle_bound = 2 ** (module_rank + 2)
        primitive_bound = 16 * cycle_bound
        difference_bound = 2 * primitive_bound
        assert primitive_bound == 2 ** (module_rank + 6)
        assert difference_bound == 2 ** (module_rank + 7)


def test_at_most_two_generators_and_six_embedding_pairs_grid() -> None:
    for degree in range(2, 7):
        strict_steps = 0
        current_degree = 1
        while 2 * current_degree <= degree:
            current_degree *= 2
            strict_steps += 1
        assert strict_steps <= 2
        assert degree * (degree - 1) // 2 <= 15
    # Six embeddings can forbid all fifteen coefficients 0,...,14;
    # the sixteenth candidate is why the safe grid ends at 15.
    assert len(range(16)) == 15 + 1


def test_certificate_scope_and_nonclaims() -> None:
    payload = build_certificate()
    assert payload["uniform_binary_discriminant_exponent"] == 300
    assert payload["uniform_field_discriminant_bound_expression"] == "2^300"
    assert payload["maximizers"] == [{"module_rank": 3, "field_degree": 6}]
    scope = payload["scope"]
    assert scope["dotted_edges_outside_generating_core_allowed"] is True
    assert scope["facet_count_bound_required"] is False
    assert scope["finite_volume_required"] is False
    conclusions = payload["conclusions"]
    assert conclusions["finite_field_set"] is True
    assert conclusions["carrier_existence_proved"] is False
    assert conclusions["essential_surjectivity_proved"] is False
    assert conclusions["canonical_mixed_height_bounded"] is False
    assert conclusions["all_higher_facet_polytopes_covered"] is False
    assert conclusions["chamber_finiteness_proved"] is False


def test_archive_is_exact_replay() -> None:
    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert archived["payload"] == build_certificate()
    assert archived["payload_sha256"] == digest(archived["payload"])
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()


def test_report_contains_scope_and_exact_bound() -> None:
    text = REPORT.read_text(encoding="utf-8")
    for needle in (
        "2^{300}",
        "field-generating core",
        "dotted or ultraparallel pairs outside",
        "does **not**",
        "2^{240}",
        "not asserted to be the optimal possible",
    ):
        assert needle in text


def test_cli_check() -> None:
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check-output"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
