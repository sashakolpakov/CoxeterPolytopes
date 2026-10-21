from __future__ import annotations

import json
import math

from scripts.certify_dessin_native_field_effectivity_h11 import (
    DEFAULT_OUTPUT,
    ROOT,
    build,
)


def test_h11_rank_locus_and_rejection_are_exact() -> None:
    certificate = build()
    assert json.loads(DEFAULT_OUTPUT.read_text()) == certificate
    locus = certificate["projective_rank_locus"]
    assert locus["geometric_points"] == [[1, 1], [1, -1]]
    assert locus["rejected"]["point"] == [1, -1]
    assert locus["accepted"]["point"] == [1, 1]
    assert locus["accepted"]["reduced_coordinate_algebra"] == "Q[b]/(b-1) = Q"


def test_h11_galois_and_vinberg_fields_are_Q() -> None:
    certificate = build()
    accepted = certificate["projective_rank_locus"]["accepted"]
    field = certificate["exact_gram_field"]
    assert accepted["galois_orbit_size"] == 1
    assert accepted["field_of_moduli"] == accepted["field_of_definition"] == "Q"
    assert field["cartan_entry_count"] == 26**2
    assert field["length_two_cyclic_product_count"] == 26 * 25 // 2
    assert field["all_cartan_entries_rational"]
    assert field["all_cyclic_products_rational"]
    assert field["recovered_Vinberg_field"] == "Q"


def test_scope_does_not_promote_one_point_to_global_finiteness() -> None:
    certificate = build()
    assert certificate["scope"]["global_field_finiteness"].startswith("open")
    report = (ROOT / "results/dessin_native_field_effectivity_theorem.md").read_text()
    assert "Prelisting finitely many" in report
    assert "zero-dimensional scheme-theoretic image" in report
    assert "positive-dimensional fibres" in report.lower()
    assert "does **not** contain a list of coefficient or Vinberg fields" in report


def test_report_has_no_corrupt_control_bytes() -> None:
    payload = (ROOT / "results/dessin_native_field_effectivity_theorem.md").read_bytes()
    forbidden = [
        (offset, byte)
        for offset, byte in enumerate(payload)
        if byte < 32 and byte not in (9, 10)
    ]
    assert forbidden == []
    assert b"\r" not in payload
    assert b"\x08" not in payload


def test_intrinsic_k3_angle_alphabet_and_scope() -> None:
    def phi(value: int) -> int:
        return sum(math.gcd(value, candidate) == 1 for candidate in range(1, value + 1))

    labels = [2] + [value for value in range(3, 100) if phi(value) <= 4]
    assert labels == [2, 3, 4, 5, 6, 8, 10, 12]
    report = (ROOT / "results/dessin_native_field_effectivity_theorem.md").read_text()
    assert r"This includes cusped \(n=3\)" in report
    assert r"It fails in \(n=2\)" in report
    assert "carrier-relative conclusions" in report
    assert "This does not bound dotted" in report
