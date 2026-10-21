from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "data/enumerations/rational_i1n_j_moduli_carriers.json"
REPORT = ROOT / "results/rational_i1n_j_moduli_carriers.md"
LEDGER = ROOT / "results/higher_dimensional_polytope_lineage_ledger.md"


def _payload() -> dict:
    return json.loads(CERTIFICATE.read_text(encoding="utf-8"))


def test_two_explicit_geometrically_integral_fibers_have_distinct_j() -> None:
    payload = _payload()
    assert payload["all_checks_pass"]
    assert all(payload["checks"].values())
    fibers = payload["explicit_fibers"]
    assert [(row["parameter"], row["j_invariant"]) for row in fibers] == [
        (0, "0"),
        (1, "6912/31"),
    ]
    for fiber in fibers:
        assert all(fiber["checks"].values())
        assert fiber["branch_data"]["C_Q_contacts"] == [1, 11]
        assert fiber["branch_data"]["Q_multiplicity_sequence_through_p15"] == [
            2,
            *([1] * 14),
        ]
        assert fiber["branch_data"]["S_multiplicity_sequence_through_p16"] == [
            *([2] * 7),
            *([1] * 9),
        ]
        assert all(
            row["absolutely_irreducible"]
            for row in fiber["absolute_irreducibility"].values()
        )


def test_complete_P12_through_P16_nef_carriers_and_cusps_are_frozen() -> None:
    payload = _payload()
    expected = {
        "12": (14, 35, 2),
        "13": (15, 43, 3),
        "14": (17, 89, 5),
        "15": (18, 139, 6),
        "16": (20, 313, 12),
    }
    assert set(payload["dimensions"]) == set(expected)
    for n, (facets, finite, cusps) in expected.items():
        record = payload["dimensions"][n]
        assert all(record["checks"].values())
        assert (record["facet_count"], record["finite_vertices"]) == (
            facets,
            finite,
        )
        assert record["cusps"]["ideal_ray_count"] == cusps
        assert len(record["cusps"]["two_universal_semiample_cusps"]) == 2
        assert record["unique_genus_one_extremal_curve"]["arithmetic_genus"] == 1


def test_second_explicit_H16_carrier_and_priority_scope_are_reported() -> None:
    payload = _payload()
    conclusion = payload["positive_conclusion"]
    assert "a=1, n=16" in conclusion["strongest_new_explicit_carrier"]
    assert conclusion["abstract_polytope_novelty"] is False
    assert conclusion["divisor_carrier_novelty_in_this_project"] is True
    report = REPORT.read_text(encoding="utf-8")
    ledger = LEDGER.read_text(encoding="utf-8")
    assert "second explicit `H16/20` carrier" in report
    assert "non-isomorphic" in report
    assert "RAT-J1-P16-H16-20" in ledger
    assert "6912/31" in ledger
