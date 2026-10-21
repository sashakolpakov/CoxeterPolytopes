from pathlib import Path

from scripts.certify_rank12_five_chamber_simultaneous_carrier import (
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results" / "rank12_five_chamber_simultaneous_carrier.md"


def test_five_chamber_carrier_is_green() -> None:
    certificate = build_certificate()
    assert certificate["all_checks_pass"]
    assert certificate["distinct_hyperbolic_polytope_count"] == 5
    assert certificate["newly_geometric_slots"] == ["H10/16", "H10/45"]


def test_all_five_vinberg_gates_and_coinvariants() -> None:
    certificate = build_certificate()
    assert [row["polytope"] for row in certificate["records"]] == [
        "H10/12",
        "H10/15",
        "H10/16",
        "H10/38",
        "H10/45",
    ]
    assert all(row["checks"]["literal_vinberg_gate"] for row in certificate["records"])
    assert all(
        row["coinvariant_gram"] == [[-4]] for row in certificate["records"]
    )


def test_finite_index_tree() -> None:
    certificate = build_certificate()
    assert [
        row["sublattice_index_in_UE8A1x2"] for row in certificate["records"]
    ] == [1, 2, 2, 4, 4]
    criterion = certificate["finite_index_divisor_carrier_criterion"]
    assert criterion["common_invariant_ample_class_exists"]
    assert criterion["all_five_divisor_chambers_can_be_oriented_around_that_class"]
    assert criterion["simultaneous_carrier_conclusion"]
    assert all(criterion["verified_hypotheses"].values())


def test_complete_rank_eleven_chamber_action_fingerprints() -> None:
    certificate = build_certificate()
    census = certificate["complete_rank_eleven_chamber_involution_census"]
    assert census["all_checks_pass"]
    assert {
        key: [row["class_index"] for row in family["classes"]]
        for key, family in census["families"].items()
    } == {
        "R12_UE8A1x2": [2],
        "R12_UD8A1x2": [2, 3],
        "R12_UD4x2A1x2": [4, 12],
    }
    assert all(
        row["checks"]["unique_complete_chamber_group_fingerprint"]
        and row["checks"]["fingerprint_forces_archived_chamber_class"]
        for row in certificate["records"]
    )


def test_source_lattice_identification_is_exact_not_histogram_only() -> None:
    certificate = build_certificate()
    for row in certificate["records"]:
        witness = row["source_lattice_discriminant_form_isometry"]
        assert witness["bijective"]
        assert witness["quadratic_form_preserved_mod_2Z"]
        assert witness["bilinear_form_preserved_mod_Z"]
        assert row["checks"]["source_lattice_identified_by_indefinite_uniqueness"]


def test_target_specific_pair_handoffs_are_only_corroboration() -> None:
    certificate = build_certificate()
    corroborated = {
        row["polytope"]: row["target_specific_exact_pair_corroboration"]
        for row in certificate["records"]
        if row["target_specific_exact_pair_corroboration"]["available"]
    }
    assert set(corroborated) == {"H10/16", "H10/45"}
    assert all(row["passes"] for row in corroborated.values())
    assert all("corroboration only" in row["role"] for row in corroborated.values())


def test_report_states_exact_simultaneous_scope() -> None:
    report = " ".join(REPORT.read_text(encoding="utf-8").lower().split())
    assert "one order-four k3 action carrying five" in report
    assert "finite-volume noncompact" in report
    assert "other four are effective finite-index" in report
    assert "divisor-subsystem" in report
    assert "not new abstract diagrams" in report
    assert "worldwide priority" in report
    assert "elliptic rank-9 edge diagrams" in report
