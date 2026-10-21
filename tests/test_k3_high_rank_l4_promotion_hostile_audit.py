import json

from scripts.certify_k3_high_rank_l4_promotion_hostile_audit import (
    OUTPUT_PATH,
    build_artifact,
)


def _artifact() -> dict:
    return json.loads(OUTPUT_PATH.read_text())


def test_frozen_hostile_audit_rebuilds_exactly() -> None:
    frozen = _artifact()
    assert build_artifact(frozen["source_check"]) == frozen


def test_carrier_and_return_pass_but_independent_ag_pushback_fails() -> None:
    frozen = _artifact()
    score = frozen["strict_score"]
    assert score["C"].startswith("PASS")
    assert score["Car"].startswith("PASS")
    assert score["Return"].startswith("PASS")
    assert score["AG"].startswith("FAIL_INDEPENDENT_EXHAUSTION")
    assert score["Pushback"].startswith("FAIL_NEW_THEOREM")
    assert frozen["verdict"]["strict_l4"] is False
    assert frozen["verdict"]["promotion"] == "rejected"
    assert frozen["upstream"]["mathematical_classification_valid"] is True


def test_raw_carrier_nonemptiness_is_not_finite_volume() -> None:
    frozen = _artifact()
    warning = frozen["correct_finite_volume_iff"]
    example = warning["counterexample"]
    assert example["lattice"] == "U+E7(-1)+A2(-1)"
    assert example["rank"] == 11
    assert example["signature"] == [1, 10]
    assert example["determinant"] == 6
    assert example["even"] is True
    assert example["primitive_k3_embedding"][
        "coordinate_sublattices_primitive"
    ] is True
    assert example["exact_picard_k3_carrier_nonempty"] is True
    assert example["full_root_chamber_returns_as_nef"] is True
    assert example["in_nineteen_finite_volume_rows"] is False
    assert example["finite_volume"] is False
    assert "finite Aut" in warning["correct_statement"]


def test_source_dependency_and_literal_audit_roles_are_separated() -> None:
    frozen = _artifact()
    dependency = frozen["classification_dependency"]
    assert "imports" in dependency["ahl_theorem_2_12_role"]
    assert "Nikulin--Vinberg" in dependency["roulleau_role"]
    assert "does not prove the list exhaustive" in dependency["literal_vinberg_role"]
    assert frozen["checks"][
        "nineteen_zero_exhaustion_is_independent_ag_native"
    ] is False
    assert frozen["checks"]["new_polytope_theorem_from_ag_pushback"] is False


def test_compact_and_noncompact_loci_are_not_conflated() -> None:
    frozen = _artifact()
    separation = frozen["compact_noncompact_separation"]
    assert separation["compact_finite_volume_locus_rank_11_20"] == "empty"
    assert separation["compact_reason_is_ag_native"] is False
    noncompact = separation["finite_volume_noncompact_locus"]
    assert noncompact == {
        "hyperbolic_dimensions": [10, 18],
        "no_H19_type": True,
        "total_ideal_vertices_in_literal_replays": 1156,
        "type_count": 19,
    }


def test_exact_missing_l4_gate_is_frozen() -> None:
    gate = _artifact()["missing_l4_gate"]
    assert "does not invoke" in gate["required"]
    assert "followed by a polytope conclusion" in gate["required"]
    assert "period existence" in gate["not_sufficient"]
    assert "literal Vinberg verification of the already listed positives" in gate[
        "not_sufficient"
    ]
