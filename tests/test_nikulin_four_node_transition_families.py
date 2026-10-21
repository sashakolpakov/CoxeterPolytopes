from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/certify_nikulin_four_node_transition_families.py"
ARCHIVE = ROOT / "data/enumerations/nikulin_four_node_transition_families.json"
REPORT = ROOT / "results/nikulin_four_node_transition_families_exact_audit.md"


def load_module():
    spec = importlib.util.spec_from_file_location("nikulin_four_node_families", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_archive() -> dict:
    return json.loads(ARCHIVE.read_text(encoding="utf-8"))


def test_frozen_archive_is_an_exact_replay() -> None:
    module = load_module()
    archived = load_archive()
    assert archived == module.build_certificate()
    assert archived["all_checks_pass"] is True


def test_support_census_has_seven_orbits_in_five_schematics() -> None:
    topology = load_archive()["topology_census"]
    assert topology["raw_masks"] == 32
    assert topology["admissible_raw_masks"] == 17
    assert topology["marked_support_orbits"] == 7
    assert topology["schematic_families"] == 5
    rows = {
        tuple(row["canonical_support"]): (row["family"], row["specialization"])
        for row in topology["orbits"]
    }
    assert rows[("13", "34")] == ("Gamma3", "s=2 (the drawn 14 edge is absent)")
    assert rows[("13", "14")] == ("Gamma4", "k=2 (the drawn 23 edge is absent)")
    assert set(family for family, _ in rows.values()) == {
        "Gamma1",
        "Gamma2",
        "Gamma3",
        "Gamma4",
        "Gamma5",
    }


def test_every_support_orbit_has_exact_hyperbolic_Q_witness() -> None:
    rows = load_archive()["topology_census"]["orbits"]
    assert [row["rational_exact_witness"]["determinant"] for row in rows] == [
        "-28",
        "-27",
        "-19",
        "-43",
        "-35",
        "-28",
        "-60",
    ]
    for row in rows:
        witness = row["rational_exact_witness"]
        assert witness["inertia"] == {"positive": 1, "negative": 3, "zero": 0}
        assert all(
            deletion["negative_definite"]
            for deletion in witness["endpoint_deletions"].values()
        )
        assert witness["v_arithmetic"]["ground_field"] == "Q"
        assert witness["v_arithmetic"]["cyclic_products_integral"] is True


def test_parameter_grammar_is_exhaustive_and_records_k2_redundancy() -> None:
    grammar = load_archive()["parameter_grammar_census"]
    assert grammar["raw_assignments"] == 7776
    assert grammar["admissible_raw_assignments"] == 249
    assert grammar["marked_label_orbits"] == 73
    assert grammar["marked_label_orbit_counts"] == {
        "Gamma1": 9,
        "Gamma2": 8,
        "Gamma3": 13,
        "Gamma4": 28,
        "Gamma5": 15,
    }
    assert grammar["coverage_failures"] == 0
    assert grammar["cross_family_ambiguities"] == 0
    assert grammar["rank_three_path_lemma"]["solutions_for_p_q_at_least_3_up_to_order"] == [
        [3, 3],
        [3, 4],
        [3, 5],
    ]
    assert all(
        row["one_unlabelled_marked_graph"]
        for row in grammar["Gamma4_k_equals_2_redundancies"]
    )


def test_symbolic_determinants_and_field_inheritance_gates() -> None:
    cert = load_archive()
    assert set(cert["symbolic_gram_matrices"]) == {
        "Gamma1",
        "Gamma2",
        "Gamma3",
        "Gamma4",
        "Gamma5",
    }
    assert all(
        row["determinant_identity_verified"]
        for row in cert["symbolic_gram_matrices"].values()
    )
    inheritance = cert["ground_field_inheritance"]
    switching = inheritance["sign_switching_finite_audit"]
    assert switching == {
        "cycle_consistent_assignments": 304,
        "edge_sign_assignments": 624,
        "four_vertex_connected_supports": 38,
        "general_lemma": "fixing every cyclic product makes the conjugate principal Gram diagonally congruent by vertex signs on each connected component",
        "switching_failures": 0,
    }
    assert inheritance["marked_u_regimes"]["u_equals_2"].endswith("hence K=Q")
    assert inheritance["quasi_lanner_and_cusp"]["compactness_used"] is False


def test_report_states_precise_verdict_and_has_clean_control_bytes() -> None:
    raw = REPORT.read_bytes()
    assert all(byte not in raw for byte in (0, 8, 11, 13))
    text = raw.decode("utf-8")
    assert "PASS after making the word “type” precise" in text
    assert "seven" in text and "five schematic" in text
    assert "Gamma3" in text and "s=2" in text
    assert "Gamma4" in text and "k=2" in text
    assert "K_I = K" in text
    assert "u=2:  the marked 2-node packet is parabolic" in text
    assert "does **not** certify the later numerical interval extrema" in text
