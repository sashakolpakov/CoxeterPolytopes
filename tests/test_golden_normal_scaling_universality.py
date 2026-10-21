import json

from coxeter_polytopes.exact.golden import Golden, PHI
from scripts.certify_golden_normal_scaling_universality import (
    OUTPUT,
    build,
    verify,
)


def test_all_frozen_m35_graphs_have_one_normal_scaling_class() -> None:
    payload = build()
    rows = payload["golden_lanner_seeds"]
    assert len(rows) == 12
    assert {row["name"] for row in rows} == {
        "(2,5,5)",
        "(3,3,5)",
        "(3,5,5)",
        "(5,5,5)",
        "[3,5,3]",
        "[5,3,5]",
        "D4(5,3,3)",
        "C4(one 5)",
        "C4(opposite 5s)",
        "[5,3,3,3]",
        "[5,3,3,5]",
        "D5(terminal 5)",
    }
    for row in rows:
        assert row["checks"]["m35_graph_connected"] is True
        assert row["checks"]["twist_audit_passes"] is True
        assert row["normal_scaling_classes_mod_oriented_vertex_units"] == 1

    bugaenko = payload["bugaenko_H8"]
    assert bugaenko["m35_graph"] == "spanning_tree"
    assert bugaenko["normal_scaling_classes_mod_oriented_vertex_units"] == 1
    assert bugaenko["checks"]["unit_twist_audit_passes"] is True


def test_bugaenko_rank_forces_the_unique_dotted_entry() -> None:
    row = build()["bugaenko_H8"]
    assert Golden.from_json(row["forced_dotted_entry"]) == -2 * PHI
    minors = row["rank_minor_factorizations"]
    assert Golden.from_json(minors["delete_facet_0"]["secondary_root"]) != (
        Golden.from_json(minors["delete_facet_1"]["secondary_root"])
    )
    assert row["checks"]["two_rank_minor_factorizations_exact"] is True
    assert row["checks"]["secondary_0_fails_minor_1"] is True
    assert row["checks"]["secondary_1_fails_minor_0"] is True
    assert row["checks"]["forced_gram_has_rank_nine"] is True
    assert row["checks"]["forced_gram_matches_archive"] is True


def test_bugaenko_wall_span_and_integral_saturation_are_unique() -> None:
    row = build()["bugaenko_H8"]
    determinant = Golden.from_json(row["first_nine_wall_determinant"])
    assert determinant.norm() == -1
    assert row["wall_lattice_discriminant_module_order"] == 4
    assert row["reflection_stable_integral_O_overclosures"]["count"] == 1
    assert row["O_stable_even_trace_overclosures"]["count"] == 1
    assert row["O_stable_even_trace_overclosures"]["nonzero_orbit"] == [
        [1, 0],
        [0, 1],
        [1, 1],
    ]
    assert row["checks"]["two_is_inert_so_discriminant_module_is_simple"]
    assert row["checks"]["full_dual_is_not_integral"]
    assert row["checks"]["trace_discriminant_action_has_no_invariant_line"]


def test_hostile_nonsquare_bridge_blocks_equal_norm() -> None:
    row = build()["hostile_boundaries"]["nonsquare_dotted_bridge"]
    assert Golden.from_json(row["field_product_p"]) == Golden(3, 1)
    assert row["field_norm"] == 11
    assert all(row["checks"].values())
    assert "equal root norm would require" in row["geometric_conclusion"]


def test_angle_only_infinity_label_has_infinite_exact_family() -> None:
    row = build()["hostile_boundaries"][
        "angle_only_infinity_label_has_infinite_realizations"
    ]
    assert row["formula"] == "p_k=3+phi^(2k), k>=1"
    terms = row["first_twelve_exact_terms"]
    assert len(terms) == 12
    values = [Golden.from_json(term["p_k"]) for term in terms]
    assert all(left < right for left, right in zip(values, values[1:]))
    assert all(term["checks"]["identity_ultraparallel"] for term in terms)
    assert all(term["checks"]["conjugate_admissible"] for term in terms)


def test_scope_is_finite_compiler_not_carrier_universality() -> None:
    scope = build()["theorem_scope"]
    assert "for each fixed wall reflection representation" in scope[
        "normal_scaling_result"
    ]
    assert "not Missing-Lemma-U carrier universality" in scope["not_proved"]
    assert "R and R^vee" in scope["saturation_result"]


def test_frozen_archive_replay() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    summary = verify(archived)
    assert summary["lanner_seed_count"] == 12
    assert summary["bugaenko_scaling_classes"] == 1
    assert Golden.from_json(summary["forced_dotted_entry"]) == -2 * PHI
