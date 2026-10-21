import json
from itertools import combinations

from coxeter_polytopes.combinatorics.facet_completion import (
    FULL_VERTEX_MASK,
    generated_face_table,
    minimal_nonfaces_from_table,
)
from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    FINITE_LABELS,
    audit_h9_14_coxeter_labels,
    five_cycle_color_orbit_codes,
    is_finite_coxeter_labeling,
    is_lanner_coxeter_labeling,
    iter_product_orthogonal_label_representatives,
    lanner_isomorphism_codes,
    lanner_placement_codes,
    minimal_nonface_component_partition,
    minimum_cross_component_nonright_edges,
    pair,
)
from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.geometric.golden_lanner import (
    golden_block_product_counterexample,
)
from scripts.certify_h9_14_product_coxeter_labels import OUTPUT, verify


def _labels_from_code(order: int, code: tuple[int, ...]) -> dict[tuple[int, int], int]:
    return dict(zip(combinations(range(order), 2), code, strict=True))


def _vertices(mask: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(14) if mask & (1 << vertex))


def _product_nonfaces(selected: list[int]) -> tuple[tuple[int, ...], ...]:
    by_order = minimal_nonfaces_from_table(generated_face_table(selected))
    return tuple(
        _vertices(mask)
        for order in sorted(by_order)
        for mask in by_order[order]
    )


def test_full_domain_lanner_tables_have_exact_classical_sizes() -> None:
    assert FINITE_LABELS == (2, 3, 4, 5)
    assert [len(lanner_placement_codes(order)) for order in (4, 5)] == [108, 420]
    assert [len(lanner_isomorphism_codes(order)) for order in (4, 5)] == [9, 5]
    for order in (4, 5):
        for code in lanner_placement_codes(order):
            labels = _labels_from_code(order, code)
            assert is_lanner_coxeter_labeling(range(order), labels)
            assert all(
                is_finite_coxeter_labeling(
                    (vertex for vertex in range(order) if vertex != deleted),
                    labels,
                )
                for deleted in range(order)
            )


def test_affine_d4_star_is_not_mistaken_for_a_lanner_diagram() -> None:
    labels = {edge: 2 for edge in combinations(range(5), 2)}
    labels.update({pair(0, leaf): 3 for leaf in range(1, 5)})
    assert all(
        is_finite_coxeter_labeling(
            (vertex for vertex in range(5) if vertex != deleted), labels
        )
        for deleted in range(5)
    )
    assert is_lanner_coxeter_labeling(range(5), labels) is False


def test_product_export_api_lists_one_representative_per_claimed_orbit() -> None:
    refinement = json.loads(
        (OUTPUT.parent / "h9_14_facet_refinement.json").read_text(encoding="utf-8")
    )
    expected = {30: (45, 88), 42: (6120, 86)}
    for index, (expected_count, expected_pair_count) in expected.items():
        selected = refinement["records"][index]["selected_facet_complement_masks"]
        nonfaces = _product_nonfaces(selected)
        assert len(minimal_nonface_component_partition(nonfaces)) == {
            30: 5,
            42: 3,
        }[index]
        assert minimum_cross_component_nonright_edges(nonfaces) == {
            30: 4,
            42: 2,
        }[index]
        representatives = list(iter_product_orthogonal_label_representatives(nonfaces))
        assert len(representatives) == expected_count
        assert len({
            (
                representative.lanner_type_indices,
                representative.cycle_color_code,
            )
            for representative in representatives
        }) == expected_count
        assert all(len(representative.labels) == expected_pair_count for representative in representatives)
        assert all(
            label in FINITE_LABELS
            for representative in (representatives[0], representatives[-1])
            for _left, _right, label in representative.labels
        )
    assert len(five_cycle_color_orbit_codes()) == 136


def test_golden_product_is_a_local_pass_global_inertia_counterexample() -> None:
    counterexample = golden_block_product_counterexample()
    dotted = {
        pair(*face) for face in counterexample.missing_faces if len(face) == 2
    }
    labels = {}
    for left, right in combinations(range(14), 2):
        edge = pair(left, right)
        if edge in dotted:
            continue
        entry = counterexample.field_gram[left][right]
        if entry == Golden(0):
            labels[edge] = 2
        elif entry == Golden(-1):
            labels[edge] = 3
        elif entry == Golden(0, -1):
            labels[edge] = 5
        else:
            raise AssertionError(f"unexpected finite Gram entry {entry}")
    audit = audit_h9_14_coxeter_labels(
        counterexample.vertices, counterexample.missing_faces, labels
    )
    assert audit.passes_discrete_gate is True
    assert counterexample.global_identity_inertia == Inertia(9, 5, 0)
    assert counterexample.passes_global_h9_gram is False


def test_frozen_product_label_certificate_replays() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    summary = verify(manifest)
    assert summary["record_indices"] == [30, 42]
    assert summary["discrete_label_statuses"] == ["sat", "sat"]
    assert summary["orthogonal_join_isomorphism_orbits"] == [45, 6120]
    assert all(
        record["orthogonal_join_subfamily"]["global_gram_obstruction"]["status"]
        == "unsat"
        for record in manifest["records"]
    )
