import json

import sympy as sp

from scripts.certify_integral_packet_dimension_bridge import (
    ARCHIVE,
    build_archive,
    exact_matrix_inertia,
    star_gram,
)


def test_star_counterfamily_is_lorentzian_of_unbounded_rank():
    for leaves in range(5, 13):
        matrix = star_gram(leaves)
        assert matrix.rank() == leaves + 1
        assert exact_matrix_inertia(matrix) == (1, leaves, 0)


def test_all_exact_four_node_principal_packets_are_elliptic():
    for leaves in range(3, 13):
        matrix = star_gram(leaves)
        no_center = matrix.extract(range(1, 5), range(1, 5)) if leaves >= 4 else None
        if no_center is not None:
            assert exact_matrix_inertia(no_center) == (0, 4, 0)
        with_center_indices = [0, 1, 2, 3]
        with_center = matrix.extract(with_center_indices, with_center_indices)
        assert exact_matrix_inertia(with_center) == (0, 4, 0)


def test_at_most_four_type_ledger_includes_disconnected_edgeless_types():
    archive = build_archive()
    types = archive["star_counterfamily"]["at_most_four_node_types"]
    descriptors = {
        (row["kind"], row["order"], row.get("star_leaves")) for row in types
    }
    assert descriptors == {
        *(('edgeless', order, None) for order in range(1, 5)),
        *(('star', leaves + 1, leaves) for leaves in range(1, 4)),
    }
    assert all(row["inertia_positive_negative_zero"][0] == 0 for row in types)
    assert all(row["inertia_positive_negative_zero"][2] == 0 for row in types)


def test_global_vertex_deletion_detects_the_unbounded_tail():
    assert exact_matrix_inertia(star_gram(4)) == (0, 4, 1)
    for leaves in range(6, 13):
        assert exact_matrix_inertia(star_gram(leaves - 1))[0] == 1


def test_characteristic_polynomial_formula():
    lam = sp.symbols("lambda")
    for leaves in range(1, 13):
        actual = star_gram(leaves).charpoly(lam).as_expr()
        expected = (lam + 2) ** (leaves - 1) * ((lam + 2) ** 2 - leaves)
        assert sp.expand(actual - expected) == 0


def test_archive_replay_and_scope_flags():
    expected = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert build_archive() == expected
    star = expected["star_counterfamily"]
    assert star["all_four_node_principal_submatrices_negative_definite"]
    assert star["lorentzian_for_all_m_at_least_5"]
    assert "only m=5" in star["scope"]
    assert "four-node type finiteness alone cannot supply this" in expected[
        "open_global_gate"
    ]
