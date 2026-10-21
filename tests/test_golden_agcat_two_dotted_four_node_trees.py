from coxeter_polytopes.classification.golden_agcat_manifest import (
    golden_finite_pair_values,
)
from coxeter_polytopes.classification.golden_four_node_packets import (
    golden_four_node_determinant,
)
from coxeter_polytopes.classification.golden_two_dotted_tree_packets import (
    conjugate_sorted_dotted_values,
    count_adjacent_path_pairs,
    count_disjoint_path_pairs,
    count_star_pairs,
    golden_two_dotted_tree_census,
)


def test_frozen_two_dotted_tree_counts() -> None:
    census = golden_two_dotted_tree_census()
    assert census.canonical_before_gram == 471_865_955
    assert census.retained_count == 97_172_716
    assert census.star_retained_count == 25_023_281
    assert census.disjoint_path_retained_count == 81
    assert census.adjacent_path_retained_count == 72_149_354


def test_two_pointer_counts_match_small_brute_force() -> None:
    values = conjugate_sorted_dotted_values()[:48]
    for row in golden_finite_pair_values()[1:]:
        finite = row.product
        star = sum(
            (first + second + finite - 4).sign(1) < 0
            for index, first in enumerate(values)
            for second in values[index:]
        )
        adjacent = sum(
            (4 * adjacent - (4 - remote) * (4 - finite)).sign(1) < 0
            for remote in values
            for adjacent in values
        )
        disjoint = sum(
            ((first - 4) * (second - 4) - 4 * finite).sign(0) < 0
            and ((4 - first) * (4 - second) - 4 * finite).sign(1) > 0
            for index, first in enumerate(values)
            for second in values[index:]
        )
        assert count_star_pairs(values, finite) == star
        assert count_disjoint_path_pairs(values, finite) == disjoint
        assert count_adjacent_path_pairs(values, finite) == adjacent


def test_tree_determinant_formulas() -> None:
    values = conjugate_sorted_dotted_values()
    x, y = values[100], values[500]
    z = golden_finite_pair_values()[2].product
    zero = (0, 0)
    raw_x = int(x.a), int(x.b)
    raw_y = int(y.a), int(y.b)
    raw_z = int(z.a), int(z.b)

    # Star edges 01,02,03.
    star = golden_four_node_determinant((raw_x, raw_y, raw_z, zero, zero, zero))
    expected_star = 16 - 4 * (x + y + z)
    assert star == (int(expected_star.a), int(expected_star.b))

    # Path edges 01=x, 12=y, 23=z; 01 and 23 are opposite.
    path = golden_four_node_determinant((raw_x, zero, zero, raw_y, zero, raw_z))
    expected_path = 16 - 4 * (x + y + z) + x * z
    assert path == (int(expected_path.a), int(expected_path.b))
