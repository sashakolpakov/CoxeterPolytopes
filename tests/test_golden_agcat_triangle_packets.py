from math import comb

from coxeter_polytopes.classification.golden_triangle_packets import (
    _golden_square_root_coefficients,
    _pair_mul,
    _pair_sign,
    _pair_square_root,
    golden_element_squareclass,
    golden_triangle_packet_census,
)
from coxeter_polytopes.exact.golden import PHI, Golden


def test_corrected_golden_square_root_solver_regressions():
    for root in (PHI, 2 + PHI, 3 + 2 * PHI, Golden(17, 23)):
        square = root * root
        assert _golden_square_root_coefficients(square) == root
        raw = (int(square.a), int(square.b))
        raw_root = (int(root.a), int(root.b))
        assert _pair_square_root(raw) == raw_root
        assert _pair_mul(raw_root, raw_root) == raw
    assert _golden_square_root_coefficients(PHI) is None


def test_element_squareclass_includes_unit_sign_and_split_ideal_data():
    assert golden_element_squareclass(Golden(1)) != golden_element_squareclass(PHI)
    assert golden_element_squareclass(Golden(1)).conjugate_sign == 1
    assert golden_element_squareclass(PHI).conjugate_sign == -1
    # Every square has the trivial element squareclass.
    for root in (PHI, 2 + PHI, 3 + 2 * PHI, Golden(17, 23)):
        squareclass = golden_element_squareclass(root * root)
        assert squareclass.ideal_parity == ()
        assert squareclass.conjugate_sign == 1


def test_raw_integer_pair_arithmetic_matches_golden_arithmetic():
    values = (Golden(1), PHI, 2 + PHI, Golden(1064, 1718))
    for left in values:
        for right in values:
            product = left * right
            assert _pair_mul(
                (int(left.a), int(left.b)), (int(right.a), int(right.b))
            ) == (int(product.a), int(product.b))
        assert _pair_sign((int(left.a), int(left.b)), 0) == left.sign(0)
        assert _pair_sign((int(left.a), int(left.b)), 1) == left.sign(1)


def test_complete_full_triangle_census_counts_strata_and_digests():
    census = golden_triangle_packet_census()
    assert census.alphabet_count == 6_875
    assert census.zero_edge_value_count == 1
    assert census.nonzero_edge_value_count == 6_874
    assert census.raw_unordered_nonzero_edge_multisets == comb(6_874 + 2, 3)
    assert census.raw_unordered_nonzero_edge_multisets == 54_158_527_500
    assert census.squareclass_bucket_count == 3_661
    assert census.maximum_squareclass_bucket_size == 109
    assert census.compatible_squareclass_bucket_triples == 27_088
    assert census.square_product_edge_multisets == 4_768_141
    assert census.conjugate_positive_packet_count == 2_633_885
    assert census.conjugate_nonpositive_packet_count == 2_134_256
    assert census.selected_positive_packet_count == 0
    assert census.selected_lorentzian_packet_count == 2_633_885
    assert census.raw_selected_zero_determinant_count == 1
    assert census.retained_selected_zero_determinant_count == 0
    assert sum(row.count for row in census.strata) == 4_768_141
    assert census.canonical_stream_sha256 == (
        "41386c1bfb64c497acc2ddb967c9a0d646f82e67dfc8907b8bc7d3b50e85bbe3"
    )
    assert census.retained_stream_sha256 == (
        "8be12b4371f49e1a0dd35f2c290c4e5735e64da56dc209655d1e0af5b49476b6"
    )
