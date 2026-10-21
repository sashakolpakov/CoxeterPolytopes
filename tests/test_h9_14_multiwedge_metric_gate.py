"""Regression certificate for the fixed-realization H9/14 metric gate.

These tests freeze exact calculations for the three displayed Delzant/
projective realizations.  They do not test, and must not be read as testing,
all realizations of the underlying combinatorial types.
"""

from sympy import zeros

from scripts.explore_h9_14_multiwedge_metric import audit, facet_covectors


EXPECTED_KERNELS = {
    43: [
        [1, 1, 0, 0],
        [-1, 0, 1, 2],
        [-1, 0, 0, 1],
        [1, -9, -7, -16],
        [1, -9, -7, -16],
        [1, -9, -7, -16],
        [1, -9, -7, -16],
        [1, -9, -7, -16],
        [-1, 0, 0, 0],
        [-1, 0, 0, 0],
        [-1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
    ],
    44: [
        [1, 1, 1, 0],
        [-1, -1, 0, 1],
        [1, -1, -2, -2],
        [1, -1, -2, -2],
        [1, -1, -2, -2],
        [1, -1, -2, -2],
        [1, -1, -2, -2],
        [-1, 0, 0, 0],
        [-1, 0, 0, 0],
        [-1, 0, 0, 0],
        [-1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
    ],
    47: [
        [1, 1, 0, 0],
        [-1, 0, 1, 0],
        [-1, 0, 0, 1],
        [1, -9, -9, -9],
        [1, -9, -9, -9],
        [1, -9, -9, -9],
        [1, -9, -9, -9],
        [1, -9, -9, -9],
        [-1, 0, 0, 0],
        [-1, 0, 0, 0],
        [-1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
    ],
}


EXPECTED_SPARSE = {
    43: (57, 40, 15, [0, 1, 2, 11, 12, 13]),
    44: (57, 39, 16, [0, 1, 11, 12, 13]),
    47: (60, 42, 13, [0, 1, 2, 11, 12, 13]),
}


def test_fixed_realizations_have_frozen_rank_kernel_and_vertex_incidence() -> None:
    for record in (43, 44, 47):
        result = audit(record)
        matrix = facet_covectors(record)
        kernel = result["kernel_basis"]

        assert matrix.shape == (14, 10)
        assert result["rank"] == 10
        assert result["kernel_dimension"] == 4
        assert kernel == EXPECTED_KERNELS[record]
        assert matrix.T * type(matrix)(kernel) == zeros(10, 4)
        assert result["vertex_count"] == 90
        assert result["vertex_complements_match_archive"] is True


def test_fixed_realization_right_angle_obstructions_are_frozen() -> None:
    for record in (43, 44, 47):
        result = audit(record)
        sparse_pair_count, sparse_rank, sparse_nullity, isotropic = EXPECTED_SPARSE[
            record
        ]

        assert result["meeting_pair_count"] == 85
        assert result["all_meeting_right_angle_rank"] == 55
        assert result["all_meeting_right_angle_solution_dimension"] == 0
        assert result["block_supported_right_pair_count"] == sparse_pair_count
        assert result["block_supported_right_angle_rank"] == sparse_rank
        assert (
            result["block_supported_right_angle_solution_dimension"]
            == sparse_nullity
        )
        assert result["block_supported_forced_isotropic_facets"] == isotropic

        assert result["cross_factor_pair_count"] == 45
        assert result["cross_factor_right_angle_rank"] == 30
        assert result["cross_factor_right_angle_solution_dimension"] == 25
        assert result["simplex_covector_span_dimension"] == 5
        assert result["residual_covector_span_dimension"] == 6
        assert result["span_intersection_dimension"] == 1
        assert result["cross_factor_orthogonality_forces_degenerate_form"] is True
