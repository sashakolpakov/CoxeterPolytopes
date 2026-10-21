from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_leech import (
    CANONICAL_MIRROR_NORMS,
    LEECH_MINIMAL_VECTOR_COUNT,
    REDUCED_TRACE_GRAM,
    first_shell_audit,
    golay_weight_distribution,
    leech_control_levels,
    minimal_leech_shell_audit,
    reduced_trace_gram_from_change_of_basis,
    trace_gram_from_norm_decomposition,
)


def test_trace_basis_change_and_positive_norm_decomposition_are_exact():
    assert reduced_trace_gram_from_change_of_basis() == REDUCED_TRACE_GRAM
    assert trace_gram_from_norm_decomposition() == REDUCED_TRACE_GRAM


def test_extended_golay_and_standard_minimal_leech_shell_are_complete():
    assert golay_weight_distribution() == {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}
    shell = minimal_leech_shell_audit()
    assert shell["type_counts"] == {
        "2^8": 97_152,
        "3,1^23": 98_304,
        "4^2": 1_104,
    }
    assert shell["total_count"] == LEECH_MINIMAL_VECTOR_COUNT
    assert shell["co0_orbit_count"] == 1
    assert shell["point_stabilizer"] == "Co2"


def test_leech_host_has_only_norm_two_and_no_first_packet_roots():
    assert CANONICAL_MIRROR_NORMS == (Cubic7(2),)
    levels = leech_control_levels(10)
    assert len(levels) == 1
    assert levels[0].control == Cubic7(-1, 1, 1)
    assert levels[0].height == Cubic7(0, 2, 1)
    assert levels[0].target == Cubic7(2, 4, 2)

    audit = first_shell_audit()
    assert audit["height_zero"]["root_count"] == 0
    assert audit["first_positive_control_window"]["candidate_root_count"] == 0
    assert audit["first_positive_control_window"]["trace_minimal_spatial_vectors_examined"] == 196_560
