from fractions import Fraction

from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_j2dim14 import (
    UNIT_CONTINUATION_HEIGHT,
    UNIT_CUTTER_DIFFERENCE,
    UNIT_CUTTER_U,
    UNIT_CUTTER_V,
    coordinate_slice_vectors,
    j2_pair,
    unit_cutter_root_audit,
    unit_h137_continuation_audit,
)


def test_exact_short_coordinate_slices():
    norm_twelve = coordinate_slice_vectors(12, 4)
    norm_fourteen = coordinate_slice_vectors(14, 4)
    assert norm_twelve.node_count == 331
    assert len(norm_twelve.vectors) == 2
    assert UNIT_CUTTER_DIFFERENCE in norm_twelve.vectors
    assert norm_fourteen.node_count == 9_647
    assert norm_fourteen.vectors == ()


def test_explicit_h137_root():
    assert unit_cutter_root_audit()
    assert j2_pair(UNIT_CUTTER_U, UNIT_CUTTER_U) == 136
    assert j2_pair(UNIT_CUTTER_V, UNIT_CUTTER_V) == 70
    assert j2_pair(UNIT_CUTTER_U, UNIT_CUTTER_V) == 97


def test_complete_h137_chronology_and_replacement():
    audit = unit_h137_continuation_audit()
    assert audit.control_count == 40
    assert audit.distinct_height_count == 27
    assert audit.nonempty_control_count == 4
    assert audit.gated_type_count == 9
    assert audit.first_old_witness_threat_height == UNIT_CONTINUATION_HEIGHT
    assert audit.cutter_first_packet_maximum == 27 + 20 * SQRT2
    assert audit.cutter_first_packet_acute_margin == 9 + 6 * SQRT2
    assert audit.cutter_old_witness_pairing == 3 + 2 * SQRT2
    assert audit.replacement_scale == Fraction(5, 4)
    assert audit.replacement_lorentz_norm == Fraction(3, 2) - SQRT2
    assert audit.replacement_first_packet_margin == Fraction(-1, 2) + 2 * SQRT2
    assert min(row.replacement_margin for row in audit.rows) == (
        Fraction(1, 2) + Fraction(1, 2) * SQRT2
    )
    assert audit.replacement_all_margins_strict
