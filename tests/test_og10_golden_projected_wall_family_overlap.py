from scripts.audit_og10_golden_projected_wall_family_overlap import audit


def test_scalar_mixed_overlap_and_union_count_are_exact() -> None:
    result = audit()
    assert result["family_indexed_wall_evaluation_count"] == 4_178_946
    assert result["scalar_mixed_overlap_eligible_direction_count"] == 9_360
    assert result["scalar_mixed_overlap_wall_count"] == 209_280
    assert result["distinct_geometric_wall_evaluation_count"] == 3_969_666
    assert result["scalar_mixed_overlap_rows"] == [
        {
            "unit_power": 1,
            "E8_vector_square": 2,
            "mixed_b_square": 2,
            "E8_vector_count": 240,
            "admissible_coefficient_count": 40,
            "overlap_wall_count": 9_600,
        },
        {
            "unit_power": 2,
            "E8_vector_square": 2,
            "mixed_b_square": 2,
            "E8_vector_count": 240,
            "admissible_coefficient_count": 24,
            "overlap_wall_count": 5_760,
        },
        {
            "unit_power": 2,
            "E8_vector_square": 4,
            "mixed_b_square": 4,
            "E8_vector_count": 2_160,
            "admissible_coefficient_count": 40,
            "overlap_wall_count": 86_400,
        },
        {
            "unit_power": 2,
            "E8_vector_square": 6,
            "mixed_b_square": 6,
            "E8_vector_count": 6_720,
            "admissible_coefficient_count": 16,
            "overlap_wall_count": 107_520,
        },
    ]
