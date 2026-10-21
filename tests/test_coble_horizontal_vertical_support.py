from scripts.certify_coble_horizontal_vertical_support import build_certificate


def test_nine_rows_and_dimension_formula() -> None:
    certificate = build_certificate()
    rows = certificate["classified_rows"]
    assert len(rows) == 9
    for row in rows:
        assert row["picard_rank"] == 10 + row["boundary_count"]
        assert row["hyperbolic_dimension"] == 9 + row["boundary_count"]


def test_horizontal_ceiling_is_strictly_lower() -> None:
    derived = build_certificate()["derived"]
    assert derived["all_hyperbolic_dimensions"] == [10, 11, 14, 17]
    assert derived["horizontal_hyperbolic_dimensions"] == [10, 11]
    assert derived["horizontal_dimension_ceiling"] == 11


def test_high_boundary_rows_are_vertical_at_three() -> None:
    certificate = build_certificate()
    assert certificate["support_by_boundary_count"] == {
        "1": "all_odd",
        "2": "all_odd",
        "5": ["3"],
        "8": ["3"],
    }
    assert certificate["derived"]["vertical_only_high_boundary_support"] == {
        "5": [3],
        "8": [3],
    }


def test_limitations_do_not_overclaim_flagship() -> None:
    limitations = build_certificate()["limitations"]
    assert limitations["does_not_prove_carrier_nonemptiness_for_raw_D_CobVin"]
    assert limitations["strict_flagship_status"].endswith("not L4")
