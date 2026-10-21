from scripts.certify_compact_k3_root_return import certify


def test_compact_k3_root_return_certificate() -> None:
    result = certify()
    assert result["all_checks_pass"] is True
    assert result["case_count"] == 8
    assert result["dimension_histogram"] == {"2": 6, "3": 2}
    assert result["pairwise_distinct_unlabelled_wall_grams"] is True
    assert result["weyl_vector_count"] == 5
    assert result["new_beyond_weyl_vector_keys"] == ["S4", "S6", "L27"]


def test_every_certified_chamber_is_compact_and_irredundant() -> None:
    result = certify()
    for case in result["cases"]:
        assert case["all_extreme_rays_finite"] is True
        assert case["facets"] == case["f_vector"][-1]
        assert case["walls_span_ns_integrally"] is True
        assert abs(case["integral_ns_change_of_basis_determinant"]) == 1


def test_exact_weyl_vector_strata_and_positive_squares() -> None:
    result = certify()
    squares = {
        case["key"]: case["rational_weyl_vector_square"]
        for case in result["cases"]
    }
    assert squares == {
        "S1": "1/2",
        "S2": "1/4",
        "S3": "1",
        "S4": None,
        "S5": "2",
        "S6": None,
        "L24": "2",
        "L27": None,
    }
