from scripts.certify_rational_pfirst_effectivity_realizability_no_go import payload


def test_exact_rank_three_canonical_escape_counterexample() -> None:
    data = payload()
    assert data["canonical_square"] == 7
    assert data["primitive_basis_determinant"] == 1
    assert data["canonical_pairings"] == [0, 2, 2]
    assert data["area_over_pi"] == [1, 4]
    assert data["proper_overgroup_area_upper_over_pi"] == [1, 8]
    assert data["minimum_noncompact_coxeter_area_over_pi"] == [1, 6]
    assert data["riemann_roch_chi_by_negative_square"] == {"1": 1, "2": 0, "4": -2}

