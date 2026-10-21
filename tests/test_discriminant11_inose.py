from scripts.discriminant11_inose_model import analyze


def test_discriminant11_inose_model() -> None:
    result = analyze()
    assert result["cm_discriminant"] == -11
    assert result["cm_j_invariant"] == "-32768"
    assert result["transcendental_lattice"] == [[2, 1], [1, 6]]
    assert result["inose_parameter_identities_verified"] is True
    assert result["c4"] == "1024*t**4"
    assert result["discriminant_factorization_verified"] is True
    assert result["fiber_euler_sum"] == 24
    assert result["shioda_tate_mordell_weil_rank"] == 1

