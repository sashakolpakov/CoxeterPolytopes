from scripts.h11_h19_symmetry_seeds import analyze


def test_h11_h19_symmetry_seeds_are_exact_loxodromic_isometries() -> None:
    result = analyze()
    assert result["status"].startswith("search seeds only")

    h11, h19 = result["seeds"]
    assert h11["lattice"] == "H11"
    assert h11["preserves_binary_gram"] is True
    assert h11["determinant"] == "1"
    assert h11["trace"] == "1 + 1*phi"
    assert h11["trace_minus_two_identity_sign"] == 1
    assert h11["pell_equation_value"] == "4"
    assert h11["infinite_projective_order"] is True
    assert h11["same_chamber_certified"] is False

    assert h19["lattice"] == "H19"
    assert h19["preserves_binary_gram"] is True
    assert h19["determinant"] == "1"
    assert h19["trace"] == "2 + 1*phi"
    assert h19["trace_minus_two_identity_sign"] == 1
    assert h19["pell_equation_value"] == "4"
    assert h19["infinite_projective_order"] is True
    assert h19["same_chamber_certified"] is False
