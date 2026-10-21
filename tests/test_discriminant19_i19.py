from scripts.discriminant19_i19_model import analyze


def test_discriminant19_i19_model() -> None:
    result = analyze()
    assert result["published_transcendental_lattice"] == [[2, 1], [1, 10]]
    assert result["published_neron_severi_lattice"] == "U + A18(-1)"
    assert result["davenport_stothers_identity_verified"] is True
    assert result["discriminant_identity_verified"] is True
    assert result["finite_discriminant_quintic_squarefree"] is True
    assert result["finite_c4_coprime_to_discriminant"] is True
    assert result["infinity_c4_order"] == 0
    assert result["infinity_discriminant_order"] == 19
    assert result["infinity_fiber"] == "I19"
    assert result["fiber_euler_sum"] == 24
    assert result["shioda_tate_mordell_weil_rank"] == 0
