from scripts.l9_seed_k3_trace import analyze


def test_l9_seed_has_k3_discriminant_length_obstruction() -> None:
    result = analyze()
    assert result["trace_scalar_is_codifferent"] is True
    assert result["rank"] == 20
    assert result["even_integral"] is True
    assert result["determinant"] == "-16"
    assert result["inertia"] == {"positive": 1, "negative": 19, "zero": 0}
    assert result["basis_change_verified"] is True
    assert result["basis_change_determinant"] in {"1", "-1"}
    assert result["discriminant_group"] == "(Z/2Z)^4"
    assert result["required_k3_complement_rank"] == 2
    assert result["primitive_k3_embedding"] is False
    assert result["obstruction_verified"] is True

