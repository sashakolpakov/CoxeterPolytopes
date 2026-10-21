from scripts.bugaenko_trace_lattice import analyze


def test_bugaenko_trace_lattice_certificate() -> None:
    result = analyze()
    assert result["normalization_status"] == "verified_against_pinned_vinberg_fixture"
    assert result["trace_scalar_is_codifferent"] is True
    assert result["rank"] == 18
    assert result["even_integral"] is True
    assert result["determinant"] == "-4"
    assert result["inertia"] == {"positive": 1, "negative": 17, "zero": 0}
    assert result["basis_change_verified"] is True
    assert result["basis_change_determinant"] in {"1", "-1"}
    assert result["primitive_k3_embedding_verified"] is True
    assert result["primitive_extension_determinant"] in {"1", "-1"}
    assert result["orthogonal_complement_verified"] is True
