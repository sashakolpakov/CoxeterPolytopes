from scripts.coupled_discriminant11_candidate import analyze


def test_coupled_discriminant11_k3_certificates() -> None:
    result = analyze()
    assert result["field_determinant"] == "-1 - 4*phi"
    assert result["field_signature_identity"] == {"positive": 9, "negative": 1, "zero": 0}
    assert result["field_signature_conjugate"] == {"positive": 10, "negative": 0, "zero": 0}
    assert result["trace_scalar_is_codifferent"] is True
    assert result["trace_even_integral"] is True
    assert result["trace_determinant"] == "-11"
    assert result["trace_signature"] == {"positive": 1, "negative": 19, "zero": 0}
    assert result["trace_basis_change_determinant"] in {"1", "-1"}
    assert result["trace_basis_change_verified"] is True
    assert result["discriminant_group"] == "Z/11Z"
    assert result["phi_action_minimal_polynomial_verified"] is True
    assert result["phi_action_self_adjoint_verified"] is True
    assert result["primitive_k3_embedding_verified"] is True
    assert result["primitive_extension_determinant"] in {"1", "-1"}
    assert result["k3_complement_gram"] == [[18, 5], [5, 2]]
    assert result["k3_complement_verified"] is True

