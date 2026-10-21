from scripts.coupled_discriminant19_candidate import analyze


def test_coupled_discriminant19_k3_certificates() -> None:
    result = analyze()
    assert result["field_determinant"] == "-1 - 5*phi"
    assert result["field_determinant_norm"] == "-19"
    assert result["field_signature_identity"] == {"positive": 9, "negative": 1, "zero": 0}
    assert result["field_signature_conjugate"] == {"positive": 10, "negative": 0, "zero": 0}
    assert result["trace_scalar"] == {"a": "-2/5", "b": "-1/5"}
    assert result["trace_scalar_is_codifferent"] is True
    assert result["trace_even_integral"] is True
    assert result["trace_determinant"] == "-19"
    assert result["trace_signature"] == {"positive": 1, "negative": 19, "zero": 0}
    assert result["trace_basis_change_determinant"] == "1"
    assert result["trace_basis_change_verified"] is True
    assert result["discriminant_group"] == "Z/19Z"
    assert result["phi_action_minimal_polynomial_verified"] is True
    assert result["phi_action_self_adjoint_verified"] is True
    assert result["primitive_k3_embedding_verified"] is True
    assert result["primitive_extension_determinant"] in {"1", "-1"}
    assert result["k3_complement_gram"] == [[2, 1], [1, 10]]
    assert result["k3_complement_verified"] is True
    assert result["eigen_positive_control"] == {
        "field_vector": [
            {"a": "2", "b": "0"},
            {"a": "0", "b": "-1"},
            *({"a": "0", "b": "0"} for _ in range(8)),
        ],
        "field_form_square": {"a": "-2", "b": "-10"},
        "identity_eigenspace_square": {"a": "14/5", "b": "32/5"},
        "restriction_of_scalars_trace_square": "12",
    }
    assert result["fixed_i19_fibration_marking"]["status"] == "not_available"
    assert "unresolved" in result["geometric_status"]
