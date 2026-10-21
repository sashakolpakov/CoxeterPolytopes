from scripts.audit_og10_golden_rejected_separator import audit


def test_malformed_separator_is_rejected_and_correction_does_not_cut() -> None:
    result = audit()
    assert result["status"] == "rejected_transcript_not_a_candidate"
    assert all(result["checks"].values())

    literal = result["literal"]
    assert literal["q_equals_2_p_square"] == [782, 1259]
    assert literal["evaluation_signs"] == [-1, 1]
    assert literal["identity_signature_positive_negative"] == [9, 2]
    assert literal["passes_required_signatures"] is False

    corrected = result["corrected"]
    assert corrected["q_equals_2_p_square"] == [2, -1]
    assert corrected["two_p_is_integral"] is True
    assert corrected["J_squared_is_identity"] is True
    assert corrected["J_is_H_isometry"] is True
    assert corrected["identity_signature_positive_negative"] == [10, 1]
    assert corrected["conjugate_signature_positive_negative"] == [11, 0]
    assert corrected["ell_on_400_times_interior"] == [-10659, -17409]
    assert corrected["ell_on_13_times_v61"] == [-259, -415]
    assert corrected["evaluation_signs"] == [-1, -1]
    assert corrected["cuts_v61"] is False
