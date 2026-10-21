import json

from scripts.certify_og10_golden_rank_two_reduction import (
    CERTIFICATE_PATH,
    _norm_six_signature_admissible,
    build_certificate,
    structural_audit,
    verify_certificate,
)


EXPECTED_CERTIFICATE_SHA256 = (
    "26293b9d366e4250ce0b940298f3ff338937b022fcd0ebe262adf4805183eebb"
)
EXPECTED_BOX_AUDIT_SHA256 = (
    "bf49a9ed42bf1760e3912fa126659546adedd9ec0395aa5b93dbf57428983d70"
)


def test_rank_two_obstruction_and_mod_three_folding_gate() -> None:
    payload = json.loads(CERTIFICATE_PATH.read_text())
    assert payload["certificate_payload_sha256"] == EXPECTED_CERTIFICATE_SHA256

    structural = structural_audit()
    assert all(structural["checks"].values())
    assert structural["base"]["trace_inertia"] == [1, 19, 0]
    assert structural["rank_two_model"]["trace_rank"] == 24
    assert structural["rank_two_model"]["trace_inertia"] == [1, 23, 0]
    assert structural["og10_beauville_lattice"]["inertia"] == [3, 21, 0]
    assert structural["norm_six_test_model"]["trace_mod_3_radical_dimension"] == 2
    assert structural["mod_3_ambient_gate"]["og10_mod_3_radical_dimension"] == 1
    assert structural["mod_3_ambient_gate"]["conclusion"] == (
        "at most one of e and (phi-1)e can have ambient divisibility 3"
    )


def test_norm_six_frozen_boxes_recompute_exactly() -> None:
    archived = json.loads(CERTIFICATE_PATH.read_text())
    fresh = build_certificate()
    assert fresh == archived
    assert archived["certificate_payload_sha256"] == EXPECTED_CERTIFICATE_SHA256

    audit = archived["norm_six_box_audit"]
    assert audit["combined_audit_sha256"] == EXPECTED_BOX_AUDIT_SHA256
    zero = audit["families"]["zero_E8"]
    scalar = audit["families"]["scalar_unit_E8"]
    mixed = audit["families"]["mixed_E8_first_b_shell"]
    assert (zero["exact_wall_test_count"], zero["cutting_covector_count"]) == (18, 0)
    assert (
        scalar[
            "signature_admissible_parameter_pair_count_summed_over_power_and_shell"
        ],
        scalar["exact_wall_test_count"],
        scalar["cutting_covector_count"],
    ) == (136, 696000, 0)
    assert (
        mixed[
            "signature_admissible_parameter_pair_count_summed_over_w_square_classes"
        ],
        mixed["exact_wall_test_count"],
        mixed["cutting_covector_count"],
    ) == (24, 5760, 0)
    assert mixed["eligible_w_count"] == 13920
    assert all(verify_certificate(archived).values())


def test_norm_six_gate_is_exact_strict_four_thirds_gate() -> None:
    assert _norm_six_signature_admissible((1, 0))
    assert not _norm_six_signature_admissible((0, 0))
    assert not _norm_six_signature_admissible((2, 0))
    # q=2-phi is small at the identity place but too large at its conjugate.
    assert not _norm_six_signature_admissible((2, -1))
