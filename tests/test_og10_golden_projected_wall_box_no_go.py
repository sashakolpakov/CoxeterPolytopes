import json
from pathlib import Path

from scripts.search_og10_golden_projected_wall import (
    CERTIFICATE_PATH,
    INV_PHI_PAIR,
    _add,
    _conjugate,
    _mul,
    _neg,
    _sign,
    audit_mixed_e8_family,
    audit_scalar_unit_family,
    audit_zero_e8,
    verify_certificate,
)


ROOT = Path(__file__).resolve().parents[1]


def test_full_frozen_boxes_recompute_and_match_digests() -> None:
    """Replay every declared wall test; do not trust frozen zero counts."""

    payload = json.loads(CERTIFICATE_PATH.read_text())
    assert payload["certificate_payload_sha256"] == (
        "48a4549d543b7f89cb69543d4656a739d5c3ca08abe8b077ffbd93af77fbfa97"
    )
    assert payload["audit"]["combined_audit_sha256"] == (
        "e13ef11274acfe60eada6beead349fe46dbc7397de27329fd2dfa0ee4cc8dcd9"
    )
    audit = payload["audit"]["families"]

    zero = audit_zero_e8(30)
    assert zero == audit["zero_E8"]
    assert zero["signature_admissible_parameter_pair_count"] == 66
    assert zero["cutting_covector_count"] == 0

    scalar = audit_scalar_unit_family(30, 8, 4)
    assert scalar == audit["scalar_unit_E8"]
    assert scalar["e8_vector_count"] == 26640
    assert scalar["e8_shell_counts"] == {
        "2": 240,
        "4": 2160,
        "6": 6720,
        "8": 17520,
    }
    assert scalar["exact_wall_test_count"] == 2706240
    assert scalar["cutting_covector_count"] == 0

    mixed = audit_mixed_e8_family(30, 1)
    assert mixed == audit["mixed_E8_first_b_shell"]
    assert mixed["b_vector_count"] == 240
    assert mixed["eligible_w_count"] == 13920
    assert mixed["eligible_w_square_class_count"] == 3
    assert mixed["exact_wall_test_count"] == 15360
    assert mixed["cutting_covector_count"] == 0

    checks = verify_certificate(payload)
    assert all(checks.values())


def test_certificate_scope_and_uncertified_larger_shell_are_explicit() -> None:
    payload = json.loads(CERTIFICATE_PATH.read_text())
    assert payload["certificate_payload_sha256"] == (
        "48a4549d543b7f89cb69543d4656a739d5c3ca08abe8b077ffbd93af77fbfa97"
    )
    assert payload["audit"]["combined_audit_sha256"] == (
        "e13ef11274acfe60eada6beead349fe46dbc7397de27329fd2dfa0ee4cc8dcd9"
    )
    assert "not a global non-existence theorem" in payload["scope"]
    exploratory = payload["exploratory_not_certified"]["mixed_b_shell_bound_2"]
    assert "excluded from theorem" in exploratory["status"]
    assert exploratory["observed_cutting_covector_count"] == 0


def test_twice_projected_square_formula_recovers_a2_case() -> None:
    # A=C=1 and w=0 gives q=-1/phi+1=2-phi, hence p^2=(2-phi)/2.
    a_square_term = _neg(_mul(_mul((1, 0), (1, 0)), INV_PHI_PAIR))
    q = _add(a_square_term, (1, 0))
    assert q == (2, -1)
    assert _sign(q) > 0
    assert _sign(_add((4, 0), _neg(q))) > 0
    assert _sign(_conjugate(_add((4, 0), _neg(q)))) > 0
