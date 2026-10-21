import hashlib
import json

from scripts.certify_og10_golden_projected_wall_mixed_shell3 import (
    CERTIFICATE_PATH,
    EXPECTED_ELIGIBLE_COUNT,
    EXPECTED_ELIGIBLE_HASH,
    EXPECTED_WALL_HASH,
    EXPECTED_WALL_TEST_COUNT,
    verify_certificate,
)
from scripts.independent_og10_golden_projected_wall_mixed_shell3 import (
    run_independent_replay,
)


EXPECTED_PAYLOAD_SHA256 = (
    "5c5f17fe36b811ed95d7aba09ea43430c3b37a2619f9cae386dea51f58b9510d"
)
EXPECTED_FILE_SHA256 = (
    "8a7c11074e5965f631b1f1c4e9d7c746275f5a2938f7e56646f7c0efb9478e08"
)


def load_certificate() -> dict[str, object]:
    return json.loads(CERTIFICATE_PATH.read_text())


def test_shell3_frozen_certificate_and_literal_comparisons() -> None:
    payload = load_certificate()
    assert payload["certificate_payload_sha256"] == EXPECTED_PAYLOAD_SHA256
    assert hashlib.sha256(CERTIFICATE_PATH.read_bytes()).hexdigest() == (
        EXPECTED_FILE_SHA256
    )
    assert all(verify_certificate(payload).values())

    production = payload["production_audit"]
    comparison = payload["comparison"]
    assert production["eligible_w_count"] == EXPECTED_ELIGIBLE_COUNT
    assert production["eligible_w_records_sha256"] == EXPECTED_ELIGIBLE_HASH
    assert production["exact_wall_test_count"] == EXPECTED_WALL_TEST_COUNT
    assert production["canonical_wall_records_sha256"] == EXPECTED_WALL_HASH
    assert production["cutting_covector_count"] == 0
    assert comparison["literal_eligible_record_count"] == EXPECTED_ELIGIBLE_COUNT
    assert comparison["literal_eligible_record_equality_before_hash"] is True
    assert comparison["literal_coefficient_record_count"] == 46899
    assert comparison["literal_coefficient_record_equality_before_hash"] is True
    assert comparison["literal_wall_record_count"] == EXPECTED_WALL_TEST_COUNT
    assert comparison["literal_wall_record_equality_before_hash"] is True


def test_standalone_shell3_weyl_transport_replay() -> None:
    payload = load_certificate()
    replay = run_independent_replay(30)
    assert replay == payload["independent_replay"]
    assert replay["imports_project_code"] is False
    assert replay["raw_a_b_pair_count_before_Weyl_quotient"] == 242_965_920
    assert replay["seed_pair_test_count_after_Weyl_quotient"] == 79_923
    assert replay["theta_pinned_shell_counts"] == {
        "2": 240,
        "4": 2160,
        "6": 6720,
        "8": 17520,
    }
    assert replay["shell_Weyl_orbit_counts"]["8"] == [17280, 240]
    assert replay["eligible_a_count_at_b_seed_by_b_square"] == {
        "2": 58,
        "4": 80,
        "6": 57,
    }
    assert replay["eligible_w_a_square_histogram"] == {
        "0": 2400,
        "2": 238800,
        "4": 321840,
        "6": 6720,
    }
    assert replay["eligible_w_b_square_histogram"] == {
        "2": 13920,
        "4": 172800,
        "6": 383040,
    }
    assert replay["separation_sign_histogram"] == {
        "-1,-1": 736320,
        "1,1": 736320,
    }
    assert replay["cutting_records"] == []


def test_shell3_scope_is_explicitly_finite_box() -> None:
    payload = load_certificate()
    assert payload["box"]["A_and_C_integral_coefficient_interval"] == [-30, 30]
    assert payload["box"]["b_square_values"] == [2, 4, 6]
    assert payload["box"]["b_shell_bound_cli_value"] == 3
    assert "Finite-box theorem only" in payload["scope"]
    assert "neither a global OG10 obstruction nor a new hyperbolic Coxeter example" in (
        payload["scope"]
    )
