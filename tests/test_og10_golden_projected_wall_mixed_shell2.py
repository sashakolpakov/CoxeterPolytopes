import hashlib
import json

from scripts.certify_og10_golden_projected_wall_mixed_shell2 import (
    CERTIFICATE_PATH,
    EXPECTED_ELIGIBLE_COUNT,
    EXPECTED_ELIGIBLE_HASH,
    EXPECTED_WALL_HASH,
    EXPECTED_WALL_TEST_COUNT,
    verify_certificate,
)
from scripts.independent_og10_golden_projected_wall_mixed_shell2 import (
    run_independent_replay,
)


EXPECTED_PAYLOAD_SHA256 = (
    "9c696c3d540fe5570b2f1becd05fb24b5bd334d614e527d5c682440b5987196f"
)
EXPECTED_FILE_SHA256 = (
    "b426e419da6b9e79a41955aba15b398e9d7063abc77b4b69f4b1f774fec7c669"
)


def load_certificate() -> dict[str, object]:
    return json.loads(CERTIFICATE_PATH.read_text())


def test_shell2_frozen_certificate_and_provenance() -> None:
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


def test_standalone_weyl_shell_replay_matches_every_frozen_summary() -> None:
    """Recompute the independent 21,890,400-pair exhaustion."""

    payload = load_certificate()
    replay = run_independent_replay(30)
    assert replay == payload["independent_replay"]
    assert replay["imports_project_code"] is False
    assert replay["raw_a_b_pair_count"] == 21_890_400
    assert replay["theta_pinned_shell_counts"] == {
        "2": 240,
        "4": 2160,
        "6": 6720,
    }
    assert replay["eligible_w_a_square_histogram"] == {
        "0": 2400,
        "2": 43920,
        "4": 140400,
    }
    assert replay["eligible_w_b_square_histogram"] == {
        "2": 13920,
        "4": 172800,
    }
    assert replay["separation_sign_histogram"] == {
        "-1,-1": 413760,
        "1,1": 413760,
    }
    assert replay["cutting_records"] == []


def test_scope_is_finite_box_and_not_a_hyperbolic_example() -> None:
    payload = load_certificate()
    assert payload["box"]["A_and_C_integral_coefficient_interval"] == [-30, 30]
    assert payload["box"]["b_square_values"] == [2, 4]
    assert payload["box"]["b_shell_bound_cli_value"] == 2
    assert "Finite-box theorem only" in payload["scope"]
    assert "neither a global OG10 obstruction nor a new hyperbolic Coxeter example" in (
        payload["scope"]
    )
