from scripts.explore_h9_14_record0_a2_source29_triple123_residual_supports import (
    EXPECTED_MATHEMATICAL_DIGEST,
    quotient,
)


def test_exact_residual_support_quotient() -> None:
    result = quotient()
    assert result["mathematical_digest"] == EXPECTED_MATHEMATICAL_DIGEST
    assert result["row_count"] == 61
    assert result["configuration_branch_instance_count"] == 144
    assert result["support_orbit_count"] == 9004
    assert result["support_branch_instance_count"] == 25407
    assert result["minimal_killer_necessary_gate_pass_count"] == 25407
    assert result["minimal_killer_necessary_gate_reject_count"] == 0
    assert result["source29_status"] == "OPEN"
    assert result["summary"] == [
        {"family": "active", "fixed_230": False, "degree": 0, "support_orbit_count": 6867},
        {"family": "active", "fixed_230": False, "degree": 1, "support_orbit_count": 729},
        {"family": "active", "fixed_230": False, "degree": 2, "support_orbit_count": 215},
        {"family": "active", "fixed_230": True, "degree": 0, "support_orbit_count": 588},
        {"family": "active", "fixed_230": True, "degree": 1, "support_orbit_count": 139},
        {"family": "active", "fixed_230": True, "degree": 2, "support_orbit_count": 54},
        {"family": "zero", "fixed_230": False, "degree": 0, "support_orbit_count": 335},
        {"family": "zero", "fixed_230": False, "degree": 1, "support_orbit_count": 55},
        {"family": "zero", "fixed_230": False, "degree": 2, "support_orbit_count": 22},
    ]
