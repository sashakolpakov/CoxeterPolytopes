import pytest

from scripts.audit_golden_a1_uniform_ray_cutter import build_audit


@pytest.fixture(scope="module")
def audit():
    return build_audit()


def test_complete_prefix_census_before_uniform_ray_cutter(audit):
    assert [control["display"] for control in audit["control_bound_proof"]["all_controls"]] == [
        "1",
        "1*phi",
        "1 + 1*phi",
        "2*phi",
        "1 + 2*phi",
    ]
    assert {
        name: (
            record["raw_pattern_count"],
            record["golay_feasible_pattern_count"],
            record["uniform_violating_pattern_count"],
        )
        for name, record in audit["pattern_census"].items()
    } == {
        "phi_squared": (43, 3, 0),
        "two_phi": (3, 1, 0),
        "cutter": (383, 33, 16),
    }
    assert audit["phi_squared_packet"]["accepted_assignment_count"] == 18216
    assert audit["phi_squared_packet"]["rejected_assignment_count"] == 6072
    assert audit["two_phi_packet"]["accepted_assignment_count"] == 0


def test_phi_cubed_root_really_cuts_the_uniform_ray(audit):
    cutter = audit["cutter"]
    assert cutter["construction_a_membership"]
    assert cutter["norm"]["display"] == "2"
    assert cutter["survives_every_earlier_wall"]
    assert {
        name: record["display"]
        for name, record in cutter["maximum_pairings_with_prior_packets"].items()
    } == {
        "height_zero": "0",
        "x_equals_1": "-2 - 4*phi",
        "x_equals_phi": "-2 - 2*phi",
        "x_equals_phi_squared_integral_orbit": "-6 - 8*phi",
        "x_equals_phi_squared_octad_orbit": "-4 - 4*phi",
    }
    assert cutter["uniform_excess_sum_d_minus_4x"]["display"] == "4 + 1*phi"
    assert cutter["uniform_ray_pairing"]["display"] == "1/2 + 5/2*phi"
    assert cutter["cuts_uniform_ray"]
    assert audit["conclusion"] == {
        "uniform_ray_persists_through_phi_squared": True,
        "uniform_ray_cut_at_phi_cubed": True,
        "cutter_is_norm_two_lattice_root": True,
        "cutter_survives_all_earlier_walls": True,
        "violating_golay_pattern_count_at_phi_cubed": 16,
    }

