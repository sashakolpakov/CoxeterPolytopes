from scripts.certify_reflection_maximal_weak_log_exact_nef_return import payload


def test_exact_weak_log_f4_calibration() -> None:
    data = payload()
    assert data["canonical_square"] == [7, 1]
    assert data["curve_squares"] == {
        "B": [-4, 1],
        "E": [-1, 1],
        "G": [-1, 1],
    }
    assert data["weak_log_square"] == [8, 1]
    assert data["weak_log_pairings_on_complete_toric_mori_generators"] == {
        "B": [0, 1],
        "E": [1, 1],
        "G": [1, 2],
    }
    assert data["facet_root_gram"] == [
        [[-2, 1], [0, 1], [1, 1]],
        [[0, 1], [-2, 1], [2, 1]],
        [[1, 1], [2, 1], [-2, 1]],
    ]
    assert data["coxeter_orders"] == {
        "B_half_E": 2,
        "B_half_G": 3,
        "E_G": "infinity",
    }
    assert data["area_over_pi"] == [1, 6]
    assert data["proper_overgroup_area_upper_over_pi"] == [1, 12]
    assert data["minimum_noncompact_coxeter_area_over_pi"] == [1, 6]
