from scripts.certify_toric_h6_congruence_smoothing import build_certificate


def test_toric_h6_congruence_and_smoothing_certificate() -> None:
    certificate = build_certificate()
    finite = certificate["finite_parabolic"]
    assert finite["type"] == "A2^3"
    assert finite["order"] == finite["mod_2_image_order"] == 216
    assert certificate["exceptional_curve_orbits"] == [9, 9, 9]
    assert certificate["orbit_union_is_exactly_all_27_cubic_line_classes"] is True
    assert certificate["all_27_are_zero_one_boundary_sums"] is True

    congruence = certificate["congruence_argument"]
    assert congruence["index_over_Gamma_2_6"] == 216
    assert congruence["P6_endpoint_chamber_count"] == 216
    assert congruence["endpoint_is_full_mod_2_preimage"] is True

    assert certificate["solomon"] == {
        "endpoint_spherical_subset_count": 322,
        "endpoint_characteristic": "-1/1728",
        "P6_right_angled_clique_vector": [1, 27, 216, 720, 1080, 648, 72],
        "P6_characteristic": "-1/8",
    }
    assert certificate["volume"]["endpoint_volume_over_pi_cubed"] == "1/3240"
    assert certificate["volume"]["P6_volume_over_pi_cubed"] == "1/15"

    smoothing = certificate["three_A2_smoothing"]
    assert smoothing["central_cubic"] == "XYZ=W^3"
    assert smoothing["marked_line_class_specialization"] == (
        "three W(A2)^3-orbits of size 9"
    )
