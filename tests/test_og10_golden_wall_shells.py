import json
from pathlib import Path

from scripts.certify_og10_golden_wall_shells import certify


def test_og10_golden_complete_wall_shell_certificate():
    result = certify()
    archived = json.loads(
        Path("data/candidates/og10_golden_complete_wall_shells.json").read_text()
    )

    assert result == archived
    assert result["verified"]
    assert result["inherited_carrier_candidate_id"] == "og10_golden_h10_carrier"
    assert result["checks"]["inherited_og10_carrier_certificate_verified"]

    mod_three = result["ambient_mod_3"]
    assert mod_three["ambient_determinant"] == -3
    assert mod_three["radical_dimension"] == 1
    assert mod_three["radical_intersection_n_dimension"] == 0
    assert mod_three["pairing_map_from_n_rank_mod_3"] == 22
    assert mod_three["golden_polynomial_values_mod_3"] == [2, 2, 1]
    assert mod_three["self_adjoint_extension_obstruction"].startswith(
        "no_integral_self_adjoint_golden_endomorphism"
    )

    wall_types = result["og10_wall_types_on_this_n"]
    assert wall_types["surviving"] == [
        {"square": -2, "ambient_divisibility": 1},
        {"square": -4, "ambient_divisibility": 1},
    ]
    assert wall_types["excluded"] == [
        {"square": -6, "ambient_divisibility": 3},
        {"square": -24, "ambient_divisibility": 3},
    ]

    assert result["complete_restricted_search_norms_text"] == [
        "2",
        "4",
        "6-2*phi",
    ]
    assert [
        orbit["field_norm"]
        for orbit in result["restricted_mirror_unit_square_orbits"]
    ] == [4, 16, 20]
    assert result["checks"]["three_representative_families_are_nonempty"]
    assert result["checks"]["three_mirror_families_are_pairwise_distinct"]
    assert result["checks"][
        "extra_kahler_mirrors_are_not_integral_field_reflections"
    ]

    factorization = result["divisorial_factorization_of_norm_2_reflection"]
    assert factorization["pair_gram"] == [[-2, 0], [0, -2]]
    assert factorization["product_on_O_root_line"] == [[-1, 0], [0, -1]]
    assert factorization["group_membership"] == (
        "restriction_to_N_of_W_Exc_intersection_Cent_O(N)(tau)"
    )
