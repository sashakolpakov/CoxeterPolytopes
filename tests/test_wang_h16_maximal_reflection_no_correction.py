import json

from scripts.certify_wang_h16_maximal_reflection_no_correction import (
    ARCHIVE_PATH,
    build_certificate,
)


def _archive() -> dict:
    return json.loads(ARCHIVE_PATH.read_text())


def test_wang_h16_exact_certificate_replays() -> None:
    archived = _archive()
    assert build_certificate() == archived
    assert archived["all_checks_pass"] is True


def test_tumarkin_identification_and_f_vector_are_frozen() -> None:
    polytope = _archive()["polytope"]
    assert (
        polytope["dimension"],
        polytope["facets"],
        polytope["finite_vertices"],
        polytope["ideal_vertices"],
        polytope["elliptic_edges"],
    ) == (16, 19, 125, 3, 1087)
    assert polytope["published_f_vector"] == [
        128,
        1087,
        4768,
        14000,
        30352,
        50960,
        67960,
        72908,
        63204,
        44200,
        24752,
        10948,
        3740,
        952,
        170,
        19,
        1,
    ]
    assert polytope["euler_characteristic"] == (
        "2499347/2360171042879569920000"
    )


def test_three_maximality_notions_and_nonreflection_overgroup_are_distinguished() -> None:
    archived = _archive()
    lattice = archived["lattice"]
    symmetry = archived["diagram_symmetry"]
    maximality = archived["maximality"]
    assert lattice["absolute_discriminant"] == 2
    assert lattice["integral_reflection_norms"] == [-2]
    assert lattice["conclusion"] == "W_all(S) = W_2(S)"
    assert symmetry["weighted_automorphism_group"] == "C2"
    assert symmetry["lorentz_quotient_eigenspace_dimensions_plus_minus"] == [9, 8]
    assert symmetry["hyperplane_reflection_required_eigenspaces_plus_minus"] == [16, 1]
    assert "maximal among discrete/arithmetic reflection" in maximality["conclusion"]
    assert maximality["nonreflection_overgroup"].startswith(
        "Theta16 = Gamma16 semidirect C2"
    )
    assert maximality["not_claimed"] == (
        "Gamma16 is not a maximal arithmetic lattice/group"
    )


def test_maximal_reflection_does_not_imply_borcherds_correction() -> None:
    archived = _archive()
    assert archived["borcherds"] == {
        "Wang_theorem_7_2_automorphic_correction": False,
        "disproved_implication": (
            "maximal arithmetic reflection group => automorphic correction"
        ),
    }
    assert archived["scope"]["new_polytope"] is False
    assert archived["scope"]["no_AG_finiteness_without_volume_input"] is True
