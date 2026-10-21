import hashlib
import json

from scripts.independent_certify_k3_rho16_h14_fold import (
    CERTIFICATE_PATH,
    build_certificate,
    verify_certificate,
)


EXPECTED_CERTIFICATE_FILE_SHA256 = (
    "9c5ddc288b39016b0ee09cd97252b2fee7a805bdb8bc79ab461423ec5ae71039"
)
EXPECTED_CERTIFICATE_PAYLOAD_SHA256 = (
    "9957719eab68a090b4f47fde2a19fa6330e1e1364a8a29d77101274c7f0fbd48"
)


def _archived() -> dict:
    return json.loads(CERTIFICATE_PATH.read_text())


def test_independent_h14_audit_recomputes_and_hashes_exactly() -> None:
    archived = _archived()
    assert hashlib.sha256(CERTIFICATE_PATH.read_bytes()).hexdigest() == (
        EXPECTED_CERTIFICATE_FILE_SHA256
    )
    assert archived["certificate_payload_sha256"] == (
        EXPECTED_CERTIFICATE_PAYLOAD_SHA256
    )
    assert archived == build_certificate()
    assert all(verify_certificate(archived).values())
    assert all(archived["checks"].values())


def test_abstract_fold_passes_every_global_hyperbolic_gate() -> None:
    payload = _archived()
    fold = payload["fixed_nef_fold"]
    literal = payload["literal_vinberg_coxiter_prop_1_10"]

    assert payload["proposed_fold"] == {
        "ambient_picard_rank": 16,
        "ambient_wall_count": 19,
        "cycles": [[index] for index in range(1, 16)]
        + [[16, 18], [17, 19]],
        "fixed_rank": 15,
        "folded_facet_count": 17,
        "hyperbolic_dimension": 14,
        "status": "abstract_lattice_and_complete_halfspace_fold_only",
    }
    assert fold["orbit_gram_rank"] == 15
    assert fold["orbit_gram_inertia"] == [1, 14, 2]
    assert fold["orbit_normal_squares"] == [-2] * 15 + [-4, -4]
    assert fold["strict_interior"] == {
        "basis_coordinates": [
            22,
            46,
            72,
            35,
            65,
            60,
            57,
            56,
            57,
            60,
            65,
            72,
            35,
            46,
            22,
        ],
        "square": 1540,
        "wall_pairings": [2] * 15 + [44, 44],
    }
    assert all(
        facet["irredundant_with_strict_relative_interior"]
        for facet in fold["facet_audits"]
    )
    assert all(
        reflection["integral"]
        and reflection["involution"]
        and reflection["isometry"]
        for reflection in fold["folded_reflections"]
    )
    assert fold["coxeter_relation_counts"] == {
        "2": 119,
        "3": 14,
        "4": 2,
        "infinity": 1,
    }
    assert (
        fold["finite_vertex_count"],
        fold["ideal_vertex_count"],
        fold["hyperideal_vertex_count"],
    ) == (89, 5, 0)

    assert literal["elliptic_edge_count"] == 704
    assert literal["extension_histogram"] == {"2": 704}
    assert literal["endpoint_type_histogram"] == {
        "0,2": 7,
        "1,1": 148,
        "2,0": 549,
    }
    assert literal["bad_edges"] == []
    assert literal["connected_affine_component_count"] == 8
    assert literal["distinct_full_zero_neighbourhood_count"] == 5
    assert literal["full_zero_neighbourhoods_match_ideal_vertices"] is True
    assert literal["finite_volume_noncompact"] is True
    assert payload["conclusion"]["abstract_hyperbolic_coxeter_polytope"] is True
    assert payload["conclusion"]["finite_volume"] is True
    assert payload["conclusion"]["compact"] is False


def test_discriminant_obstruction_rejects_the_k3_realization() -> None:
    payload = _archived()
    ambient = payload["ambient_lattice"]
    algebraic = payload["algebraic_realization_audit"]

    assert ambient["basis_determinant"] == -4
    assert ambient["basis_smith_diagonal"] == [1] * 14 + [2, 2]
    assert algebraic["weighted_graph_automorphism_count"] == 4
    assert algebraic["weighted_graph_group"] == "V4"
    assert algebraic["proposed_action_trivial_on_A_NS"] is False
    assert algebraic["published_deck_action_trivial_on_A_NS"] is True
    assert algebraic["dual_discriminant_witness"]["difference_is_integral"] is False
    assert payload["conclusion"]["algebraic_K3_quotient_nef_realization"] is False
    assert payload["conclusion"]["new_hyperbolic_example_claim"] is False
