import json

from scripts.independent_certify_k3_rho15_h12_cyclic_fold import (
    CERTIFICATE_PATH,
    EXPECTED_F_VECTOR,
    _canonical_sha256,
)


def _archived() -> dict:
    return json.loads(CERTIFICATE_PATH.read_text())


def test_independent_h12_archive_is_self_authenticating() -> None:
    payload = _archived()
    body = dict(payload)
    digest = body.pop("certificate_payload_sha256")
    assert digest == _canonical_sha256(body)
    assert all(payload["checks"].values())


def test_order_three_fold_passes_the_exact_volume_gate() -> None:
    payload = _archived()
    fold = payload["order_three_H12_fold"]
    literal = payload["literal_vinberg_coxiter_prop_1_10"]

    assert fold["hyperbolic_dimension"] == 12
    assert fold["facet_count"] == 17
    assert fold["orbit_gram_rank"] == 13
    assert fold["orbit_gram_inertia"] == [1, 12, 4]
    assert fold["orbit_normal_squares"] == [-2] * 13 + [-6, -6, -2, -2]
    assert (
        fold["finite_vertex_count"],
        fold["ideal_vertex_count"],
        fold["hyperideal_vertex_count"],
    ) == (93, 3, 0)
    assert fold["coxeter_relation_counts"] == {
        "2": 116,
        "3": 12,
        "6": 2,
        "infinity": 3,
        "ultraparallel": 3,
    }

    assert literal["elliptic_edge_count"] == 623
    assert literal["extension_histogram"] == {"2": 623}
    assert literal["endpoint_type_histogram"] == {
        "0,2": 3,
        "1,1": 124,
        "2,0": 496,
    }
    assert literal["bad_edges"] == []
    assert literal["finite_volume_noncompact"] is True


def test_full_face_lattice_and_order_six_rejection() -> None:
    payload = _archived()
    faces = payload["complete_face_lattice"]
    census = payload["weighted_wall_graph_cyclic_census"]

    assert faces["facet_subset_count"] == 2**17
    assert faces["nonempty_face_count"] == 34733
    assert faces["f_vector_vertices_through_polytope"] == EXPECTED_F_VECTOR
    assert faces["alternating_euler_sum"] == 1

    assert census["group_identification"] == "S3 x C2"
    assert census["element_order_histogram"] == {
        "1": 1,
        "2": 7,
        "3": 2,
        "6": 2,
    }
    assert census["order_three_elements_are_inverse_with_same_partition"]
    assert all(
        item["valid_spacelike_Coxeter_wall_partition"] is False
        and item["orbit_gram_inertia"] == [1, 6, 2]
        and item["orbit_normal_squares"][-2:] == [0, 8]
        for item in census["order_six_elements"]
    )
