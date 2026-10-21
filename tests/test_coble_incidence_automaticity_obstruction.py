from scripts.certify_coble_incidence_automaticity_obstruction import (
    build_certificate,
)


def test_plane_sextic_partition_identities_and_kondo_pruning() -> None:
    certificate = build_certificate()
    assert certificate["all_checks_pass"]
    plane = certificate["plane_rational_sextic_factory"]
    assert len(plane["rows"]) == 11
    assert plane["Kondo_finite_Aut_surviving_degree_partitions"] == [
        [6],
        [5, 1],
        [4, 2],
        [3, 3],
        [2, 1, 1, 1, 1],
    ]
    for row in plane["rows"]:
        assert row["blowup_count"] == 9 + row["boundary_count"]
        assert row["K_square"] == -row["boundary_count"]
        assert set(row["component_self_intersections"]) == {-4}


def test_complete_grid_has_one_finite_volume_stabilizer() -> None:
    certificate = build_certificate()
    rows = certificate["P1xP1_complete_grid"][
        "four_point_projectivity_stabilizers_and_exact_cones"
    ]
    assert [rows[key]["graph_wall_count"] for key in rows] == [4, 8, 12, 24]
    assert [key for key, row in rows.items() if row["finite_volume"]] == [
        "characteristic_3"
    ]
    assert rows["generic"]["hyperideal"] == 2500
    assert rows["harmonic"]["hyperideal"] == 3344
    assert rows["equianharmonic"]["hyperideal"] == 5340
    assert rows["characteristic_3"]["hyperideal"] == 0


def test_tumarkin_route_misses_both_calibrations() -> None:
    certificate = build_certificate()
    rows = certificate["Tumarkin_n_plus_3_route"]["known_D1_split_CM_chambers"]
    assert rows["I"]["non_simple_ideal_vertex_count"] == 3
    assert rows["II"]["non_simple_ideal_vertex_count"] == 3

