from coxeter_polytopes.search.niemeier_host_triage import (
    UNSEARCHED_ROOTED_NIEMEIER_HOSTS,
    highest_root_marks,
    niemeier_triage,
)


def test_inventory_is_exactly_the_nineteen_requested_unscreened_hosts():
    names = {host.name for host in UNSEARCHED_ROOTED_NIEMEIER_HOSTS}
    assert len(names) == 19
    assert not names & {"Leech", "A1^24", "D24", "E8^3", "A24"}
    assert all(sum(rank for _family, rank in host.components) == 24 for host in UNSEARCHED_ROOTED_NIEMEIER_HOSTS)


def test_affine_marks_and_ranked_front_are_frozen():
    assert highest_root_marks(("E", 6)) == (1, 2, 3, 2, 1, 2)
    assert highest_root_marks(("E", 7)) == (2, 3, 4, 3, 2, 1, 2)
    assert highest_root_marks(("E", 8)) == (2, 4, 6, 5, 4, 3, 2, 3)
    rows = niemeier_triage()
    assert [row.host for row in rows[:4]] == ["A17 E7", "A15 D9", "D12^2", "A12^2"]
    assert [(row.affine_product_ray_count, row.affine_finite_count, row.affine_hyperideal_count) for row in rows[:4]] == [
        (144, 2, 142),
        (160, 2, 158),
        (169, 3, 166),
        (169, 1, 168),
    ]
    assert all(row.affine_cusp_count == 0 for row in rows)
    assert all(row.glue_index * row.glue_index == row.root_lattice_determinant for row in rows)
    assert rows[-1].host == "D16 E8"


def test_e6_tetracode_triage_row_matches_full_backend():
    row = next(row for row in niemeier_triage() if row.host == "E6^4")
    assert row.glue_index == 9
    assert row.minimum_coupling_support == 3
    assert (row.affine_product_ray_count, row.affine_finite_count, row.affine_hyperideal_count) == (2401, 5, 2396)

