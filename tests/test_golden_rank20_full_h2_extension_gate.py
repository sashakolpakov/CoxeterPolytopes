from scripts.certify_golden_rank20_full_h2_extension_gate import analyze


def test_golden_rank20_extension_gate() -> None:
    result = analyze()
    assert len(result["packages"]) == 2

    packages = {row["candidate_id"]: row for row in result["packages"]}
    d11 = packages["golden_discriminant_11_rank20"]
    d19 = packages["golden_discriminant_19_rank20"]

    assert d11["complement_gram"] == [[18, 5], [5, 2]]
    assert d11["gluing_basis_determinant"] == -11
    assert d19["complement_gram"] == [[2, 1], [1, 10]]
    assert d19["gluing_basis_determinant"] == -19

    for package in packages.values():
        assert package["integral_complement_solution_count"] == 4
        assert package["z_sqrt5_full_lattice_extension_count"] == 2
        assert package["z_phi_full_lattice_extension_count"] == 0
        assert package["transcendental_rank_over_q_sqrt5"] == 1
        assert not package["hodge_coherent_extension_exists"]


def test_exact_complement_solution_lists() -> None:
    packages = {row["candidate_id"]: row for row in analyze()["packages"]}
    d11_actions = [
        row["complement_sqrt5_action"]
        for row in packages["golden_discriminant_11_rank20"][
            "integral_complement_solutions"
        ]
    ]
    d19_actions = [
        row["complement_sqrt5_action"]
        for row in packages["golden_discriminant_19_rank20"][
            "integral_complement_solutions"
        ]
    ]

    assert d11_actions == [
        [[-4, -1], [11, 4]],
        [[-1, -1], [-4, 1]],
        [[1, 1], [4, -1]],
        [[4, 1], [-11, -4]],
    ]
    assert d19_actions == [
        [[-1, 4], [1, 1]],
        [[0, -5], [-1, 0]],
        [[0, 5], [1, 0]],
        [[1, -4], [-1, -1]],
    ]
