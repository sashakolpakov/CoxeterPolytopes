from scripts.verify_k3_high_rank_nef_chambers import analyze


def test_high_rank_k3_wall_grams() -> None:
    report = analyze()
    assert report["artifact_status"] == "source_audited_gram_verified"
    assert report["ahl_rank_15_20_boundary_recorded"] is True
    assert report["all_checks_pass"] is True

    cases = {case["rho"]: case for case in report["cases"]}
    assert set(cases) == {15, 16, 17, 18, 19}
    assert [cases[rho]["wall_count"] for rho in range(15, 20)] == [21, 19, 19, 19, 24]
    assert [cases[rho]["rank"] for rho in range(15, 20)] == [15, 16, 17, 18, 19]
    assert [cases[rho]["inertia"] for rho in range(15, 20)] == [
        [1, 14, 6],
        [1, 15, 3],
        [1, 16, 2],
        [1, 17, 1],
        [1, 18, 5],
    ]
    assert [cases[rho]["basis_determinant"] for rho in range(15, 20)] == [
        8,
        -4,
        2,
        -1,
        2,
    ]


def test_rank15_redraw_error_is_detected_exactly() -> None:
    report = analyze()
    rank15 = next(case for case in report["cases"] if case["rho"] == 15)
    literal = rank15["literal_rank15_redraw"]
    assert literal["rank"] == 19
    assert literal["inertia"] == [3, 16, 2]
    assert literal["checks"]["cannot_embed_in_claimed_ns"] is True


def test_derived_rank19_inose_equation() -> None:
    report = analyze()
    model = report["derived_rank19_inose_model"]
    assert model["status"] == "derived_and_symbolically_verified_not_roulleau_quoted"
    assert model["local_valuations_at_zero"] == {
        "a4": 4,
        "a6": 5,
        "c4": 4,
        "discriminant": 10,
    }
    assert model["local_valuations_at_infinity"] == {
        "a4": 4,
        "a6": 5,
        "c4": 4,
        "discriminant": 10,
    }
    assert model["generic_fibers"] == ["II*", "II*", "I2", "I1", "I1"]
    assert model["fiber_euler_sum"] == 24
    assert model["trivial_lattice_rank"] == 19
    assert model["all_checks_pass"] is True
