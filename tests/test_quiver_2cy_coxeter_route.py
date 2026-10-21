from scripts.certify_quiver_2cy_coxeter_route import build_certificate


def test_quiver_2cy_control_replay() -> None:
    certificate = build_certificate()
    assert certificate["all_checks_pass"] is True

    h21, x3 = certificate["controls"]
    assert (h21["facet_count"], h21["gram_rank"], h21["gram_corank"]) == (
        210,
        22,
        188,
    )
    assert (x3["facet_count"], x3["gram_rank"], x3["gram_corank"]) == (
        36,
        20,
        16,
    )

    assert h21["coxeter_parameter_histogram"]["2"] == 2016
    assert x3["coxeter_parameter_histogram"]["3"] == 36
