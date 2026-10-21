from scripts.certify_global_k3_root_automorphic_return import analyze


def test_global_rank_and_compactness_counts() -> None:
    report = analyze()
    assert report["marked_finite_volume_types"] == 118
    assert report["compact_types"] == 8
    assert report["finite_volume_noncompact_types"] == 110
    assert report["maximum_hyperbolic_dimension"] == 18


def test_wang_automorphic_correction_split() -> None:
    report = analyze()
    correction = report["automorphic_correction_range"]
    assert correction["types"] == 69
    assert correction["corrected"] == 51
    assert correction["uncorrected"] == 18


def test_high_rank_intersection() -> None:
    report = analyze()
    assert report["high_rank_corrected"] == 7
    assert report["high_rank_uncorrected"] == 12
    corrected = {
        row["lattice"]
        for row in report["high_rank_rows"]
        if row["automorphic_correction"]
    }
    assert corrected == {
        "U+E8+A1",
        "U+D8+A1",
        "U+2D4+A1",
        "U+D4+5A1",
        "U+E8+2A1",
        "U+2E8",
        "U+2E8+A1",
    }

