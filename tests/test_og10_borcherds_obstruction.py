from scripts.certify_og10_borcherds_obstruction import certify


def test_og10_borcherds_lattice_certificate():
    result = certify()
    assert result["verified"]
    assert result["ordinary_type_iv"]["signature"] == {
        "positive": 22,
        "negative": 2,
        "zero": 0,
    }
    assert result["ordinary_type_iv"]["smith_nonunits"] == [2, 22]
    assert result["ordinary_type_iv"]["level"] == 44
    assert result["totally_real_extension"]["field_signatures"] == [
        [10, 2],
        [12, 0],
    ]
    assert result["totally_real_extension"]["trace_smith_nonunits"] == [
        2,
        2,
        2,
        22,
    ]
    assert result["totally_real_extension"]["whittaker_input_weight"] == [
        -4,
        6,
    ]
    assert result["totally_real_extension"]["obstruction_cusp_weight"] == [
        6,
        6,
    ]
