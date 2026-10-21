from scripts.certify_positive_example_lineage_header_controls import (
    ROOT,
    build_certificate,
)


def test_pinned_corpus_is_not_a_worldwide_priority_denominator() -> None:
    certificate = build_certificate()
    assert certificate["all_checks_pass"]
    assert len(certificate["known_zero_header_controls"]) == 6
    assert all(
        control["ordinary_header_count"] == 0
        and control["doubling_header_count"] == 0
        for control in certificate["known_zero_header_controls"]
    )
    assert all(
        (ROOT / control["exact_local_lineage"]).is_file()
        for control in certificate["known_zero_header_controls"]
    )


def test_ledger_header_census_and_x2_scope() -> None:
    certificate = build_certificate()
    ledger = certificate["ledger"]
    assert ledger["positive_row_count"] == 56
    assert ledger["distinct_header_count"] == 43
    assert ledger["zero_header_distinct_count"] == 30
    assert ledger["zero_header_row_count"] == 35
    assert ledger["dimension_histogram"] == {
        "H10": 29,
        "H11": 14,
        "H12": 6,
        "H13": 2,
        "H14": 3,
        "H19": 1,
        "H21": 1,
    }
    x2 = certificate["x2_strengthened_bounded_audit"]
    assert x2["same_node_records_checked"] == 45
    assert x2["histogram_match_count"] == 0
    assert x2["worldwide_priority"] == "open"
