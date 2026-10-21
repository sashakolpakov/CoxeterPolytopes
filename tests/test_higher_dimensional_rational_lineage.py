from scripts.certify_higher_dimensional_rational_lineage import build_certificate


def test_rational_lineage_certificate() -> None:
    certificate = build_certificate()
    assert certificate["all_checks_pass"]


def test_jacobian_kummer_commensurability_conclusions() -> None:
    comparisons = build_certificate()["comparisons"]
    assert comparisons["JK-H10-0007_vs_odd_unimodular_H10"]["rationally_similar"]
    assert comparisons["JK-H12_vs_odd_unimodular_H12"]["rationally_similar"]
    assert comparisons["JK-H12_vs_McLeod_f3_H12"]["rationally_similar"]
    assert not comparisons["K3-D4-H12_vs_odd_unimodular_H12"]["rationally_similar"]
