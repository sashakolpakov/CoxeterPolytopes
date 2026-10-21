from scripts.certify_coble_prime_binomial_k3_rank_pilot import (
    build_certificate,
    ceiling,
    certify,
    k3_rank,
    polytope_dimension,
)


def test_full_hostile_audit_certificate():
    certificate = build_certificate()
    certify(certificate)
    assert certificate["all_checks_pass"]
    assert not certificate["missing_lemma"]["proved"]
    assert not certificate["verdict"]["strict_L4"]


def test_canonical_cover_rank_pushback():
    assert k3_rank(10) == 20
    assert polytope_dimension(10) == 19
    assert k3_rank(12) == 22
    assert polytope_dimension(12) == 21
    assert k3_rank(11) > ceiling(0)
    assert k3_rank(13) > ceiling(3)


def test_four_controls_are_separately_routed_not_one_class():
    certificate = build_certificate()
    rows = certificate["separately_routed_controls"]
    assert {row["id"] for row in rows} == {
        "H10_30",
        "H11_26",
        "H14_60",
        "H17_48",
    }
    assert not certificate["verdict"][
        "four_rows_pass_one_invariant_raw_prime_binomial_recognizer"
    ]
    assert not certificate["operational_audit"][
        "candidate_has_arbitrary_input_option"
    ]
    assert certificate["operational_audit"][
        "candidate_hardwired_source_constants"
    ] == ["CI4", "H10", "H11", "HV"]


def test_prime_binomial_alone_has_exact_false_positives():
    certificate = build_certificate()
    failures = certificate["prime_binomial_false_positives"]
    assert {row["component"] for row in failures} == {"c=1", "a=-b"}
    assert all(row["prime_binomial"] for row in failures)


def test_no_synthetic_pre_ag_delta_is_counted():
    certificate = build_certificate()
    assert not certificate["verdict"][
        "pre_AG_allowed_b_values_wider_than_K3_bound_proved"
    ]
    generic = certificate["only_current_generic_raw_class"]
    assert max(generic["boundary_counts"]) == 8
    assert max(generic["dimensions"]) == 17
    assert generic["K3_rank_pushback_removes_recognized_members"] == 0
