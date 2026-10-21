from scripts.certify_supersingular_golden_h10_obstruction import certify


def test_supersingular_golden_h10_obstruction_certificate():
    result = certify()
    assert result["verified"]
    assert result["conclusion"] == (
        "no_supersingular_K3_NS_has_a_self_adjoint_golden_action"
    )
    assert result["two_trace_functionals_independent"]
    assert not any(
        row["compatible"] for row in result["characteristic_two_parity_table"]
    )
    assert result["naive_rank11_trace_lattice"]["nikulin_invariants"] == [22, 6, 1]
