"""Regression tests for the terminal AGCat composition audit."""

from __future__ import annotations

from fractions import Fraction

import pytest

from scripts.certify_agcat_global_maximal_arithmetic_coxeter_finiteness_independent_audit import (
    FROZEN_INPUTS,
    build_certificate,
    field_branch,
    fork_identity,
    sha256_file,
)


INDEPENDENT_CERTIFIER = (
    FROZEN_INPUTS["composite_owner"][0].parents[1]
    / "scripts/certify_agcat_global_maximal_arithmetic_coxeter_finiteness_independent_audit.py"
)
INDEPENDENT_CERTIFIER_SHA256 = (
    "8acee124164a974157ea19a32b5a8d765ad2252f15078fd9882909089bf189c7"
)


def test_every_frozen_mathematical_input_matches_its_terminal_pin():
    for path, expected in FROZEN_INPUTS.values():
        assert sha256_file(path) == expected
    assert sha256_file(INDEPENDENT_CERTIFIER) == INDEPENDENT_CERTIFIER_SHA256


def test_field_branch_partition_is_exhaustive_for_n_at_least_two():
    with pytest.raises(ValueError):
        field_branch(1, compact=True)
    for dimension in range(2, 1501):
        assert field_branch(dimension, compact=False) == "rational_cusp"
        expected = (
            "real_dm_orbicurve"
            if dimension == 2
            else "compact_h3_general_fork"
            if dimension == 3
            else "codimension_descent"
            if dimension in (4, 5)
            else "minimality14_n_ge_6"
        )
        assert field_branch(dimension, compact=True) == expected


def test_general_fork_beta_is_the_rescaling_invariant_cycle_formula():
    for u, a, b in (
        (Fraction(5, 2), Fraction(1), Fraction(3, 2)),
        (Fraction(3), Fraction(4, 3), Fraction(5, 4)),
        (Fraction(13), Fraction(7, 5), Fraction(8, 7)),
    ):
        assert fork_identity(u, a, b)["equal"] is True


def test_dimension_to_rank_quantifier_has_no_fixed_dimension_leak():
    certificate = build_certificate()
    assert certificate["dimension_quantifier"] == {
        "hyperbolic_dimensions": [2, 995],
        "native_ranks": [3, 996],
        "rank_count": 994,
    }


def test_composite_interface_links_tex_and_dependency_firewalls():
    certificate = build_certificate()
    assert certificate["all_checks_pass"] is True
    assert certificate["verdict"] == "PASS"
    assert all(certificate["categorical_and_scope_interface"].values())
    assert all(certificate["field_branch_partition"].values())
    assert certificate["local_links"]["missing"] == []
    assert all(certificate["tex"].values())
