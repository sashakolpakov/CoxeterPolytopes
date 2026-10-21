from __future__ import annotations

import json

from scripts.certify_rational_i1n_symbolic_cusp_family import (
    CASES,
    OUTPUT,
    _actual_multiplicity_gate,
    _coprimality_gates,
    _polynomial_to_string,
    _section_expression,
    _verify_symbolic_kernel,
    affine_monomials,
    integer_jet_condition_matrix,
    symbolic_jet_condition_matrix,
    validate_archive,
    verify_symbolic_cluster_jet,
)
from scripts.certify_rational_i1n_cusp_fibrations import (
    _jet_condition_matrix as original_jet_condition_matrix,
)


def _archive() -> dict[str, object]:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def _section_from_terms(
    degree: int, terms: list[list[object]]
) -> tuple[tuple[int, ...], ...]:
    monomials = affine_monomials(degree)
    lookup = {monomial: index for index, monomial in enumerate(monomials)}
    answer: list[tuple[int, ...]] = [()] * len(monomials)
    for i, j, coefficients in terms:
        answer[lookup[(int(i), int(j))]] = tuple(map(int, coefficients))
    return tuple(answer)


def test_symbolic_family_archive_and_zero_exceptional_polynomial() -> None:
    archive = validate_archive()
    assert archive["scope"]["additional_cusp_exceptional_polynomial"] == "1"
    assert archive["totals"] == {
        "additional_parameter_exceptions": 0,
        "cusps": 28,
        "exact_index_four_pencils": 1,
        "primitive_base_point_free_pencils": 27,
    }
    assert all(
        record["chart_gate_gcd"] == "1"
        for record in archive["moving_systems"].values()
    )


def test_fast_integer_builder_matches_the_original_two_fibres() -> None:
    cluster = verify_symbolic_cluster_jet()
    assert cluster["post_p15_exceptional_restrictions"] == {
        "Q": "z",
        "S": "z+1",
    }
    assert cluster["selected_S_jet_coefficients"] == {
        "3": [1],
        "7": [0, 1],
        "9": [1],
        "11": [-1, 0, 2],
        "13": [1, 5],
        "15": [1, -6, 0, 5],
    }
    cases = (
        (4, (1,) * 16),
        (11, (*((3,) * 13), 2)),
        (16, (4,) * 16),
    )
    for parameter in (0, 1):
        for degree, multiplicities in cases:
            new_monomials, new_rows, _new_sparse = integer_jet_condition_matrix(
                degree, multiplicities, parameter
            )
            old_monomials, old_rows, _old_sparse = original_jet_condition_matrix(
                degree, multiplicities, parameter
            )
            assert new_monomials == old_monomials
            assert new_rows == old_rows


def test_exact_rank_minors_and_two_chart_gates() -> None:
    archive = _archive()
    expected = {
        "genus_5": [
            ("-a - 3", "a + 3"),
            ("-3", "3*a + 1"),
        ],
        "genus_10": [("1", "1")],
        "genus_11": [
            ("a - 1", "a - 1"),
            ("-7", "7*a - 1"),
        ],
        "fourth_multiple": [
            ("-12*(a + 1)", "a + 1"),
            (
                "2*(a**2 - 10*a - 13)",
                "(2*a + 3)*(a**2 - 10*a - 13)",
            ),
        ],
    }
    for name, pairs in expected.items():
        charts = archive["moving_systems"][name]["charts"]
        assert [
            (
                chart["rank_minor_determinant_factorization"],
                chart["chart_gate_factorization"],
            )
            for chart in charts
        ] == pairs
    assert {
        name: record["rank_minor_determinant"]
        for name, record in archive["rigid_lower_multiples"].items()
    } == {"F": "1", "2F": "-2", "3F": "2"}


def test_archived_second_sections_are_symbolic_kernel_members() -> None:
    archive = _archive()
    for name, case in CASES.items():
        degree = int(case["degree"])
        multiplicities = tuple(map(int, case["multiplicities"]))
        monomials, rows, _sparse = symbolic_jet_condition_matrix(
            degree, multiplicities
        )
        for chart in archive["moving_systems"][name]["charts"]:
            section = _section_from_terms(
                degree, chart["second_section_sparse_terms"]
            )
            _verify_symbolic_kernel(rows, section)
            actual, gates = _actual_multiplicity_gate(
                section, monomials, multiplicities
            )
            assert actual == multiplicities
            assert [_polynomial_to_string(value) for value in gates] == chart[
                "actual_multiplicity_leading_gcds"
            ]
            coprimality = _coprimality_gates(
                _section_expression(section, monomials),
                case["wall_factors"].keys(),
            )
            assert {
                label: _polynomial_to_string(value)
                for label, value in coprimality.items()
            } == chart["coprimality_remainder_coefficient_gcds"]
