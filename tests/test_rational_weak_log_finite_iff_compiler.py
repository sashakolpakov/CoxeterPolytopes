from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import pytest

from scripts.certify_rational_weak_log_finite_iff_compiler import (
    ARCHIVE,
    _standard_p9_cusp,
    build,
    ceil_fraction,
    delta246_mod3_descent_certificate,
    delta246_scaled_form,
    determinant_3x3,
    fujino_multiplier,
    lorentz_dot,
    physical_square_assignment_bound,
    positive_cusp_relation,
    sequential_program_bound,
    sha256_json,
    weak_log_big_cartier_multiplier,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_weak_log_finite_iff_compiler.md"


def frozen() -> dict[str, object]:
    return json.loads(ARCHIVE.read_text())


def test_live_certificate_matches_frozen_archive() -> None:
    assert build() == frozen()
    payload = frozen().copy()
    digest = payload.pop("content_sha256")
    assert digest == sha256_json(payload)


def test_effective_big_branch_has_one_fixed_multiplier() -> None:
    assert ceil_fraction(Fraction(0)) == 0
    assert ceil_fraction(Fraction(1, 3)) == 1
    assert fujino_multiplier(2, 0) == 96
    assert weak_log_big_cartier_multiplier() == 192
    row = frozen()["nef_routes"]["big"]  # type: ignore[index]
    assert row["cartier_divisor"] == "L_D=2D"
    assert row["fixed_basepoint_free_system"] == "|96(2D)|=|192D|"
    assert row["multiplier_of_D"] == 192


def test_effective_multiplier_rejects_hostile_numeric_coercions() -> None:
    with pytest.raises(ValueError):
        fujino_multiplier(True, 0)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        fujino_multiplier(-1, 0)
    with pytest.raises(ValueError):
        fujino_multiplier(2, Fraction(-1, 2))


def test_standard_anticanonical_cusp_is_positive_affine_divisor() -> None:
    target, roots, names = _standard_p9_cusp()
    coefficients = positive_cusp_relation(target, roots)
    assert dict(zip(names, coefficients)) == {
        "d1": 2,
        "d2": 4,
        "d3": 6,
        "d4": 5,
        "d5": 4,
        "d6": 3,
        "d7": 2,
        "d8": 1,
        "c": 3,
    }
    assert lorentz_dot(target, target) == 0
    assert all(lorentz_dot(target, root) == 0 for root in roots)
    assert tuple(
        sum(coefficients[j] * roots[j][row] for j in range(len(roots)))
        for row in range(len(target))
    ) == target


def test_cusp_relation_is_invariant_under_exceptional_coordinate_permutation() -> None:
    target, roots, _ = _standard_p9_cusp()
    permutation = [0, 3, 2, 1, 4, 5, 6, 7, 8, 9]
    moved_target = tuple(target[index] for index in permutation)
    moved_roots = [tuple(root[index] for index in permutation) for root in roots]
    assert positive_cusp_relation(moved_target, moved_roots) == (
        2,
        4,
        6,
        5,
        4,
        3,
        2,
        1,
        3,
    )


def test_cusp_gate_rejects_wrong_orientation_and_nonincident_roots() -> None:
    target, roots, _ = _standard_p9_cusp()
    with pytest.raises(ValueError, match="strictly positive"):
        positive_cusp_relation(tuple(-entry for entry in target), roots)
    broken = list(roots)
    broken[0] = tuple(
        entry + int(index == 0) for index, entry in enumerate(broken[0])
    )
    with pytest.raises(ValueError, match="incident"):
        positive_cusp_relation(target, broken)


def test_finite_host_and_mmp_bounds_are_literal() -> None:
    assert physical_square_assignment_bound(12) == 531441
    assert sequential_program_bound(12, 11) == 344058145
    for bad in (True, -1, 2.5):
        with pytest.raises(ValueError):
            physical_square_assignment_bound(bad)  # type: ignore[arg-type]


def test_p10_and_p11_geometric_meanings_are_not_conflated() -> None:
    routes = frozen()["p10_p11_routing"]  # type: ignore[index]
    assert routes["decorated_P10"]["weak_log_compiler_expected_answer_over_C"] == (
        "EMPTY"
    )
    assert routes["decorated_P10"]["maximality_forces_prime_effectivity"] is False
    assert routes["uncoloured_P10"]["known_exact_flex_carrier"] is True
    assert routes["uncoloured_P11"]["known_exact_flex_carrier"] is True
    assert routes["uncoloured_P10"]["weak_log_smooth_rational_species"] is False
    assert routes["uncoloured_P11"]["weak_log_smooth_rational_species"] is False


def test_raw_maximal_nonrational_triangle_kills_every_rational_log_return() -> None:
    row = frozen()["raw_maximal_rational_log_no_go"]  # type: ignore[index]
    assert row["coxeter_triangle"] == [2, 5, 6]
    assert row["compact"] is True
    assert row["absolute_reflection_maximality"] is True
    assert row["minimal_vinberg_field"] == "Q(sqrt(5))"
    assert row["literal_rational_full_nef_field"] == "Q"
    assert row["every_rational_or_rational_log_full_nef_carrier_empty"] is True
    assert row["other_carrier_species_excluded"] is False
    assert (ROOT / row["companion_report"]).is_file()


def test_raw_maximal_rational_246_is_rejected_by_the_host_gate() -> None:
    row = frozen()["raw_maximal_rational_field_no_go"]  # type: ignore[index]
    assert row["coxeter_triangle"] == [2, 4, 6]
    assert row["minimal_vinberg_field"] == "Q"
    assert row["compact"] is True
    assert row["absolute_reflection_maximality"] is True
    assert determinant_3x3(row["gram"]) == row["gram_determinant"] == 36
    assert row["anisotropy_certificate"] == delta246_mod3_descent_certificate()
    assert row["anisotropy_certificate"]["divided_binary_zero_residues"] == [[0, 0]]
    assert row["rank3_rational_picard_null_vector"] == [1, 1, 0]
    assert row["host_gate"].startswith("EMPTY")
    assert row["geometric_gate"].startswith("EMPTY")
    assert row["every_rational_or_rational_log_full_nef_carrier_empty"] is True
    assert (ROOT / row["companion_report"]).is_file()


def test_246_mod3_descent_is_exhaustive_and_public_form_is_strict() -> None:
    residues = [
        (x, y, z)
        for x in range(3)
        for y in range(3)
        for z in range(3)
        if delta246_scaled_form(x, y, z) % 3 == 0
    ]
    assert all(z == 0 for _, _, z in residues)
    assert [
        (x, y)
        for x in range(3)
        for y in range(3)
        if (2 * x * x + x * y + y * y) % 3 == 0
    ] == [(0, 0)]
    with pytest.raises(ValueError):
        delta246_scaled_form(True, 0, 0)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        determinant_3x3([[1, 0], [0, 1]])


def test_report_states_iff_but_not_universal_property4() -> None:
    text = REPORT.read_text()
    for required in (
        "finite iff/`EMPTY` compiler",
        "positive cusp divisor",
        "|192D|",
        "RLogComp` iff/`EMPTY` theorem",
        "Gr\\\"obner elimination",
        "Decorated adapted `P10`",
        "Raw-maximal chambers with every rational/log return empty",
        "rational-field example",
        "does **not** assert",
    ):
        assert required in text
    scope = frozen()["scope"]  # type: ignore[index]
    assert scope["finite_iff_criterion_proved"] is True
    assert scope["universal_property4_claimed"] is False
    assert scope["universal_nonemptiness_claimed"] is False
    assert scope["hostile_direct_full_ns_lane_closed"] is False
    assert scope["generic_groebner_elimination_implemented"] is False
    assert scope["maximality_used_to_create_sections"] is False
    assert scope["maximality_used_only_after_prime_and_nef_gates"] is True
