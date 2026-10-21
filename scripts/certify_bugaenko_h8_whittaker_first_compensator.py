#!/usr/bin/env python3
"""Certify the first complete Bugaenko H8 Whittaker compensator shell.

For the determinant-line extension

    Hhat = <-2*phi> + E8 + <-4*phi>,
    Q = phi*Hhat/(2*sqrt(5)),

the discriminant module is A=O/(2)+O/(4).  A class is represented by

    mu(x,y) = (x/2)u + (y/4)e.

The auxiliary coordinate y is decisive on the canonical hyperbolic slice:
if y is nonzero modulo 4, every vector in mu+L has nonzero e-coordinate,
so the complete Heegner divisor is disjoint from that slice.

This script proves two bounded statements.

1. Among all 64 discriminant classes, the globally smallest totally
   positive admissible index has norm 1/80.  Its complete positive-unit
   shell consists of one three-component orbit, represented by

       m1=phi/(4*sqrt(5)),  mu1=(phi,2*phi^2).

   Every component is slice-disjoint, and an explicit vector proves that
   the special divisor is nonempty away from the slice.
2. The weight-(5,5) Hilbert--Weil Poincare series indexed by the target and
   this complete minimal shell are linearly independent.  Hence no
   weakly-holomorphic principal part supported on these two divisor orbits
   can retain a nonzero target coefficient.

The second assertion is a certified 2x2 diagonal-dominance calculation.
Finite unit Bessel sums use Arb.  Infinite unit and ideal tails are bounded
by exact rational inequalities.  Nothing here bounds larger compensator
support.
"""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from fractions import Fraction
from math import factorial
from pathlib import Path

from flint import arb, ctx

from coxeter_polytopes.exact.golden import Golden, PHI
from scripts.certify_bugaenko_whittaker_gate import (
    IDEAL_NORM_CUTOFF,
    PI_UPPER,
    R,
    SQRT5_LOWER,
    ideal_count,
    nonunit_ideal_sum_bound,
    unit_bessel_orbit,
)


ZERO = Golden(0)
ONE = Golden(1)
SQRT5 = 2 * PHI - 1

DISCRIMINANT_ORDER = 64
INDEX_RESIDUE_DENOMINATOR = 8

M_TARGET = PHI / SQRT5
W_MINIMAL = 2 * PHI
M_MINIMAL = W_MINIMAL / (INDEX_RESIDUE_DENOMINATOR * SQRT5)

MOD2_REPRESENTATIVES = tuple(
    Golden(a, b) for a in range(2) for b in range(2)
)
MOD4_REPRESENTATIVES = tuple(
    Golden(a, b) for a in range(4) for b in range(4)
)

# mu=(x/2)u+(y/4)e.  This representative satisfies the congruence for
# w=8*sqrt(5)*m=2*phi.
MINIMAL_COSET = (PHI, 2 * PHI**2)

UNIT_SUM_THRESHOLDS = {
    "minimal_diagonal": Fraction(37, 100_000),
    "target_minimal_cross": Fraction(27, 1000),
}

# For Bessel order four the target small-argument prefactor is <13/5.
# If r=N(m_minimal)/N(m_target)=1/16, the minimal diagonal and cross
# prefactors scale respectively by r^2 and r.
NONUNIT_BESSEL_PREFACTOR_UPPER = {
    "minimal_diagonal": Fraction(1, 96),
    "target_minimal_cross": Fraction(1, 6),
}

# This is deliberately retained even though every component of the first
# shell is self-opposite.  It makes the estimate independent of the choice
# between vector-valued rows and rows quotiented by the metaplectic center.
NONZERO_COMPONENT_CENTER_MULTIPLICITY = 2

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
BASE_CERTIFIER = REPOSITORY_ROOT / "scripts/certify_bugaenko_whittaker_gate.py"
BASE_CERTIFICATE = (
    REPOSITORY_ROOT / "data/candidates/bugaenko_h8_whittaker_gate.json"
)


def _decimal(value: Fraction, digits: int = 15) -> str:
    scale = 10**digits
    integer = value.numerator * scale // value.denominator
    return f"{integer // scale}.{integer % scale:0{digits}d}"


def _arb_upper_string(value: arb) -> str:
    return str(value.upper())


def _fraction_sha256(value: Fraction) -> str:
    return hashlib.sha256(str(value).encode("ascii")).hexdigest()


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _residue(value: Golden, modulus: int) -> tuple[int, int]:
    if not value.is_integral:
        raise ValueError(f"expected an integral golden value, got {value}")
    return (int(value.a) % modulus, int(value.b) % modulus)


def _reduce(value: Golden, modulus: int) -> Golden:
    a, b = _residue(value, modulus)
    return Golden(a, b)


def _q_residue_numerator(x: Golden, y: Golden) -> Golden:
    """Return 8*sqrt(5)*Q(mu(x,y)) modulo 8O."""

    return -(PHI**2) * (2 * x**2 + y**2)


def _all_discriminant_residues() -> dict[
    tuple[int, int], list[tuple[Golden, Golden]]
]:
    residues: dict[
        tuple[int, int], list[tuple[Golden, Golden]]
    ] = defaultdict(list)
    for x in MOD2_REPRESENTATIVES:
        for y in MOD4_REPRESENTATIVES:
            residues[_residue(_q_residue_numerator(x, y), 8)].append((x, y))
    return dict(residues)


def _unit_action(
    coset: tuple[Golden, Golden],
) -> tuple[Golden, Golden]:
    """Multiply mu by phi, while the corresponding index scales by phi^2."""

    x, y = coset
    return (_reduce(PHI * x, 2), _reduce(PHI * y, 4))


def _negate_coset(
    coset: tuple[Golden, Golden],
) -> tuple[Golden, Golden]:
    x, y = coset
    return (_reduce(-x, 2), _reduce(-y, 4))


def _unit_coset_orbit(
    seed: tuple[Golden, Golden],
) -> tuple[tuple[Golden, Golden], ...]:
    orbit: list[tuple[Golden, Golden]] = []
    current = seed
    while current not in orbit:
        orbit.append(current)
        current = _unit_action(current)
    return tuple(orbit)


def _minimal_shell_data() -> dict[str, object]:
    residues = _all_discriminant_residues()
    residue_keys = set(residues)
    norm_residues = {
        (a * a + a * b - b * b) % 8 for a, b in residue_keys
    }

    # Every totally positive m gives w=8*sqrt(5)*m of signs (+,-), hence
    # N(w)<0.  The residue set excludes N(w)=-1,-2,-3.  The first possible
    # magnitude is four and w=2*phi realizes it.
    smaller_negative_norm_residues = {
        magnitude: (-magnitude) % 8 for magnitude in range(1, 4)
    }

    w_residue_orbit: list[tuple[int, int]] = []
    current_w = W_MINIMAL
    while _residue(current_w, 8) not in w_residue_orbit:
        w_residue_orbit.append(_residue(current_w, 8))
        current_w *= PHI**2

    minimal_cosets = {
        coset
        for residue in w_residue_orbit
        for coset in residues.get(residue, ())
    }
    coset_orbit = set(_unit_coset_orbit(MINIMAL_COSET))

    return {
        "residues": residues,
        "residue_keys": residue_keys,
        "norm_residues": norm_residues,
        "smaller_negative_norm_residues": smaller_negative_norm_residues,
        "w_residue_orbit": tuple(w_residue_orbit),
        "minimal_cosets": minimal_cosets,
        "coset_orbit": coset_orbit,
    }


def _unit_bessel_sum(
    left_embeddings: tuple[arb, arb],
    right_embeddings: tuple[arb, arb],
    cutoff: int = 10,
) -> tuple[arb, arb, arb]:
    """Enclose a complete positive-unit orbit of J_4-products."""

    pi = arb.pi()
    phi = (1 + arb(5).sqrt()) / 2
    base_left = 4 * pi * (left_embeddings[0] * right_embeddings[0]).sqrt()
    base_right = 4 * pi * (left_embeddings[1] * right_embeddings[1]).sqrt()

    finite = arb(0)
    for exponent in range(-cutoff, cutoff + 1):
        first = (base_left * phi**exponent).bessel_j(4)
        second = (base_right * phi ** (-exponent)).bessel_j(4)
        finite += abs(first * second)

    left_small = (base_left / 2) ** 4 / factorial(4)
    right_small = (base_right / 2) ** 4 / factorial(4)
    tail = (
        (left_small + right_small)
        * phi ** (-4 * (cutoff + 1))
        / (1 - phi**-4)
    )
    return finite, tail, finite + tail


def _generic_nonunit_orbit_bound(norm: int, a_upper: Fraction) -> Fraction:
    """Bound a Bessel unit orbit for a lower-left ideal of norm >=4."""

    if norm < 4:
        raise ValueError("the first nonunit ideal norm is four")
    a_value = a_upper / norm**2
    return min(
        (2 * cutoff + 1) * a_value**2
        + 2
        * a_value
        * R ** (4 * cutoff + 2)
        / (1 - R**4)
        for cutoff in range(0, 2 * norm.bit_length() + 3)
    )


def _generic_nonunit_ideal_sum(
    a_upper: Fraction,
) -> tuple[Fraction, Fraction, Fraction]:
    """Bound all nonunit rows using only unitarity of the Weil action."""

    finite = sum(
        (
            ideal_count(norm)
            * _generic_nonunit_orbit_bound(norm, a_upper)
            for norm in range(4, IDEAL_NORM_CUTOFF + 1)
        ),
        Fraction(0),
    )

    # This is the order-four tail used in the frozen target certificate,
    # now with a variable archimedean prefactor and no Weil saving.
    tail = Fraction(1, IDEAL_NORM_CUTOFF**2) * (
        Fraction(49, 4) * a_upper**2
        + a_upper * R**2 / (1 - R**4)
    )
    return finite, tail, finite + tail


def _coset_json(coset: tuple[Golden, Golden]) -> dict[str, str]:
    x, y = coset
    return {"x_mod_2": str(x), "auxiliary_y_mod_4": str(y)}


def certify() -> dict[str, object]:
    shell = _minimal_shell_data()
    ctx.prec = 220

    sqrt5 = arb(5).sqrt()
    phi = (1 + sqrt5) / 2
    target_embeddings = (phi / sqrt5, 1 / (phi * sqrt5))
    minimal_embeddings = (
        target_embeddings[0] / 4,
        target_embeddings[1] / 4,
    )

    bessel_sums = {
        "minimal_diagonal": _unit_bessel_sum(
            minimal_embeddings, minimal_embeddings
        ),
        "target_minimal_cross": _unit_bessel_sum(
            target_embeddings, minimal_embeddings
        ),
    }
    bessel_checks = {
        name: total < arb(threshold.numerator) / threshold.denominator
        for name, (_, _, total) in bessel_sums.items()
        for threshold in [UNIT_SUM_THRESHOLDS[name]]
    }

    # Reproduce the frozen sharp target-diagonal estimate, including the
    # exact cancellation of its ideal-norm-four Kloosterman row.
    _, _, target_unit_sum = unit_bessel_orbit()
    _, _, target_nonunit_sum = nonunit_ideal_sum_bound()
    scalar_prefactor_arb = 4 * arb.pi() ** 2 / sqrt5
    scalar_prefactor_upper = 4 * PI_UPPER**2 / SQRT5_LOWER
    target_diagonal_error = (
        scalar_prefactor_arb * target_unit_sum / 8
        + arb(
            (scalar_prefactor_upper * target_nonunit_sum).numerator
        )
        / (scalar_prefactor_upper * target_nonunit_sum).denominator
    )

    nonunit_data = {
        name: _generic_nonunit_ideal_sum(prefactor)
        for name, prefactor in NONUNIT_BESSEL_PREFACTOR_UPPER.items()
    }
    minimal_diagonal_error = NONZERO_COMPONENT_CENTER_MULTIPLICITY * (
        scalar_prefactor_upper
        / 8
        * UNIT_SUM_THRESHOLDS["minimal_diagonal"]
        + scalar_prefactor_upper * nonunit_data["minimal_diagonal"][2]
    )
    target_minimal_cross = NONZERO_COMPONENT_CENTER_MULTIPLICITY * (
        scalar_prefactor_upper
        / 8
        * UNIT_SUM_THRESHOLDS["target_minimal_cross"]
        + scalar_prefactor_upper * nonunit_data["target_minimal_cross"][2]
    )

    # If C_ij is the Fourier coefficient matrix, conjugation by
    # diag(N(m_i)^2) removes the weight-five norm-ratio prefactor.
    target_row_margin = 1 - Fraction(3, 4) - Fraction(1, 6)
    minimal_row_margin = 1 - Fraction(1, 500) - Fraction(1, 6)

    x_minimal, y_minimal = MINIMAL_COSET
    explicit_vector_h_value = (
        -2 * PHI * (x_minimal / 2) ** 2
        + 2 * (PHI**2) ** 2
        - 4 * PHI * (y_minimal / 4) ** 2
    )
    minimal_orbit = tuple(_unit_coset_orbit(MINIMAL_COSET))
    q_difference = (
        W_MINIMAL - _q_residue_numerator(*MINIMAL_COSET)
    ) / 8
    ratio = M_MINIMAL.norm() / M_TARGET.norm()

    checks = {
        "all_64_discriminant_classes_enumerated": sum(
            len(cosets) for cosets in shell["residues"].values()
        )
        == DISCRIMINANT_ORDER,
        "forty_distinct_quadratic_residues": len(shell["residue_keys"])
        == 40,
        "allowed_norm_residues_mod_8_are_0_1_3_4": shell["norm_residues"]
        == {0, 1, 3, 4},
        "negative_norm_magnitudes_1_2_3_excluded": all(
            residue not in shell["norm_residues"]
            for residue in shell["smaller_negative_norm_residues"].values()
        ),
        "minimal_w_has_norm_minus_4_and_correct_signs": W_MINIMAL.norm()
        == -4
        and W_MINIMAL.sign(0) > 0
        and W_MINIMAL.sign(1) < 0,
        "unique_ideal_of_norm_4": ideal_count(4) == 1,
        "positive_norm_minus_4_elements_form_one_unit_orbit": tuple(
            shell["w_residue_orbit"]
        )
        == ((0, 2), (2, 4), (6, 2)),
        "minimal_index_is_admissible": q_difference.is_integral,
        "minimal_index_is_totally_positive": M_MINIMAL.sign(0) > 0
        and M_MINIMAL.sign(1) > 0,
        "minimal_index_norm_is_1_over_80": M_MINIMAL.norm()
        == Fraction(1, 80),
        "minimal_shell_is_one_three_component_unit_orbit": shell[
            "minimal_cosets"
        ]
        == shell["coset_orbit"]
        and len(shell["coset_orbit"]) == 3,
        "minimal_components_are_self_opposite": all(
            coset == _negate_coset(coset) for coset in minimal_orbit
        ),
        "minimal_shell_has_nonzero_auxiliary_coordinate": all(
            _residue(coset[1], 4) != (0, 0) for coset in minimal_orbit
        ),
        "minimal_special_divisor_is_nonempty": explicit_vector_h_value
        == Fraction(1, 2),
        "minimal_vector_has_Q_equal_to_m1": (
            PHI * explicit_vector_h_value / (2 * SQRT5)
        )
        == M_MINIMAL,
        "minimal_to_target_norm_ratio_is_1_over_16": ratio
        == Fraction(1, 16),
        "order_four_archimedean_prefactor_majorants": Fraction(13, 5)
        * ratio**2
        < NONUNIT_BESSEL_PREFACTOR_UPPER["minimal_diagonal"]
        and Fraction(13, 5) * ratio
        < NONUNIT_BESSEL_PREFACTOR_UPPER["target_minimal_cross"],
        **{f"arb_unit_sum_{name}": value for name, value in bessel_checks.items()},
        "frozen_target_diagonal_error_below_3_over_4": target_diagonal_error
        < arb(3) / 4,
        "minimal_diagonal_error_below_1_over_500": minimal_diagonal_error
        < Fraction(1, 500),
        "both_cross_errors_below_1_over_6": target_minimal_cross
        < Fraction(1, 6),
        "scaled_target_row_strictly_diagonally_dominant": target_row_margin
        > 0,
        "scaled_minimal_row_strictly_diagonally_dominant": minimal_row_margin
        > 0,
    }

    return {
        "analysis_id": "bugaenko_h8_whittaker_first_compensator_v1",
        "verified": all(checks.values()),
        "theorem": (
            "the target and the complete globally minimal compensator-shell "
            "Poincare series are linearly independent"
        ),
        "checks": checks,
        "discriminant_module": {
            "decomposition": "O/(2) + O/(4)_aux",
            "order": DISCRIMINANT_ORDER,
            "coset_coordinates": "mu(x,y)=(x/2)u+(y/4)e",
            "sqrt5_Q_mod_O": "-phi^2*x^2/4-phi^2*y^2/8",
            "slice_disjointness_criterion": "y != 0 mod 4",
        },
        "complete_minimal_shell": {
            "w_equals_8_sqrt5_m": str(W_MINIMAL),
            "w_norm": str(W_MINIMAL.norm()),
            "index": "phi/(4*sqrt(5))",
            "index_norm": "1/80",
            "target_index_norm": "1/5",
            "representative_coset": _coset_json(MINIMAL_COSET),
            "positive_unit_component_orbit": [
                _coset_json(coset) for coset in minimal_orbit
            ],
            "component_orbit_count_modulo_units_and_sign": 1,
            "components_in_full_orbit": 3,
            "every_component_is_self_opposite": True,
            "represented_vector": (
                "lambda=(phi/2)u+phi^2*r+(phi^2/2)e, "
                "where (r,r)_E8=2"
            ),
            "represented_Hhat_value": str(explicit_vector_h_value),
        },
        "real_slice_classification": {
            "facet_sheets": 0,
            "lower_face_contact_sheets": 0,
            "interior_cut_sheets": 0,
            "disposition": "entire nonempty Heegner divisor is disjoint",
            "proof": (
                "every lambda in mu+L has auxiliary coefficient y/4+a; "
                "it cannot vanish when y is nonzero modulo 4"
            ),
        },
        "scaled_fourier_matrix_bounds": {
            "matrix_size": 2,
            "indices": ["target_(phi/sqrt(5),0)", "minimal_shell_orbit"],
            "similarity_scaling": "s_i=N(m_i)^2",
            "nonzero_component_center_multiplicity_upper": (
                NONZERO_COMPONENT_CENTER_MULTIPLICITY
            ),
            "target_diagonal_error_upper": "3/4",
            "target_diagonal_certified_arb_upper": _arb_upper_string(
                target_diagonal_error
            ),
            "minimal_diagonal_error_exact_sha256": _fraction_sha256(
                minimal_diagonal_error
            ),
            "minimal_diagonal_error_decimal": _decimal(
                minimal_diagonal_error
            ),
            "minimal_diagonal_error_upper": "1/500",
            "cross_error_exact_sha256": _fraction_sha256(
                target_minimal_cross
            ),
            "cross_error_decimal": _decimal(target_minimal_cross),
            "each_cross_error_upper": "1/6",
            "target_row_dominance_margin_lower": str(target_row_margin),
            "minimal_row_dominance_margin_lower": str(minimal_row_margin),
        },
        "unit_bessel_sums": {
            name: {
                "finite_interval": str(finite),
                "tail_interval": str(tail),
                "total_upper": _arb_upper_string(total),
                "certified_rational_threshold": str(
                    UNIT_SUM_THRESHOLDS[name]
                ),
            }
            for name, (finite, tail, total) in bessel_sums.items()
        },
        "nonunit_bounds": {
            name: {
                "finite_sha256": _fraction_sha256(finite),
                "tail": str(tail),
                "total_sha256": _fraction_sha256(total),
                "total_decimal": _decimal(total),
            }
            for name, (finite, tail, total) in nonunit_data.items()
        },
        "bounded_conclusion": (
            "modulo positive-unit and sign relations, no weakly holomorphic "
            "principal part supported on the target and the complete "
            "globally minimal norm-1/80 shell can have a nonzero target "
            "coefficient"
        ),
        "nonclaims": [
            "no support bound is proved beyond the globally minimal shell",
            "larger-index or multi-shell compensator kernels are not excluded",
            "finite-index level refinements are not decided",
            "the independently certified compact Coxeter geometry is unchanged",
        ],
        "dependencies": {
            "frozen_target_certifier_sha256": _file_sha256(BASE_CERTIFIER),
            "frozen_target_certificate_sha256": _file_sha256(BASE_CERTIFICATE),
        },
    }


if __name__ == "__main__":
    print(json.dumps(certify(), indent=2, sort_keys=True))
