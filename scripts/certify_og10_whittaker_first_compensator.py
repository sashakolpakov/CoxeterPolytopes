#!/usr/bin/env python3
"""Certify the first complete OG10 off-slice compensator shell.

The OG10 Hilbert discriminant module decomposes as

    O/(1+4*phi) + O/(2) + O/(2),

where the last factor is contributed by the added sign-changing line.  A
nonzero class in that last factor forces every vector in its coset to have
nonzero added-line coordinate.  Its complete special divisor therefore
misses the embedded real H^10 slice.

This script proves two bounded statements.

1. Among *all* such discriminant classes, the smallest possible totally
   positive Whittaker index has norm 9/880.  Up to multiplication by units
   and the tautological mu -> -mu symmetry, the complete minimal shell has
   exactly two discriminant-component orbits.
2. The three Poincare series indexed by the target and these two minimal
   off-slice orbits are linearly independent.  Hence no weakly holomorphic
   principal part supported on this complete three-orbit shell can have a
   nonzero target coefficient.

Finite Bessel sums are enclosed by Arb through python-flint.  Infinite unit
and ideal tails are bounded separately by exact rational inequalities.
"""

from __future__ import annotations

import json
from fractions import Fraction

from flint import arb, ctx

from coxeter_polytopes.exact.golden import Golden, PHI
from scripts.certify_og10_whittaker_nonvanishing import (
    IDEAL_NORM_CUTOFF,
    PI_UPPER,
    R,
    SQRT176_LOWER,
    SQRT5_LOWER,
    SQRT5_UPPER,
    ideal_count,
)


ZERO = Golden(0)
ONE = Golden(1)
SQRT5 = 2 * PHI - 1

# The field block has determinant -(1+4*phi), of absolute norm 11.  The two
# O/(2) factors come from the positive <2> line in H and the auxiliary line.
SPLIT_ELEVEN_GENERATOR = 1 + 4 * PHI
AUXILIARY_LINE_COEFFICIENT = 2 * (1 - PHI)
DISCRIMINANT_ORDER = 11 * 4 * 4
INDEX_RESIDUE_DENOMINATOR = 44

M0 = PHI / SQRT5
W0 = -ONE + 3 * PHI
W_MINIMAL = 3 * W0
M_MINIMAL = W_MINIMAL / (INDEX_RESIDUE_DENOMINATOR * SQRT5)

MOD2_REPRESENTATIVES = (ZERO, ONE, PHI, PHI**2)
NONZERO_MOD2_REPRESENTATIVES = MOD2_REPRESENTATIVES[1:]

# A discriminant class is recorded as (t,x,y), with
#
#   t in O/(1+4phi) = F_11,  x in O/(2),  y in O/(2).
#
# The corresponding dual vector is t*v + (x/2)u + (phi*y/2)e.  The two
# minimal orbits below are representatives modulo units and mu -> -mu.
MINIMAL_COSET_A = (4, ZERO, ONE)
MINIMAL_COSET_B = (4, PHI, PHI)

UNIT_SUM_THRESHOLDS = {
    "target_diagonal": Fraction(127, 1000),
    "minimal_diagonal": Fraction(1, 300_000),
    "target_minimal_cross": Fraction(1, 500),
}

# If A_ij denotes the fifth-order small-Bessel prefactor for the (i,j)
# coefficient, these are exact convenient majorants.  The target value is
# <3/2.  With r=N(m_min)/N(m_0)=9/176, the other two follow from
# A_11/A_00=r^(5/2)<r^2 and A_01/A_00=r^(5/4)<r.
NONUNIT_BESSEL_PREFACTOR_UPPER = {
    "target_diagonal": Fraction(3, 2),
    "minimal_diagonal": Fraction(1, 200),
    "target_minimal_cross": Fraction(2, 25),
}

# Bruinier's translation subgroup contains the metaplectic kernel but not the
# lift of -I.  For a nonzero discriminant component, passing from scalar
# center-quotiented rows to the exact vector Poincare series can pair the
# seeds mu and -mu.  Bounding both representatives separately costs at most
# this factor.  It makes the certificate independent of the center/sign
# normalization convention; for the zero target component the factor cancels
# against its doubled identity term and is not needed.
NONZERO_COMPONENT_CENTER_MULTIPLICITY = 2


def _decimal(value: Fraction, digits: int = 15) -> str:
    scale = 10**digits
    integer = value.numerator * scale // value.denominator
    return f"{integer // scale}.{integer % scale:0{digits}d}"


def _arb_upper_string(value: arb) -> str:
    return str(value.upper())


def _residue(value: Golden, modulus: int) -> tuple[int, int]:
    """Return the coefficient pair of an integral field element modulo n."""

    if not value.is_integral:
        raise ValueError(f"residue requested for nonintegral element {value}")
    return (int(value.a) % modulus, int(value.b) % modulus)


def _representative(pair: tuple[int, int]) -> Golden:
    return Golden(pair[0], pair[1])


def _reduce(value: Golden, modulus: int) -> Golden:
    return _representative(_residue(value, modulus))


def _q_residue_numerator(t: int, x: Golden, y: Golden) -> Golden:
    """Return 44*sqrt(5)*Q(mu) modulo 44 O for mu=(t,x,y).

    If v=A^{-1}e_1 for A=[[-2phi,1],[1,2]], then

        sqrt(5) Q(v)       = (phi-4)/11,
        sqrt(5) Q(xu/2)    = phi*x^2/4,
        sqrt(5) Q(phi*y*e/2) = -phi^2*y^2/4.

    Multiplication by 44 gives the integral expression below.
    """

    return (
        4 * t * t * (PHI - 4)
        + 11 * PHI * x**2
        - 11 * PHI**2 * y**2
    )


def _all_off_slice_residues() -> dict[
    tuple[int, int], list[tuple[int, Golden, Golden]]
]:
    residues: dict[tuple[int, int], list[tuple[int, Golden, Golden]]] = {}
    for t in range(11):
        for x in MOD2_REPRESENTATIVES:
            for y in NONZERO_MOD2_REPRESENTATIVES:
                key = _residue(_q_residue_numerator(t, x, y), 44)
                residues.setdefault(key, []).append((t, x, y))
    return residues


def _unit_action(
    coset: tuple[int, Golden, Golden],
) -> tuple[int, Golden, Golden]:
    """Multiply a discriminant class by phi.

    Modulo p=(1+4phi), phi is 8 in F_11.  On both O/(2) factors it acts by
    ordinary multiplication.  The associated index is multiplied by phi^2.
    """

    t, x, y = coset
    return (8 * t % 11, _reduce(PHI * x, 2), _reduce(PHI * y, 2))


def _negate_coset(
    coset: tuple[int, Golden, Golden],
) -> tuple[int, Golden, Golden]:
    t, x, y = coset
    return (-t % 11, _reduce(-x, 2), _reduce(-y, 2))


def _unit_coset_orbit(
    seed: tuple[int, Golden, Golden],
) -> set[tuple[int, Golden, Golden]]:
    orbit: set[tuple[int, Golden, Golden]] = set()
    current = seed
    while current not in orbit:
        orbit.add(current)
        current = _unit_action(current)
    return orbit


def _positive_unit_residue_orbit(value: Golden) -> set[tuple[int, int]]:
    """Residues modulo 44 of the sign-(+,-) norm-one unit orbit."""

    if value.sign(0) > 0 and value.sign(1) < 0:
        oriented = value
    elif value.sign(0) < 0 and value.sign(1) > 0:
        oriented = -value
    else:
        raise ValueError("expected an element of negative field norm")
    return {_residue(oriented * PHI ** (2 * exponent), 44) for exponent in range(15)}


def _minimal_shell_data() -> dict[str, object]:
    allowed = _all_off_slice_residues()
    allowed_keys = set(allowed)
    norm_residues = {
        (a * a + a * b - b * b) % 44 for a, b in allowed_keys
    }
    smaller_absolute_norms = tuple(
        value for value in range(1, 99) if (-value) % 44 in norm_residues
    )

    # Class number one follows from the real-quadratic Minkowski bound
    # sqrt(5)/2<2.  Consequently the possible element orbits at norms 11,
    # 44, and 55 are generated by the two factors of 11, with the forced
    # factors 2 or phi*sqrt(5).  Norms 33, 77, and 88 are impossible because
    # an inert prime occurs to odd valuation.
    possible_below = {
        11: (W0, W0.conjugate()),
        44: (2 * W0, 2 * W0.conjugate()),
        55: (
            PHI * SQRT5 * W0,
            PHI * SQRT5 * W0.conjugate(),
        ),
    }
    possible_intersections = {
        norm: [
            sorted(_positive_unit_residue_orbit(base) & allowed_keys)
            for base in bases
        ]
        for norm, bases in possible_below.items()
    }

    minimal_prime_orbits = (
        3 * W0,
        3 * W0.conjugate(),
    )
    minimal_intersections = [
        _positive_unit_residue_orbit(base) & allowed_keys
        for base in minimal_prime_orbits
    ]
    minimal_allowed_cosets = {
        coset
        for key in minimal_intersections[0]
        for coset in allowed[key]
    }
    orbit_a = _unit_coset_orbit(MINIMAL_COSET_A)
    orbit_b = _unit_coset_orbit(MINIMAL_COSET_B)

    return {
        "allowed": allowed,
        "allowed_keys": allowed_keys,
        "norm_residues": norm_residues,
        "smaller_absolute_norms": smaller_absolute_norms,
        "possible_intersections": possible_intersections,
        "minimal_intersections": minimal_intersections,
        "minimal_allowed_cosets": minimal_allowed_cosets,
        "orbit_a": orbit_a,
        "orbit_b": orbit_b,
    }


def _unit_bessel_sum(
    left_embeddings: tuple[arb, arb],
    right_embeddings: tuple[arb, arb],
    cutoff: int = 8,
) -> tuple[arb, arb, arb]:
    """Rigorous unit-orbit Bessel sum, split into finite part and tail."""

    pi = arb.pi()
    phi = (1 + arb(5).sqrt()) / 2
    base_left = 4 * pi * (left_embeddings[0] * right_embeddings[0]).sqrt()
    base_right = 4 * pi * (left_embeddings[1] * right_embeddings[1]).sqrt()

    finite = arb(0)
    for exponent in range(-cutoff, cutoff + 1):
        first = (base_left * phi**exponent).bessel_j(5)
        second = (base_right * phi ** (-exponent)).bessel_j(5)
        finite += abs(first * second)

    # For t>cutoff the second argument is small; for t<-cutoff the first is
    # small.  Use |J_5(x)| <= (x/2)^5/5! for the small factor and <=1 for
    # the other one.
    left_small_constant = (base_left / 2) ** 5 / 120
    right_small_constant = (base_right / 2) ** 5 / 120
    tail = (
        (left_small_constant + right_small_constant)
        * phi ** (-5 * (cutoff + 1))
        / (1 - phi**-5)
    )
    return finite, tail, finite + tail


def _generic_nonunit_orbit_bound(norm: int, a_upper: Fraction) -> Fraction:
    """Bessel unit-orbit bound for every lower-left ideal of norm >=4."""

    cutoff = norm.bit_length()
    central = Fraction(2 * cutoff + 1, norm**5) * a_upper**2
    tail = (
        Fraction(2, norm**2)
        * a_upper
        * R ** (5 * cutoff + 2)
        / (1 - R**5)
    )
    return central + tail


def _generic_nonunit_ideal_sum(a_upper: Fraction) -> Fraction:
    finite = sum(
        (
            ideal_count(norm) * _generic_nonunit_orbit_bound(norm, a_upper)
            for norm in range(4, IDEAL_NORM_CUTOFF + 1)
        ),
        Fraction(0),
    )

    # For n>B, bit_length(n)<=log_2(n)+1, (R^5)^K<n^-3,
    # and a_F(n)<=n.  Integrating the resulting decreasing n^-4 majorant
    # at B=2^10 and using log(2)>2/3 gives this rational tail.
    constant = 3 * a_upper**2 + 2 * a_upper * R**2 / (1 - R**5)
    tail = Fraction(1, IDEAL_NORM_CUTOFF**3) * (
        7 * a_upper**2 + constant / 3
    )
    return finite + tail


def _coset_json(coset: tuple[int, Golden, Golden]) -> dict[str, object]:
    t, x, y = coset
    return {
        "field_block_t_mod_11": t,
        "positive_line_x_mod_2": str(x),
        "auxiliary_line_y_mod_2": str(y),
    }


def certify() -> dict[str, object]:
    shell = _minimal_shell_data()

    ctx.prec = 180
    sqrt5 = arb(5).sqrt()
    phi = (1 + sqrt5) / 2

    target_embeddings = (phi / sqrt5, 1 / (phi * sqrt5))
    minimal_embeddings = (
        (-3 + 9 * phi) / (44 * sqrt5),
        (9 * phi - 6) / (44 * sqrt5),
    )

    bessel_sums: dict[str, tuple[arb, arb, arb]] = {
        "target_diagonal": _unit_bessel_sum(
            target_embeddings, target_embeddings
        ),
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

    scalar_prefactor_upper = 4 * PI_UPPER**2 / SQRT5_LOWER
    unit_weil_factor_upper = scalar_prefactor_upper / SQRT176_LOWER

    nonunit_bounds = {
        name: _generic_nonunit_ideal_sum(
            NONUNIT_BESSEL_PREFACTOR_UPPER[name]
        )
        for name in NONUNIT_BESSEL_PREFACTOR_UPPER
    }
    target_diagonal_error = (
        unit_weil_factor_upper * UNIT_SUM_THRESHOLDS["target_diagonal"]
        + scalar_prefactor_upper * nonunit_bounds["target_diagonal"]
    )
    minimal_diagonal_error = NONZERO_COMPONENT_CENTER_MULTIPLICITY * (
        unit_weil_factor_upper * UNIT_SUM_THRESHOLDS["minimal_diagonal"]
        + scalar_prefactor_upper * nonunit_bounds["minimal_diagonal"]
    )
    target_minimal_cross = NONZERO_COMPONENT_CENTER_MULTIPLICITY * (
        unit_weil_factor_upper
        * UNIT_SUM_THRESHOLDS["target_minimal_cross"]
        + scalar_prefactor_upper
        * nonunit_bounds["target_minimal_cross"]
    )

    # If C_ij is the normalized coefficient matrix and
    # s_i=N(m_i)^(5/2), then A=diag(s_i)^(-1) C diag(s_i) removes the
    # outer Fourier ratio.  The following are strict row-diagonal-dominance
    # margins for A.  The two minimal-minimal off-diagonal entries have the
    # same analytic bound as a minimal diagonal error, but no identity term.
    target_row_margin = 1 - target_diagonal_error - 2 * target_minimal_cross
    minimal_row_margin = (
        1 - 2 * minimal_diagonal_error - target_minimal_cross
    )

    q_block = -PHI / SPLIT_ELEVEN_GENERATOR
    q_positive_line = PHI / 4
    q_auxiliary_line = -(PHI**2) / 4
    ratio = M_MINIMAL.norm() / M0.norm()

    unit_iterate = MINIMAL_COSET_A
    for _ in range(15):
        unit_iterate = _unit_action(unit_iterate)

    checks = {
        "field_block_determinant_norm_minus_11": (
            -SPLIT_ELEVEN_GENERATOR
        ).norm()
        == -11,
        "field_block_dual_q_residue": q_block == (PHI - 4) / 11,
        "positive_line_q_residue": q_positive_line == PHI / 4,
        "auxiliary_line_q_residue": q_auxiliary_line == -(PHI**2) / 4,
        "auxiliary_discriminant_factor_is_O_mod_2": (
            AUXILIARY_LINE_COEFFICIENT * PHI / 2
        ).norm()
        in (1, -1),
        "discriminant_module_order_176": DISCRIMINANT_ORDER == 176,
        "all_132_off_slice_classes_enumerated": sum(
            len(cosets) for cosets in shell["allowed"].values()
        )
        == 11 * 4 * 3,
        "off_slice_index_residue_count_54": len(shell["allowed_keys"]) == 54,
        "allowed_norm_residues_mod_44": shell["norm_residues"]
        == {0, 11, 33},
        "only_six_smaller_norm_candidates": shell["smaller_absolute_norms"]
        == (11, 33, 44, 55, 77, 88),
        "minkowski_proves_class_number_one": SQRT5_UPPER < 4,
        "inert_primes_exclude_33_77_88": ideal_count(2) == 0
        and ideal_count(3) == 0
        and ideal_count(7) == 0,
        "five_is_ramified": ideal_count(5) == 1 and SQRT5.norm() == -5,
        "eleven_splits_into_two_prime_ideals": ideal_count(11) == 2
        and W0.norm() == W0.conjugate().norm() == -11,
        "positive_unit_order_mod_44_is_15": _residue(PHI**30, 44)
        == (1, 0)
        and all(
            _residue(PHI ** (2 * exponent), 44) != (1, 0)
            for exponent in range(1, 15)
        ),
        "all_possible_11_44_55_orbits_miss_allowed_residues": all(
            not intersection
            for intersections in shell["possible_intersections"].values()
            for intersection in intersections
        ),
        "minimal_w_has_norm_minus_99": W_MINIMAL.norm() == -99,
        "both_norm_99_prime_orbits_classified": (3 * W0).norm()
        == (3 * W0.conjugate()).norm()
        == -99,
        "minimal_index_totally_positive": M_MINIMAL.sign(0) > 0
        and M_MINIMAL.sign(1) > 0,
        "minimal_index_norm_9_over_880": M_MINIMAL.norm()
        == Fraction(9, 880),
        "minimal_coset_a_has_candidate_residue": _residue(
            _q_residue_numerator(*MINIMAL_COSET_A), 44
        )
        == _residue(W_MINIMAL, 44),
        "minimal_coset_b_has_candidate_residue": _residue(
            _q_residue_numerator(*MINIMAL_COSET_B), 44
        )
        == _residue(W_MINIMAL, 44),
        "only_one_norm_99_prime_orbit_is_admissible": len(
            shell["minimal_intersections"][0]
        )
        == 15
        and not shell["minimal_intersections"][1],
        "complete_minimal_shell_is_two_unit_orbits": shell[
            "minimal_allowed_cosets"
        ]
        == shell["orbit_a"] | shell["orbit_b"]
        and not (shell["orbit_a"] & shell["orbit_b"])
        and len(shell["orbit_a"]) == len(shell["orbit_b"]) == 30,
        "unit_fifteenth_power_is_coset_negation": unit_iterate
        == _negate_coset(MINIMAL_COSET_A),
        "minimal_to_target_norm_ratio_9_over_176": ratio
        == Fraction(9, 176),
        "nonunit_prefactor_majorants": Fraction(3, 2) * ratio**2
        < NONUNIT_BESSEL_PREFACTOR_UPPER["minimal_diagonal"]
        and Fraction(3, 2) * ratio
        < NONUNIT_BESSEL_PREFACTOR_UPPER["target_minimal_cross"],
        **{f"arb_unit_sum_{name}": value for name, value in bessel_checks.items()},
        "scaled_target_row_strictly_diagonally_dominant": target_row_margin
        > 0,
        "scaled_minimal_rows_strictly_diagonally_dominant": minimal_row_margin
        > 0,
        "published_compact_matrix_bounds": target_diagonal_error
        < Fraction(68, 125)
        and minimal_diagonal_error < Fraction(1, 31_250)
        and target_minimal_cross < Fraction(1, 125)
        and target_row_margin > Fraction(441, 1000)
        and minimal_row_margin > Fraction(124, 125),
    }

    return {
        "verified": all(checks.values()),
        "theorem": (
            "the target and both globally minimal slice-invisible "
            "compensator Poincare series are linearly independent"
        ),
        "checks": checks,
        "slice_invisibility_criterion": (
            "nonzero projection y to the auxiliary O/(2) discriminant factor"
        ),
        "discriminant_module": {
            "decomposition": "O/(1+4*phi) + O/(2) + O/(2)_aux",
            "order": DISCRIMINANT_ORDER,
            "coset_coordinates": "(t,x,y), with y nonzero for off-slice support",
            "sqrt5_Q_mod_O": (
                "t^2*(phi-4)/11 + phi*x^2/4 - phi^2*y^2/4"
            ),
        },
        "minimal_shell": {
            "index": str(M_MINIMAL),
            "sqrt5_times_index": str(W_MINIMAL / 44),
            "index_norm": "9/880",
            "representative_cosets": [
                _coset_json(MINIMAL_COSET_A),
                _coset_json(MINIMAL_COSET_B),
            ],
            "unit_orbit_count_modulo_mu_sign": 2,
            "classes_in_each_full_phi_orbit": 30,
            "excluded_smaller_negative_element_norms": list(
                shell["smaller_absolute_norms"]
            ),
        },
        "unit_bessel_sums": {
            name: {
                "finite_interval": str(finite),
                "tail_interval": str(tail),
                "total_upper": _arb_upper_string(total),
                "certified_rational_threshold": str(UNIT_SUM_THRESHOLDS[name]),
            }
            for name, (finite, tail, total) in bessel_sums.items()
        },
        "scaled_fourier_matrix_bounds": {
            "matrix_size": 3,
            "similarity_scaling": "s_i=N(m_i)^(5/2)",
            "nonzero_component_center_multiplicity_upper": (
                NONZERO_COMPONENT_CENTER_MULTIPLICITY
            ),
            "target_diagonal_error_upper": str(target_diagonal_error),
            "target_diagonal_error_upper_decimal": _decimal(
                target_diagonal_error
            ),
            "minimal_diagonal_error_upper": str(minimal_diagonal_error),
            "minimal_diagonal_error_upper_decimal": _decimal(
                minimal_diagonal_error
            ),
            "target_minimal_cross_upper": str(target_minimal_cross),
            "target_minimal_cross_upper_decimal": _decimal(
                target_minimal_cross
            ),
            "minimal_minimal_cross_upper": str(minimal_diagonal_error),
            "minimal_minimal_cross_upper_decimal": _decimal(
                minimal_diagonal_error
            ),
            "target_row_dominance_margin_lower": str(target_row_margin),
            "target_row_dominance_margin_lower_decimal": _decimal(
                target_row_margin
            ),
            "minimal_row_dominance_margin_lower": str(minimal_row_margin),
            "minimal_row_dominance_margin_lower_decimal": _decimal(
                minimal_row_margin
            ),
        },
        "bounded_conclusion": (
            "modulo tautological unit/sign relations, no weakly holomorphic "
            "principal part supported on the target orbit and the complete "
            "globally minimal norm-9/880 off-slice shell can have a nonzero "
            "target coefficient"
        ),
        "nonclaims": [
            "no completeness beyond the globally minimal off-slice shell",
            "larger-index compensators and multi-shell kernels are not excluded",
            "the reflection index and chamber covolume remain undecided",
        ],
    }


if __name__ == "__main__":
    print(json.dumps(certify(), indent=2, sort_keys=True))
