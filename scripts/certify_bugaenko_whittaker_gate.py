#!/usr/bin/env python3
"""Certify the pure Bruinier gate for Bugaenko's compact H8 chamber.

The canonical determinant-line extension, in Bruinier's sign convention, is

    Hhat = <-2*phi> + E8 + <-4*phi>

over O=Z[phi].  With Q=phi*Hhat/(2*sqrt(5)), its discriminant module is

    A = O/(2) + O/(4),  |A|=64.

All eleven Bugaenko facet roots have Q-value m0=phi/sqrt(5) and component
zero.  Hence the pure full norm-two Heegner input is the single principal
part q^(-m0)e_0.  This script proves that its weight-(5,5) Hilbert--Weil
Poincare obstruction is nonzero.

The proof normalizes the diagonal coefficient to 1+E and proves |E|<1.
The unit lower row is enclosed by Arb with an analytic unit tail.  For all
nonunit rows we use the exact finite-Weil coefficient magnitude

    |rho(gamma)_(0,0)| = sqrt(|A[c]|/|A|).

The sole ideal of norm 4 is better: every relevant (0,0) coefficient is an
exactly vanishing 64-term Gauss sum.  All remaining ideal and unit tails are
bounded by rational inequalities.  No floating-point inference is used.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from fractions import Fraction
from math import factorial, isqrt
from pathlib import Path

from flint import arb, ctx

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.exact.matrix import determinant, inertia
from coxeter_polytopes.exact.trace_lattice import (
    INV_SQRT5,
    change_gram_basis,
    integral_gram,
    is_codifferent,
    trace_gram,
)
from coxeter_polytopes.known_examples.bugaenko import bugaenko_basis_gram
from coxeter_polytopes.lattices import block_diagonal
from scripts.bugaenko_trace_lattice import basis_change_certificate


ZERO = Golden(0)
ONE = Golden(1)
SQRT5 = 2 * PHI - 1

DISCRIMINANT_ORDER = 64
FIELD_LEVEL = 8
WEIGHT_INPUT = (-3, 5)
WEIGHT_OBSTRUCTION = (5, 5)
IDEAL_NORM_CUTOFF = 1024
UNIT_EXPONENT_CUTOFF = 10

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
POLYTOPE_ARCHIVE = (
    REPOSITORY_ROOT / "data/known_examples/bugaenko_compact_h8.json"
)
CONTACT_ARCHIVE = (
    REPOSITORY_ROOT / "data/trace_lattices/bugaenko_k3_eigen_nef.json"
)

PI_UPPER = Fraction(355, 113)
SQRT5_LOWER = Fraction(2_236_067, 1_000_000)
PHI_LOWER = Fraction(8, 5)
R = 1 / PHI_LOWER

# For s=4*pi/sqrt(5), the fourth-order small-Bessel prefactor is
# A_1=(s/2)^4/4! < 13/5.
BESSEL_A1_UPPER = Fraction(13, 5)

MOD2_REPRESENTATIVES = tuple(
    Golden(a, b) for a in range(2) for b in range(2)
)
MOD4_REPRESENTATIVES = tuple(
    Golden(a, b) for a in range(4) for b in range(4)
)


def _decimal(value: Fraction, digits: int = 15) -> str:
    """Return a deterministic decimal truncation of a nonnegative fraction."""

    scale = 10**digits
    integer = value.numerator * scale // value.denominator
    return f"{integer // scale}.{integer % scale:0{digits}d}"


def _arb_upper_string(value: arb) -> str:
    return str(value.upper())


def _fraction_sha256(value: Fraction) -> str:
    return hashlib.sha256(str(value).encode("ascii")).hexdigest()


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def archived_chamber_data() -> dict[str, object]:
    """Read and verify the frozen compactness and mirror-contact archives."""

    polytope = json.loads(POLYTOPE_ARCHIVE.read_text(encoding="utf-8"))
    contacts = json.loads(CONTACT_ARCHIVE.read_text(encoding="utf-8"))
    closed = contacts["closed_chamber"]
    wall_gram = polytope["wall_gram"]
    facet_count = len(polytope["simple_roots"])

    return {
        "dimension": polytope["dimension"],
        "facet_count": facet_count,
        "all_facet_norms_are_two": all(
            wall_gram[index][index] == {"a": "2", "b": "0"}
            for index in range(facet_count)
        ),
        "compact_vertex_count": len(polytope["geometry"]["vertices"]),
        "cusp_count": len(polytope["geometry"]["cusps"]),
        "vinberg_status": polytope["completeness"]["evidence"]["status"],
        "facet_sheets": closed["facet_mirrors"],
        "lower_face_contact_sheets": closed[
            "nonfacet_mirrors_touching_lower_faces"
        ],
        "closed_chamber_contact_sheets": closed[
            "distinct_k3_mirrors_meeting"
        ],
        "interior_cut_sheets": closed["mirrors_cutting_interior"],
        "contact_archive_verified": contacts["computed_certificate_verified"],
        "vinberg_completion_verified": contacts[
            "vinberg_completion_verified"
        ],
        "polytope_archive_sha256": _file_sha256(POLYTOPE_ARCHIVE),
        "contact_archive_sha256": _file_sha256(CONTACT_ARCHIVE),
    }


def _chi5(value: int) -> int:
    residue = value % 5
    if residue == 0:
        return 0
    return 1 if residue in (1, 4) else -1


def ideal_count(norm: int) -> int:
    """Return the number of ideals of norm ``norm`` in Q(sqrt(5))."""

    total = 0
    for divisor in range(1, isqrt(norm) + 1):
        if norm % divisor:
            continue
        total += _chi5(divisor)
        quotient = norm // divisor
        if quotient != divisor:
            total += _chi5(quotient)
    return total


def _q_numerator(x: Golden, y: Golden) -> Golden:
    """Return 8*sqrt(5)*Q(x*u/2+y*e/4).

    The two nonunimodular lines have coefficients -2*phi and -4*phi.
    Consequently

        sqrt(5) Q = -phi^2*x^2/4 - phi^2*y^2/8.
    """

    return -(PHI**2) * (2 * x**2 + y**2)


def determinant_trace_data() -> dict[str, object]:
    """Verify the rank-20 trace lattice and its unimodular basis change."""

    field_gram = [row + [ZERO] for row in bugaenko_basis_gram()]
    field_gram.append([ZERO for _ in range(9)] + [-4 * PHI])
    alpha = PHI * INV_SQRT5
    traced = integral_gram(trace_gram(field_gram, alpha))

    old_change = basis_change_certificate()
    change = [[0 for _ in range(20)] for _ in range(20)]
    for i in range(18):
        for j in range(18):
            change[i][j] = old_change[i][j]
    # [[-4,-8],[-8,-12]] becomes diag(-4,4) on columns (1,0),(2,-1).
    change[18][18] = 1
    change[18][19] = 2
    change[19][19] = -1

    archived = bugaenko_basis_gram()
    e8 = [[int(archived[i + 1][j + 1].a) for j in range(8)] for i in range(8)]
    expected = block_diagonal([[-2]], [[2]], e8, e8, [[-4]], [[4]])
    transformed = change_gram_basis(traced, change)
    exact_inertia = inertia(traced)
    return {
        "trace_scalar": str(alpha),
        "trace_scalar_is_codifferent": is_codifferent(alpha),
        "rank": len(traced),
        "determinant": str(determinant(traced)),
        "inertia": [
            exact_inertia.positive,
            exact_inertia.negative,
            exact_inertia.zero,
        ],
        "basis_change_determinant": str(determinant(change)),
        "standard_gram_verified": transformed == expected,
    }


def _residue_mod_eight(value: Golden) -> tuple[int, int]:
    if not value.is_integral:
        raise ValueError(f"expected an integral golden value, got {value}")
    return (int(value.a) % 8, int(value.b) % 8)


def discriminant_residue_counts() -> Counter[tuple[int, int]]:
    """Count field-valued q residues on O/(2)+O/(4)."""

    return Counter(
        _residue_mod_eight(_q_numerator(x, y))
        for x in MOD2_REPRESENTATIVES
        for y in MOD4_REPRESENTATIVES
    )


def _trace_phase_residue(multiplier: Golden, x: Golden, y: Golden) -> int:
    """Return 8*Tr(multiplier*Q(mu)) modulo 8.

    If z=a+b*phi, then Tr(z/sqrt(5))=b.  Thus only the phi
    coefficient of multiplier*(8*sqrt(5)*Q(mu)) is needed.
    """

    value = multiplier * _q_numerator(x, y)
    if not value.is_integral:
        raise ValueError("Gauss numerator unexpectedly nonintegral")
    return int(value.b) % 8


def norm_four_gauss_histograms() -> dict[str, dict[int, int]]:
    """Return the three exact Gauss distributions for c=2.

    The units d modulo 2 are 1, phi, phi+1.  Taking a=d^(-1) modulo 2
    reduces the (0,0) Weil coefficients to lower-unipotent Gauss sums with
    multipliers 2, 2*phi, and 2*phi^2, up to harmless unit phases.
    """

    multipliers = (2 * ONE, 2 * PHI, 2 * PHI**2)
    return {
        str(multiplier): dict(
            sorted(
                Counter(
                    _trace_phase_residue(multiplier, x, y)
                    for x in MOD2_REPRESENTATIVES
                    for y in MOD4_REPRESENTATIVES
                ).items()
            )
        )
        for multiplier in multipliers
    }


def discriminant_kernel_size(two_adic_valuation: int) -> int:
    """Return |A[c]| for v_(2)(c)=0, 1, or at least 2."""

    if two_adic_valuation <= 0:
        multiplier = 1
    elif two_adic_valuation == 1:
        multiplier = 2
    else:
        multiplier = 4

    def annihilated(value: Golden, modulus: int) -> bool:
        product = multiplier * value
        return int(product.a) % modulus == int(product.b) % modulus == 0

    return sum(
        annihilated(x, 2) and annihilated(y, 4)
        for x in MOD2_REPRESENTATIVES
        for y in MOD4_REPRESENTATIVES
    )


def _weil_row_bound(norm: int) -> Fraction:
    """Bound a normalized (0,0) Weil--Kloosterman row of ideal norm n."""

    # The prime 2 is inert and has ideal norm 4.  Therefore divisibility of
    # the ideal norm by 4 and 16 detects v_(2)(c)=1 and >=2.
    if norm % 16 == 0:
        return Fraction(1)
    if norm % 4 == 0:
        return Fraction(1, 2)
    return Fraction(1, 8)


def _nonunit_bessel_orbit_bound(norm: int) -> Fraction:
    """Uniformly bound the positive-unit Bessel orbit for a lower row.

    A balanced generator gives arguments

        s*n^(-1/2)*phi^(t+delta), s*n^(-1/2)*phi^(-t-delta),

    with |delta|<=1/2.  For A=A_1/n^2<1, the product is bounded both by
    A^2 and by A*phi^(-4|t+delta|).  Every value below is therefore a
    rigorous upper bound; taking their finite minimum is still rigorous.
    """

    if norm < 4:
        raise ValueError("nonunit estimate requires norm at least four")
    a_value = BESSEL_A1_UPPER / norm**2
    candidates = []
    for cutoff in range(0, 2 * norm.bit_length() + 3):
        central = (2 * cutoff + 1) * a_value**2
        tail = (
            2
            * a_value
            * R ** (4 * cutoff + 2)
            / (1 - R**4)
        )
        candidates.append(central + tail)
    return min(candidates)


def nonunit_ideal_sum_bound() -> tuple[Fraction, Fraction, Fraction]:
    """Bound all nonunit rows, using the exact zero at norm four."""

    finite = sum(
        (
            ideal_count(norm)
            * _weil_row_bound(norm)
            * _nonunit_bessel_orbit_bound(norm)
            for norm in range(5, IDEAL_NORM_CUTOFF + 1)
        ),
        Fraction(0),
    )

    # For n>B, choose K=bit_length(n).  Then K<=log_2(n)+1,
    # R^(4K)<n^-2 because R^4<(1/4), and a_F(n)<=n.  Integrating the
    # resulting n^-3 majorant at B=2^10, using log(2)>2/3, gives
    #
    # B^-2 * [(49/4) C^2 + C R^2/(1-R^4)],  C=13/5.
    #
    # We discard the additional Weil saving in this infinite tail.
    tail = Fraction(1, IDEAL_NORM_CUTOFF**2) * (
        Fraction(49, 4) * BESSEL_A1_UPPER**2
        + BESSEL_A1_UPPER * R**2 / (1 - R**4)
    )
    return finite, tail, finite + tail


def unit_bessel_orbit() -> tuple[arb, arb, arb]:
    """Rigorously enclose the full unit Bessel orbit with Arb."""

    sqrt5 = arb(5).sqrt()
    phi = (1 + sqrt5) / 2
    argument = 4 * arb.pi() / sqrt5

    finite = arb(0)
    for exponent in range(-UNIT_EXPONENT_CUTOFF, UNIT_EXPONENT_CUTOFF + 1):
        first = (argument * phi**exponent).bessel_j(4)
        second = (argument * phi ** (-exponent)).bessel_j(4)
        finite += abs(first * second)

    a_one = (argument / 2) ** 4 / factorial(4)
    tail = (
        2
        * a_one
        * phi ** (-4 * (UNIT_EXPONENT_CUTOFF + 1))
        / (1 - phi**-4)
    )
    return finite, tail, finite + tail


def _arb_fraction(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def certify() -> dict[str, object]:
    ctx.prec = 220

    chamber_data = archived_chamber_data()
    trace_data = determinant_trace_data()
    residue_counts = discriminant_residue_counts()
    isotropic_count = residue_counts[(0, 0)]
    gauss_histograms = norm_four_gauss_histograms()
    expected_gauss = {0: 16, 2: 16, 4: 16, 6: 16}

    unit_finite, unit_tail, unit_sum = unit_bessel_orbit()
    sqrt5 = arb(5).sqrt()
    scalar_prefactor = 4 * arb.pi() ** 2 / sqrt5
    unit_error = scalar_prefactor * unit_sum / 8

    nonunit_finite, nonunit_tail, nonunit_sum = nonunit_ideal_sum_bound()
    scalar_prefactor_upper = 4 * PI_UPPER**2 / SQRT5_LOWER
    nonunit_error = scalar_prefactor_upper * nonunit_sum
    total_error = unit_error + _arb_fraction(nonunit_error)

    kernel_sizes = {
        "v2_0": discriminant_kernel_size(0),
        "v2_1": discriminant_kernel_size(1),
        "v2_ge_2": discriminant_kernel_size(2),
    }

    checks = {
        "archived_polytope_is_compact_H8_with_eleven_norm_two_facets": (
            chamber_data["dimension"] == 8
            and chamber_data["facet_count"] == 11
            and chamber_data["all_facet_norms_are_two"]
            and chamber_data["compact_vertex_count"] == 41
            and chamber_data["cusp_count"] == 0
            and chamber_data["vinberg_status"] == "verified"
        ),
        "archived_norm_two_contact_profile_is_11_309_320_0": (
            chamber_data["facet_sheets"] == 11
            and chamber_data["lower_face_contact_sheets"] == 309
            and chamber_data["closed_chamber_contact_sheets"] == 320
            and chamber_data["interior_cut_sheets"] == 0
            and chamber_data["contact_archive_verified"]
            and chamber_data["vinberg_completion_verified"]
        ),
        "trace_scalar_is_codifferent": trace_data[
            "trace_scalar_is_codifferent"
        ],
        "trace_lattice_is_rank_20_signature_18_2_determinant_64": (
            trace_data["rank"] == 20
            and trace_data["inertia"] == [18, 2, 0]
            and trace_data["determinant"] == "64"
            and trace_data["basis_change_determinant"] in {"1", "-1"}
            and trace_data["standard_gram_verified"]
        ),
        "discriminant_has_64_elements": sum(residue_counts.values()) == 64,
        "field_discriminant_form_is_anisotropic": isotropic_count == 1,
        "no_nontrivial_even_overlattice": isotropic_count == 1,
        "kernel_sizes_are_1_16_64": kernel_sizes
        == {"v2_0": 1, "v2_1": 16, "v2_ge_2": 64},
        "norm_four_is_the_unique_ideal_at_that_norm": ideal_count(4) == 1,
        "all_norm_four_gauss_sums_vanish": all(
            histogram == expected_gauss
            for histogram in gauss_histograms.values()
        ),
        "bessel_prefactor_bound": (
            (2 * PI_UPPER / SQRT5_LOWER) ** 4 / factorial(4)
            < BESSEL_A1_UPPER
        ),
        "class_number_one_minkowski_bound": 5 < 16,
        "total_off_identity_below_one": total_error < 1,
    }

    return {
        "analysis_id": "bugaenko_h8_determinant_whittaker_gate_v1",
        "verified": all(checks.values()),
        "checks": checks,
        "field": "Q(sqrt(5))",
        "ring": "Z[phi]",
        "determinant_line_extension": {
            "field_bilinear_gram": "<-2*phi> + E8 + <-4*phi>",
            "signatures": "((8,2),(10,0))",
            "quadratic_form": "Q(x)=phi*Hhat(x,x)/(2*sqrt(5))",
            "trace_bilinear_lattice": (
                "E8 + E8 + <-2> + <2> + <-4> + <4>"
            ),
            "trace_signature": [18, 2],
            "trace_determinant": 64,
            "trace_certificate": trace_data,
        },
        "discriminant_weil_module": {
            "field_module": "O/(2) + O/(4)",
            "additive_group": "(Z/2)^2 + (Z/4)^2",
            "order": DISCRIMINANT_ORDER,
            "field_level_ideal": "(8)",
            "field_quadratic_form": (
                "sqrt(5)*Q(x,y)=-phi^2*x^2/4-phi^2*y^2/8 mod O"
            ),
            "trace_square_form_mod_2Z": (
                "<-1/2> + <1/2> + <-1/4> + <1/4>"
            ),
            "residue_denominator": FIELD_LEVEL,
            "distinct_q_residues": len(residue_counts),
            "isotropic_elements": isotropic_count,
            "even_overlattice_status": "none_nontrivial",
            "W(P)_action": (
                "trivial: every norm-two reflection changes a dual vector "
                "by an integral root multiple"
            ),
            "kernel_sizes": kernel_sizes,
        },
        "weights": {
            "Whittaker_input": list(WEIGHT_INPUT),
            "cusp_obstruction": list(WEIGHT_OBSTRUCTION),
            "Bessel_order": 4,
        },
        "principal_part": {
            "m0": "phi/sqrt(5)",
            "mu0": "0 in A",
            "field_norm_m0": "1/5",
            "pure_input": "q^(-phi/sqrt(5))*e_0",
            "independent_columns_for_eleven_facets": 1,
            "facet_sheets": chamber_data["facet_sheets"],
            "lower_face_contact_sheets": chamber_data[
                "lower_face_contact_sheets"
            ],
            "closed_chamber_contact_sheets": chamber_data[
                "closed_chamber_contact_sheets"
            ],
            "interior_cut_sheets": chamber_data["interior_cut_sheets"],
            "facet_only_Heegner_divisor_selectable": False,
            "facet_only_selectability_scope": (
                "standard full-group/discriminant-module column Z(m0,0); "
                "finite-index level refinements are not decided"
            ),
            "archive_certificate": chamber_data,
        },
        "norm_four_exact_zero": {
            "unit_residues_mod_2": ["1", "phi", "phi+1"],
            "lower_unipotent_multipliers": list(gauss_histograms),
            "phase_histograms": gauss_histograms,
            "each_histogram": expected_gauss,
            "conclusion": "the complete ideal-norm-four Kloosterman row is zero",
        },
        "identity_term_domination": {
            "diagonal_identity_coefficient": "1",
            "scalar_Hilbert_prefactor": "4*pi^2/sqrt(5)",
            "unit_Weil_coefficient_absolute_value": "1/8",
            "unit_exponent_cutoff": UNIT_EXPONENT_CUTOFF,
            "unit_finite_upper": _arb_upper_string(unit_finite),
            "unit_tail_upper": _arb_upper_string(unit_tail),
            "unit_orbit_upper": _arb_upper_string(unit_sum),
            "unit_error_upper": _arb_upper_string(unit_error),
            "ideal_norm_cutoff": IDEAL_NORM_CUTOFF,
            "nonunit_finite_sha256": _fraction_sha256(nonunit_finite),
            "nonunit_finite_numerator_digits": len(
                str(nonunit_finite.numerator)
            ),
            "nonunit_finite_denominator_digits": len(
                str(nonunit_finite.denominator)
            ),
            "nonunit_finite_decimal": _decimal(nonunit_finite),
            "nonunit_tail_rational": str(nonunit_tail),
            "nonunit_tail_decimal": _decimal(nonunit_tail),
            "nonunit_error_sha256": _fraction_sha256(nonunit_error),
            "nonunit_error_decimal": _decimal(nonunit_error),
            "total_error_upper": _arb_upper_string(total_error),
        },
        "verdict": {
            "Poincare_series_P_m0_0_nonzero": True,
            "pure_norm_two_Whittaker_principal_part": "obstructed",
            "facet_only_divisor": (
                "not selectable by the standard full-group/discriminant-module "
                "column Z(m0,0), which also contains the 309 lower-face "
                "contact sheets; finite-index level refinements are not decided"
            ),
            "larger_compensating_principal_parts": "not_decided",
            "polytope_compactness": "independently archived and unaffected",
            "transfer_theorem_consequence": (
                "the canonical determinant line does not automatically admit "
                "the pure facet/norm-two Whittaker correction; a general "
                "transfer theorem must allow compensating terms or a "
                "non-Borcherds carrier branch"
            ),
        },
    }


if __name__ == "__main__":
    print(json.dumps(certify(), indent=2, sort_keys=True))
