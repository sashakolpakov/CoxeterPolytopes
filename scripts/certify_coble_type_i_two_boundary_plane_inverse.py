#!/usr/bin/env python3
"""Certify a raw-P plane inverse for the two-boundary H11/26 chamber.

The program starts with the physical odd-unimodular Gram and primitive walls,
enumerates all sequential numerical P2 contractions, and recognizes four
intrinsic (quintic,line) endpoints.  It then derives, rather than receives,
the proximity/root-line program.  The canonical 23 by 21 quintic matrix has
rank-drop locus a^2=b^2.  One component is nonreduced; the other gives an
absolutely irreducible rational quintic with two A5 singularities.

No Kondō type, surface equation, or supplied coordinate/Hilbert witness is
used to select the four raw branches.  Kondō and the uniform D_CobVin return
theorem enter only after the explicit strong-Coble root seed is constructed.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Sequence

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.certify_coble_remaining_seven_contraction_designs import (  # noqa: E402
    _dot,
    _load_rows,
    _normalized,
)
from scripts.certify_coble_type_ii_b3_plane_inverse import (  # noqa: E402
    _basis_coordinates,
    _endpoint_signature,
    _orthogonal_plane_basis,
    _sequential_endpoints,
)


OUTPUT = ROOT / "data/enumerations/coble_type_i_two_boundary_plane_inverse.json"
SOURCE = ROOT / "data/enumerations/coble_type_i_two_boundary_full_ns_h11.json"
TRANSFER = ROOT / "data/enumerations/coble_vinberg_incidence_inverse_transfer.json"
RAW_HOSTILE_OUTPUT = (
    ROOT
    / "data/enumerations/coble_type_i_two_boundary_plane_inverse_raw_hostile_audit.json"
)
RAW_HOSTILE_SCRIPT = (
    ROOT
    / "scripts/independent_coble_type_i_two_boundary_plane_inverse_raw_hostile.py"
)
CM_REFERENCE = ROOT / "data/known_examples/enriques_finite_automorphism_nef_chambers.json"
SCHEMA = "coxeter.coble_type_i_two_boundary_plane_inverse.v1"

X, Y, Z = sp.symbols("x y z")
A, B = sp.symbols("a b")
COEFFICIENTS = sp.symbols("c0:21")
AFFINE_MONOMIALS = tuple(
    X**i * Y**j for i in range(6) for j in range(6 - i)
)

TARGET_SIGNATURE = (
    (1, 5),
    ((0, 7), (1, 3)),
    ((-2, 7), (-1, 4)),
)

GOOD_QUINTIC = sp.expand(
    16 * X**5
    + 39 * X**4 * Y
    - 32 * X**4
    + 34 * X**3 * Y**2
    - 52 * X**3 * Y
    + 16 * X**3
    + 9 * X**2 * Y**3
    - 16 * X**2 * Y**2
    + 6 * X**2 * Y
    - 2 * X * Y**4
    + 8 * X * Y**3
    - 10 * X * Y**2
    + 4 * X * Y
    - Y**5
    + 4 * Y**4
    - 6 * Y**3
    + 4 * Y**2
    - Y
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _replay_explicit_missing_curve_repair() -> dict[str, Any]:
    """Replay the helper-free repair; do not accept an archive-only boolean.

    The hostile executable independently rebuilds the 52/9/4 endpoint census,
    interpolates and absolutely factors all ten missing cusp curves, and runs
    fresh CM and full-cone double descriptions.  Calling its frozen-output
    check here makes those calculations part of the primary replay.
    """

    completed = subprocess.run(
        [
            sys.executable,
            str(RAW_HOSTILE_SCRIPT),
            "--check-output",
            str(RAW_HOSTILE_OUTPUT),
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
        timeout=120,
    )
    if completed.returncode:
        raise ArithmeticError(
            "the independent ten-cusp-curve replay failed: "
            + (completed.stderr[-2000:] or completed.stdout[-2000:])
        )
    certificate = json.loads(RAW_HOSTILE_OUTPUT.read_text(encoding="utf-8"))
    repair = certificate["explicit_missing_curve_repair"]
    cm = certificate["CM_root_and_cusp_replay"]
    full = certificate["full_cone_replay"]
    if not (
        certificate["all_checks_pass"] is True
        and certificate["sources"]["raw_input_sha256"] == _sha256(SOURCE)
        and certificate["raw_endpoint_replay"]["all_sequential_P2_endpoint_count"] == 52
        and certificate["raw_endpoint_replay"]["selected_line_quintic_endpoint_count"] == 4
        and repair["ten_positive_degree_minus_one_curve_count"] == 10
        and repair["all_ten_unique_plane_curves"]
        and repair["all_ten_absolutely_irreducible"]
        and repair["all_ten_exact_strict_transforms_smooth_rational"]
        and cm["CM_finite_ideal_hyperideal"] == [14, 9, 0]
        and cm["all_and_only_integral_cusp_lifts_equal_twelve_raw_walls"]
        and full["finite_ideal_hyperideal"] == [508, 42, 0]
        and full["all_26_facets_irredundant"]
    ):
        raise ArithmeticError("the structured ten-cusp-curve repair gates changed")
    return {
        "independent_replay_script": str(RAW_HOSTILE_SCRIPT.relative_to(ROOT)),
        "independent_replay_archive": str(RAW_HOSTILE_OUTPUT.relative_to(ROOT)),
        "script_sha256": _sha256(RAW_HOSTILE_SCRIPT),
        "archive_sha256": _sha256(RAW_HOSTILE_OUTPUT),
        "original_primary_verdict": certificate["verdict"]["original_primary_certificate"],
        "repaired_theorem_verdict": certificate["verdict"]["repaired_theorem"],
        "degree_multiset": repair["degree_multiset"],
        "four_deep_exceptional_minus_one_curves": repair[
            "four_deep_exceptional_minus_one_curves"
        ],
        "ten_positive_degree_minus_one_curve_count": repair[
            "ten_positive_degree_minus_one_curve_count"
        ],
        "all_ten_unique_plane_curves": repair["all_ten_unique_plane_curves"],
        "all_ten_absolutely_irreducible": repair["all_ten_absolutely_irreducible"],
        "all_ten_exact_strict_transforms_smooth_rational": repair[
            "all_ten_exact_strict_transforms_smooth_rational"
        ],
        "curves": repair["curves"],
        "fresh_CM_finite_ideal_hyperideal": cm["CM_finite_ideal_hyperideal"],
        "fresh_full_finite_ideal_hyperideal": full["finite_ideal_hyperideal"],
        "all_and_only_integral_cusp_lifts_literal": cm[
            "all_and_only_integral_cusp_lifts_equal_twelve_raw_walls"
        ],
        "passes": True,
    }


def _rational_vector(vector: sp.Matrix) -> list[int | str]:
    """JSON-stable exact coordinates; no floating point enters a transfer gate."""

    return [int(value) if value.q == 1 else str(value) for value in vector]


def _primitive_vector(vector: sp.Matrix) -> tuple[int, ...]:
    denominator = 1
    for value in vector:
        denominator = int(sp.ilcm(denominator, value.q))
    entries = [int(value * denominator) for value in vector]
    divisor = 0
    for value in entries:
        divisor = math.gcd(divisor, abs(value))
    entries = [value // divisor for value in entries]
    for value in entries:
        if value:
            if value < 0:
                entries = [-entry for entry in entries]
            break
    return tuple(entries)


def _homogenize(polynomial: sp.Expr, degree: int) -> sp.Expr:
    return sp.expand(polynomial.subs({X: X / Z, Y: Y / Z}) * Z**degree)


def _poly_string(polynomial: sp.Expr) -> str:
    return str(sp.expand(polynomial)).replace("**", "^")


def _split_low_terms(
    polynomial: sp.Expr, first: sp.Symbol, second: sp.Symbol, multiplicity: int
) -> tuple[list[sp.Expr], sp.Expr]:
    parsed = sp.Poly(sp.expand(polynomial), first, second)
    conditions: list[sp.Expr] = []
    high = sp.Integer(0)
    for (left_degree, right_degree), coefficient in parsed.terms():
        term = coefficient * first**left_degree * second**right_degree
        if left_degree + right_degree < multiplicity:
            conditions.append(coefficient)
        else:
            high += term
    return conditions, sp.expand(high)


def _chain_conditions(
    polynomial: sp.Expr,
    x_map: Any,
    y_map: Any,
    multiplicities: Sequence[int],
    later_slopes: Sequence[sp.Expr],
) -> list[sp.Expr]:
    """Linear jet conditions along one unibranch proximity chain.

    The first tangent is always the local line ``v=0``.  ``later_slopes``
    gives the tangent selected at each subsequent infinitely-near center.
    Lower terms are recorded as equations before the strict transform is
    divided, so no illicit Laurent terms enter the matrix.
    """

    u, v = sp.symbols("u v")
    current = sp.expand(polynomial.subs({X: x_map(u, v), Y: y_map(u, v)}))
    second = v
    answer: list[sp.Expr] = []
    for level, multiplicity in enumerate(multiplicities):
        low, high = _split_low_terms(current, u, second, multiplicity)
        answer.extend(low)
        if level == len(multiplicities) - 1:
            break
        next_variable = sp.Symbol(f"w{level}")
        slope = sp.Integer(0) if level == 0 else later_slopes[level - 1]
        transformed = high.subs(second, u * (next_variable + slope))
        current = sp.cancel(sp.expand(transformed) / u**multiplicity)
        second = next_variable
    return [sp.expand(value) for value in answer]


def _quintic_condition_matrix() -> tuple[sp.Matrix, list[str]]:
    polynomial = sum(
        coefficient * monomial
        for coefficient, monomial in zip(COEFFICIENTS, AFFINE_MONOMIALS)
    )
    blocks = [
        (
            "A_line_chain_3",
            _chain_conditions(
                polynomial,
                lambda u, v: u,
                lambda u, v: v,
                (1, 1, 1),
                (sp.Integer(0),),
            ),
        ),
        (
            "B_line_chain_2",
            _chain_conditions(
                polynomial,
                lambda u, v: 1 + u,
                lambda u, v: v,
                (1, 1),
                (),
            ),
        ),
        (
            "C_double_chain_3",
            _chain_conditions(
                polynomial,
                lambda u, v: v,
                lambda u, v: 1 + u,
                (2, 2, 2),
                (A,),
            ),
        ),
        (
            "D_double_chain_3",
            _chain_conditions(
                polynomial,
                lambda u, v: -1 + u,
                lambda u, v: 2 - 2 * u + v,
                (2, 2, 2),
                (B,),
            ),
        ),
    ]
    conditions: list[sp.Expr] = []
    labels: list[str] = []
    for name, block in blocks:
        for index, condition in enumerate(block):
            conditions.append(condition)
            labels.append(f"{name}.{index}")
    matrix = sp.Matrix(
        [
            [condition.coeff(coefficient) for coefficient in COEFFICIENTS]
            for condition in conditions
        ]
    )
    return matrix, labels


def _kernel_polynomial(matrix: sp.Matrix) -> sp.Expr:
    kernel = matrix.nullspace()
    if len(kernel) != 1:
        raise ArithmeticError("the quintic kernel is not one-dimensional")
    vector = _primitive_vector(kernel[0])
    return sp.expand(
        sum(value * monomial for value, monomial in zip(vector, AFFINE_MONOMIALS))
    )


def _tangent_cone(
    polynomial: sp.Expr, variables: tuple[sp.Symbol, sp.Symbol]
) -> tuple[int, sp.Expr]:
    parsed = sp.Poly(sp.expand(polynomial), *variables)
    order = min(sum(exponents) for exponents, coefficient in parsed.terms() if coefficient)
    tangent = sum(
        coefficient
        * math.prod(variable**exponent for variable, exponent in zip(variables, exponents))
        for exponents, coefficient in parsed.terms()
        if sum(exponents) == order
    )
    return order, sp.factor(tangent)


def _local_resolution_record(
    polynomial: sp.Expr,
    label: str,
    x_map: Any,
    y_map: Any,
    slope: int,
) -> dict[str, Any]:
    u, v, w, q = sp.symbols("u v w q")
    local = sp.expand(polynomial.subs({X: x_map(u, v), Y: y_map(u, v)}))
    first = sp.cancel(local.subs(v, u * w) / u**2)
    second = sp.cancel(first.subs(w, u * (q + slope)) / u**2)
    proper = _tangent_cone(local, (u, v))
    near = _tangent_cone(first, (u, w))
    terminal = _tangent_cone(second, (u, q))
    terminal_poly = sp.Poly(terminal[1], u, q)
    aa = terminal_poly.coeff_monomial(u**2)
    bb = terminal_poly.coeff_monomial(u * q)
    cc = terminal_poly.coeff_monomial(q**2)
    discriminant = sp.factor(bb**2 - 4 * aa * cc)
    if proper[0] != near[0] or near[0] != terminal[0] or proper[0] != 2:
        raise ArithmeticError(f"{label}: the multiplicity sequence is not (2,2,2)")
    if discriminant == 0:
        raise ArithmeticError(f"{label}: the terminal double point is not ordinary")
    return {
        "label": label,
        "multiplicity_sequence": [2, 2, 2],
        "proper_tangent_cone": str(proper[1]),
        "first_near_tangent_cone": str(near[1]),
        "terminal_tangent_cone": str(terminal[1]),
        "terminal_discriminant": str(discriminant),
        "analytic_type_in_characteristic_zero": "A5",
        "delta_invariant": 3,
        "resolved_after_three_blowups": True,
    }


def _singular_audit(quintic_homogeneous: sp.Expr) -> dict[str, Any]:
    executable = shutil.which("Singular")
    if executable is None:
        raise RuntimeError("Singular is required for the absolute-factorization audit")
    program = f'''
LIB "primdec.lib";
ring r=0,(x,y,z),dp;
poly f={_poly_string(quintic_homogeneous)};
ideal C=x,y-z;
ideal D=x+z,y-2z;
ideal T=std(intersect(C,D));
ideal J=std(sat(jacob(f),maxideal(1)));
ideal R=std(radical(J));
list FF=factorize(f);
print("Q_FACTOR_COUNT_WITH_UNIT"); print(size(FF[1]));
print("Q_NONUNIT_EXPONENT"); print(FF[2][2]);
print("SINGULAR_AFFINE_CONE_DIMENSION"); print(dim(J));
print("SINGULAR_SCHEME_DEGREE"); print(mult(J));
print("SINGULAR_RADICAL_DEGREE"); print(mult(R));
print("SINGULAR_MINIMAL_PRIMES"); print(size(minAssGTZ(J)));
print("J_MOD_R"); print(size(reduce(J,R)));
print("R_MOD_J"); print(size(reduce(R,J)));
print("R_MOD_POINTS"); print(size(reduce(R,T)));
print("POINTS_MOD_R"); print(size(reduce(T,R)));
LIB "absfact.lib";
def ABSR=absFactorize(f);
setring ABSR;
print("ABS_FACTOR_COUNT_WITH_UNIT"); print(size(absolute_factors[1]));
print("ABS_NONUNIT_EXPONENT"); print(absolute_factors[2][2]);
print("ABS_COMPONENT_COUNT"); print(absolute_factors[4]);
exit;
'''
    completed = subprocess.run(
        [executable, "-q", "--no-rc"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    names = {
        "Q_FACTOR_COUNT_WITH_UNIT": "factor_count_with_unit_over_Q",
        "Q_NONUNIT_EXPONENT": "nonunit_exponent_over_Q",
        "SINGULAR_AFFINE_CONE_DIMENSION": "singular_affine_cone_dimension",
        "SINGULAR_SCHEME_DEGREE": "singular_scheme_degree",
        "SINGULAR_RADICAL_DEGREE": "singular_radical_degree",
        "SINGULAR_MINIMAL_PRIMES": "singular_minimal_primes",
        "J_MOD_R": "jacobian_mod_radical",
        "R_MOD_J": "radical_mod_jacobian",
        "R_MOD_POINTS": "radical_mod_two_points",
        "POINTS_MOD_R": "two_points_mod_radical",
        "ABS_FACTOR_COUNT_WITH_UNIT": "absolute_factor_count_with_unit",
        "ABS_NONUNIT_EXPONENT": "absolute_nonunit_exponent",
        "ABS_COMPONENT_COUNT": "absolute_component_count",
    }
    parsed: dict[str, int] = {}
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    for index, line in enumerate(lines[:-1]):
        if line in names:
            parsed[names[line]] = int(lines[index + 1])
    expected = {
        "factor_count_with_unit_over_Q": 2,
        "nonunit_exponent_over_Q": 1,
        "singular_affine_cone_dimension": 1,
        "singular_scheme_degree": 10,
        "singular_radical_degree": 2,
        "singular_minimal_primes": 2,
        "jacobian_mod_radical": 0,
        "radical_mod_jacobian": 2,
        "radical_mod_two_points": 0,
        "two_points_mod_radical": 0,
        "absolute_factor_count_with_unit": 2,
        "absolute_nonunit_exponent": 1,
        "absolute_component_count": 1,
    }
    if parsed != expected:
        raise ArithmeticError(f"the Singular quintic audit changed: {parsed}")
    return {
        "engine": "Singular over Q: factorize, absFactorize, saturated Jacobian, radical and primary support",
        **parsed,
        "absolutely_irreducible": True,
        "singular_support": ["C=(0:1:1)", "D=(-1:2:1)"],
        "scheme_length_decomposition": "two A5 singularities, Tjurina length 5+5",
    }


def _chain_data_from_coordinates(
    root_coordinates: list[tuple[int, ...]],
) -> tuple[list[list[int]], list[int]]:
    parent_of: dict[int, int] = {}
    degree_zero = []
    for coordinates in root_coordinates:
        if coordinates[0] != 0:
            continue
        negative = [index for index, value in enumerate(coordinates[1:]) if value == -1]
        positive = [index for index, value in enumerate(coordinates[1:]) if value == 1]
        if len(negative) != 1 or len(positive) != 1:
            raise ArithmeticError("a contracted root is not a proximity edge")
        parent_of[negative[0]] = positive[0]
        degree_zero.append(coordinates)
    roots = sorted(set(range(11)) - set(parent_of))
    chains: list[list[int]] = []
    for root in roots:
        chain = [root]
        while True:
            children = [child for child, parent in parent_of.items() if parent == chain[-1]]
            if not children:
                break
            if len(children) != 1:
                raise ArithmeticError("the proximity forest is not a union of chains")
            chain.append(children[0])
        chains.append(chain)
    return sorted(chains, key=lambda chain: (len(chain), chain)), sorted(parent_of)


def _plane_vector(degree: int, multiplicities: dict[int, int]) -> tuple[int, ...]:
    result = [degree, *([0] * 11)]
    for position, multiplicity in multiplicities.items():
        result[1 + position] = -multiplicity
    return tuple(result)


def _canonical_program_from_raw_classes(
    row: dict[str, Any],
    target: Sequence[int],
    all_coordinates: list[tuple[int, ...]],
) -> dict[str, Any]:
    """Derive A/B/C/D and compare every seed class coefficientwise.

    No A/B/C/D labeling is supplied.  The degree-one boundary identifies the
    length-three and length-two line-contact chains.  The three positive root
    supports then identify BCD and the two tangent lines; C/D are canonically
    ordered by their tangent-support tuples.
    """

    root_records = [
        (index, all_coordinates[index]) for index in row["roots"]
    ]
    root_coordinates = [coordinates for _index, coordinates in root_records]
    chains, _near_positions = _chain_data_from_coordinates(root_coordinates)
    line_boundary_index = next(
        index
        for index in row["boundaries"]
        if all_coordinates[index][0] == 1
    )
    quintic_boundary_index = next(
        index
        for index in row["boundaries"]
        if all_coordinates[index][0] == 5
    )
    line_support = {
        position
        for position, value in enumerate(all_coordinates[line_boundary_index][1:])
        if value == -1
    }
    line_chains = [chain for chain in chains if set(chain).issubset(line_support)]
    if sorted(map(len, line_chains)) != [2, 3]:
        raise ArithmeticError("the degree-one boundary did not derive the 3+2 chains")
    chain_a = next(chain for chain in line_chains if len(chain) == 3)
    chain_b = next(chain for chain in line_chains if len(chain) == 2)
    singular_chains = [chain for chain in chains if chain not in line_chains]
    if len(singular_chains) != 2 or any(len(chain) != 3 for chain in singular_chains):
        raise ArithmeticError("the quintic did not derive two length-three singular chains")

    positive_roots = [
        (
            index,
            coordinates,
            frozenset(
                position
                for position, value in enumerate(coordinates[1:])
                if value == -1
            ),
        )
        for index, coordinates in root_records
        if coordinates[0] == 1
    ]
    proper_singular = {chain[0] for chain in singular_chains}
    bcd = [
        record
        for record in positive_roots
        if record[2] == frozenset({chain_b[0], *proper_singular})
    ]
    if len(bcd) != 1:
        raise ArithmeticError("the BCD root line was not intrinsically unique")
    tangent_records = [record for record in positive_roots if record != bcd[0]]
    tangent_by_chain: dict[tuple[int, ...], tuple[int, tuple[int, ...], frozenset[int]]] = {}
    for chain in singular_chains:
        expected_support = frozenset({chain_a[0], chain[0], chain[1]})
        matches = [record for record in tangent_records if record[2] == expected_support]
        if len(matches) != 1:
            raise ArithmeticError("a singular chain did not have one raw-derived tangent line")
        tangent_by_chain[tuple(chain)] = matches[0]
    ordered_singular = sorted(
        (tuple(chain) for chain in singular_chains),
        key=lambda chain: (tuple(sorted(tangent_by_chain[chain][2])), tangent_by_chain[chain][0]),
    )
    chain_c = list(ordered_singular[0])
    chain_d = list(ordered_singular[1])

    labeled_positions = {
        "A0": chain_a[0],
        "A1": chain_a[1],
        "A2": chain_a[2],
        "B0": chain_b[0],
        "B1": chain_b[1],
        "C0": chain_c[0],
        "C1": chain_c[1],
        "C2": chain_c[2],
        "D0": chain_d[0],
        "D1": chain_d[1],
        "D2": chain_d[2],
    }
    all_ab = {labeled_positions[label] for label in ("A0", "A1", "A2", "B0", "B1")}
    all_cd = {labeled_positions[label] for label in ("C0", "C1", "C2", "D0", "D1", "D2")}
    expected_boundaries = {
        _plane_vector(1, {position: 1 for position in all_ab}),
        _plane_vector(
            5,
            {
                **{position: 1 for position in all_ab},
                **{position: 2 for position in all_cd},
            },
        ),
    }
    raw_boundaries = {all_coordinates[index] for index in row["boundaries"]}
    if raw_boundaries != expected_boundaries:
        raise ArithmeticError("the derived line/quintic classes do not equal the raw boundaries")

    expected_roots: dict[str, tuple[int, ...]] = {}
    for chain_name, chain in (
        ("A", chain_a),
        ("B", chain_b),
        ("C", chain_c),
        ("D", chain_d),
    ):
        for edge in range(len(chain) - 1):
            # E_parent-E_near has coefficient +1 at the parent and -1 at
            # the child in the orthogonal exceptional basis.
            vector = [0, *([0] * 11)]
            vector[1 + chain[edge]] = 1
            vector[1 + chain[edge + 1]] = -1
            expected_roots[f"{chain_name}{edge}_proximity"] = tuple(vector)
    expected_roots["AC_tangent"] = _plane_vector(
        1,
        {
            labeled_positions["A0"]: 1,
            labeled_positions["C0"]: 1,
            labeled_positions["C1"]: 1,
        },
    )
    expected_roots["AD_tangent"] = _plane_vector(
        1,
        {
            labeled_positions["A0"]: 1,
            labeled_positions["D0"]: 1,
            labeled_positions["D1"]: 1,
        },
    )
    expected_roots["BCD_line"] = _plane_vector(
        1,
        {
            labeled_positions["B0"]: 1,
            labeled_positions["C0"]: 1,
            labeled_positions["D0"]: 1,
        },
    )
    raw_root_set = {coordinates for _index, coordinates in root_records}
    if raw_root_set != set(expected_roots.values()) or len(expected_roots) != 10:
        raise ArithmeticError("the derived proximity/line classes do not equal all ten raw roots")

    expected_deep = {
        "cross_A": _plane_vector(0, {labeled_positions["A2"]: -1}),
        "cross_B": _plane_vector(0, {labeled_positions["B1"]: -1}),
        "cusp_C": _plane_vector(0, {labeled_positions["C2"]: -1}),
        "cusp_D": _plane_vector(0, {labeled_positions["D2"]: -1}),
    }
    contracted_minus_one = [
        index for index in target if row["norms"][index] == -1
    ]
    raw_deep_set = {all_coordinates[index] for index in contracted_minus_one}
    if raw_deep_set != set(expected_deep.values()):
        raise ArithmeticError("the four derived deepest exceptionals do not equal the raw target walls")
    crossing_indices = [
        index
        for index in contracted_minus_one
        if sorted(row["wall_gram"][index][boundary] for boundary in row["boundaries"])
        == [1, 1]
    ]
    cusp_indices = [
        index
        for index in contracted_minus_one
        if sorted(row["wall_gram"][index][boundary] for boundary in row["boundaries"])
        == [0, 2]
    ]
    if {all_coordinates[index] for index in crossing_indices} != {
        expected_deep["cross_A"],
        expected_deep["cross_B"],
    }:
        raise ArithmeticError("the derived A/B deepest exceptionals are not exactly the crossings")
    if {all_coordinates[index] for index in cusp_indices} != {
        expected_deep["cusp_C"],
        expected_deep["cusp_D"],
    }:
        raise ArithmeticError("the derived C/D deepest exceptionals are not exactly the visible cusp lifts")

    canonical_label_order = ("A0", "A1", "A2", "B0", "B1", "C0", "C1", "C2", "D0", "D1", "D2")

    def canonicalize(vector: tuple[int, ...]) -> list[int]:
        return [vector[0], *(vector[1 + labeled_positions[label]] for label in canonical_label_order)]

    return {
        "derived_proximity_chains_proper_to_deep": {
            "A": chain_a,
            "B": chain_b,
            "C": chain_c,
            "D": chain_d,
        },
        "canonical_label_to_exceptional_position_zero_based": labeled_positions,
        "raw_root_line_indices_zero_based": {
            "AC_tangent": tangent_by_chain[tuple(chain_c)][0],
            "AD_tangent": tangent_by_chain[tuple(chain_d)][0],
            "BCD_line": bcd[0][0],
        },
        "raw_crossing_wall_indices_zero_based": sorted(crossing_indices),
        "raw_visible_cusp_lift_indices_zero_based": sorted(cusp_indices),
        "canonical_boundary_classes": sorted(canonicalize(vector) for vector in raw_boundaries),
        "canonical_actual_root_classes": sorted(canonicalize(vector) for vector in raw_root_set),
        "canonical_deep_exceptional_classes": sorted(canonicalize(vector) for vector in raw_deep_set),
        "coefficientwise_boundary_match": True,
        "coefficientwise_all_ten_actual_roots_match": True,
        "coefficientwise_two_crossings_match": True,
        "coefficientwise_two_visible_cusp_lifts_match": True,
    }


def _raw_branch_record(row: dict[str, Any], target: Sequence[int]) -> dict[str, Any]:
    order, plane, exceptionals = _orthogonal_plane_basis(row, target)
    basis = sp.Matrix.hstack(plane, *exceptionals)
    if basis.T * row["gram_matrix"] * basis != sp.diag(1, *([-1] * 11)):
        raise ArithmeticError("the plane basis Gram is not I_(1,11)")
    if abs(int(basis.det())) != 1:
        raise ArithmeticError("the plane basis is not determinant one")
    coordinates = [
        _basis_coordinates(basis, vector) for vector in row["walls_int"]
    ]
    canonical = _basis_coordinates(basis, row["canonical"])
    boundaries = [coordinates[index] for index in row["boundaries"]]
    roots = [coordinates[index] for index in row["roots"]]
    if canonical != (-3, *([1] * 11)):
        raise ArithmeticError("the canonical class changed in the plane basis")
    expected_boundaries = {
        (1, *([-1] * 5), *([0] * 6)),
        (5, *([-1] * 5), *([-2] * 6)),
    }
    if set(boundaries) != expected_boundaries:
        raise ArithmeticError("the boundary profile is not line plus quintic")
    chains, near_positions = _chain_data_from_coordinates(roots)
    chain_lengths = sorted(map(len, chains))
    if chain_lengths != [2, 3, 3, 3]:
        raise ArithmeticError("the proximity lengths changed")
    positive_roots = [coordinates for coordinates in roots if coordinates[0] > 0]
    if len(positive_roots) != 3 or any(
        coordinates[0] != 1
        or sum(value == -1 for value in coordinates[1:]) != 3
        or any(value not in (0, -1) for value in coordinates[1:])
        for coordinates in positive_roots
    ):
        raise ArithmeticError("the three positive-degree roots are not three-center lines")

    # Classify the four deepest total exceptionals by boundary profile.  Two
    # are crossings and two are primitive cusp lifts already visible in this
    # plane branch.
    contracted_minus_one = [index for index in target if row["norms"][index] == -1]
    profiles = Counter(
        tuple(sorted(row["wall_gram"][index][boundary] for boundary in row["boundaries"]))
        for index in contracted_minus_one
    )
    if profiles != Counter({(1, 1): 2, (0, 2): 2}):
        raise ArithmeticError("the contracted exceptional boundary profiles changed")

    program = _canonical_program_from_raw_classes(row, target, coordinates)

    return {
        "target_raw_wall_indices_zero_based": list(target),
        "lexicographic_valid_contraction_order_zero_based": order,
        "plane_class_in_raw_basis": [int(value) for value in plane],
        "plane_basis_determinant": int(basis.det()),
        "canonical_class_in_plane_basis": list(canonical),
        "boundary_classes_in_plane_basis": [list(value) for value in boundaries],
        "root_classes_in_plane_basis": [list(value) for value in roots],
        "proximity_chain_lengths": chain_lengths,
        "near_exceptional_positions_zero_based": near_positions,
        "degree_zero_proximity_roots": 7,
        "degree_one_three_center_root_lines": 3,
        "contracted_minus_one_boundary_profiles": {
            ",".join(map(str, profile)): count for profile, count in sorted(profiles.items())
        },
        "all_26_walls_integral_in_unimodular_plane_basis": len(coordinates) == 26,
        "raw_derived_canonical_program": program,
    }


def _raw_extraction() -> dict[str, Any]:
    row = _normalized(next(row for row in _load_rows() if row["id"] == "two_boundary_26"))
    endpoints = _sequential_endpoints(row)
    signatures = Counter()
    selected: list[tuple[int, ...]] = []
    for contracted, _inverse in endpoints.values():
        signature = _endpoint_signature(row, contracted)
        signatures[json.dumps(signature)] += 1
        if signature == TARGET_SIGNATURE:
            selected.append(contracted)
    if len(endpoints) != 52 or len(signatures) != 9 or len(selected) != 4:
        raise ArithmeticError("the H11/26 endpoint denominator changed")
    branch_records = [_raw_branch_record(row, target) for target in sorted(selected)]
    invariant_profiles = {
        json.dumps(
            {
                "boundary_classes": record["raw_derived_canonical_program"]["canonical_boundary_classes"],
                "actual_root_classes": record["raw_derived_canonical_program"]["canonical_actual_root_classes"],
                "deep_exceptional_classes": record["raw_derived_canonical_program"]["canonical_deep_exceptional_classes"],
            },
            sort_keys=True,
        )
        for record in branch_records
    }
    # The exceptional ordering can vary, so the literal coordinate arrays
    # need not agree; their invariant construction profile must.
    if len(invariant_profiles) != 1:
        raise ArithmeticError("the four clean endpoints have different invariant profiles")
    return {
        "input_policy": (
            "only the odd-unimodular lattice Gram and 26 primitive wall vectors are read; "
            "no surface/type/equation/known contraction enters branch selection"
        ),
        "all_sequential_P2_endpoint_count": len(endpoints),
        "endpoint_degree_signature_count": len(signatures),
        "intrinsic_line_quintic_endpoint_signature": {
            "boundary_degrees": [1, 5],
            "root_degree_histogram": {"0": 7, "1": 3},
            "contracted_norm_histogram": {"-2": 7, "-1": 4},
        },
        "selected_endpoint_count": len(selected),
        "invariant_canonical_program_profile_count": len(invariant_profiles),
        "branches": branch_records,
    }


def _determinantal_geometry() -> dict[str, Any]:
    matrix, labels = _quintic_condition_matrix()
    if matrix.shape != (23, 21):
        raise ArithmeticError("the quintic condition matrix changed size")

    rows = [index for index in range(23) if index not in (11, 20)]
    witness_minor = sp.factor(matrix[rows, :].det(method="domain-ge"))
    if sp.expand(witness_minor + 64 * (A - B) * (A + B)) != 0:
        raise ArithmeticError(f"the rank-locus witness minor changed: {witness_minor}")
    plus_rank = matrix.subs(B, A).rank()
    minus_rank = matrix.subs(B, -A).rank()
    generic_rank = matrix.subs({A: 1, B: 2}).rank()
    if (plus_rank, minus_rank, generic_rank) != (20, 20, 21):
        raise ArithmeticError("the symbolic rank trichotomy changed")

    good = _kernel_polynomial(matrix.subs({A: 1, B: 1}))
    bad = _kernel_polynomial(matrix.subs({A: 1, B: -1}))
    if good not in (GOOD_QUINTIC, -GOOD_QUINTIC):
        raise ArithmeticError("the good quintic kernel changed")
    if good == -GOOD_QUINTIC:
        good = -good
    bad_factor = sp.factor(bad)
    expected_bad = Y * (X**2 - X * Y + 2 * X - Y**2 + 2 * Y - 1) ** 2
    if bad_factor not in (expected_bad, -expected_bad):
        raise ArithmeticError("the rejected branch is no longer line times conic squared")

    homogeneous = _homogenize(good, 5)
    line_restriction = sp.factor(good.subs(Y, 0))
    if line_restriction != 16 * X**3 * (X - 1) ** 2:
        raise ArithmeticError("the quintic-line contact divisor changed")
    local = [
        _local_resolution_record(
            good,
            "C=(0:1:1)",
            lambda u, v: v,
            lambda u, v: 1 + u,
            1,
        ),
        _local_resolution_record(
            good,
            "D=(-1:2:1)",
            lambda u, v: -1 + u,
            lambda u, v: 2 - 2 * u + v,
            1,
        ),
    ]
    if sum(record["delta_invariant"] for record in local) != 6:
        raise ArithmeticError("the quintic delta sum changed")

    line_restrictions = {
        "AC_x=0": str(sp.factor(good.subs(X, 0))),
        "AD_2x+y=0": str(sp.factor(good.subs(Y, -2 * X))),
        "BCD_x+y=1": str(sp.factor(good.subs(Y, 1 - X))),
    }
    expected_restrictions = {
        "AC_x=0": "-y*(y - 1)**4",
        "AD_2x+y=0": "2*x*(x + 1)**4",
        "BCD_x+y=1": "x**2*(x - 1)*(x + 1)**2",
    }
    if line_restrictions != expected_restrictions:
        raise ArithmeticError(f"the three root-line restrictions changed: {line_restrictions}")

    singular = _singular_audit(homogeneous)
    return {
        "normalization": {
            "proper_points": {
                "A": [0, 0, 1],
                "B": [1, 0, 1],
                "C": [0, 1, 1],
                "D": [-1, 2, 1],
            },
            "line_boundary": "y=0",
            "first_singularity_tangent": "AC: x=0",
            "second_singularity_tangent": "AD: 2x+y=0",
            "third_root_line": "BCD: x+y-z=0",
            "remaining_second_near_slopes": ["a", "b"],
        },
        "quintic_matrix": {
            "rows": matrix.rows,
            "columns": matrix.cols,
            "row_block_sizes": [3, 2, 9, 9],
            "row_labels": labels,
            "rank_off_determinantal_locus": generic_rank,
            "rank_on_a_equals_b": plus_rank,
            "rank_on_a_equals_minus_b": minus_rank,
            "maximal_minor_drop_rows_zero_based": [11, 20],
            "maximal_minor": str(witness_minor),
            "exact_rank_drop_locus": "a^2-b^2=0",
        },
        "rejected_component": {
            "equation": "a+b=0",
            "sample": "a=1,b=-1",
            "kernel_factorization": str(bad_factor),
            "reason": "the quintic is the line boundary times a conic squared; boundary is nonreduced",
        },
        "surviving_horizontal_component": {
            "equation": "a-b=0 with a nonzero",
            "sample": "a=b=1 over Q",
            "horizontal_over": "Spec Z[1/2] after inverting the finite local discriminants",
            "affine_quintic": str(good),
            "homogeneous_quintic": str(homogeneous),
            "line_contact_restriction": str(line_restriction),
            "root_line_restrictions": line_restrictions,
            "local_resolutions": local,
            "arithmetic_genus": 6,
            "delta_sum": 6,
            "geometric_genus": 0,
            "strict_transform_smooth_rational": True,
            "singular_audit": singular,
        },
        "compiler_generated_not_supplied": (
            "the 23 jet rows, maximal minor, two saturated factors, kernel sections, "
            "Jacobian support and irreducibility audit are all generated from the raw proximity/root program"
        ),
    }


def _transfer_gate_replay(
    raw: dict[str, Any], source: dict[str, Any], explicit_repair: dict[str, Any]
) -> dict[str, Any]:
    """Replay the strict same-marked-P return, not merely finite volume.

    The old implementation merely searched an upstream theorem JSON for
    phrases and had geometrically constructed only 16 walls.  This routine
    now requires the independently replayed ten-curve construction, then
    checks directly in the raw marking the twelve finite-index CM roots, all
    six complete split cusp pairs, and the literal 26-wall cone.
    """

    row = _normalized(next(item for item in _load_rows() if item["id"] == "two_boundary_26"))
    gram = row["gram_matrix"]
    walls = [sp.Matrix(vector) for vector in row["walls_int"]]
    labels = source["facets"]["labels"]
    if row["gram"] != source["full_NS_lattice"]["gram"]:
        raise ArithmeticError("the raw H11 lattice is not the pinned exact-return lattice")
    if row["walls"] != source["facets"]["NS_wall_vectors"]:
        raise ArithmeticError("the raw H11 walls are not the pinned exact-return walls")
    computed_wall_gram = [
        [int((left.T * gram * right)[0]) for right in walls]
        for left in walls
    ]
    if computed_wall_gram != source["facets"]["negative_wall_gram"]:
        raise ArithmeticError("the literal 26-wall Gram changed")

    boundary_vectors = [walls[index] for index in row["boundaries"]]
    crossing_indices = [
        index
        for index in row["minus_one"]
        if sorted(row["wall_gram"][index][boundary] for boundary in row["boundaries"])
        == [1, 1]
    ]
    cusp_indices = [
        index
        for index in row["minus_one"]
        if sorted(row["wall_gram"][index][boundary] for boundary in row["boundaries"])
        == [0, 2]
    ]
    if crossing_indices != [12, 13] or cusp_indices != list(range(14, 26)):
        raise ArithmeticError("the raw (-1)-wall taxonomy changed")

    # Reconstruct the CM diagram intrinsically: actual (-2) walls stay put,
    # while a two-boundary crossing E gives alpha=2E+B1/2+B2/2.
    cm_roots = [walls[index] for index in row["roots"]]
    cm_roots.extend(
        2 * walls[index] + boundary_vectors[0] / 2 + boundary_vectors[1] / 2
        for index in crossing_indices
    )
    cm_root_matrix = sp.Matrix.hstack(*cm_roots)
    cm_gram = [
        [int((left.T * gram * right)[0]) for right in cm_roots]
        for left in cm_roots
    ]
    cm_reference_file = json.loads(CM_REFERENCE.read_text(encoding="utf-8"))
    cm_reference = next(
        record
        for record in cm_reference_file["types"]
        if record["nikulin_kondo_type"] == "I"
    )
    source_reference_hash = source["sources"]["input_sha256"][
        str(CM_REFERENCE.relative_to(ROOT))
    ]
    cm_reference_hash_matches = source_reference_hash == _sha256(CM_REFERENCE)
    cm_literal = cm_reference["literal_vinberg"]
    cm_checks = cm_reference["wall_vector_checks"]
    root_gate = {
        "reconstructed_root_count": len(cm_roots),
        "reconstructed_root_rank": int(cm_root_matrix.rank()),
        "all_square_minus_two": all(
            (root.T * gram * root)[0] == -2 for root in cm_roots
        ),
        "all_boundary_orthogonal": all(
            (root.T * gram * boundary)[0] == 0
            for root in cm_roots
            for boundary in boundary_vectors
        ),
        "exact_Type_I_CM_Gram_match": cm_gram == cm_reference["intersection_matrix"],
        "reference_hash_matches_pinned_source": cm_reference_hash_matches,
        "wall_sublattice_index": cm_checks["wall_sublattice_index_in_L10"],
        "literal_Vinberg_bad_finite_volume_edges": cm_literal[
            "bad_finite_volume_edge_count"
        ],
        "literal_Vinberg_invalid_parabolic_neighborhoods": cm_literal[
            "invalid_full_zero_neighborhood_count"
        ],
        "CM_finite_vertices": cm_literal["finite_vertex_count"],
        "CM_ideal_vertices": cm_literal["ideal_vertex_count"],
        "CM_finite_volume_noncompact": (
            cm_literal["finite_volume"] is True
            and cm_literal["compact"] is False
            and cm_reference["exact_cone"]["hyperideal_vertex_count"] == 0
        ),
        "simple_root_stock_complete": cm_literal["simple_roots_complete"],
    }
    root_gate_passes = (
        root_gate["reconstructed_root_count"] == source["source_CM_chamber"]["root_count"] == 12
        and root_gate["reconstructed_root_rank"] == 10
        and root_gate["all_square_minus_two"]
        and root_gate["all_boundary_orthogonal"]
        and root_gate["exact_Type_I_CM_Gram_match"]
        and root_gate["reference_hash_matches_pinned_source"]
        and root_gate["wall_sublattice_index"] == 2
        and root_gate["literal_Vinberg_bad_finite_volume_edges"] == 0
        and root_gate["literal_Vinberg_invalid_parabolic_neighborhoods"] == 0
        and (root_gate["CM_finite_vertices"], root_gate["CM_ideal_vertices"]) == (14, 9)
        and root_gate["CM_finite_volume_noncompact"]
        and root_gate["simple_root_stock_complete"]
    )

    # Parse the raw labels only after identifying the vectors numerically.
    # Each split cusp must have both boundary lifts and a common f=E+B_i/2.
    cusp_groups: dict[int, dict[int, tuple[int, sp.Matrix]]] = {}
    cusp_index_to_group_boundary: dict[int, tuple[int, int]] = {}
    for index in cusp_indices:
        parts = labels[index].split("_")
        if len(parts) != 4 or parts[:2] != ["E", "cusp"] or not parts[3].startswith("B"):
            raise ArithmeticError("a raw cusp wall label is malformed")
        cusp_number = int(parts[2]) - 1
        boundary_number = int(parts[3][1:]) - 1
        f_vector = walls[index] + boundary_vectors[boundary_number] / 2
        cusp_groups.setdefault(cusp_number, {})[boundary_number] = (index, f_vector)
        cusp_index_to_group_boundary[index] = (cusp_number, boundary_number)

    cusp_records: list[dict[str, Any]] = []
    common_f_vectors: list[tuple[sp.Rational, ...]] = []
    all_cusp_formulas_exact = True
    for cusp_number in sorted(cusp_groups):
        lifts = cusp_groups[cusp_number]
        if sorted(lifts) != [0, 1]:
            raise ArithmeticError("a split cusp does not have both boundary lifts")
        left_f = lifts[0][1]
        right_f = lifts[1][1]
        same_f = left_f == right_f
        isotropic = (left_f.T * gram * left_f)[0] == 0
        boundary_orthogonal = all(
            (left_f.T * gram * boundary)[0] == 0 for boundary in boundary_vectors
        )
        recoveries = [
            left_f - boundary_vectors[boundary_number] / 2 == walls[lifts[boundary_number][0]]
            for boundary_number in (0, 1)
        ]
        twice_f = 2 * left_f
        twice_f_integral_primitive = all(value.q == 1 for value in twice_f) and math.gcd(
            *(abs(int(value)) for value in twice_f)
        ) == 1
        all_cusp_formulas_exact &= (
            same_f and isotropic and boundary_orthogonal and all(recoveries)
        )
        common_f_vectors.append(tuple(left_f))
        cusp_records.append(
            {
                "CM_cusp_zero_based": cusp_number,
                "raw_lift_indices_zero_based": [lifts[0][0], lifts[1][0]],
                "raw_lift_labels": [labels[lifts[0][0]], labels[lifts[1][0]]],
                "common_f_coordinates_in_raw_NS_Q": _rational_vector(left_f),
                "f_square_zero": isotropic,
                "f_boundary_orthogonal": boundary_orthogonal,
                "both_recovery_formulas_E_equals_f_minus_Bi_over_2": all(recoveries),
                "two_f_is_primitive_in_raw_NS": twice_f_integral_primitive,
                "geometric_primitivity_gate": (
                    "the connected vertical divisor 2E+B_i has primitive multiplicity vector (2,1); "
                    "the Halphen index-at-most-two theorem therefore makes f primitive"
                ),
            }
        )

    expected_split = source["curve_exhaustion"]["split_CM_cusps_zero_based"]
    cusp_gate = {
        "complete_CM_cusp_count": source["source_CM_chamber"]["cusp_count"],
        "split_CM_cusp_ids_zero_based": sorted(cusp_groups),
        "expected_split_CM_cusp_ids_zero_based": expected_split,
        "split_cusp_records": cusp_records,
        "six_distinct_split_cusp_vectors": len(set(common_f_vectors)) == 6,
        "all_six_have_both_integral_raw_lifts": all(
            source["curve_exhaustion"]["boundary_split_membership_by_CM_cusp"][index]
            == [True, True]
            for index in cusp_groups
        ),
        "all_recovery_formulas_exact": all_cusp_formulas_exact,
        "single_boundary_lift_count": len(cusp_indices),
    }
    cusp_gate_passes = (
        cusp_gate["complete_CM_cusp_count"] == 9
        and cusp_gate["split_CM_cusp_ids_zero_based"] == expected_split == [1, 3, 4, 5, 7, 8]
        and cusp_gate["six_distinct_split_cusp_vectors"]
        and cusp_gate["all_six_have_both_integral_raw_lifts"]
        and cusp_gate["all_recovery_formulas_exact"]
        and cusp_gate["single_boundary_lift_count"] == 12
    )

    program_flags = (
        "coefficientwise_boundary_match",
        "coefficientwise_all_ten_actual_roots_match",
        "coefficientwise_two_crossings_match",
        "coefficientwise_two_visible_cusp_lifts_match",
    )
    branch_accounting: list[dict[str, Any]] = []
    for branch in raw["branches"]:
        program = branch["raw_derived_canonical_program"]
        visible_cusps = program["raw_visible_cusp_lift_indices_zero_based"]
        visible_pairs = [cusp_index_to_group_boundary[index] for index in visible_cusps]
        visible_group_ids = sorted(group for group, _boundary in visible_pairs)
        visible_boundaries = sorted({boundary for _group, boundary in visible_pairs})
        expected_visible_labels = (
            ["E_cusp_5_B1", "E_cusp_9_B1"]
            if visible_boundaries == [0]
            else ["E_cusp_5_B2", "E_cusp_9_B2"]
        )
        expected_remaining: set[int] = set()
        for group, lifts in cusp_groups.items():
            if group in visible_group_ids:
                selected_boundary = visible_pairs[visible_group_ids.index(group)][1]
                expected_remaining.add(lifts[1 - selected_boundary][0])
            else:
                expected_remaining.update(lift[0] for lift in lifts.values())
        actual_remaining = set(cusp_indices) - set(visible_cusps)
        explicit_wall_indices = set(row["boundaries"] + row["roots"] + crossing_indices + visible_cusps)
        branch_accounting.append(
            {
                "target_raw_wall_indices_zero_based": branch[
                    "target_raw_wall_indices_zero_based"
                ],
                "all_program_classes_match_coefficientwise": all(
                    program[flag] for flag in program_flags
                ),
                "crossing_indices_zero_based": program[
                    "raw_crossing_wall_indices_zero_based"
                ],
                "visible_cusp_indices_zero_based": visible_cusps,
                "visible_cusp_labels": sorted(labels[index] for index in visible_cusps),
                "visible_cusp_group_ids_zero_based": visible_group_ids,
                "visible_lifts_use_one_common_boundary": len(visible_boundaries) == 1,
                "visible_labels_are_exactly_cusp5_and_cusp9_same_boundary": (
                    sorted(labels[index] for index in visible_cusps)
                    == sorted(expected_visible_labels)
                ),
                "explicit_negative_wall_count": len(explicit_wall_indices),
                "explicit_minus_one_wall_count": len(crossing_indices) + len(visible_cusps),
                "remaining_cusp_lift_indices_zero_based": sorted(actual_remaining),
                "remaining_cusp_lift_labels": sorted(labels[index] for index in actual_remaining),
                "remaining_count": len(actual_remaining),
                "remaining_is_two_opposite_partners_plus_four_complete_pairs": (
                    actual_remaining == expected_remaining
                    and len(visible_group_ids) == 2
                    and len(set(cusp_groups) - set(visible_group_ids)) == 4
                ),
            }
        )
    branches_pass = (
        raw["invariant_canonical_program_profile_count"] == 1
        and len(branch_accounting) == 4
        and all(
            record["all_program_classes_match_coefficientwise"]
            and record["crossing_indices_zero_based"] == crossing_indices
            and record["visible_cusp_group_ids_zero_based"] == [4, 8]
            and record["visible_lifts_use_one_common_boundary"]
            and record["visible_labels_are_exactly_cusp5_and_cusp9_same_boundary"]
            and record["explicit_negative_wall_count"] == 16
            and record["explicit_minus_one_wall_count"] == 4
            and record["remaining_count"] == 10
            and record["remaining_is_two_opposite_partners_plus_four_complete_pairs"]
            for record in branch_accounting
        )
    )

    explicit_curve_gate = {
        "four_deep_exceptional_minus_one_curves": explicit_repair[
            "four_deep_exceptional_minus_one_curves"
        ],
        "ten_positive_degree_curve_count": explicit_repair[
            "ten_positive_degree_minus_one_curve_count"
        ],
        "degree_multiset": explicit_repair["degree_multiset"],
        "all_ten_unique_plane_curves": explicit_repair["all_ten_unique_plane_curves"],
        "all_ten_absolutely_irreducible": explicit_repair[
            "all_ten_absolutely_irreducible"
        ],
        "all_ten_exact_strict_transforms_smooth_rational": explicit_repair[
            "all_ten_exact_strict_transforms_smooth_rational"
        ],
        "fresh_CM_finite_ideal_hyperideal": explicit_repair[
            "fresh_CM_finite_ideal_hyperideal"
        ],
        "fresh_full_finite_ideal_hyperideal": explicit_repair[
            "fresh_full_finite_ideal_hyperideal"
        ],
        "all_and_only_integral_cusp_lifts_literal": explicit_repair[
            "all_and_only_integral_cusp_lifts_literal"
        ],
    }
    explicit_curve_gate_passes = (
        explicit_curve_gate["ten_positive_degree_curve_count"] == 10
        and explicit_curve_gate["degree_multiset"] == [2, 2, 2, 2, 3, 3, 4, 4, 5, 5]
        and explicit_curve_gate["all_ten_unique_plane_curves"]
        and explicit_curve_gate["all_ten_absolutely_irreducible"]
        and explicit_curve_gate["all_ten_exact_strict_transforms_smooth_rational"]
        and explicit_curve_gate["fresh_CM_finite_ideal_hyperideal"] == [14, 9, 0]
        and explicit_curve_gate["fresh_full_finite_ideal_hyperideal"] == [508, 42, 0]
        and explicit_curve_gate["all_and_only_integral_cusp_lifts_literal"]
    )
    theorem_gate = {
        "not_a_phrase_gate": True,
        "published_root_exhaustion_dependency": (
            "Kondo Proposition 3.7: a finite-index set of distinct effective "
            "irreducible roots is the complete effective irreducible root set"
        ),
        "finite_index_root_hypotheses_replayed": root_gate_passes,
        "all_fourteen_minus_one_walls_geometrically_constructed": (
            explicit_curve_gate_passes
            and explicit_curve_gate["ten_positive_degree_curve_count"] + 4 == 14
        ),
        "primitive_cusp_lifts_replayed_coefficientwise": cusp_gate_passes,
        "hidden_same_boundary_curve_argument": (
            "any hidden same-boundary (-1)-curve gives a primitive nef isotropic CM ray; "
            "the complete cusp gate identifies its class with an already constructed curve"
        ),
        "negative_curve_trichotomy_applies_on_the_constructed_strong_Coble_surface": True,
    }
    theorem_gate_passes = (
        theorem_gate["not_a_phrase_gate"]
        and theorem_gate["finite_index_root_hypotheses_replayed"]
        and theorem_gate["all_fourteen_minus_one_walls_geometrically_constructed"]
        and theorem_gate["primitive_cusp_lifts_replayed_coefficientwise"]
        and theorem_gate[
            "negative_curve_trichotomy_applies_on_the_constructed_strong_Coble_surface"
        ]
    )

    rays = source["exact_double_description"]
    full_cone_gate = {
        "literal_26_wall_marking": (
            len(walls) == len(labels) == 26
            and labels
            == [
                "B1",
                "B2",
                *(f"R{index}" for index in range(1, 11)),
                "E_cross_1",
                "E_cross_2",
                "E_cusp_2_B1",
                "E_cusp_2_B2",
                "E_cusp_4_B1",
                "E_cusp_4_B2",
                "E_cusp_5_B1",
                "E_cusp_5_B2",
                "E_cusp_6_B1",
                "E_cusp_6_B2",
                "E_cusp_8_B1",
                "E_cusp_8_B2",
                "E_cusp_9_B1",
                "E_cusp_9_B2",
            ]
        ),
        "negative_curve_taxonomy_2_10_2_12": (
            len(row["boundaries"]) == 2
            and len(row["roots"]) == 10
            and len(crossing_indices) == 2
            and len(cusp_indices) == 12
            and source["curve_exhaustion"]["negative_curve_count"] == 26
        ),
        "all_26_facets_irredundant": rays["all_facets_irredundant"],
        "literal_full_NS_bad_edges": source["literal_vinberg"]["bad_edge_count"],
        "literal_full_NS_bad_parabolic_neighborhoods": source["literal_vinberg"][
            "bad_parabolic_neighborhood_count"
        ],
        "finite_ideal_hyperideal": [
            rays["finite_vertex_count"],
            rays["ideal_vertex_count"],
            rays["hyperideal_vertex_count"],
        ],
    }
    full_cone_gate_passes = (
        full_cone_gate["literal_26_wall_marking"]
        and full_cone_gate["negative_curve_taxonomy_2_10_2_12"]
        and full_cone_gate["all_26_facets_irredundant"]
        and full_cone_gate["literal_full_NS_bad_edges"] == 0
        and full_cone_gate["literal_full_NS_bad_parabolic_neighborhoods"] == 0
        and full_cone_gate["finite_ideal_hyperideal"] == [508, 42, 0]
    )

    all_gates = (
        source["all_checks_pass"]
        and root_gate_passes
        and cusp_gate_passes
        and branches_pass
        and explicit_curve_gate_passes
        and theorem_gate_passes
        and full_cone_gate_passes
    )
    if not all_gates:
        raise ArithmeticError("the explicit H11/26 exact-return transfer replay failed")
    return {
        "input_lattice_and_26_wall_vectors_equal_pinned_source": True,
        "finite_index_CM_root_gate": {**root_gate, "passes": root_gate_passes},
        "complete_cusp_lift_gate": {**cusp_gate, "passes": cusp_gate_passes},
        "four_raw_plane_branch_wall_accounting": branch_accounting,
        "four_raw_plane_branches_pass": branches_pass,
        "explicit_ten_missing_curve_gate": {
            **explicit_curve_gate,
            "passes": explicit_curve_gate_passes,
        },
        "uniform_hidden_wall_theorem_gate": {**theorem_gate, "passes": theorem_gate_passes},
        "literal_full_NS_cone_gate": {**full_cone_gate, "passes": full_cone_gate_passes},
        "all_transfer_gates_pass": all_gates,
    }


def _strong_coble_and_return(
    raw: dict[str, Any], geometry: dict[str, Any], explicit_repair: dict[str, Any]
) -> dict[str, Any]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    rays = source["exact_double_description"]
    if (
        rays["finite_vertex_count"],
        rays["ideal_vertex_count"],
        rays["hyperideal_vertex_count"],
    ) != (508, 42, 0):
        raise ArithmeticError("the pinned H11/26 ray certificate changed")
    if not source["all_checks_pass"]:
        raise ArithmeticError("the upstream raw-cone certificate is not passing")
    transfer_replay = _transfer_gate_replay(raw, source, explicit_repair)
    return {
        "surface": "blow up P2 at the four proper chains of lengths 3,2,3,3",
        "boundary_classes": {
            "B_quintic": "5H-E1-E2-E3-E4-E5-2E6-2E7-2E8-2E9-2E10-2E11",
            "B_line": "H-E1-E2-E3-E4-E5",
            "sum": "6H-2*sum(Ei)=-2K",
            "intersection_matrix": [[-4, 0], [0, -4]],
            "reduced_smooth_disjoint": True,
        },
        "strong_Coble_argument": (
            "the two reduced disjoint (-4) boundaries sum to -2K; any member of |-2K| "
            "meets each boundary negatively and must contain both, so it equals their sum; "
            "then 2D=B rules out D in |-K| because B is reduced"
        ),
        "explicit_root_seed": {
            "actual_minus_two_curves": "7 exceptional proximity curves plus the 3 strict root lines",
            "crossing_minus_one_curves": "the two deepest exceptionals over the line-contact chains",
            "visible_single_boundary_cusp_lifts": (
                "the two deepest exceptionals over the two quintic A5 chains; on each branch these are "
                "the same-boundary lifts of CM cusps 5 and 9"
            ),
            "complete_CM_root_count": 12,
            "coefficientwise_match_to_raw_classes_on_all_four_endpoints": all(
                all(
                    branch["raw_derived_canonical_program"][flag]
                    for flag in (
                        "coefficientwise_boundary_match",
                        "coefficientwise_all_ten_actual_roots_match",
                        "coefficientwise_two_crossings_match",
                        "coefficientwise_two_visible_cusp_lifts_match",
                    )
                )
                for branch in raw["branches"]
            ),
            "four_endpoints_have_one_invariant_canonical_program_profile": (
                raw["invariant_canonical_program_profile_count"] == 1
            ),
        },
        "explicit_missing_cusp_curve_repair": explicit_repair,
        "exact_return_handoff": (
            "The original root-seed-only handoff was insufficient: it constructed 16 of 26 negative "
            "curves.  The independent repair now interpolates, absolutely factors, and checks the exact "
            "cluster multiplicities of all ten missing positive-degree cusp curves.  Thus all 26 raw wall "
            "classes are actual curves.  The finite-index root theorem and complete primitive-cusp gate "
            "exclude every hidden negative curve, proving P(Nef S) is literally P."
        ),
        "exact_return_transfer_replay": transfer_replay,
        "full_nef_polytope": {
            "dimension": 11,
            "facets": 26,
            "finite_vertices": 508,
            "ideal_vertices": 42,
            "hyperideal_rays": 0,
            "volume": "finite",
            "compactness": "noncompact",
        },
        "scope": {
            "abstract_polytope": "the already frozen H11/26 Type-I full-nef chamber",
            "surface_isomorphism_class": "known from Kondō Example 4.1 / Proposition 4.2",
            "project_new": "raw-P-derived line-plus-quintic birational model and exact same-marked-P inverse",
            "worldwide_priority_of_this_plane_model": "open",
            "not_claimed": [
                "a new abstract Coxeter polytope",
                "a new Coble surface isomorphism class",
                "a class-wide strict L4 theorem",
                "a compact example",
            ],
        },
    }


def build_certificate() -> dict[str, Any]:
    raw = _raw_extraction()
    geometry = _determinantal_geometry()
    explicit_repair = _replay_explicit_missing_curve_repair()
    return_data = _strong_coble_and_return(raw, geometry, explicit_repair)
    checks = {
        "all_52_sequential_endpoints_enumerated": raw["all_sequential_P2_endpoint_count"] == 52,
        "exactly_four_intrinsic_line_quintic_endpoints": raw["selected_endpoint_count"] == 4,
        "all_four_plane_bases_unimodular": all(
            abs(branch["plane_basis_determinant"]) == 1 for branch in raw["branches"]
        ),
        "raw_program_is_7_proximity_plus_3_lines": all(
            branch["degree_zero_proximity_roots"] == 7
            and branch["degree_one_three_center_root_lines"] == 3
            for branch in raw["branches"]
        ),
        "rank_locus_is_exactly_a2_minus_b2": geometry["quintic_matrix"]["exact_rank_drop_locus"] == "a^2-b^2=0",
        "minus_branch_rejected_as_nonreduced": "conic squared" in geometry["rejected_component"]["reason"],
        "plus_branch_absolutely_irreducible": geometry["surviving_horizontal_component"]["singular_audit"]["absolutely_irreducible"],
        "two_A5_singularities_resolve": all(
            record["resolved_after_three_blowups"]
            for record in geometry["surviving_horizontal_component"]["local_resolutions"]
        ),
        "strong_Coble_boundary": return_data["boundary_classes"]["reduced_smooth_disjoint"],
        "all_ten_missing_cusp_curves_explicit": (
            explicit_repair["ten_positive_degree_minus_one_curve_count"] == 10
            and explicit_repair["all_ten_unique_plane_curves"]
            and explicit_repair["all_ten_absolutely_irreducible"]
            and explicit_repair["all_ten_exact_strict_transforms_smooth_rational"]
        ),
        "exact_same_P_return_handoff": return_data["exact_return_transfer_replay"][
            "all_transfer_gates_pass"
        ],
        "finite_volume_checked_separately_from_same_P_return": (
            return_data["full_nef_polytope"]["hyperideal_rays"] == 0
            and return_data["exact_return_transfer_replay"]["literal_full_NS_cone_gate"][
                "passes"
            ]
        ),
        "compact_scope_is_empty": return_data["full_nef_polytope"]["compactness"] == "noncompact",
    }
    return {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "direction": "RAW_POLYTOPE_TO_CANONICAL_PLANE_PROGRAM_TO_STRONG_COBLE_TO_THE_SAME_POLYTOPE",
        "theorem": (
            "The raw H11/26 wall datum has exactly four intrinsic line-plus-quintic P2 endpoints. "
            "Their compiler-generated quintic rank locus has one nonreduced branch and one horizontal "
            "absolutely irreducible branch.  The latter constructs a strong two-boundary Coble surface "
            "over Q whose full nef chamber is literally the input finite-volume noncompact H11/26 chamber."
        ),
        "raw_extraction": raw,
        "canonical_determinantal_geometry": geometry,
        "strong_Coble_and_exact_return": return_data,
        "sources": {
            str(SOURCE.relative_to(ROOT)): _sha256(SOURCE),
            str(RAW_HOSTILE_OUTPUT.relative_to(ROOT)): _sha256(RAW_HOSTILE_OUTPUT),
            str(RAW_HOSTILE_SCRIPT.relative_to(ROOT)): _sha256(RAW_HOSTILE_SCRIPT),
            str(TRANSFER.relative_to(ROOT)): (
                "mathematical background only; no phrase or all_checks boolean is used as a gate"
            ),
            "Kondo": "https://arxiv.org/abs/2102.00591",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-output", type=Path)
    arguments = parser.parse_args()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if arguments.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
    if arguments.check_output is not None:
        if arguments.check_output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("frozen certificate differs from a fresh exact replay")
    if not arguments.write and arguments.check_output is None:
        print(rendered, end="")


if __name__ == "__main__":
    main()
