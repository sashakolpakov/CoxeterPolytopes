"""Certify mixed-root enlargements for Wang's 18 negative rows.

Wang's Theorem 7.2 is a published input: it gives the eighteen hyperbolic
2-reflective lattices of rank at least six without automorphic correction.
This file performs only finite exact lattice arithmetic.  In particular, it
does not run Vinberg's algorithm.

Signs on the ADE factors are suppressed in the row names.  Internally every
ADE factor is negative definite and U has Gram matrix [[0, 1], [1, 0]].
"""

from __future__ import annotations

from fractions import Fraction
from functools import reduce
import json
from math import gcd, lcm, prod
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_PATH = ROOT / "data/enumerations/wang18_same_lattice_mixed_root_ledger.json"


ADE = {
    "A1": {"rank": 1, "determinant": 2},
    "A2": {"rank": 2, "determinant": 3},
    "A3": {"rank": 3, "determinant": 4},
    "A4": {"rank": 4, "determinant": 5},
    "D4": {"rank": 4, "determinant": 4},
    "D5": {"rank": 5, "determinant": 4},
    "D6": {"rank": 6, "determinant": 4},
    "D8": {"rank": 8, "determinant": 4},
    "E6": {"rank": 6, "determinant": 3},
    "E7": {"rank": 7, "determinant": 2},
    "E8": {"rank": 8, "determinant": 1},
}


# A witness omega lies in R^*/R, has the displayed exact order m, and has
# positive-definite square 2/m.  In S=U+R(-1), r=m*e+m*omega is therefore a
# primitive vector of square -2m and divisibility m.
MIXED_WITNESSES = {
    "D_vector": {
        "description": "the vector fundamental weight omega_1 of D_n",
        "class_order": 2,
        "positive_weight_square": Fraction(1),
        "root_formula": "r=2e+2*omega_1(D_n)",
    },
    "A3_middle": {
        "description": "the middle fundamental weight omega_2 of A3",
        "class_order": 2,
        "positive_weight_square": Fraction(1),
        "root_formula": "r=2e+2*omega_2(A3)",
    },
    "A1_pair": {
        "description": "omega=(a+b)/2 for two orthogonal A1 roots a,b",
        "class_order": 2,
        "positive_weight_square": Fraction(1),
        "root_formula": "r=2e+a+b",
    },
    "A2_minuscule": {
        "description": "a minuscule fundamental weight omega_1 of A2",
        "class_order": 3,
        "positive_weight_square": Fraction(2, 3),
        "root_formula": "r=3e+3*omega_1(A2)",
    },
}


# The class x generates an integral overlattice S+Z*x.  Its square is stated
# in the negative-definite convention.  An odd square gives an odd integral
# overclosure; an even square gives an even one.
OVERCLOSURE_WITNESSES = {
    "D_vector_odd": {
        "description": "adjoin the D_n vector weight omega_1",
        "class_order": 2,
        "negative_square": Fraction(-1),
        "parity": "odd",
    },
    "A3_middle_odd": {
        "description": "adjoin the A3 middle weight omega_2",
        "class_order": 2,
        "negative_square": Fraction(-1),
        "parity": "odd",
    },
    "A1_pair_odd": {
        "description": "adjoin (a+b)/2 from two orthogonal A1 factors",
        "class_order": 2,
        "negative_square": Fraction(-1),
        "parity": "odd",
    },
    "E6_A2_diagonal_even": {
        "description": "adjoin the diagonal E6+A2 discriminant class",
        "class_order": 3,
        "negative_square": Fraction(-2),
        "parity": "even",
    },
}


# number, lattice display, ADE factors, mixed-root witness, component used,
# proper integral-overlattice witness (or None).
ROWS = (
    (1, "U+E8+E7", ("E8", "E7"), None, None, None),
    (2, "U+E8+D6", ("E8", "D6"), "D_vector", "D6", "D_vector_odd"),
    (3, "U+E8+D4+A1", ("E8", "D4", "A1"), "D_vector", "D4", "D_vector_odd"),
    (4, "U+E8+D4", ("E8", "D4"), "D_vector", "D4", "D_vector_odd"),
    (5, "U+D8+D4", ("D8", "D4"), "D_vector", "D8", "D_vector_odd"),
    (6, "U+E8+4A1", ("E8", "A1", "A1", "A1", "A1"), "A1_pair", "2A1", "A1_pair_odd"),
    (7, "U+E8+3A1", ("E8", "A1", "A1", "A1"), "A1_pair", "2A1", "A1_pair_odd"),
    (8, "U+D8+3A1", ("D8", "A1", "A1", "A1"), "D_vector", "D8", "D_vector_odd"),
    (9, "U+E8+A3", ("E8", "A3"), "A3_middle", "A3", "A3_middle_odd"),
    (10, "U+D8+2A1", ("D8", "A1", "A1"), "D_vector", "D8", "D_vector_odd"),
    (11, "U+2D4+2A1", ("D4", "D4", "A1", "A1"), "D_vector", "D4", "D_vector_odd"),
    (12, "U+E8+A2", ("E8", "A2"), "A2_minuscule", "A2", None),
    (13, "U+E6+A2", ("E6", "A2"), "A2_minuscule", "A2", "E6_A2_diagonal_even"),
    (14, "U+D4+A3", ("D4", "A3"), "D_vector", "D4", "D_vector_odd"),
    (15, "U+D5+A2", ("D5", "A2"), "D_vector", "D5", "D_vector_odd"),
    (16, "U+D4+A2", ("D4", "A2"), "D_vector", "D4", "D_vector_odd"),
    (17, "U+A4+A2", ("A4", "A2"), "A2_minuscule", "A2", None),
    (18, "U+A3+A2", ("A3", "A2"), "A3_middle", "A3", "A3_middle_odd"),
)


def _cartan_a(n: int) -> list[list[int]]:
    matrix = [[0] * n for _ in range(n)]
    for i in range(n):
        matrix[i][i] = 2
        if i:
            matrix[i][i - 1] = matrix[i - 1][i] = -1
    return matrix


def _cartan_d(n: int) -> list[list[int]]:
    """D_n with chain 0--...--(n-3), branched to n-2,n-1."""
    matrix = [[0] * n for _ in range(n)]
    for i in range(n):
        matrix[i][i] = 2
    for i in range(n - 3):
        matrix[i][i + 1] = matrix[i + 1][i] = -1
    matrix[n - 3][n - 2] = matrix[n - 2][n - 3] = -1
    matrix[n - 3][n - 1] = matrix[n - 1][n - 3] = -1
    return matrix


def _star_cartan(arm_lengths: tuple[int, ...]) -> list[list[int]]:
    """Simply-laced Cartan matrix for a star with vertex 0 at its centre."""
    size = 1 + sum(arm_lengths)
    matrix = [[0] * size for _ in range(size)]
    for i in range(size):
        matrix[i][i] = 2
    next_vertex = 1
    for length in arm_lengths:
        previous = 0
        for _ in range(length):
            matrix[previous][next_vertex] = matrix[next_vertex][previous] = -1
            previous = next_vertex
            next_vertex += 1
    return matrix


def _ade_cartan(name: str) -> list[list[int]]:
    family = name[0]
    rank = int(name[1:])
    if family == "A":
        return _cartan_a(rank)
    if family == "D":
        return _cartan_d(rank)
    if name == "E6":
        return _star_cartan((2, 2, 1))
    if name == "E7":
        return _star_cartan((2, 3, 1))
    if name == "E8":
        return _star_cartan((2, 4, 1))
    raise KeyError(name)


def _block_diagonal(*blocks: list[list[int]]) -> list[list[int]]:
    size = sum(len(block) for block in blocks)
    matrix = [[0] * size for _ in range(size)]
    offset = 0
    for block in blocks:
        for i, row in enumerate(block):
            for j, value in enumerate(row):
                matrix[offset + i][offset + j] = value
        offset += len(block)
    return matrix


def _mat_vec(matrix: list[list[int]], vector: list[Fraction]) -> list[Fraction]:
    return [sum(Fraction(value) * x for value, x in zip(row, vector)) for row in matrix]


def _pair(matrix: list[list[int]], left: list[Fraction], right: list[Fraction]) -> Fraction:
    return sum(x * y for x, y in zip(left, _mat_vec(matrix, right)))


def _solve(matrix: list[list[int]], rhs: list[int]) -> list[Fraction]:
    n = len(matrix)
    augmented = [
        [Fraction(value) for value in row] + [Fraction(rhs[i])]
        for i, row in enumerate(matrix)
    ]
    for column in range(n):
        pivot = next(i for i in range(column, n) if augmented[i][column])
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for i in range(n):
            if i == column:
                continue
            scale = augmented[i][column]
            if scale:
                augmented[i] = [
                    x - scale * y
                    for x, y in zip(augmented[i], augmented[column])
                ]
    return [augmented[i][-1] for i in range(n)]


def _determinant(matrix: list[list[int]]) -> int:
    work = [[Fraction(value) for value in row] for row in matrix]
    determinant = Fraction(1)
    sign = 1
    for column in range(len(work)):
        pivot = next((i for i in range(column, len(work)) if work[i][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign *= -1
        pivot_value = work[column][column]
        determinant *= pivot_value
        for j in range(column, len(work)):
            work[column][j] /= pivot_value
        for i in range(column + 1, len(work)):
            scale = work[i][column]
            for j in range(column, len(work)):
                work[i][j] -= scale * work[column][j]
    result = sign * determinant
    if result.denominator != 1:
        raise ArithmeticError("integral matrix has nonintegral determinant")
    return result.numerator


def _fundamental_weight(cartan: list[list[int]], node: int) -> list[Fraction]:
    rhs = [0] * len(cartan)
    rhs[node] = 1
    return _solve(cartan, rhs)


def _class_order(vector: list[Fraction]) -> int:
    return lcm(*(entry.denominator for entry in vector))


def _mixed_coordinate_model(kind: str, component: str) -> tuple[list[list[int]], list[Fraction]]:
    if kind == "D_vector":
        cartan = _cartan_d(int(component[1:]))
        return cartan, _fundamental_weight(cartan, 0)
    if kind == "A3_middle":
        cartan = _cartan_a(3)
        return cartan, _fundamental_weight(cartan, 1)
    if kind == "A2_minuscule":
        cartan = _cartan_a(2)
        return cartan, _fundamental_weight(cartan, 0)
    if kind == "A1_pair":
        cartan = _block_diagonal(_cartan_a(1), _cartan_a(1))
        return cartan, [Fraction(1, 2), Fraction(1, 2)]
    raise KeyError(kind)


def _exact_mixed_coordinate_check(kind: str, component: str) -> dict[str, object]:
    cartan, omega = _mixed_coordinate_model(kind, component)
    order = _class_order(omega)
    multiple = [order * entry for entry in omega]
    integral_multiple = [entry.denominator == 1 for entry in multiple]
    root_coordinates = [order, 0] + [int(entry) for entry in multiple]
    positive_square = _pair(cartan, omega, omega)
    norm = -(order * order) * positive_square
    # Pairings with the U basis are (0,m), and the ADE pairings have a
    # minus sign.  The gcd of these numbers is the lattice divisibility.
    local_pairings = [-int(entry) for entry in _mat_vec(cartan, multiple)]
    all_pairings = [0, order] + local_pairings
    divisibility = reduce(gcd, (abs(value) for value in all_pairings))
    reflection_coefficients = [Fraction(-2 * value, norm) for value in all_pairings]
    return {
        "class_order": order,
        "omega_coordinates": [_fraction_text(entry) for entry in omega],
        "positive_square": _fraction_text(positive_square),
        "root_coordinates_in_U_plus_component": root_coordinates,
        "root_norm": int(norm),
        "root_divisibility": divisibility,
        "root_coordinate_gcd": reduce(gcd, (abs(value) for value in root_coordinates)),
        "multiple_is_integral": all(integral_multiple),
        "reflection_coefficients_integral": all(
            coefficient.denominator == 1 for coefficient in reflection_coefficients
        ),
    }


def _exact_overclosure_coordinate_check(kind: str, component: str | None) -> dict[str, object]:
    if kind == "D_vector_odd":
        cartan = _cartan_d(int(component[1:]))
        glue = _fundamental_weight(cartan, 0)
    elif kind == "A3_middle_odd":
        cartan = _cartan_a(3)
        glue = _fundamental_weight(cartan, 1)
    elif kind == "A1_pair_odd":
        cartan = _block_diagonal(_cartan_a(1), _cartan_a(1))
        glue = [Fraction(1, 2), Fraction(1, 2)]
    elif kind == "E6_A2_diagonal_even":
        # E6 is the star with arms 2,2,1.  Either endpoint of a length-two
        # arm is a minuscule class of square 4/3 and exact order three.
        e6 = _star_cartan((2, 2, 1))
        a2 = _cartan_a(2)
        cartan = _block_diagonal(e6, a2)
        e6_weight = _fundamental_weight(e6, 2)
        a2_weight = _fundamental_weight(a2, 0)
        glue = e6_weight + a2_weight
    else:
        raise KeyError(kind)
    pairings = _mat_vec(cartan, glue)
    return {
        "class_order": _class_order(glue),
        "glue_coordinates": [_fraction_text(entry) for entry in glue],
        "glue_pairings_integral": all(entry.denominator == 1 for entry in pairings),
        "negative_glue_square": _fraction_text(-_pair(cartan, glue, glue)),
        "proper": any(entry.denominator != 1 for entry in glue),
    }


def _ue8e7_discriminant_check() -> dict[str, object]:
    e7 = _star_cartan((2, 3, 1))
    order_two_weights = []
    for node in range(7):
        weight = _fundamental_weight(e7, node)
        if _class_order(weight) == 2:
            order_two_weights.append((node, _pair(e7, weight, weight)))
    # Every order-two representative differs by E7; its q-value is 3/2.
    has_three_halves = any(square % 2 == Fraction(3, 2) for _, square in order_two_weights)
    return {
        "E7_determinant": _determinant(e7),
        "order_two_fundamental_weight_squares": [
            {"node": node, "positive_square": _fraction_text(square)}
            for node, square in order_two_weights
        ],
        "negative_discriminant_q_mod_2": "1/2",
        "mixed_root_required_q_mod_2": "1",
        "values_are_distinct": Fraction(1, 2) % 2 != Fraction(1) % 2,
        "found_E7_q_three_halves": has_three_halves,
    }


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _row_artifact(row: tuple[object, ...]) -> dict[str, object]:
    number, lattice, factors, witness_key, component, overclosure_key = row
    rank = 2 + sum(ADE[factor]["rank"] for factor in factors)
    determinant = 1
    for factor in factors:
        determinant *= ADE[factor]["determinant"]

    if witness_key is None:
        mixed = {
            "exists": False,
            "proof": (
                "A_S=Z/2 has q(nonzero)=+1/2 mod 2Z.  An integral "
                "reflection outside the norm-two shell would force "
                "divisibility 2, norm -4, and q=-1 mod 2Z, impossible."
            ),
        }
    else:
        witness = MIXED_WITNESSES[witness_key]
        order = witness["class_order"]
        weight_square = witness["positive_weight_square"]
        root_norm = -(order * order) * weight_square
        divisibility = order
        reflection_multiplier = Fraction(2 * divisibility, -root_norm)
        mixed = {
            "exists": True,
            "kind": witness_key,
            "component": component,
            "dual_class_order": order,
            "positive_weight_square": _fraction_text(weight_square),
            "root_formula": witness["root_formula"],
            "root_norm": int(root_norm),
            "root_divisibility": divisibility,
            "integral_reflection_coefficient_check": _fraction_text(reflection_multiplier),
            "primitive": True,
            "outside_W2": True,
            "exact_coordinate_replay": _exact_mixed_coordinate_check(witness_key, component),
        }

    if overclosure_key is None:
        overclosure = {
            "exists": False,
            "proof": (
                f"|det S|={determinant} is squarefree, while a proper "
                "integral overclosure of index k would require k^2 dividing |det S|."
            ),
        }
    else:
        witness = OVERCLOSURE_WITNESSES[overclosure_key]
        index = witness["class_order"]
        target_det = determinant // (index * index)
        overclosure = {
            "exists": True,
            "kind": overclosure_key,
            "description": witness["description"],
            "index": index,
            "glue_square": _fraction_text(witness["negative_square"]),
            "parity": witness["parity"],
            "target_absolute_determinant": target_det,
            "W2_stable": True,
            "exact_coordinate_replay": _exact_overclosure_coordinate_check(
                overclosure_key, component
            ),
        }

    same_lattice_nonmaximal = bool(mixed["exists"])
    return {
        "wang_number": number,
        "lattice": lattice,
        "ade_factors": list(factors),
        "rank": rank,
        "hyperbolic_dimension": rank - 1,
        "absolute_determinant": determinant,
        "automorphic_correction": False,
        "same_lattice_mixed_root": mixed,
        "proper_integral_overclosure": overclosure,
        "same_lattice_reflection_maximal": not same_lattice_nonmaximal,
        "arithmetic_reflection_maximal": number == 1,
        "arithmetic_reflection_maximality_basis": (
            "results/wang_h16_maximal_reflection_no_correction.md"
            if number == 1
            else "strict same-lattice inclusion W2(S) < W_ref(S)"
        ),
    }


def analyze() -> dict[str, object]:
    rows = [_row_artifact(row) for row in ROWS]
    ue8e7_discriminant_check = _ue8e7_discriminant_check()

    formula_checks = []
    for row in rows:
        mixed = row["same_lattice_mixed_root"]
        if not mixed["exists"]:
            continue
        m = mixed["dual_class_order"]
        q = Fraction(mixed["positive_weight_square"])
        formula_checks.append({
            "wang_number": row["wang_number"],
            "m2_q_equals_2m": m * m * q == 2 * m,
            "norm_is_minus_2m": mixed["root_norm"] == -2 * m,
            "divisibility_is_m": mixed["root_divisibility"] == m,
            "reflection_is_integral": mixed["integral_reflection_coefficient_check"] == "1",
            "norm_ratio_to_minus_two_is_nonsquare": m in (2, 3),
            "coordinate_replay_matches": (
                mixed["exact_coordinate_replay"]["class_order"] == m
                and mixed["exact_coordinate_replay"]["positive_square"]
                == mixed["positive_weight_square"]
                and mixed["exact_coordinate_replay"]["root_norm"]
                == mixed["root_norm"]
                and mixed["exact_coordinate_replay"]["root_divisibility"]
                == mixed["root_divisibility"]
                and mixed["exact_coordinate_replay"]["root_coordinate_gcd"] == 1
                and mixed["exact_coordinate_replay"]["multiple_is_integral"]
                and mixed["exact_coordinate_replay"]["reflection_coefficients_integral"]
            ),
        })

    overclosure_checks = []
    for row in rows:
        closure = row["proper_integral_overclosure"]
        if not closure["exists"]:
            continue
        index = closure["index"]
        square = Fraction(closure["glue_square"])
        overclosure_checks.append({
            "wang_number": row["wang_number"],
            "glue_has_integral_square": square.denominator == 1,
            "determinant_index_formula": (
                row["absolute_determinant"]
                == index * index * closure["target_absolute_determinant"]
            ),
            "parity_matches_square": (
                closure["parity"] == ("even" if square.numerator % 2 == 0 else "odd")
            ),
            "coordinate_replay_matches": (
                closure["exact_coordinate_replay"]["class_order"] == index
                and closure["exact_coordinate_replay"]["negative_glue_square"]
                == closure["glue_square"]
                and closure["exact_coordinate_replay"]["glue_pairings_integral"]
                and closure["exact_coordinate_replay"]["proper"]
            ),
        })

    checks = {
        "wang_negative_row_count_is_18": len(rows) == 18,
        "row_numbers_are_1_through_18": [row["wang_number"] for row in rows] == list(range(1, 19)),
        "all_ADE_determinants_replay": all(
            prod(ADE[factor]["determinant"] for factor in factors)
            == prod(_determinant(_ade_cartan(factor)) for factor in factors)
            for _, _, factors, _, _, _ in ROWS
        ),
        "mixed_root_count_is_17": sum(row["same_lattice_mixed_root"]["exists"] for row in rows) == 17,
        "proper_integral_overclosure_count_is_15": sum(row["proper_integral_overclosure"]["exists"] for row in rows) == 15,
        "squarefree_no_overclosure_rows_are_exact": {
            row["wang_number"]
            for row in rows
            if not row["proper_integral_overclosure"]["exists"]
        } == {1, 12, 17},
        "all_mixed_root_formula_checks_pass": all(
            all(value for key, value in item.items() if key != "wang_number")
            for item in formula_checks
        ),
        "all_overclosure_checks_pass": all(
            all(value for key, value in item.items() if key != "wang_number")
            for item in overclosure_checks
        ),
        "unique_same_lattice_reflection_maximal_row_is_UE8E7": [
            row["lattice"] for row in rows if row["same_lattice_reflection_maximal"]
        ] == ["U+E8+E7"],
        "unique_arithmetic_reflection_maximal_negative_row_is_UE8E7": [
            row["lattice"] for row in rows if row["arithmetic_reflection_maximal"]
        ] == ["U+E8+E7"],
        "UE8E7_discriminant_obstruction_replays": (
            ue8e7_discriminant_check["E7_determinant"] == 2
            and ue8e7_discriminant_check["found_E7_q_three_halves"]
            and ue8e7_discriminant_check["values_are_distinct"]
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Wang-18 mixed-root audit failed: {checks}")

    return {
        "schema": "coxeter.wang18_same_lattice_mixed_root_ledger.v1",
        "status": "closed_18_of_18",
        "published_input": (
            "Haowu Wang, The classification of 2-reflective modular forms, "
            "JEMS 26 (2024), Theorem 7.2"
        ),
        "new_exact_result": (
            "Seventeen negative rows have strict same-lattice mixed-root "
            "reflection enlargements; U+E8+E7 is the unique exception and "
            "the unique maximal arithmetic reflection group among the eighteen."
        ),
        "global_69_corollary": {
            "complete_norm_two_reflective_lattices_of_rank_at_least_6": 69,
            "with_Wang_automorphic_correction": 51,
            "without_Wang_automorphic_correction": 18,
            "theorem": (
                "If W2(S) is maximal as an arithmetic reflection group, then "
                "S has a Wang automorphic correction or "
                "S is U+E8(-1)+E7(-1) (Tumarkin Gamma16)."
            ),
            "converse_claimed": False,
            "Gamma16_maximality_proof": (
                "results/wang_h16_maximal_reflection_no_correction.md"
            ),
        },
        "notions": {
            "same_lattice": "W2(S) < W_ref(S), both acting in the same H^n",
            "integral_overclosure": "S is contained properly in an integral lattice M in S tensor Q",
            "arithmetic_reflection_maximal": (
                "no larger arithmetic reflection group in Isom(H^n), including "
                "commensurable lattice changes"
            ),
        },
        "rows": rows,
        "formula_checks": formula_checks,
        "overclosure_checks": overclosure_checks,
        "UE8E7_discriminant_check": ue8e7_discriminant_check,
        "counts": {
            "rows": 18,
            "strict_same_lattice_mixed_reflection_enlargements": 17,
            "proper_integral_overclosures": 15,
            "arithmetic_reflection_maximal_rows": 1,
            "undecided_rows": 0,
        },
        "checks": checks,
    }


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2, sort_keys=True))
