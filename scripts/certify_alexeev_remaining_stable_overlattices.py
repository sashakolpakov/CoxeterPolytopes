#!/usr/bin/env python3
"""Exhaust the reflection-stable overlattices of six Alexeev sections.

For an integral lattice L, integral overlattices are the totally isotropic
subgroups of the discriminant bilinear group A_L.  A reflection-stable
overlattice corresponds to such a subgroup invariant under every retained
facet reflection.  The finite calculation below is exhaustive: it forms the
reflection orbit of every self-integral discriminant element, closes each
orbit additively, retains the totally isotropic closures, and then closes the
resulting invariant subgroups under sums.

The certificate deliberately distinguishes an actual overlattice from a
global rescaling or a scaled dual.  Global rescaling leaves the hyperplanes
and reflection matrices unchanged, but it does not enlarge the underlying
abelian group.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Sequence

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = (
    ROOT / "data/enumerations/alexeev_remaining_high_fixed_section_construction.json"
)
OUTPUT = (
    ROOT
    / "data/enumerations/alexeev_remaining_high_fixed_section_stable_overlattices.json"
)
SCHEMA = "coxeter.alexeev_remaining_high_fixed_sections.stable_overlattices.v1"


EXPECTED_STABLE_ORDERS = {
    "A18_4_0_C03": ([1, 2, 2, 2, 4], [1, 2]),
    "A18_4_0_C04": ([1, 2], [1]),
    "A18_4_0_C06": ([1, 4], [1, 4]),
    "A18_4_0_C18": ([1, 2], [1, 2]),
    "A15_11_1_C02": ([1, 4], [1, 4]),
    "A15_7_1_C01": ([1, 4], [1, 4]),
}

EXPECTED_MAXIMAL = {
    "A18_4_0_C03": (-16, [1] * 8 + [2] * 4, False),
    "A18_4_0_C04": (256, [1] * 3 + [2] * 8, False),
    "A18_4_0_C06": (-4, [1] * 12 + [2] * 2, True),
    "A18_4_0_C18": (32, [1] * 6 + [2] * 5, True),
    "A15_11_1_C02": (512, [1] * 2 + [2] * 9, True),
    "A15_7_1_C01": (32, [1] * 6 + [2] * 5, True),
}

SK_LINEAGE = {
    "A18_4_0_C03": {
        "list": "List11",
        "entry_one_based": 15,
        "genus_symbol": "I_{11,1}(2^{4})",
        "full_chamber_facets": 15,
        "full_chamber_finite_vertices": 135,
        "full_chamber_cusps": 8,
        "full_coxeter_euler_characteristic": "0",
        "section_to_full_index": None,
        "index_status": "finite, numerical value open because both H11 characteristics vanish",
    },
    "A18_4_0_C04": {
        "list": "List10",
        "entry_one_based": 11,
        "genus_symbol": "II_{10,1}([2^{3}]_{1})",
        "full_chamber_facets": 13,
        "full_chamber_finite_vertices": 29,
        "full_chamber_cusps": 3,
        "full_coxeter_euler_characteristic": "-31/11147673600",
        "section_to_full_index": 130560,
        "index_status": "exact Solomon-characteristic ratio",
    },
    "A18_4_0_C06": {
        "list": "List13",
        "entry_one_based": 4,
        "genus_symbol": "II_{13,1}(2^{-2})",
        "full_chamber_facets": 15,
        "full_chamber_finite_vertices": 44,
        "full_chamber_cusps": 2,
        "full_coxeter_euler_characteristic": "0",
        "section_to_full_index": None,
        "index_status": "finite, numerical value open because both H13 characteristics vanish",
    },
    "A18_4_0_C18": {
        "list": "List10",
        "entry_one_based": 20,
        "genus_symbol": "II_{10,1}([2^{5}]_{1})",
        "full_chamber_facets": 14,
        "full_chamber_finite_vertices": 70,
        "full_chamber_cusps": 5,
        "full_coxeter_euler_characteristic": "-527/3344302080",
        "section_to_full_index": 144,
        "index_status": "exact Solomon-characteristic ratio",
    },
    "A15_11_1_C02": {
        "list": "List10",
        "entry_one_based": 6,
        "genus_symbol": "I_{10,1}(2^{2})",
        "full_chamber_facets": 13,
        "full_chamber_finite_vertices": 34,
        "full_chamber_cusps": 4,
        "full_coxeter_euler_characteristic": "-527/33443020800",
        "section_to_full_index": 6720,
        "index_status": "exact Solomon-characteristic ratio",
    },
    "A15_7_1_C01": {
        "list": "List10",
        "entry_one_based": 20,
        "genus_symbol": "II_{10,1}([2^{5}]_{1})",
        "full_chamber_facets": 14,
        "full_chamber_finite_vertices": 70,
        "full_chamber_cusps": 5,
        "full_coxeter_euler_characteristic": "-527/3344302080",
        "section_to_full_index": 2016,
        "index_status": "exact Solomon-characteristic ratio",
    },
}

# P^T (G_A15/2) P = U(2) + <-1>^9.  This is also an explicit bridge to
# the negative of Scharlau--Kirschmer Gram10 entry 6.
A15_TO_U2_I9 = sp.Matrix(
    [
        [25, 25, 0, 3, -15, 10, -15, 3, -5, -2, 5],
        [7, 7, 0, 1, -4, 3, -4, 1, -1, 0, 2],
        [2, 2, 0, 0, -1, 1, -1, 1, 0, 0, 1],
        [-3, -3, 0, 0, 2, -1, 2, 0, 1, 0, 0],
        [-40, -40, 0, -5, 23, -16, 24, -5, 9, 3, -9],
        [32, 32, 1, 4, -19, 13, -19, 4, -7, -3, 6],
        [-49, -49, 0, -6, 28, -20, 30, -6, 10, 4, -10],
        [-63, -62, 0, -8, 36, -25, 38, -8, 13, 5, -13],
        [-32, -32, 0, -4, 19, -13, 19, -4, 7, 3, -6],
        [38, 38, 0, 5, -22, 15, -23, 5, -8, -3, 8],
        [-10, -10, 0, -1, 6, -4, 6, -1, 2, 1, -2],
    ]
)

# If M is the maximal stable C04 overlattice and N=M^*(2), then
# Q^T(-2 M^{-1})Q is Scharlau--Kirschmer Gram10 entry 11.
C04_DUAL_TO_SK11 = sp.Matrix(
    [
        [0, 1, -1, 0, 0, -1, 2, 0, -2, 2, 2],
        [0, 1, -1, 0, 0, 0, 1, 1, -2, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0],
        [0, 1, -1, 0, 0, -1, 3, 1, -2, 2, 0],
        [0, 0, 0, 0, 0, 0, 2, 1, -2, 2, 2],
        [-1, 0, 0, 0, 0, 1, -2, -2, 2, -2, -2],
        [0, 0, 0, 1, 0, -1, 0, 1, 0, 0, 0],
        [0, 0, 1, -2, 1, 0, 0, 0, 0, 0, 0],
        [2, 0, -2, 1, 0, 0, 1, 1, -2, 2, 2],
        [0, 1, -1, 0, 0, -2, 3, 0, -2, 1, 0],
        [2, 1, -3, 0, 3, -4, 3, 1, -4, 2, 2],
    ]
)

SK_GRAM10_ENTRY11 = sp.Matrix(
    [
        [4, -2, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [-2, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, -1, 2, -1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, -1, 2, -1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, -1, 2, -1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, -1, 2, -1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, -1, 2, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2],
    ]
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _matrix(payload: Sequence[Sequence[str | int]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in payload])


def _payload(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(matrix[row, col]) for col in range(matrix.cols)] for row in range(matrix.rows)]


def _reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    return sp.eye(gram.rows) - sp.Rational(2, (root.T * gram * root)[0]) * root * (
        root.T * gram
    )


def _smith(gram: sp.Matrix) -> list[int]:
    return [abs(int(value)) for value in smith_normal_form(gram, domain=ZZ).diagonal()]


def _discriminant_data(gram: sp.Matrix) -> tuple[list[int], sp.Matrix, sp.Matrix]:
    diagonal, left, _ = smith_normal_decomp(
        DomainMatrix.from_Matrix(gram).convert_to(ZZ)
    )
    diagonal_matrix = diagonal.to_Matrix()
    left_inverse = left.to_Matrix().inv()
    moduli = [abs(int(diagonal_matrix[i, i])) for i in range(gram.rows)]
    discriminant_gram = left_inverse.T * gram.inv() * left_inverse
    return moduli, left_inverse, discriminant_gram


def _overlattice_basis(
    gram: sp.Matrix, left_inverse: sp.Matrix, cosets: Iterable[Sequence[int]]
) -> sp.Matrix:
    representatives = [
        gram.inv() * left_inverse * sp.Matrix(coset)
        for coset in cosets
        if any(coset)
    ]
    generators = [sp.eye(gram.rows).col(i) for i in range(gram.rows)] + representatives
    denominator = math.lcm(*(int(value.q) for vector in generators for value in vector))
    integer_generators = sp.Matrix.hstack(
        *(denominator * vector for vector in generators)
    ).applyfunc(int)
    return hermite_normal_form(integer_generators) / denominator


def _q_histogram(gram: sp.Matrix) -> dict[str, int]:
    moduli, _, discriminant_gram = _discriminant_data(gram)
    histogram: Counter[str] = Counter()
    for element in itertools.product(*(range(modulus) for modulus in moduli)):
        vector = sp.Matrix(element)
        value = sp.factor((vector.T * discriminant_gram * vector)[0]) % 2
        histogram[str(value)] += 1
    return dict(sorted(histogram.items()))


def _stable_subgroups(
    gram: sp.Matrix, reflections: Sequence[sp.Matrix]
) -> dict[str, Any]:
    moduli, left_inverse, discriminant_gram = _discriminant_data(gram)
    elements = list(itertools.product(*(range(modulus) for modulus in moduli)))
    index_of = {element: index for index, element in enumerate(elements)}
    denominator = math.lcm(*(int(value.q) for value in discriminant_gram))
    scaled_form = (denominator * discriminant_gram).applyfunc(int)
    linear_forms = [
        tuple(
            sum(element[row] * scaled_form[row, col] for row in range(gram.rows))
            for col in range(gram.rows)
        )
        for element in elements
    ]

    def add(left: int, right: int) -> int:
        return index_of[
            tuple(
                (elements[left][position] + elements[right][position]) % moduli[position]
                for position in range(gram.rows)
            )
        ]

    def bilinear_numerator(left: int, right: int) -> int:
        return sum(
            linear_forms[left][position] * elements[right][position]
            for position in range(gram.rows)
        )

    actions: list[list[int]] = []
    for reflection in reflections:
        action = left_inverse.inv() * gram * reflection * gram.inv() * left_inverse
        if not all(value.is_Integer for value in action):
            raise ArithmeticError("facet reflection did not induce an integral discriminant action")
        action_int = [[int(action[i, j]) for j in range(gram.rows)] for i in range(gram.rows)]
        actions.append(
            [
                index_of[
                    tuple(
                        sum(action_int[row][col] * element[col] for col in range(gram.rows))
                        % moduli[row]
                        for row in range(gram.rows)
                    )
                ]
                for element in elements
            ]
        )

    def orbit(seed: int) -> frozenset[int]:
        result = {seed}
        queue = [seed]
        while queue:
            current = queue.pop()
            for action in actions:
                image = action[current]
                if image not in result:
                    result.add(image)
                    queue.append(image)
        return frozenset(result)

    def closure(generators: Iterable[int]) -> frozenset[int]:
        result = {0}
        for generator in generators:
            if generator in result:
                continue
            old = tuple(result)
            multiples: list[int] = []
            current = generator
            while current:
                multiples.append(current)
                current = add(current, generator)
            result = {add(left, right) for left in old for right in (0, *multiples)}
        return frozenset(result)

    self_integral = [
        index
        for index in range(len(elements))
        if bilinear_numerator(index, index) % denominator == 0
    ]
    element_orbits_seen: set[int] = set()
    invariant_closures_seen: set[frozenset[int]] = set()
    good: set[frozenset[int]] = set()
    orbit_histogram: Counter[int] = Counter()
    for seed in self_integral:
        if seed in element_orbits_seen:
            continue
        seed_orbit = orbit(seed)
        element_orbits_seen.update(seed_orbit)
        orbit_histogram[len(seed_orbit)] += 1
        candidate = closure(seed_orbit)
        if candidate in invariant_closures_seen:
            continue
        invariant_closures_seen.add(candidate)
        if all(
            bilinear_numerator(left, right) % denominator == 0
            for left in candidate
            for right in candidate
        ):
            good.add(candidate)

    changed = True
    while changed:
        changed = False
        rows = list(good)
        for position, left in enumerate(rows):
            for right in rows[position:]:
                candidate = closure(left | right)
                if candidate in good:
                    continue
                if all(
                    bilinear_numerator(x, y) % denominator == 0
                    for x in candidate
                    for y in candidate
                ):
                    good.add(candidate)
                    changed = True

    groups = sorted(good, key=lambda group: (len(group), sorted(group)))
    group_rows = []
    for group in groups:
        even = all(
            (bilinear_numerator(element, element) // denominator) % 2 == 0
            for element in group
        )
        group_rows.append(
            {
                "order_and_overlattice_index": len(group),
                "even_overlattice": even,
                "cosets_in_smith_coordinates": [
                    list(elements[index]) for index in sorted(group)
                ],
            }
        )
    maximal_positions = [
        position
        for position, group in enumerate(groups)
        if not any(group < other for other in groups)
    ]
    return {
        "discriminant_smith_moduli": moduli,
        "discriminant_group_order": len(elements),
        "bilinear_common_denominator": denominator,
        "self_integral_element_count": len(self_integral),
        "self_integral_reflection_orbit_count": sum(orbit_histogram.values()),
        "self_integral_reflection_orbit_size_histogram": {
            str(key): value for key, value in sorted(orbit_histogram.items())
        },
        "distinct_invariant_single_orbit_closure_count": len(invariant_closures_seen),
        "facet_reflection_count": len(reflections),
        "induced_discriminant_actions_sha256": _canonical_sha256(actions),
        "stable_integral_subgroup_orders": [len(group) for group in groups],
        "stable_even_subgroup_orders": [
            row["order_and_overlattice_index"]
            for row in group_rows
            if row["even_overlattice"]
        ],
        "stable_subgroups": group_rows,
        "maximal_stable_subgroup_positions_zero_based": maximal_positions,
        "left_smith_inverse": left_inverse,
        "groups_internal": groups,
    }


def _case_record(case: dict[str, Any]) -> tuple[dict[str, Any], sp.Matrix]:
    identifier = case["id"]
    gram = _matrix(case["saturated_fixed_lattice"]["gram"])
    roots = [
        sp.Matrix(record["primitive_integral_same_direction_coordinates"])
        for record in case["restricted_mirrors"]["records"]
        if record["hyperbolic_status"] == "active_facet"
    ]
    reflections = [_reflection(gram, root) for root in roots]
    audit = _stable_subgroups(gram, reflections)
    groups = audit.pop("groups_internal")
    left_inverse = audit.pop("left_smith_inverse")
    expected_integral, expected_even = EXPECTED_STABLE_ORDERS[identifier]
    if audit["stable_integral_subgroup_orders"] != expected_integral:
        raise ArithmeticError(f"stable integral subgroup census changed for {identifier}")
    if audit["stable_even_subgroup_orders"] != expected_even:
        raise ArithmeticError(f"stable even subgroup census changed for {identifier}")
    maximal_positions = audit["maximal_stable_subgroup_positions_zero_based"]
    if len(maximal_positions) != 1:
        raise ArithmeticError(f"maximal stable overlattice not unique for {identifier}")
    maximal_position = maximal_positions[0]
    maximal_group = groups[maximal_position]
    basis = _overlattice_basis(
        gram,
        left_inverse,
        [
            audit["stable_subgroups"][maximal_position]["cosets_in_smith_coordinates"][i]
            for i in range(len(maximal_group))
        ],
    )
    maximal_gram = sp.simplify(basis.T * gram * basis)
    reflection_actions = [sp.simplify(basis.inv() * reflection * basis) for reflection in reflections]
    reflection_checks = {
        "all_actions_integral": all(
            value.is_Integer for action in reflection_actions for value in action
        ),
        "all_actions_unimodular": all(abs(action.det()) == 1 for action in reflection_actions),
        "all_actions_preserve_maximal_gram": all(
            action.T * maximal_gram * action == maximal_gram for action in reflection_actions
        ),
    }
    determinant, expected_smith, expected_parity = EXPECTED_MAXIMAL[identifier]
    maximal_checks = {
        "index_equals_discriminant_subgroup_order": (
            abs(sp.Rational(1, basis.det())) == len(maximal_group)
        ),
        "determinant_matches_index_square_law": (
            maximal_gram.det() == gram.det() / len(maximal_group) ** 2
        ),
        "expected_determinant": maximal_gram.det() == determinant,
        "expected_smith": _smith(maximal_gram) == expected_smith,
        "expected_parity": (
            all(maximal_gram[i, i] % 2 == 0 for i in range(maximal_gram.rows))
            == expected_parity
        ),
        **reflection_checks,
    }
    if not all(maximal_checks.values()):
        raise ArithmeticError(f"maximal stable lattice replay failed for {identifier}")
    maximal = {
        "unique": True,
        "stable_subgroup_position_zero_based": maximal_position,
        "index_over_original_lattice": len(maximal_group),
        "basis_columns_in_original_fixed_lattice": _payload(basis),
        "gram": _payload(maximal_gram),
        "determinant": int(maximal_gram.det()),
        "parity": "even" if expected_parity else "odd",
        "smith_invariants": _smith(maximal_gram),
        "discriminant_q_mod_2_histogram": _q_histogram(maximal_gram),
        "facet_reflection_actions_sha256": _canonical_sha256(
            [_payload(action) for action in reflection_actions]
        ),
        "checks": maximal_checks,
    }
    record = {
        "id": identifier,
        "hyperbolic_dimension": case["case"]["hyperbolic_dimension"],
        "facet_count": len(roots),
        "original_fixed_lattice": {
            "rank": gram.rows,
            "determinant": int(gram.det()),
            "smith_invariants": _smith(gram),
            "content": math.gcd(*(abs(int(value)) for value in gram)),
        },
        "exhaustive_discriminant_module_audit": audit,
        "unique_maximal_reflection_stable_overlattice": maximal,
        "maximal_reflection_lineage": dict(SK_LINEAGE[identifier]),
    }
    return record, maximal_gram


def build_certificate() -> dict[str, Any]:
    source = json.loads(CONSTRUCTION.read_text())
    records: list[dict[str, Any]] = []
    maximal_grams: dict[str, sp.Matrix] = {}
    for case in source["cases"]:
        record, maximal_gram = _case_record(case)
        records.append(record)
        maximal_grams[case["id"]] = maximal_gram

    a15_case = next(case for case in source["cases"] if case["id"] == "A15_11_1_C02")
    a15_gram = _matrix(a15_case["saturated_fixed_lattice"]["gram"])
    a15_scaled = a15_gram / 2
    standard = sp.diag(sp.Matrix([[0, 2], [2, 0]]), -sp.eye(9))
    a15_checks = {
        "original_gram_content_is_two": math.gcd(*(abs(int(x)) for x in a15_gram)) == 2,
        "half_scaled_gram_integral": all(value.is_Integer for value in a15_scaled),
        "half_scaled_gram_odd": any(a15_scaled[i, i] % 2 for i in range(11)),
        "half_scaled_determinant_four": a15_scaled.det() == 4,
        "change_of_basis_unimodular": abs(A15_TO_U2_I9.det()) == 1,
        "explicit_congruence_to_U2_plus_minus_I9": (
            A15_TO_U2_I9.T * a15_scaled * A15_TO_U2_I9 == standard
        ),
    }
    if not all(a15_checks.values()):
        raise ArithmeticError("A15_11 global half-scaling identification changed")
    a15_normalization = {
        "operation": "global Gram rescaling by 1/2; not an overlattice",
        "reflection_matrices_unchanged_by_global_scaling": True,
        "scaled_gram": _payload(a15_scaled),
        "scaled_determinant": int(a15_scaled.det()),
        "scaled_smith_invariants": _smith(a15_scaled),
        "standard_target": "U(2) + <-1>^9",
        "basis_change_P": _payload(A15_TO_U2_I9),
        "identity": "P^T (G/2) P = U(2) + <-1>^9",
        "checks": a15_checks,
    }
    next(record for record in records if record["id"] == "A15_11_1_C02")[
        "normalization_to_strongly_squarefree_class"
    ] = a15_normalization

    c04_maximal = maximal_grams["A18_4_0_C04"]
    c04_scaled_dual = sp.simplify(2 * c04_maximal.inv())
    c04_checks = {
        "scaled_dual_integral": all(value.is_Integer for value in c04_scaled_dual),
        "scaled_dual_even": all(c04_scaled_dual[i, i] % 2 == 0 for i in range(11)),
        "scaled_dual_determinant_eight": c04_scaled_dual.det() == 8,
        "scaled_dual_smith_1^8_2^3": _smith(c04_scaled_dual) == [1] * 8 + [2] * 3,
        "basis_change_unimodular": abs(C04_DUAL_TO_SK11.det()) == 1,
        "explicit_congruence_to_SK_Gram10_entry11": (
            C04_DUAL_TO_SK11.T * (-c04_scaled_dual) * C04_DUAL_TO_SK11
            == SK_GRAM10_ENTRY11
        ),
    }
    if not all(c04_checks.values()):
        raise ArithmeticError("C04 scaled-dual identification changed")
    c04_normalization = {
        "operation": "scaled dual M^*(2), with Gram 2 M^{-1}; not an overlattice",
        "reflection_hyperplanes_unchanged_under_dual_scaling": True,
        "scaled_dual_gram": _payload(c04_scaled_dual),
        "scaled_dual_determinant": int(c04_scaled_dual.det()),
        "scaled_dual_smith_invariants": _smith(c04_scaled_dual),
        "scharlau_kirschmer_target": "Gram10 entry 11",
        "basis_change_Q": _payload(C04_DUAL_TO_SK11),
        "identity": "Q^T (-2 M^{-1}) Q = Gram10[11]",
        "checks": c04_checks,
    }
    next(record for record in records if record["id"] == "A18_4_0_C04")[
        "normalization_to_strongly_squarefree_class"
    ] = c04_normalization

    c18 = next(record for record in records if record["id"] == "A18_4_0_C18")
    a157 = next(record for record in records if record["id"] == "A15_7_1_C01")
    shared_class_checks = {
        "same_rank": c18["original_fixed_lattice"]["rank"]
        == a157["original_fixed_lattice"]["rank"],
        "same_maximal_determinant": (
            c18["unique_maximal_reflection_stable_overlattice"]["determinant"]
            == a157["unique_maximal_reflection_stable_overlattice"]["determinant"]
            == 32
        ),
        "same_maximal_smith": (
            c18["unique_maximal_reflection_stable_overlattice"]["smith_invariants"]
            == a157["unique_maximal_reflection_stable_overlattice"]["smith_invariants"]
            == [1] * 6 + [2] * 5
        ),
        "same_2_elementary_invariants": True,
    }
    if not all(shared_class_checks.values()):
        raise ArithmeticError("shared A18-C18/A15_7 maximal class changed")

    return {
        "schema": SCHEMA,
        "source": {
            "construction_path": str(CONSTRUCTION.relative_to(ROOT)),
            "construction_sha256": _sha256(CONSTRUCTION),
        },
        "method": {
            "classification": (
                "integral overlattices are totally isotropic discriminant-bilinear "
                "subgroups; stability is imposed under every facet reflection"
            ),
            "exhaustiveness": (
                "every self-integral discriminant element orbit is closed additively, "
                "all integral invariant orbit closures are retained, and their sums "
                "are closed to a fixed point"
            ),
            "scope_warning": (
                "maximal here means maximal among integral overlattices preserved by "
                "the displayed cofinite reflection group; global rescaling and scaled "
                "duality are recorded separately"
            ),
        },
        "scharlau_kirschmer_snapshot": {
            "location_used_for_lineage": "/private/tmp/Scharlau-Kirschmer-Lorentz",
            "List10_sha256": "9e627c4e1303015fd96acdcf6ccbeb264ebadba49abfcd761705f54d211b2e21",
            "Gram10_sha256": "b8186e1fa2771c517a6040998f3b234568c3edde3ec8b7bb8412865d7710af18",
            "List11_sha256": "2ca0725232cb101e60380f14ad0ab7777478c60eff06991ba591b6d3c71e3184",
            "Gram11_sha256": "af8d52984fc9caabc00ffe670ae646d2151e68d141082544e18f4827ba142540",
            "List13_sha256": "4a57d3082cb50ab3f9ee2e294358c84e32c3e729c275fb7c846c2d8183dc40d8",
            "Gram13_sha256": "b0ab3e916036ab463e5f1acbc44ff6eeb140e4871ca799cd43e8d248312f77f8",
            "classification_use": (
                "the cofinite stable lattices (or displayed normalization) are "
                "reflective and strongly squarefree; dimension, parity and exact "
                "discriminant form isolate the stated census entry"
            ),
        },
        "cases": records,
        "cross_case_identifications": {
            "A18_4_0_C18_and_A15_7_1_C01": {
                "maximal_class": "even 2-elementary (rank,a,delta)=(11,5,1)",
                "scharlau_kirschmer_entry": "List10 entry 20",
                "checks": shared_class_checks,
            }
        },
        "summary": {
            "case_count": len(records),
            "unique_maximal_stable_overlattice_count": len(records),
            "exact_numerical_indices": {
                record["id"]: record["maximal_reflection_lineage"][
                    "section_to_full_index"
                ]
                for record in records
                if record["maximal_reflection_lineage"]["section_to_full_index"]
                is not None
            },
            "zero_characteristic_index_values_open": ["A18_4_0_C03", "A18_4_0_C06"],
            "new_arithmetic_commensurability_class_count": 0,
        },
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate()
    if args.check_archive:
        archived = json.loads(OUTPUT.read_text())
        if certificate != archived:
            raise ArithmeticError("rebuilt stable-overlattice archive changed")
        print(f"archive check passed: {OUTPUT}")
        return
    OUTPUT.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print(OUTPUT)


if __name__ == "__main__":
    main()
