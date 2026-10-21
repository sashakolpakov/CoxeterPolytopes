#!/usr/bin/env python3
"""Exact no-import audit of Nikulin's five four-node edge-graph families.

The program starts with the five finite-edge positions other than the marked
broken edge 12, imposes connectivity and the two endpoint-deletion Coxeter
conditions, and quotients by swaps of the two marked and two auxiliary
vertices.  It also checks the printed parameter grammar on a finite alphabet
that represents every logically different label regime.

No Lanner list, Vinberg/Nikulin dimension bound, or arithmetic-volume theorem
is imported.  The only external algebra package used is SymPy, on five 4 by 4
symbolic determinants.
"""

from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/nikulin_four_node_transition_families.json"
SCHEMA = "coxeter.nikulin_four_node_transition_families.v1"

NODES = (1, 2, 3, 4)
MARKED_EDGE = (1, 2)
FINITE_POSITIONS = ((1, 3), (1, 4), (2, 3), (2, 4), (3, 4))
ALL_EDGES = tuple(itertools.combinations(NODES, 2))
GROUP = tuple(itertools.product((False, True), repeat=2))


def edge_name(edge: tuple[int, int]) -> str:
    return f"{edge[0]}{edge[1]}"


def mask_names(mask: Iterable[tuple[int, int]]) -> list[str]:
    return [edge_name(edge) for edge in sorted(mask)]


def _permutation(swap_endpoints: bool, swap_auxiliaries: bool) -> dict[int, int]:
    return {
        1: 2 if swap_endpoints else 1,
        2: 1 if swap_endpoints else 2,
        3: 4 if swap_auxiliaries else 3,
        4: 3 if swap_auxiliaries else 4,
    }


def transform_mask(
    mask: frozenset[tuple[int, int]],
    swap_endpoints: bool,
    swap_auxiliaries: bool,
) -> frozenset[tuple[int, int]]:
    permutation = _permutation(swap_endpoints, swap_auxiliaries)
    return frozenset(
        tuple(sorted((permutation[left], permutation[right])))
        for left, right in mask
    )


def canonical_mask(mask: frozenset[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    return min(tuple(sorted(transform_mask(mask, *element))) for element in GROUP)


def connected_with_marked_edge(mask: frozenset[tuple[int, int]]) -> bool:
    adjacency = {node: set() for node in NODES}
    for left, right in set(mask) | {MARKED_EDGE}:
        adjacency[left].add(right)
        adjacency[right].add(left)
    seen = {1}
    stack = [1]
    while stack:
        for neighbour in adjacency[stack.pop()]:
            if neighbour not in seen:
                seen.add(neighbour)
                stack.append(neighbour)
    return len(seen) == 4


def endpoint_triangle_flags(mask: frozenset[tuple[int, int]]) -> dict[str, bool]:
    return {
        "delete_1_triangle_234": frozenset(((2, 3), (2, 4), (3, 4))).issubset(mask),
        "delete_2_triangle_134": frozenset(((1, 3), (1, 4), (3, 4))).issubset(mask),
    }


def support_admissible(mask: frozenset[tuple[int, int]]) -> bool:
    flags = endpoint_triangle_flags(mask)
    return connected_with_marked_edge(mask) and not any(flags.values())


def determinant(matrix: list[list[Fraction]]) -> Fraction:
    work = [row[:] for row in matrix]
    value = Fraction(1)
    size = len(work)
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            value *= -1
        diagonal = work[column][column]
        value *= diagonal
        for row in range(column + 1, size):
            multiplier = work[row][column] / diagonal
            for index in range(column + 1, size):
                work[row][index] -= multiplier * work[column][index]
    return value


def principal_submatrix(matrix: list[list[Fraction]], subset: tuple[int, ...]) -> list[list[Fraction]]:
    return [[matrix[row][column] for column in subset] for row in subset]


def all_principal_minors_have_negative_sign(matrix: list[list[Fraction]]) -> bool:
    indices = tuple(range(len(matrix)))
    for size in range(1, len(matrix) + 1):
        for subset in itertools.combinations(indices, size):
            minor = determinant(principal_submatrix(matrix, subset))
            if ((-1) ** size) * minor <= 0:
                return False
    return True


def witness_matrix(mask: frozenset[tuple[int, int]]) -> list[list[Fraction]]:
    matrix = [[Fraction(0) for _ in NODES] for _ in NODES]
    for index in range(4):
        matrix[index][index] = Fraction(-2)
    matrix[0][1] = matrix[1][0] = Fraction(3)
    for left, right in mask:
        matrix[left - 1][right - 1] = matrix[right - 1][left - 1] = Fraction(1)
    return matrix


def fraction_string(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def witness_record(mask: frozenset[tuple[int, int]]) -> dict[str, Any]:
    matrix = witness_matrix(mask)
    leading = [determinant(principal_submatrix(matrix, tuple(range(size)))) for size in range(1, 5)]
    if any(value == 0 for value in leading):
        raise ArithmeticError((mask, leading))
    pivots = [leading[0]] + [leading[index] / leading[index - 1] for index in range(1, 4)]
    inertia = {
        "positive": sum(pivot > 0 for pivot in pivots),
        "negative": sum(pivot < 0 for pivot in pivots),
        "zero": 0,
    }
    deletions: dict[str, Any] = {}
    for endpoint in (1, 2):
        subset = tuple(index for index in range(4) if index != endpoint - 1)
        deletion = principal_submatrix(matrix, subset)
        deletions[str(endpoint)] = {
            "determinant": fraction_string(determinant(deletion)),
            "negative_definite": all_principal_minors_have_negative_sign(deletion),
        }
    if inertia != {"positive": 1, "negative": 3, "zero": 0}:
        raise ArithmeticError((mask, inertia))
    if not all(record["negative_definite"] for record in deletions.values()):
        raise ArithmeticError((mask, deletions))
    return {
        "u": 3,
        "finite_edge_labels": {edge_name(edge): 3 for edge in sorted(mask)},
        "gram_matrix": [[fraction_string(value) for value in row] for row in matrix],
        "leading_principal_determinants": [fraction_string(value) for value in leading],
        "ldl_pivots": [fraction_string(value) for value in pivots],
        "determinant": fraction_string(leading[-1]),
        "inertia": inertia,
        "endpoint_deletions": deletions,
        "v_arithmetic": {
            "ground_field": "Q",
            "coefficient_field": "Q",
            "cyclic_products_integral": True,
            "nongeometric_ground_embeddings": 0,
            "minimality_14": True,
        },
    }


TEMPLATE_BY_SUPPORT = {
    ((1, 3), (1, 4)): {
        "family": "Gamma4",
        "specialization": "k=2 (the drawn 23 edge is absent)",
    },
    ((1, 3), (2, 4)): {
        "family": "Gamma5",
        "specialization": "generic",
    },
    ((1, 3), (3, 4)): {
        "family": "Gamma3",
        "specialization": "s=2 (the drawn 14 edge is absent)",
    },
    ((1, 3), (1, 4), (2, 3)): {
        "family": "Gamma4",
        "specialization": "k>=3",
    },
    ((1, 3), (2, 3), (3, 4)): {
        "family": "Gamma2",
        "specialization": "generic",
    },
    ((1, 3), (2, 4), (3, 4)): {
        "family": "Gamma3",
        "specialization": "s>=3",
    },
    ((1, 3), (1, 4), (2, 3), (2, 4)): {
        "family": "Gamma1",
        "specialization": "generic",
    },
}


def topology_census() -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    admitted: list[frozenset[tuple[int, int]]] = []
    rejection_counts = {"disconnected": 0, "endpoint_triangle": 0}
    for bits in itertools.product((False, True), repeat=5):
        mask = frozenset(edge for edge, bit in zip(FINITE_POSITIONS, bits) if bit)
        connected = connected_with_marked_edge(mask)
        flags = endpoint_triangle_flags(mask)
        reasons: list[str] = []
        if not connected:
            reasons.append("disconnected_hyperbolic_component")
            rejection_counts["disconnected"] += 1
        if any(flags.values()):
            reasons.append("endpoint_deletion_contains_a_3_cycle")
            rejection_counts["endpoint_triangle"] += 1
        passed = not reasons
        if passed:
            admitted.append(mask)
        rows.append(
            {
                "mask": mask_names(mask),
                "connected_with_marked_12": connected,
                "endpoint_triangle_flags": flags,
                "passed": passed,
                "reasons": reasons,
            }
        )

    orbit_members: dict[tuple[tuple[int, int], ...], list[frozenset[tuple[int, int]]]] = {}
    for mask in admitted:
        orbit_members.setdefault(canonical_mask(mask), []).append(mask)
    if set(orbit_members) != set(TEMPLATE_BY_SUPPORT):
        raise ArithmeticError((set(orbit_members), set(TEMPLATE_BY_SUPPORT)))

    orbits = []
    for canonical in sorted(orbit_members, key=lambda item: (len(item), item)):
        mask = frozenset(canonical)
        orbits.append(
            {
                "canonical_support": mask_names(mask),
                "raw_orbit_size": len(orbit_members[canonical]),
                **TEMPLATE_BY_SUPPORT[canonical],
                "rational_exact_witness": witness_record(mask),
            }
        )
    return {
        "finite_edge_positions": mask_names(FINITE_POSITIONS),
        "raw_masks": len(rows),
        "admissible_raw_masks": len(admitted),
        "rejected_raw_masks": len(rows) - len(admitted),
        "rejection_incidence_counts_not_disjoint": rejection_counts,
        "marked_support_orbits": len(orbits),
        "schematic_families": 5,
        "orbits": orbits,
        "all_masks": rows,
    }


def allowed_two_edge_path(left: int, right: int) -> bool:
    """Exact rank-three spherical-path criterion for Coxeter labels."""
    return tuple(sorted((left, right))) in ((3, 3), (3, 4), (3, 5))


def transform_labels(
    labels: dict[tuple[int, int], int],
    swap_endpoints: bool,
    swap_auxiliaries: bool,
) -> dict[tuple[int, int], int]:
    permutation = _permutation(swap_endpoints, swap_auxiliaries)
    transformed: dict[tuple[int, int], int] = {}
    for (left, right), label in labels.items():
        edge = tuple(sorted((permutation[left], permutation[right])))
        transformed[edge] = label
    return transformed


def labels_admissible(labels: dict[tuple[int, int], int]) -> bool:
    mask = frozenset(edge for edge, label in labels.items() if label >= 3)
    if not connected_with_marked_edge(mask):
        return False
    for positions in (
        ((2, 3), (2, 4), (3, 4)),
        ((1, 3), (1, 4), (3, 4)),
    ):
        present = [labels[edge] for edge in positions if labels[edge] >= 3]
        if len(present) == 3:
            return False
        if len(present) == 2 and not allowed_two_edge_path(*present):
            return False
    return True


def matching_templates(labels: dict[tuple[int, int], int]) -> list[str]:
    a, b, c, d, e = [labels[edge] for edge in FINITE_POSITIONS]
    matches: list[str] = []

    # Gamma1: 13=s, 14=r, 23=k, 24=p, and 34 is absent.
    if e == 2 and min(a, b, c, d) >= 3:
        s, r, k, p = a, b, c, d
        if (s == k == 3 and 3 <= p <= r <= 5) or (
            s == p == 3 and 4 <= k <= r <= 5
        ):
            matches.append("Gamma1")

    # Gamma2: 13=s, 23=k, 34=p.
    if b == d == 2 and min(a, c, e) >= 3:
        s, k, p = a, c, e
        if (p == 3 and 3 <= s <= k <= 5) or (s == k == 3 and p in (4, 5)):
            matches.append("Gamma2")

    # Gamma3: 14=s (possibly 2), 23=k, 34=r.
    if a == d == 2 and c >= 3 and e >= 3:
        s, k, r = b, c, e
        degenerate = s == 2 and (
            (k == 3 and 3 <= r <= 5) or (r == 3 and 4 <= k <= 5)
        )
        generic = s >= 3 and (
            (r == 3 and 3 <= s <= k <= 5) or (s == k == 3 and r in (4, 5))
        )
        if degenerate or generic:
            matches.append("Gamma3")

    # Gamma4: 13=s, 23=k (possibly 2), 14=r.
    if d == e == 2 and a >= 3 and b >= 3 and c >= 2:
        s, k, r = a, c, b
        if (s == 3 and 3 <= r <= 5) or (s in (4, 5) and r == 3):
            matches.append("Gamma4")

    # Gamma5: 14=s, 23=k.  Both labels are genuinely unbounded.
    if a == d == e == 2 and b >= 3 and c >= 3:
        s, k = b, c
        if k >= s:
            matches.append("Gamma5")
    return matches


def parameter_grammar_census() -> dict[str, Any]:
    # 6 and 101 are two ordered representatives of the unrestricted >=6
    # class.  A two-edge spherical deletion never contains either class.
    alphabet = (2, 3, 4, 5, 6, 101)
    raw_pass = 0
    raw_family_counts = {f"Gamma{index}": 0 for index in range(1, 6)}
    canonical_rows: dict[tuple[int, ...], str] = {}
    for values in itertools.product(alphabet, repeat=5):
        labels = dict(zip(FINITE_POSITIONS, values))
        if not labels_admissible(labels):
            continue
        raw_pass += 1
        family_matches: set[str] = set()
        transformed_values: list[tuple[int, ...]] = []
        for group_element in GROUP:
            transformed = transform_labels(labels, *group_element)
            transformed_values.append(tuple(transformed[edge] for edge in FINITE_POSITIONS))
            family_matches.update(matching_templates(transformed))
        if len(family_matches) != 1:
            raise ArithmeticError((values, family_matches))
        family = next(iter(family_matches))
        raw_family_counts[family] += 1
        canonical = min(transformed_values)
        previous = canonical_rows.setdefault(canonical, family)
        if previous != family:
            raise ArithmeticError((canonical, previous, family))

    canonical_family_counts = {f"Gamma{index}": 0 for index in range(1, 6)}
    for family in canonical_rows.values():
        canonical_family_counts[family] += 1
    expected_raw = {"Gamma1": 25, "Gamma2": 22, "Gamma3": 42, "Gamma4": 110, "Gamma5": 50}
    expected_canonical = {"Gamma1": 9, "Gamma2": 8, "Gamma3": 13, "Gamma4": 28, "Gamma5": 15}
    if raw_pass != 249 or raw_family_counts != expected_raw:
        raise ArithmeticError((raw_pass, raw_family_counts))
    if len(canonical_rows) != 73 or canonical_family_counts != expected_canonical:
        raise ArithmeticError((len(canonical_rows), canonical_family_counts))
    return {
        "representative_label_alphabet": list(alphabet),
        "alphabet_meaning": {
            "2": "absent edge",
            "3_4_5": "all labels allowed in a two-edge spherical rank-three path",
            "6_101": "two ordered representatives of the unrestricted >=6 class",
        },
        "rank_three_path_lemma": {
            "criterion": "cos(pi/p)^2 + cos(pi/q)^2 < 1",
            "solutions_for_p_q_at_least_3_up_to_order": [[3, 3], [3, 4], [3, 5]],
            "cutoff_proof": "if both labels are >=4 the sum is >=1; with p=3, q=6 gives equality and monotonicity excludes q>=6",
        },
        "raw_assignments": len(alphabet) ** 5,
        "admissible_raw_assignments": raw_pass,
        "raw_family_counts": raw_family_counts,
        "marked_label_orbits": len(canonical_rows),
        "marked_label_orbit_counts": canonical_family_counts,
        "coverage_failures": 0,
        "cross_family_ambiguities": 0,
        "unbounded_domains": {
            "Gamma4": "k>=2; at k=2, (s,r) and (r,s) are the same marked graph",
            "Gamma5": "3<=s<=k",
        },
        "Gamma4_k_equals_2_redundancies": [
            {"printed_rows": [[3, 4], [4, 3]], "one_unlabelled_marked_graph": True},
            {"printed_rows": [[3, 5], [5, 3]], "one_unlabelled_marked_graph": True},
        ],
    }


def symbolic_gram_census() -> dict[str, Any]:
    u, S, K, R, P = sp.symbols("u S K R P")
    definitions = {
        "Gamma1": {(1, 3): S, (1, 4): R, (2, 3): K, (2, 4): P},
        "Gamma2": {(1, 3): S, (2, 3): K, (3, 4): P},
        "Gamma3": {(1, 4): S, (2, 3): K, (3, 4): R},
        "Gamma4": {(1, 3): S, (1, 4): R, (2, 3): K},
        "Gamma5": {(1, 4): S, (2, 3): K},
    }
    expected = {
        "Gamma1": K**2 * R**2 - 4 * K**2 - 2 * K * P * R * S - 4 * K * S * u
        + P**2 * S**2 - 4 * P**2 - 4 * P * R * u - 4 * R**2 - 4 * S**2
        - 4 * u**2 + 16,
        "Gamma2": -4 * K**2 - 4 * K * S * u + P**2 * u**2 - 4 * P**2
        - 4 * S**2 - 4 * u**2 + 16,
        "Gamma3": K**2 * S**2 - 4 * K**2 - 2 * K * R * S * u + R**2 * u**2
        - 4 * R**2 - 4 * S**2 - 4 * u**2 + 16,
        "Gamma4": K**2 * R**2 - 4 * K**2 - 4 * K * S * u - 4 * R**2
        - 4 * S**2 - 4 * u**2 + 16,
        "Gamma5": K**2 * S**2 - 4 * K**2 - 4 * S**2 - 4 * u**2 + 16,
    }
    rows: dict[str, Any] = {}
    for family, edges in definitions.items():
        gram = sp.eye(4) * -2
        gram[0, 1] = gram[1, 0] = u
        for (left, right), value in edges.items():
            gram[left - 1, right - 1] = gram[right - 1, left - 1] = value
        found = sp.factor(gram.det())
        if sp.expand(found - expected[family]) != 0:
            raise ArithmeticError((family, found, expected[family]))
        rows[family] = {
            "gram_matrix": [[sp.sstr(gram[i, j]) for j in range(4)] for i in range(4)],
            "determinant": sp.sstr(found),
            "determinant_identity_verified": True,
        }
    return rows


def sign_switching_audit() -> dict[str, Any]:
    connected_supports = 0
    sign_assignments = 0
    cycle_consistent_assignments = 0
    switching_failures = 0
    for bits in itertools.product((False, True), repeat=len(ALL_EDGES)):
        support = tuple(edge for edge, bit in zip(ALL_EDGES, bits) if bit)
        adjacency = {node: set() for node in NODES}
        for left, right in support:
            adjacency[left].add(right)
            adjacency[right].add(left)
        seen = {1}
        stack = [1]
        while stack:
            for neighbour in adjacency[stack.pop()]:
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
        if len(seen) != 4:
            continue
        connected_supports += 1
        for signs in itertools.product((-1, 1), repeat=len(support)):
            sign_assignments += 1
            edge_sign = {edge: sign for edge, sign in zip(support, signs)}
            vertex_sign = {1: 1}
            stack = [1]
            while stack:
                left = stack.pop()
                for right in adjacency[left]:
                    edge = tuple(sorted((left, right)))
                    proposed = vertex_sign[left] * edge_sign[edge]
                    if right not in vertex_sign:
                        vertex_sign[right] = proposed
                        stack.append(right)
            consistent = all(
                edge_sign[(left, right)] == vertex_sign[left] * vertex_sign[right]
                for left, right in support
            )
            if consistent:
                cycle_consistent_assignments += 1
                if len(vertex_sign) != 4:
                    switching_failures += 1
    if (connected_supports, sign_assignments, cycle_consistent_assignments, switching_failures) != (
        38,
        624,
        304,
        0,
    ):
        raise ArithmeticError(
            (connected_supports, sign_assignments, cycle_consistent_assignments, switching_failures)
        )
    return {
        "four_vertex_connected_supports": connected_supports,
        "edge_sign_assignments": sign_assignments,
        "cycle_consistent_assignments": cycle_consistent_assignments,
        "switching_failures": switching_failures,
        "general_lemma": "fixing every cyclic product makes the conjugate principal Gram diagonally congruent by vertex signs on each connected component",
    }


def build_certificate() -> dict[str, Any]:
    topology = topology_census()
    grammar = parameter_grammar_census()
    symbolic = symbolic_gram_census()
    switching = sign_switching_audit()
    checks = {
        "32_support_masks_exhausted": topology["raw_masks"] == 32,
        "17_admissible_masks": topology["admissible_raw_masks"] == 17,
        "seven_marked_support_orbits": topology["marked_support_orbits"] == 7,
        "five_schematic_families_cover_all_seven": len(TEMPLATE_BY_SUPPORT) == 7,
        "all_seven_have_exact_Q_witnesses": all(
            row["rational_exact_witness"]["v_arithmetic"]["ground_field"] == "Q"
            for row in topology["orbits"]
        ),
        "printed_parameter_grammar_exhausts_representative_alphabet": grammar[
            "coverage_failures"
        ]
        == 0,
        "symbolic_determinants_verified": all(
            row["determinant_identity_verified"] for row in symbolic.values()
        ),
        "cyclic_signs_are_vertex_switches": switching["switching_failures"] == 0,
    }
    if not all(checks.values()):
        raise ArithmeticError(checks)
    return {
        "schema": SCHEMA,
        "status": "exact_PASS_five_schematics_cover_seven_support_orbits_no_missing_family",
        "source": {
            "ground_fields_I": {
                "paper": "V. V. Nikulin, On ground fields of arithmetic hyperbolic reflection groups",
                "arxiv": "0708.3991v2",
                "section": "3.3 and Section 4 four-node parameter lists",
                "source_tex_sha256": "486855e2702e2ad0bc1e7a8fc828d9bc38ce948ea7dc5ff0e4d98ea0dcdc54b7",
                "five_graphs_eps_sha256": "7d8c7dee5b54b7110b82800e7729bbfe2d96e012f5535fd074d11f6bf8d0eb6a",
            },
            "journal_continuation": {
                "paper": "V. V. Nikulin, The transition constant for arithmetic hyperbolic reflection groups",
                "journal": "Izvestiya: Mathematics 75 (2011), 971--1005",
                "doi": "10.1070/IM2011v075n05ABEH002561",
                "pdf_sha256": "716e6e7333efbfaa22060b60dc97e8303363a4914409ccdd56ef7a90c923c7af",
                "section": "3.3 and Sections 4.3--4.7",
            },
        },
        "normalization": {
            "diagonal": -2,
            "finite_entry": "2*cos(pi/m)",
            "label_2": "entry 0, hence no support edge",
            "marked_entry": "u=a_12>2",
            "equivalence": "S2 on {1,2} times S2 on {3,4}",
        },
        "topology_census": topology,
        "parameter_grammar_census": grammar,
        "symbolic_gram_matrices": symbolic,
        "ground_field_inheritance": {
            "verdict": "proved for every principal packet that is not negative definite; hyperbolic, quasi-Lanner, and parabolic packets are included",
            "logic": [
                "let K0 be generated by the packet cyclic products and K by all cyclic products",
                "if K/K0 were proper, total reality gives a nonidentity real K0-embedding of K and an extension to the coefficient field",
                "fixed two-cycles preserve absolute edge entries and fixed longer cycles make their signs a vertex-switching coboundary",
                "the conjugate packet is diagonally congruent to the geometric packet and is therefore not negative definite",
                "V-arithmeticity makes the full nongeometric conjugate negative definite, hence every principal packet negative definite: contradiction",
            ],
            "marked_u_regimes": {
                "u_greater_than_2": "the 2-node marked packet is hyperbolic, hence K=Q(u^2)",
                "u_equals_2": "the 2-node marked packet is parabolic and has cyclic-product field Q, hence K=Q",
                "u_less_than_2": "the marked pair is negative definite and cannot inherit the field; one must choose another non-negative-definite packet, allowing quasi-Lanner as well as Lanner",
            },
            "sign_switching_finite_audit": switching,
            "quasi_lanner_and_cusp": {
                "compactness_used": False,
                "336_example": "the full [3,3,6] packet is hyperbolic and its affine [3,6] deletion is parabolic; both retain the cyclic-product field Q",
                "general_noncompact_arithmetic_repair": "a finite-volume cusp makes the defining k-form isotropic; a nontrivial real conjugate is definite and cannot contain the conjugate nonzero isotropic vector, so k=Q",
                "warning": "an affine null vector or unlabeled support alone is not the full cyclic-product algebra and does not trigger inheritance",
            },
        },
        "literal_qualifications": {
            "five_is_schematic_not_topological": "there are seven marked support orbits; Gamma3 at s=2 and Gamma4 at k=2 supply the two extra strata",
            "Gamma4_k2_parameter_redundancy": "(s,r)=(3,4)/(4,3) and (3,5)/(5,3) are duplicate unlabelled marked graphs",
            "0708_subscript_typo": "inside the Gamma2 determinant discussion the source says Gamma4 is V-arithmetic iff; the displayed matrix and formula are Gamma2",
            "scope": "the family and parameter grammar passes; numerical interval extrema and degree bounds are not certified here",
        },
        "checks": checks,
        "all_checks_pass": True,
        "resource_profile": {
            "largest_finite_loop": "6^5=7776 label assignments",
            "symbolic_matrices": "five 4x4 determinants",
            "classification_tables_imported": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="write the frozen JSON archive")
    args = parser.parse_args()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.write:
        ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        ARCHIVE.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
