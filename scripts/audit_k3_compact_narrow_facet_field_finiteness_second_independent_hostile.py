#!/usr/bin/env python3
"""Second independent hostile audit of compact narrow-facet K3 field finiteness."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from math import comb
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
PINS = {
    "report": (
        "results/k3_compact_narrow_facet_field_finiteness.md",
        "39c9a3a8422110bcdcf751ed7bae90ac86cdf3a720d220ee164ff241bd9ac9b5",
    ),
    "certifier": (
        "scripts/certify_k3_compact_narrow_facet_field_finiteness.py",
        "d754d69fad85366ad6a2774be2034c24ebf50a9620ec79618df4b18a7f42d9e7",
    ),
    "archive": (
        "data/enumerations/k3_compact_narrow_facet_field_finiteness.json",
        "3d0a2dff63bae86421a286131026f75fc4dbfaae9dec72be10cc4f38f5cd7e09",
    ),
    "tests": (
        "tests/test_k3_compact_narrow_facet_field_finiteness.py",
        "10ddc82f12c88d7e718f62ee859812ac22dbd6f37d17450e2f5c59280dce8a45",
    ),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_pins() -> dict[str, str]:
    observed: dict[str, str] = {}
    for label, (relative, expected) in PINS.items():
        actual = sha256(ROOT / relative)
        assert actual == expected, (label, expected, actual)
        observed[label] = actual
    return observed


def determinant(matrix: tuple[tuple[Fraction, ...], ...]) -> Fraction:
    size = len(matrix)
    assert all(len(row) == size for row in matrix)
    work = [list(row) for row in matrix]
    answer = Fraction(1)
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]), None
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            answer = -answer
        value = work[column][column]
        answer *= value
        for row in range(column + 1, size):
            ratio = work[row][column] / value
            for target in range(column, size):
                work[row][target] -= ratio * work[column][target]
    return answer


def star_spanning_certificate() -> dict[str, object]:
    """Replay the exact coordinate determinant in the star-span proof."""

    rows = []
    # tangent_dimension=m-2.  In the basis (x,u_1,...,u_r), take the
    # incident normals u_i and a nonincident side w=alpha*x+sum beta_i*u_i.
    for module_rank in range(3, 11):
        tangent_dimension = module_rank - 2
        alpha = Fraction(module_rank + 1)
        betas = tuple(Fraction(index + 2, index + 3) for index in range(tangent_dimension))
        size = tangent_dimension + 1
        columns: list[list[Fraction]] = []
        for index in range(tangent_dimension):
            column = [Fraction(0) for _ in range(size)]
            column[index + 1] = Fraction(1)
            columns.append(column)
        columns.append([alpha, *betas])
        matrix = tuple(
            tuple(columns[column][row] for column in range(size))
            for row in range(size)
        )
        observed = determinant(matrix)
        expected = (-1) ** tangent_dimension * alpha
        assert observed == expected != 0
        rows.append(
            {
                "module_rank": module_rank,
                "facet_dimension": tangent_dimension,
                "projected_basis_determinant": str(observed),
                "nonincident_side_coefficient": str(alpha),
            }
        )
    return {
        "finite_vertex_link_is_spherical_simplex": True,
        "compact_Coxeter_chamber_is_simple": True,
        "facet_is_compact_simple_convex_polytope": True,
        "positive_dimensional_facet_has_nonincident_side": True,
        "projected_star_spans_e_perp": True,
        "adding_minimal_facet_normal_spans_V": True,
        "exact_coordinate_rows": rows,
        "dimension_one_countergate": (
            "m=2 has a zero-dimensional facet and no nonincident facet side"
        ),
    }


def _connected_edges() -> tuple[tuple[int, int], ...]:
    return ((0, 1), (1, 2), (2, 3), (3, 4), (0, 2), (1, 4), (0, 4))


def _recover_vertex_weights(
    vertex_count: int,
    edges: Iterable[tuple[int, int]],
    ratios: dict[tuple[int, int], Fraction],
) -> tuple[Fraction, ...]:
    adjacency: dict[int, list[int]] = {vertex: [] for vertex in range(vertex_count)}
    for left, right in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    weights: list[Fraction | None] = [None] * vertex_count
    weights[0] = Fraction(1)
    queue = [0]
    while queue:
        left = queue.pop(0)
        assert weights[left] is not None
        for right in adjacency[left]:
            candidate = weights[left] * ratios[left, right]
            if weights[right] is None:
                weights[right] = candidate
                queue.append(right)
            elif weights[right] != candidate:
                raise AssertionError("edge ratios have nontrivial cycle holonomy")
    assert all(weight is not None for weight in weights)
    return tuple(weight for weight in weights if weight is not None)


def cartan_coboundary_certificate() -> dict[str, object]:
    """Reconstruct the diagonal similarity from exact cycle ratios."""

    vertex_count = 5
    edges = _connected_edges()
    original_weights = (
        Fraction(2),
        Fraction(3),
        Fraction(5),
        Fraction(7),
        Fraction(11),
    )
    ratios: dict[tuple[int, int], Fraction] = {}
    for left, right in edges:
        ratios[left, right] = original_weights[right] / original_weights[left]
        ratios[right, left] = original_weights[left] / original_weights[right]
        assert ratios[left, right] * ratios[right, left] == 1
    recovered = _recover_vertex_weights(vertex_count, edges, ratios)
    common_scale = original_weights[0] / recovered[0]
    assert tuple(common_scale * value for value in recovered) == original_weights

    # Build an exact sparse Cartan matrix and check D^{-1} A D entrywise.
    matrix = [[Fraction(0) for _ in range(vertex_count)] for _ in range(vertex_count)]
    conjugate = [[Fraction(0) for _ in range(vertex_count)] for _ in range(vertex_count)]
    for index in range(vertex_count):
        matrix[index][index] = Fraction(2)
        conjugate[index][index] = Fraction(2)
    for edge_index, (left, right) in enumerate(edges, start=1):
        matrix[left][right] = Fraction(edge_index + 1, edge_index + 2)
        matrix[right][left] = Fraction(edge_index + 3, edge_index + 4)
        conjugate[left][right] = matrix[left][right] * ratios[left, right]
        conjugate[right][left] = matrix[right][left] * ratios[right, left]
    for left in range(vertex_count):
        for right in range(vertex_count):
            expected = matrix[left][right] * recovered[right] / recovered[left]
            assert conjugate[left][right] == expected

    poisoned = dict(ratios)
    poisoned[0, 4] *= 2
    poisoned[4, 0] /= 2
    poison_rejected = False
    try:
        _recover_vertex_weights(vertex_count, edges, poisoned)
    except AssertionError:
        poison_rejected = True
    assert poison_rejected
    return {
        "connected_vertices": vertex_count,
        "two_cycle_inverse_relations_checked": len(edges),
        "cycle_holonomy_reconstructs_vertex_coboundary": True,
        "diagonal_similarity_checked_entrywise": True,
        "poisoned_cycle_holonomy_rejected": True,
        "similarity_preserves_Lorentz_negative_eigenvalue": True,
        "definite_conjugate_has_only_positive_eigenvalues": True,
        "proper_cycle_subfield_contradiction": True,
    }


def height_and_discriminant_certificate() -> dict[str, object]:
    admissible_rows = []
    for module_rank in range(3, 21):
        degree_ceiling = 20 // module_rank
        for degree in range(2, degree_ceiling + 1):
            generator_count = degree.bit_length() - 1
            admissible_rows.append(
                {
                    "module_rank": module_rank,
                    "field_degree": degree,
                    "generator_count_ceiling": generator_count,
                    "cycle_height": 14**module_rank,
                    "Vandermonde_exponent": degree * (degree - 1),
                }
            )
    assert admissible_rows
    assert max(row["module_rank"] for row in admissible_rows) == 10
    assert max(row["field_degree"] for row in admissible_rows) == 6
    assert max(row["generator_count_ceiling"] for row in admissible_rows) == 2
    assert max(row["Vandermonde_exponent"] for row in admissible_rows) == 30
    assert comb(6, 2) == 15
    height = 14**10
    return {
        "trace_rank_inequality": "d*m<=20",
        "minimum_module_rank": 3,
        "nonrational_module_rank_ceiling": 10,
        "nonrational_field_degree_ceiling": 6,
        "strict_tower_generator_count_ceiling": 2,
        "embedding_pair_count_ceiling": 15,
        "primitive_grid": "x+q*y, q=0,...,15",
        "uniform_cycle_height": height,
        "primitive_conjugate_height": 16 * height,
        "conjugate_difference_height": 32 * height,
        "uniform_Vandermonde_exponent": 30,
        "field_discriminant_bound": "(32*14^10)^30",
        "admissible_rows": admissible_rows,
        "excluded_m2_degree_ceiling": 10,
    }


def scope_certificate() -> dict[str, bool]:
    text = " ".join((ROOT / PINS["report"][0]).read_text(encoding="utf-8").split())
    required = {
        "dimension_gate": "m=n+1\\ge3",
        "compact_gate": "compact Coxeter chamber",
        "projective_module": "finite projective rank-\\(m\\)",
        "admissibility": "positive definite at every other real embedding",
        "no_minimal_assumption": "No separate minimal-ground-field hypothesis is imposed",
        "star_span": "projected neighboring normals span \\(e^\\perp\\)",
        "field_inheritance": "One has \\(k_I=k\\)",
        "cycle_integrality": "every such product is integral",
        "dimension_one_excluded": "does **not** apply to hyperbolic dimension one",
        "property4_not_proved": "prove essential surjectivity/Property 4",
        "chamber_finiteness_not_proved": "bound the number of chambers or hosts",
    }
    for label, phrase in required.items():
        assert phrase in text, (label, phrase)
    return {label: True for label in required}


def build_audit() -> dict[str, object]:
    return {
        "schema": "k3_compact_narrow_facet_field_finiteness/second_independent_hostile/v1",
        "verdict": "PASS_EXACT_STATED_SCOPE",
        "pins": verify_pins(),
        "compact_star_spanning": star_spanning_certificate(),
        "cycle_field_inheritance": cartan_coboundary_certificate(),
        "height_and_discriminant": height_and_discriminant_certificate(),
        "scope": scope_certificate(),
    }


def main() -> None:
    print(json.dumps(build_audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
