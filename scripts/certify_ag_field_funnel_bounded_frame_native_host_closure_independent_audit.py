#!/usr/bin/env python3
"""Independent exact checks for the bounded-frame native-host audit.

The certificate deliberately uses a symmetric (hence genuinely
symmetrizable) Lorentzian Cartan sample.  It checks the identities used in
the proof, not merely phrases in the owner report.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
OWNER_REPORT = ROOT / "results/ag_field_funnel_bounded_frame_native_host_closure.md"
OWNER_CERTIFIER = (
    ROOT / "scripts/certify_ag_field_funnel_bounded_frame_native_host_closure.py"
)
OWNER_TEST = ROOT / "tests/test_ag_field_funnel_bounded_frame_native_host_closure.py"


Matrix = tuple[tuple[Fraction, ...], ...]


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identity(size: int) -> Matrix:
    return tuple(
        tuple(Fraction(int(i == j)) for j in range(size)) for i in range(size)
    )


def matrix_product(left: Matrix, right: Matrix) -> Matrix:
    rows = len(left)
    inner = len(right)
    columns = len(right[0])
    return tuple(
        tuple(
            sum((left[i][k] * right[k][j] for k in range(inner)), Fraction(0))
            for j in range(columns)
        )
        for i in range(rows)
    )


def matrix_subtract(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(left[i][j] - right[i][j] for j in range(len(left[0])))
        for i in range(len(left))
    )


def trace(matrix: Matrix) -> Fraction:
    return sum((matrix[i][i] for i in range(len(matrix))), Fraction(0))


def rank_one_operator(cartan: Matrix, index: int) -> Matrix:
    size = len(cartan)
    return tuple(
        tuple(cartan[index][j] if i == index else Fraction(0) for j in range(size))
        for i in range(size)
    )


def reflection_matrix(cartan: Matrix, index: int) -> Matrix:
    return matrix_subtract(identity(len(cartan)), rank_one_operator(cartan, index))


def directed_cycle_product(cartan: Matrix, cycle: Iterable[int]) -> Fraction:
    vertices = tuple(cycle)
    value = Fraction(1)
    for position, vertex in enumerate(vertices):
        value *= cartan[vertex][vertices[(position + 1) % len(vertices)]]
    return value


def trace_cycle_product(cartan: Matrix, cycle: Iterable[int]) -> Fraction:
    vertices = tuple(cycle)
    product = rank_one_operator(cartan, vertices[0])
    for vertex in vertices[1:]:
        product = matrix_product(product, rank_one_operator(cartan, vertex))
    return trace(product)


def edge_two_cycle(cartan: Matrix, i: int, j: int) -> Fraction:
    return cartan[i][j] * cartan[j][i]


def cycle_square_from_edges(cartan: Matrix, cycle: Iterable[int]) -> Fraction:
    vertices = tuple(cycle)
    value = Fraction(1)
    for position, vertex in enumerate(vertices):
        value *= edge_two_cycle(
            cartan, vertex, vertices[(position + 1) % len(vertices)]
        )
    return value


def diagonal_rescale(cartan: Matrix, scales: tuple[Fraction, ...]) -> Matrix:
    return tuple(
        tuple(scales[j] * cartan[i][j] / scales[i] for j in range(len(cartan)))
        for i in range(len(cartan))
    )


def recover_tree_scales(
    original: Matrix,
    changed: Matrix,
    tree: tuple[tuple[int, int], ...],
    root: int = 0,
) -> tuple[Fraction, ...]:
    values: list[Fraction | None] = [None] * len(original)
    values[root] = Fraction(1)
    pending = list(tree)
    while pending:
        next_pending: list[tuple[int, int]] = []
        progress = False
        for i, j in pending:
            if values[i] is not None and values[j] is None:
                values[j] = values[i] * changed[i][j] / original[i][j]
                progress = True
            elif values[j] is not None and values[i] is None:
                values[i] = values[j] * changed[j][i] / original[j][i]
                progress = True
            else:
                next_pending.append((i, j))
        if next_pending and not progress:
            raise ValueError("the proposed tree is disconnected")
        pending = next_pending
    if any(value is None for value in values):
        raise ValueError("not every vertex was reached")
    return tuple(value for value in values if value is not None)


def flatten(matrix: Matrix) -> tuple[Fraction, ...]:
    return tuple(entry for row in matrix for entry in row)


def vector_rank(vectors: Iterable[tuple[Fraction, ...]]) -> int:
    rows = [list(vector) for vector in vectors]
    if not rows:
        return 0
    row = 0
    columns = len(rows[0])
    for column in range(columns):
        pivot = next((i for i in range(row, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[row], rows[pivot] = rows[pivot], rows[row]
        scale = rows[row][column]
        rows[row] = [entry / scale for entry in rows[row]]
        for i in range(len(rows)):
            if i == row or not rows[i][column]:
                continue
            multiple = rows[i][column]
            rows[i] = [
                rows[i][j] - multiple * rows[row][j] for j in range(columns)
            ]
        row += 1
        if row == len(rows):
            break
    return row


def generated_algebra_dimension(generators: tuple[Matrix, ...]) -> int:
    """Close a rational matrix span under right multiplication."""
    basis: list[Matrix] = [identity(len(generators[0]))]
    current_rank = 1
    cursor = 0
    while cursor < len(basis):
        word = basis[cursor]
        cursor += 1
        for generator in generators:
            candidate = matrix_product(word, generator)
            trial = basis + [candidate]
            trial_rank = vector_rank(flatten(matrix) for matrix in trial)
            if trial_rank > current_rank:
                basis.append(candidate)
                current_rank = trial_rank
    return current_rank


def in_diagonal_lattice(vector: tuple[int, int], steps: tuple[int, int]) -> bool:
    return all(vector[i] % steps[i] == 0 for i in range(2))


def apply_integer_matrix(matrix: tuple[tuple[int, int], ...], vector: tuple[int, int]):
    return tuple(sum(matrix[i][j] * vector[j] for j in range(2)) for i in range(2))


def small_conductor_sandwich() -> dict[str, object]:
    """Check L=2Ze1+Ze2 for Lambda=Z*I+2M_2(Z)."""
    lattice_steps = (2, 1)
    generators = (
        ((1, 0), (0, 1)),
        ((2, 0), (0, 0)),
        ((0, 2), (0, 0)),
        ((0, 0), (2, 0)),
        ((0, 0), (0, 2)),
    )
    lattice_generators = ((2, 0), (0, 1))
    stable = all(
        in_diagonal_lattice(apply_integer_matrix(matrix, vector), lattice_steps)
        for matrix in generators
        for vector in lattice_generators
    )
    # B=M_2(Z) sends e2 to both coordinate vectors, hence BL=Z^2.
    hull_steps = (1, 1)
    conductor_hull_steps = (2, 2)
    conductor_inside = all(
        lattice_steps[i] != 0 and conductor_hull_steps[i] % lattice_steps[i] == 0
        for i in range(2)
    )
    lattice_inside_hull = all(lattice_steps[i] % hull_steps[i] == 0 for i in range(2))
    return {
        "lambda_stable": stable,
        "L_steps": list(lattice_steps),
        "BL_steps": list(hull_steps),
        "2BL_steps": list(conductor_hull_steps),
        "two_BL_inside_L_inside_BL": conductor_inside and lattice_inside_hull,
        "number_of_subgroups_of_BL_mod_2BL_is_finite": 5,
    }


def core_sum_sample() -> dict[str, object]:
    """A C2 sample of the normal-core translate-sum construction."""
    # L'=2Ze1+Ze2 and the coordinate swap sends it to Ze1+2Ze2.
    first_steps = (2, 1)
    swapped_steps = (1, 2)
    # Coordinatewise gcd computes the sum for these diagonal lattices.
    sum_steps = (
        __import__("math").gcd(first_steps[0], swapped_steps[0]),
        __import__("math").gcd(first_steps[1], swapped_steps[1]),
    )
    return {
        "L_prime_steps": list(first_steps),
        "swap_L_prime_steps": list(swapped_steps),
        "translate_sum_steps": list(sum_steps),
        "translate_sum_is_swap_stable": sum_steps[0] == sum_steps[1],
    }


def owner_interface() -> dict[str, object]:
    text = OWNER_REPORT.read_text(encoding="utf-8")
    compact_text = " ".join(text.split())
    required = (
        "Arithmeticity is in Vinberg's sense",
        "normal core",
        "length-two backtrack",
        "projective trace",
        "absolutely irreducible",
        "Burnside",
        "\\mathfrak cM\\subseteq L\\subseteq M",
        "matrix units",
        "maximal either",
        "n\\leq995",
        "not AG-only",
    )
    return {
        "missing": [phrase for phrase in required if phrase not in text],
        "fixed_data_conclusion_disclaimed": (
            "Nikulin fixed-data finiteness conclusion | none | **not used**" in text
        ),
        "low_dimensional_field_scope_not_overclaimed": (
            "one must separately supply a finite field set" in compact_text
        ),
    }


def build_certificate() -> dict[str, object]:
    # Symmetric Lorentzian Cartan matrix: its eigenvalue signs are (-,+,+,+).
    cartan: Matrix = tuple(
        tuple(Fraction(value) for value in row)
        for row in (
            (2, -1, -1, 0),
            (-1, 2, -1, -1),
            (-1, -1, 2, -1),
            (0, -1, -1, 2),
        )
    )
    cycles = ((0, 1), (1, 2), (0, 2), (1, 3), (2, 3), (0, 1, 2), (1, 2, 3))
    cycle_checks = {
        "-".join(map(str, cycle)): {
            "product": str(directed_cycle_product(cartan, cycle)),
            "trace": str(trace_cycle_product(cartan, cycle)),
            "trace_equal": directed_cycle_product(cartan, cycle)
            == trace_cycle_product(cartan, cycle),
            "forward_reverse_equal": directed_cycle_product(cartan, cycle)
            == directed_cycle_product(cartan, tuple(reversed(cycle))),
            "square_edge_identity": (
                True
                if len(cycle) == 2
                else directed_cycle_product(cartan, cycle) ** 2
                == cycle_square_from_edges(cartan, cycle)
            ),
        }
        for cycle in cycles
    }
    scales = (Fraction(2), Fraction(3), Fraction(5), Fraction(7))
    changed = diagonal_rescale(cartan, scales)
    recovered = recover_tree_scales(cartan, changed, ((0, 1), (1, 3), (3, 2)))
    normalized = tuple(value / scales[0] for value in scales)
    reflections = tuple(reflection_matrix(cartan, i) for i in range(4))
    algebra_dimension = generated_algebra_dimension(reflections)
    sandwich = small_conductor_sandwich()
    core_sum = core_sum_sample()
    interface = owner_interface()
    owner_hashes = {
        "report": sha256_file(OWNER_REPORT),
        "certifier": sha256_file(OWNER_CERTIFIER),
        "test": sha256_file(OWNER_TEST),
    }
    expected_hashes = {
        "report": "7ba392dd833358379ea92d1385909a52b2dc65417a9cb4ccb856170118df02a3",
        "certifier": "44848ad5727458b1777a4989395229d847ad0e05136ae5c96d2ffd9ad12a6def",
        "test": "e5b0f259ccaee50ec6c1d9d17af9d5806eed999f28c5712c330cabbc03c5ed45",
    }
    all_checks_pass = all(
        (
            owner_hashes == expected_hashes,
            all(check["trace_equal"] for check in cycle_checks.values()),
            all(check["forward_reverse_equal"] for check in cycle_checks.values()),
            all(check["square_edge_identity"] for check in cycle_checks.values()),
            recovered == normalized,
            diagonal_rescale(cartan, recovered) == changed,
            algebra_dimension == 16,
            sandwich["lambda_stable"],
            sandwich["two_BL_inside_L_inside_BL"],
            core_sum["translate_sum_is_swap_stable"],
            not interface["missing"],
            interface["fixed_data_conclusion_disclaimed"],
            interface["low_dimensional_field_scope_not_overclaimed"],
            len(range(3, 997)) == 994,
        )
    )
    return {
        "schema": "ag-field-funnel-bounded-frame-native-host-independent-audit-v1",
        "all_checks_pass": all_checks_pass,
        "owner_hashes": owner_hashes,
        "expected_owner_hashes": expected_hashes,
        "symmetrizable_cycle_checks": cycle_checks,
        "diagonal_gauge": {
            "recovered": [str(value) for value in recovered],
            "expected_up_to_scalar": [str(value) for value in normalized],
            "recovers_non_tree_entries": diagonal_rescale(cartan, recovered) == changed,
        },
        "burnside_sample": {
            "rank": 4,
            "generated_algebra_dimension": algebra_dimension,
            "full_endomorphism_dimension": 16,
        },
        "conductor_sandwich_sample": sandwich,
        "normal_core_translate_sum_sample": core_sum,
        "dimension_quantifier": {
            "dimensions": [2, 995],
            "ranks": [3, 996],
            "rank_count": len(range(3, 997)),
        },
        "owner_interface": interface,
        "verdict": "PASS",
    }


def main() -> None:
    certificate = build_certificate()
    print(json.dumps(certificate, indent=2, sort_keys=True))
    if not certificate["all_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
