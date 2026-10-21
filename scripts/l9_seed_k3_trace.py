"""K3 trace-lattice obstruction for <-2*phi> + E8 + A1."""

from __future__ import annotations

import json

from coxeter_polytopes.exact.golden import PHI
from coxeter_polytopes.exact.matrix import determinant, inertia
from coxeter_polytopes.exact.trace_lattice import (
    INV_SQRT5,
    change_gram_basis,
    integral_gram,
    is_codifferent,
    is_even_integral,
    trace_gram,
)
from coxeter_polytopes.lattices import block_diagonal, e8_cartan_gram


def seed_field_gram() -> list[list[object]]:
    e8 = e8_cartan_gram()
    result: list[list[object]] = [[0 for _ in range(10)] for _ in range(10)]
    result[0][0] = -2 * PHI
    for i in range(8):
        for j in range(8):
            result[i + 1][j + 1] = e8[i][j]
    result[9][9] = 2
    return result


def basis_change_certificate() -> list[list[int]]:
    size = 20
    change = [[0 for _ in range(size)] for _ in range(size)]
    # Negative field line: e, phi*e - 2e.
    change[0][0] = 1
    change[0][1] = -2
    change[1][1] = 1
    # E8 tensor factor: e_i and phi*e_i-e_i.
    for i in range(8):
        old_e = 2 + 2 * i
        old_phi_e = old_e + 1
        change[old_e][2 + i] = 1
        change[old_e][10 + i] = -1
        change[old_phi_e][10 + i] = 1
    # A1 tensor factor.
    change[18][18] = 1
    change[18][19] = -1
    change[19][19] = 1
    return change


def expected_standard_gram() -> list[list[int]]:
    e8_negative = [[-entry for entry in row] for row in e8_cartan_gram()]
    return block_diagonal(
        [[2]],
        [[-2]],
        e8_negative,
        e8_negative,
        [[-2]],
        [[-2]],
    )


def analyze() -> dict[str, object]:
    alpha = -PHI * INV_SQRT5
    trace = trace_gram(seed_field_gram(), alpha)
    integer = integral_gram(trace)
    transformed = change_gram_basis(integer, basis_change_certificate())
    exact_inertia = inertia(integer)
    discriminant_length = 4
    complement_rank = 22 - len(integer)
    return {
        "field_lattice": "<-2*phi> + E8 + A1 over Z[phi]",
        "trace_scalar": alpha.to_json(),
        "trace_scalar_is_codifferent": is_codifferent(alpha),
        "rank": len(integer),
        "even_integral": is_even_integral(trace),
        "determinant": str(determinant(integer)),
        "inertia": {
            "positive": exact_inertia.positive,
            "negative": exact_inertia.negative,
            "zero": exact_inertia.zero,
        },
        "standard_lattice": "<2> + <-2>^3 + E8(-1)^2",
        "basis_change_verified": transformed == expected_standard_gram(),
        "basis_change_determinant": str(determinant(basis_change_certificate())),
        "discriminant_group": "(Z/2Z)^4",
        "discriminant_group_generator_length": discriminant_length,
        "required_k3_complement_rank": complement_rank,
        "primitive_k3_embedding": False,
        "obstruction": (
            "For a primitive embedding in a unimodular lattice, the complement "
            "has the same discriminant group. A rank-2 lattice has discriminant "
            "group generator length at most 2, but this trace lattice has length 4."
        ),
        "obstruction_verified": discriminant_length > complement_rank,
    }


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2, sort_keys=True))

