"""Verify the K3 trace lattice attached to Bugaenko's audited L8."""

from __future__ import annotations

import json

from coxeter_polytopes.exact.golden import PHI
from coxeter_polytopes.exact.matrix import determinant, inertia, rectangular_pairing
from coxeter_polytopes.exact.trace_lattice import (
    INV_SQRT5,
    change_gram_basis,
    integral_gram,
    is_codifferent,
    is_even_integral,
    trace_gram,
)
from coxeter_polytopes.lattices import (
    block_diagonal,
    e8_cartan_gram,
    hyperbolic_plane,
    k3_lattice_gram,
)


def bugaenko_field_gram() -> list[list[object]]:
    e8 = e8_cartan_gram()
    result: list[list[object]] = [[0 for _ in range(9)] for _ in range(9)]
    result[0][0] = -2 * PHI
    for i in range(8):
        for j in range(8):
            result[i + 1][j + 1] = e8[i][j]
    return result


def basis_change_certificate() -> list[list[int]]:
    """Columns take the trace basis to <2,-2> + E8(-1) + E8(-1)."""

    size = 18
    change = [[0 for _ in range(size)] for _ in range(size)]
    # Negative O_F-line: e, phi*e - 2e.
    change[0][0] = 1
    change[0][1] = -2
    change[1][1] = 1
    # First E8 copy: e_i.  Second copy: phi*e_i - e_i.
    for i in range(8):
        old_e = 2 + 2 * i
        old_phi_e = old_e + 1
        change[old_e][2 + i] = 1
        change[old_e][10 + i] = -1
        change[old_phi_e][10 + i] = 1
    return change


def expected_standard_gram() -> list[list[int]]:
    e8_negative = [[-entry for entry in row] for row in e8_cartan_gram()]
    return block_diagonal([[2]], [[-2]], e8_negative, e8_negative)


def k3_embedding() -> list[list[int]]:
    """Embed the standard rank-18 trace lattice in U^3 + E8(-1)^2."""

    embedding = [[0 for _ in range(18)] for _ in range(22)]
    # <2>: e_1 + f_1.  <-2>: e_2 - f_2.
    embedding[0][0] = embedding[1][0] = 1
    embedding[2][1] = 1
    embedding[3][1] = -1
    # The two negative E8 summands map identically.
    for i in range(16):
        embedding[6 + i][2 + i] = 1
    return embedding


def primitive_extension() -> list[list[int]]:
    """Extend the embedding columns to a Z-basis, proving primitivity."""

    embedding = k3_embedding()
    extension = [row + [0, 0, 0, 0] for row in embedding]
    extension[1][18] = 1  # f_1
    extension[3][19] = 1  # f_2
    extension[4][20] = 1  # e_3
    extension[5][21] = 1  # f_3
    return extension


def orthogonal_complement() -> list[list[int]]:
    """Columns span <-2> + <2> + U, the embedding's orthogonal complement."""

    complement = [[0 for _ in range(4)] for _ in range(22)]
    complement[0][0] = 1
    complement[1][0] = -1
    complement[2][1] = complement[3][1] = 1
    complement[4][2] = 1
    complement[5][3] = 1
    return complement


def analyze() -> dict[str, object]:
    alpha = -PHI * INV_SQRT5
    trace = trace_gram(bugaenko_field_gram(), alpha)
    integer = integral_gram(trace)
    transformed = change_gram_basis(integer, basis_change_certificate())
    expected = expected_standard_gram()
    k3 = k3_lattice_gram()
    embedding = k3_embedding()
    complement = orthogonal_complement()
    embedded_gram = rectangular_pairing(embedding, k3, embedding)
    complement_gram = rectangular_pairing(complement, k3, complement)
    cross_pairing = rectangular_pairing(embedding, k3, complement)
    exact_inertia = inertia(integer)
    result = {
        "normalization_status": "verified_against_pinned_vinberg_fixture",
        "field_lattice": "<-2*phi> + E8 over Z[phi]",
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
        "standard_lattice": "<2> + <-2> + E8(-1) + E8(-1)",
        "basis_change_verified": transformed == expected,
        "basis_change_determinant": str(determinant(basis_change_certificate())),
        "primitive_k3_embedding_verified": (
            embedded_gram == expected
            and abs(determinant(primitive_extension()).a) == 1
        ),
        "primitive_extension_determinant": str(determinant(primitive_extension())),
        "orthogonal_complement": "<-2> + <2> + U",
        "orthogonal_complement_verified": (
            complement_gram == block_diagonal([[-2]], [[2]], hyperbolic_plane())
            and not any(entry for row in cross_pairing for entry in row)
        ),
    }
    return result


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2, sort_keys=True))
