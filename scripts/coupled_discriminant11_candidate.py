"""Exact K3 certificates for the first nonorthogonal rank-10 field form."""

from __future__ import annotations

import json

from coxeter_polytopes.exact.golden import Golden, ONE, PHI
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
T11 = [[18, 5], [5, 2]]
T11_REDUCED = [[2, 1], [1, 6]]


def field_gram() -> list[list[Golden]]:
    """The form [[-2phi,1],[1,2]] + E8 over Z[phi]."""

    e8 = e8_cartan_gram()
    result = [[Golden(0) for _ in range(10)] for _ in range(10)]
    result[0][0] = -2 * PHI
    result[0][1] = result[1][0] = ONE
    result[1][1] = Golden(2)
    for i in range(8):
        for j in range(8):
            result[i + 2][j + 2] = Golden(e8[i][j])
    return result


def trace_basis_change() -> list[list[int]]:
    """Take the trace lattice to U + (-T11) + E8(-1)^2."""

    size = 20
    change = [[0 for _ in range(size)] for _ in range(size)]
    # In the old trace basis (e0,phi*e0,e1,phi*e1), these vectors give
    # U followed by the negative binary lattice -T11.
    coupled_vectors = [
        (-3, -1, 1, -3),
        (-2, 0, 0, -1),
        (-4, -4, -3, -4),
        (-4, -1, 0, -3),
    ]
    for column, vector in enumerate(coupled_vectors):
        for row, entry in enumerate(vector):
            change[row][column] = entry

    # Split the E8 trace factor using e_i and phi*e_i-e_i.
    for i in range(8):
        old_e = 4 + 2 * i
        old_phi_e = old_e + 1
        change[old_e][4 + i] = 1
        change[old_e][12 + i] = -1
        change[old_phi_e][12 + i] = 1
    return change


def expected_trace_gram() -> list[list[int]]:
    e8_negative = [[-entry for entry in row] for row in e8_cartan_gram()]
    negative_t11 = [[-entry for entry in row] for row in T11]
    return block_diagonal(hyperbolic_plane(), negative_t11, e8_negative, e8_negative)


def k3_embedding() -> list[list[int]]:
    """Embed U + (-T11) + E8(-1)^2 into U^3 + E8(-1)^2."""

    embedding = [[0 for _ in range(20)] for _ in range(22)]
    # U maps identically to U_1.
    embedding[0][0] = 1
    embedding[1][1] = 1
    # For T11=[[2a,b],[b,2c]] with (a,b,c)=(9,5,1), embed -T11 in U_2+U_3.
    embedding[2][2] = 1
    embedding[3][2] = -9
    embedding[3][3] = -5
    embedding[4][3] = 1
    embedding[5][3] = -1
    # Negative E8 summands map identically.
    for i in range(16):
        embedding[6 + i][4 + i] = 1
    return embedding


def primitive_extension() -> list[list[int]]:
    embedding = k3_embedding()
    extension = [row + [0, 0] for row in embedding]
    extension[3][20] = 1  # f_2
    extension[5][21] = 1  # f_3
    return extension


def orthogonal_complement() -> list[list[int]]:
    """The positive T11 complement inside U_2 + U_3."""

    complement = [[0 for _ in range(2)] for _ in range(22)]
    complement[2][0] = 1
    complement[3][0] = 9
    complement[5][0] = 5
    complement[4][1] = 1
    complement[5][1] = 1
    return complement


def multiplication_by_phi(rank: int) -> list[list[int]]:
    result = [[0 for _ in range(2 * rank)] for _ in range(2 * rank)]
    for i in range(rank):
        first = 2 * i
        result[first + 1][first] = 1
        result[first][first + 1] = 1
        result[first + 1][first + 1] = 1
    return result


def matrix_product(left, right):
    return [
        [sum(left[i][k] * right[k][j] for k in range(len(right))) for j in range(len(right[0]))]
        for i in range(len(left))
    ]


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def golden_integer_matrix(matrix):
    result = []
    for row in matrix:
        integer_row = []
        for entry in row:
            entry = Golden.coerce(entry)
            if entry.b or entry.a.denominator != 1:
                raise ValueError("expected an integer-valued matrix")
            integer_row.append(entry.a.numerator)
        result.append(integer_row)
    return result


def analyze() -> dict[str, object]:
    h = field_gram()
    h_conjugate = [[entry.conjugate() for entry in row] for row in h]
    alpha = -PHI * INV_SQRT5
    trace = trace_gram(h, alpha)
    integer = integral_gram(trace)
    transform = trace_basis_change()
    transformed = change_gram_basis(integer, transform)
    expected = expected_trace_gram()

    k3 = k3_lattice_gram()
    embedding = k3_embedding()
    complement = orthogonal_complement()
    embedded_gram = rectangular_pairing(embedding, k3, embedding)
    complement_gram = rectangular_pairing(complement, k3, complement)
    cross_pairing = rectangular_pairing(embedding, k3, complement)

    phi_action = multiplication_by_phi(10)
    self_adjoint = matrix_product(transpose(phi_action), integer) == matrix_product(integer, phi_action)
    phi_polynomial = matrix_product(phi_action, phi_action)
    phi_plus_identity = [
        [phi_action[i][j] + (1 if i == j else 0) for j in range(20)]
        for i in range(20)
    ]

    signature_identity = inertia(h)
    signature_conjugate = inertia(h_conjugate)
    trace_signature = inertia(integer)
    return {
        "candidate_id": "coupled_discriminant11",
        "field": "Q(sqrt(5))",
        "field_gram": [[entry.to_json() for entry in row] for row in h],
        "field_determinant": str(determinant(h)),
        "field_signature_identity": signature_identity.__dict__ if hasattr(signature_identity, "__dict__") else {
            "positive": signature_identity.positive,
            "negative": signature_identity.negative,
            "zero": signature_identity.zero,
        },
        "field_signature_conjugate": {
            "positive": signature_conjugate.positive,
            "negative": signature_conjugate.negative,
            "zero": signature_conjugate.zero,
        },
        "trace_scalar": alpha.to_json(),
        "trace_scalar_is_codifferent": is_codifferent(alpha),
        "trace_even_integral": is_even_integral(trace),
        "trace_determinant": str(determinant(integer)),
        "trace_signature": {
            "positive": trace_signature.positive,
            "negative": trace_signature.negative,
            "zero": trace_signature.zero,
        },
        "trace_standard_lattice": "U + T11(-1) + E8(-1)^2",
        "T11": T11,
        "T11_reduced": T11_REDUCED,
        "trace_basis_change_determinant": str(determinant(transform)),
        "trace_basis_change_verified": transformed == expected,
        "discriminant_group": "Z/11Z",
        "phi_action_minimal_polynomial_verified": phi_polynomial == phi_plus_identity,
        "phi_action_self_adjoint_verified": self_adjoint,
        "primitive_k3_embedding_verified": (
            embedded_gram == expected and abs(determinant(primitive_extension()).a) == 1
        ),
        "primitive_extension_determinant": str(determinant(primitive_extension())),
        "k3_complement_gram": golden_integer_matrix(complement_gram),
        "k3_complement_verified": (
            complement_gram == T11
            and not any(entry for row in cross_pairing for entry in row)
        ),
        "geometric_status": "lattice_realizable_nef_marking_checked_separately",
    }


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2, sort_keys=True))
