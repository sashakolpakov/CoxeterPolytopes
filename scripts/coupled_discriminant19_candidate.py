"""Exact K3 certificates for the coupled discriminant-19 field form."""

from __future__ import annotations

import json

from coxeter_polytopes.exact.golden import Golden, PHI
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


T19 = [[2, 1], [1, 10]]


def field_gram() -> list[list[Golden]]:
    """Return ``[[-2*phi,phi],[phi,2]] + E8`` over ``Z[phi]``."""

    e8 = e8_cartan_gram()
    result = [[Golden(0) for _ in range(10)] for _ in range(10)]
    result[0][0] = -2 * PHI
    result[0][1] = result[1][0] = PHI
    result[1][1] = Golden(2)
    for i in range(8):
        for j in range(8):
            result[i + 2][j + 2] = Golden(e8[i][j])
    return result


def trace_basis_change() -> list[list[int]]:
    """Take the trace lattice to ``U + T19(-1) + E8(-1)^2``.

    On the coupled rank-four block, the columns are

    ``(-1,1,-1,0), (-2,2,-2,1), (2,-2,3,-1), (7,-6,4,-1)``.

    Their determinant is one.  The first two span a hyperbolic plane and the
    last two have negative Gram matrix ``-T19``.
    """

    size = 20
    change = [[0 for _ in range(size)] for _ in range(size)]
    coupled_vectors = [
        (-1, 1, -1, 0),
        (-2, 2, -2, 1),
        (2, -2, 3, -1),
        (7, -6, 4, -1),
    ]
    for column, vector in enumerate(coupled_vectors):
        for row, entry in enumerate(vector):
            change[row][column] = entry

    # In each E8 coordinate pair, replace (e_i, phi*e_i) by
    # (e_i, phi*e_i-e_i), then group the two E8 copies.
    for i in range(8):
        old_e = 4 + 2 * i
        old_phi_e = old_e + 1
        change[old_e][4 + i] = 1
        change[old_e][12 + i] = -1
        change[old_phi_e][12 + i] = 1
    return change


def expected_trace_gram() -> list[list[int]]:
    e8_negative = [[-entry for entry in row] for row in e8_cartan_gram()]
    negative_t19 = [[-entry for entry in row] for row in T19]
    return block_diagonal(hyperbolic_plane(), negative_t19, e8_negative, e8_negative)


def k3_embedding() -> list[list[int]]:
    """Embed ``U + T19(-1) + E8(-1)^2`` in the K3 lattice."""

    embedding = [[0 for _ in range(20)] for _ in range(22)]
    # U maps identically to U_1.
    embedding[0][0] = 1
    embedding[1][1] = 1
    # For T19=[[2a,b],[b,2c]], (a,b,c)=(1,1,5).  In U_2+U_3,
    # z_1=e_2-a*f_2 and z_2=-b*f_2+e_3-c*f_3 have Gram -T19.
    embedding[2][2] = 1
    embedding[3][2] = -1
    embedding[3][3] = -1
    embedding[4][3] = 1
    embedding[5][3] = -5
    # Negative E8 summands map identically.
    for i in range(16):
        embedding[6 + i][4 + i] = 1
    return embedding


def primitive_extension() -> list[list[int]]:
    """Extend the embedding columns to an integral basis of the K3 lattice."""

    embedding = k3_embedding()
    extension = [row + [0, 0] for row in embedding]
    extension[3][20] = 1  # f_2
    extension[5][21] = 1  # f_3
    return extension


def orthogonal_complement() -> list[list[int]]:
    """Return a positive ``T19`` complement inside ``U_2 + U_3``."""

    complement = [[0 for _ in range(2)] for _ in range(22)]
    # w_1=e_2+f_2+f_3, w_2=e_3+5*f_3.
    complement[2][0] = 1
    complement[3][0] = 1
    complement[5][0] = 1
    complement[4][1] = 1
    complement[5][1] = 5
    return complement


def multiplication_by_phi(rank: int) -> list[list[int]]:
    """Multiplication by phi in each interleaved ``(e,phi*e)`` pair."""

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


def _signature(inertia_result) -> dict[str, int]:
    return {
        "positive": inertia_result.positive,
        "negative": inertia_result.negative,
        "zero": inertia_result.zero,
    }


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
    phi_polynomial = matrix_product(phi_action, phi_action)
    phi_plus_identity = [
        [phi_action[i][j] + (1 if i == j else 0) for j in range(20)]
        for i in range(20)
    ]
    self_adjoint = matrix_product(transpose(phi_action), integer) == matrix_product(
        integer, phi_action
    )

    field_determinant = determinant(h)
    signature_identity = inertia(h)
    signature_conjugate = inertia(h_conjugate)
    trace_signature = inertia(integer)
    control = [Golden(2), -PHI, *(Golden(0) for _ in range(8))]
    control_square = sum(
        (control[i] * h[i][j] * control[j] for i in range(10) for j in range(10)),
        start=Golden(0),
    )
    identity_eigenspace_square = alpha * control_square
    return {
        "candidate_id": "coupled_discriminant19",
        "field": "Q(sqrt(5))",
        "field_gram": [[entry.to_json() for entry in row] for row in h],
        "field_determinant": str(field_determinant),
        "field_determinant_norm": str(field_determinant.norm()),
        "field_signature_identity": _signature(signature_identity),
        "field_signature_conjugate": _signature(signature_conjugate),
        "trace_scalar": alpha.to_json(),
        "trace_scalar_is_codifferent": is_codifferent(alpha),
        "trace_even_integral": is_even_integral(trace),
        "trace_determinant": str(determinant(integer)),
        "trace_signature": _signature(trace_signature),
        "trace_standard_lattice": "U + T19(-1) + E8(-1)^2",
        "T19": T19,
        "trace_basis_change_determinant": str(determinant(transform)),
        "trace_basis_change_verified": transformed == expected,
        "discriminant_group": "Z/19Z",
        "phi_action_minimal_polynomial_verified": phi_polynomial == phi_plus_identity,
        "phi_action_self_adjoint_verified": self_adjoint,
        "primitive_k3_embedding_verified": (
            embedded_gram == expected and abs(determinant(primitive_extension()).a) == 1
        ),
        "primitive_extension_determinant": str(determinant(primitive_extension())),
        "k3_complement_gram": golden_integer_matrix(complement_gram),
        "k3_complement_verified": (
            complement_gram == T19
            and not any(entry for row in cross_pairing for entry in row)
        ),
        "published_k3_match": {
            "surface": "Shioda's unique elliptic K3 surface S19 with an I19 fiber",
            "transcendental_lattice": T19,
            "neron_severi_lattice": "U + A18(-1)",
            "doi": "10.1016/j.crma.2003.07.007",
        },
        "eigen_positive_control": {
            "field_vector": [entry.to_json() for entry in control],
            "field_form_square": control_square.to_json(),
            "identity_eigenspace_square": identity_eigenspace_square.to_json(),
            "restriction_of_scalars_trace_square": str(
                identity_eigenspace_square.trace()
            ),
        },
        "fixed_i19_fibration_marking": {
            "status": "not_available",
            "missing_datum": (
                "oriented integral isometry to Shioda's "
                "fiber-section-component basis"
            ),
        },
        "geometric_status": (
            "singular_k3_carrier_and_eigen_positive_cone_constructed_"
            "fixed_eigen_nef_and_reflectivity_unresolved"
        ),
    }


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2, sort_keys=True))
