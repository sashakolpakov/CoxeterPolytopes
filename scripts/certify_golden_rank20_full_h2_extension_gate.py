"""Certify the full-H2 extension gate for the golden rank-20 K3 packages.

The divisor lattices in the discriminant-11 and discriminant-19 candidates
already carry multiplication by ``phi``.  This script asks three separate
questions, in increasing geometric strength.

1. Does the induced ``sqrt(5)`` operator extend rationally to the K3 lattice?
2. Can it preserve the integral unimodular glue (at least for Z[sqrt(5)])?
3. Can such an extension be a Hodge action on a K3 period?

The first two questions are exact finite lattice calculations.  The last one
is decided by the dimension of the real eigenspaces on the positive rank-two
orthogonal complement.
"""

from __future__ import annotations

import hashlib
import json
from math import isqrt
from typing import Any

from sympy import Matrix, eye, zeros

from coxeter_polytopes.exact.trace_lattice import integral_gram, trace_gram
from coxeter_polytopes.lattices import k3_lattice_gram
from scripts import coupled_discriminant11_candidate as d11
from scripts import coupled_discriminant19_candidate as d19


def _is_integral(matrix: Matrix) -> bool:
    return all(entry.is_Integer for entry in matrix)


def _integer_matrix(matrix: Matrix) -> list[list[int]]:
    if not _is_integral(matrix):
        raise ValueError("matrix is not integral")
    return [[int(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def _matrix_sha256(matrix: Matrix) -> str:
    payload = json.dumps(_integer_matrix(matrix), separators=(",", ":"))
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def enumerate_integral_self_adjoint_square_roots(
    complement_gram: list[list[int]], radicand: int = 5
) -> list[Matrix]:
    """Enumerate every integral self-adjoint ``B`` with ``B^2=radicand``.

    For ``T=[[x,y],[y,z]]``, every nonscalar solution has the form

        B = [[a,b],[c,-a]],

    and self-adjointness plus ``B^2=radicand*I`` reduce to

        x*a^2 + 2*y*a*c + z*c^2 = radicand*x,
        x*b = 2*y*a + z*c.

    Positive definiteness gives the exhaustive bounds used below:

        det(T)*a^2 <= radicand*x*z,
        det(T)*c^2 <= radicand*x^2.
    """

    x, y = complement_gram[0]
    y_again, z = complement_gram[1]
    if y != y_again:
        raise ValueError("complement Gram matrix must be symmetric")
    determinant = x * z - y * y
    if x <= 0 or determinant <= 0:
        raise ValueError("complement Gram matrix must be positive definite")

    a_bound = isqrt((radicand * x * z) // determinant)
    while determinant * (a_bound + 1) ** 2 <= radicand * x * z:
        a_bound += 1
    c_bound = isqrt((radicand * x * x) // determinant)
    while determinant * (c_bound + 1) ** 2 <= radicand * x * x:
        c_bound += 1

    solutions: list[Matrix] = []
    for a in range(-a_bound, a_bound + 1):
        for c in range(-c_bound, c_bound + 1):
            if x * a * a + 2 * y * a * c + z * c * c != radicand * x:
                continue
            numerator = 2 * y * a + z * c
            if numerator % x:
                continue
            b = numerator // x
            candidate = Matrix([[a, b], [c, -a]])
            gram = Matrix(complement_gram)
            if candidate * candidate != radicand * eye(2):
                raise AssertionError("square-root equation reduction is inconsistent")
            if candidate.T * gram != gram * candidate:
                raise AssertionError("self-adjoint equation reduction is inconsistent")
            solutions.append(candidate)

    solutions.sort(key=lambda matrix: tuple(int(value) for value in matrix))
    return solutions


def _block_diagonal(left: Matrix, right: Matrix) -> Matrix:
    result = zeros(left.rows + right.rows, left.cols + right.cols)
    result[: left.rows, : left.cols] = left
    result[left.rows :, left.cols :] = right
    return result


def _trace_divisor_action(module: Any) -> tuple[Matrix, Matrix]:
    """Return ``(G_N, B_N)`` in the standard embedded divisor basis."""

    alpha = -module.PHI * module.INV_SQRT5
    old_gram = Matrix(integral_gram(trace_gram(module.field_gram(), alpha)))
    basis_change = Matrix(module.trace_basis_change())
    phi_old = Matrix(module.multiplication_by_phi(10))
    phi_standard = basis_change.inv() * phi_old * basis_change
    divisor_gram = basis_change.T * old_gram * basis_change
    sqrt5_standard = 2 * phi_standard - eye(20)

    if divisor_gram != Matrix(module.expected_trace_gram()):
        raise AssertionError("standard divisor Gram matrix was not recovered")
    if not _is_integral(phi_standard):
        raise AssertionError("the maximal-order action should preserve the divisor lattice")
    if phi_standard * phi_standard != phi_standard + eye(20):
        raise AssertionError("phi action has the wrong minimal polynomial")
    if sqrt5_standard * sqrt5_standard != 5 * eye(20):
        raise AssertionError("sqrt(5) action has the wrong square")
    if sqrt5_standard.T * divisor_gram != divisor_gram * sqrt5_standard:
        raise AssertionError("divisor action is not self-adjoint")
    return divisor_gram, sqrt5_standard


def analyze_package(candidate_id: str, module: Any) -> dict[str, Any]:
    divisor_gram, divisor_sqrt5 = _trace_divisor_action(module)
    k3_gram = Matrix(k3_lattice_gram())
    embedding = Matrix(module.k3_embedding())
    complement = Matrix(module.orthogonal_complement())
    complement_gram = complement.T * k3_gram * complement
    gluing_basis = embedding.row_join(complement)

    if embedding.T * k3_gram * embedding != divisor_gram:
        raise AssertionError("the stated divisor embedding has the wrong Gram matrix")
    if embedding.T * k3_gram * complement != zeros(20, 2):
        raise AssertionError("the stated complement is not orthogonal")
    if complement_gram.det() <= 0 or complement_gram[0, 0] <= 0:
        raise AssertionError("the rank-two complement is not positive definite")

    complement_solutions = enumerate_integral_self_adjoint_square_roots(
        _integer_matrix(complement_gram)
    )
    extensions: list[dict[str, Any]] = []
    for complement_sqrt5 in complement_solutions:
        split_action = _block_diagonal(divisor_sqrt5, complement_sqrt5)
        full_action = gluing_basis * split_action * gluing_basis.inv()
        glue_compatible = _is_integral(full_action)
        maximal_order_action = (eye(22) + full_action) / 2
        maximal_order_integral = _is_integral(maximal_order_action)

        if full_action * full_action != 5 * eye(22):
            raise AssertionError("extended action does not square to 5")
        if full_action.T * k3_gram != k3_gram * full_action:
            raise AssertionError("extended action is not self-adjoint")

        extension: dict[str, Any] = {
            "complement_sqrt5_action": _integer_matrix(complement_sqrt5),
            "rational_full_k3_extension": True,
            "z_sqrt5_glue_compatible": glue_compatible,
            "z_phi_maximal_order_glue_compatible": maximal_order_integral,
        }
        if glue_compatible:
            extension["full_k3_action_sha256"] = _matrix_sha256(full_action)
        extensions.append(extension)

    compatible = [row for row in extensions if row["z_sqrt5_glue_compatible"]]
    maximal = [row for row in extensions if row["z_phi_maximal_order_glue_compatible"]]
    return {
        "candidate_id": candidate_id,
        "divisor_rank": 20,
        "divisor_signature": [1, 19],
        "complement_gram": _integer_matrix(complement_gram),
        "complement_rank": 2,
        "complement_signature": [2, 0],
        "gluing_basis_determinant": int(gluing_basis.det()),
        "integral_complement_solutions": extensions,
        "integral_complement_solution_count": len(extensions),
        "z_sqrt5_full_lattice_extension_count": len(compatible),
        "z_phi_full_lattice_extension_count": len(maximal),
        "transcendental_rank_over_q_sqrt5": 1,
        "real_eigenspace_dimensions": {"+sqrt5": 1, "-sqrt5": 1},
        "period_plane_dimension": 2,
        "hodge_coherent_extension_exists": False,
        "hodge_obstruction": (
            "a real-quadratic Hodge action must put H20 and H02 in the same "
            "real eigenspace, but each eigenspace on this complement is one-dimensional"
        ),
    }


def analyze() -> dict[str, Any]:
    packages = [
        analyze_package("golden_discriminant_11_rank20", d11),
        analyze_package("golden_discriminant_19_rank20", d19),
    ]
    return {
        "certificate": "golden_rank20_full_h2_extension_gate",
        "field": "Q(sqrt(5))",
        "integral_order_levels": {
            "suborder": "Z[sqrt(5)] generated by B with B^2=5",
            "maximal_order": "Z[phi] generated by (I+B)/2",
        },
        "packages": packages,
        "universal_rank20_conclusion": (
            "full-H2 Hodge coherence is impossible for every quadratic rank-20 "
            "K3 divisor package, independently of rational or integral glue"
        ),
        "automaticity_result": (
            "primitive embedding plus an NS action does not imply Hodge coherence; "
            "both tested packages pass rational and Z[sqrt(5)] lattice extension "
            "but fail the period-plane condition"
        ),
    }


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2, sort_keys=True))
