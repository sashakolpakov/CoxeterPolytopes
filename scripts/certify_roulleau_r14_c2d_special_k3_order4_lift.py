#!/usr/bin/env python3
"""Certify a genuine order-four K3 lift of the R14 C2D section.

The Roulleau rank-14 lattice is

    N = U + E8(-1) + A1(-1)^4.

Its double-transposition wall symmetry is glued equivariantly to an
order-four isometry of

    T = U + U + A1(-1)^4.

The graph glue is the even unimodular K3 lattice.  An explicit algebraic
period in the i-eigenspace has no rational orthogonal vector in T, so its
Picard lattice is exactly N.  The archived symmetry preserves a strict
interior point of Roulleau's complete nef chamber.  Strong Torelli therefore
turns the glued Hodge isometry into a genuine purely non-symplectic order-four
K3 automorphism.  Its fixed nef-chamber section is the already certified
finite-volume noncompact H11/19 Coxeter polytope.

Optionally, ``--brandhorst-hofmann-order4`` replays the full published
Brandhorst--Hofmann order-four table and identifies the construction with the
unique compatible row 0.4.3.4.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import re
from typing import Any, Sequence

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.domains import ZZ as POLY_ZZ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT / "data" / "enumerations" / "roulleau_r14_ue8a1x4_all_subgroup_folds.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "roulleau_r14_c2d_special_k3_order4_lift.json"
)

EXPECTED_SOURCE_SHA256 = (
    "5c39aa26fc4e2c7f2af5bb08815532ec9cce7b2222279c65f257e383f356f4ef"
)

BH_DOI = "10.5281/zenodo.10410278"
BH_VERSION = "2.1"
BH_ORDER4_SHA256 = "d5a61603c1539278e8411e188656c9417662336ff4282b714e939c40d38c94e8"
BH_ORDER4_ROW_COUNT = 79
BH_MATCHING_ROWS = (
    {
        "id": "0.4.3.4",
        "row_one_based": 28,
        "row_sha256_with_newline": (
            "303fd69804142bc57a27799551fa96de8dc10bfe764449e48e1ffe5d55cd7cf3"
        ),
        "fixed_determinant": -16,
        "fixed_smith": [*[1] * 10, 4, 4],
    },
    {
        "id": "0.4.3.12",
        "row_one_based": 36,
        "row_sha256_with_newline": (
            "130fb1503a43c1de96e627517d66849011731e62bf72a49873ac3f3255d328ad"
        ),
        "fixed_determinant": -64,
        "fixed_smith": [*[1] * 8, 2, 2, 4, 4],
    },
)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def matrix_payload(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(value) for value in matrix.row(index)] for index in range(matrix.rows)]


def matrix_digest(matrix: sp.Matrix) -> str:
    payload = json.dumps(matrix_payload(matrix), separators=(",", ":")).encode()
    return sha256_bytes(payload)


def _domain_matrix(matrix: sp.Matrix) -> DomainMatrix:
    if any(sp.Rational(value).q != 1 for value in matrix):
        raise AssertionError("an integral Smith computation received fractions")
    return DomainMatrix.from_Matrix(sp.Matrix(matrix)).convert_to(POLY_ZZ)


def saturated_row_kernel(equations: sp.Matrix) -> sp.Matrix:
    """Return a saturated Z-basis for row vectors x with x*equations=0."""

    transposed = equations.T
    denominators = [sp.Rational(value).q for value in transposed]
    denominator = int(sp.ilcm(*denominators)) if denominators else 1
    integral = sp.Matrix(denominator * transposed)
    diagonal, _, right = smith_normal_decomp(_domain_matrix(integral))
    diagonal_matrix = diagonal.to_Matrix()
    right_matrix = right.to_Matrix()
    zero_columns = [
        column
        for column in range(integral.cols)
        if all(
            diagonal_matrix[row, column] == 0
            for row in range(diagonal_matrix.rows)
        )
    ]
    basis = right_matrix[:, zero_columns].T
    if basis * equations != sp.zeros(basis.rows, equations.cols):
        raise AssertionError("Smith-kernel replay failed")
    return basis


def induced_row_action(basis: sp.Matrix, ambient_action: sp.Matrix) -> sp.Matrix:
    pivots = basis.rref()[1]
    square = basis[:, list(pivots)]
    action = (basis * ambient_action)[:, list(pivots)] * square.inv()
    if not (
        all(sp.Rational(value).q == 1 for value in action)
        and action * basis == basis * ambient_action
    ):
        raise AssertionError("induced row action is not integral")
    return action


def same_row_lattice(left: sp.Matrix, right: sp.Matrix) -> bool:
    if left.shape != right.shape or left.rank() != left.rows:
        return False
    pivots = right.rref()[1]
    change = left[:, list(pivots)] * right[:, list(pivots)].inv()
    return (
        all(sp.Rational(value).q == 1 for value in change)
        and abs(int(change.det())) == 1
        and change * right == left
    )


def smith_diagonal(matrix: sp.Matrix) -> list[int]:
    normal = smith_normal_form(matrix, domain=ZZ)
    return [abs(int(normal[index, index])) for index in range(min(normal.shape))]


def gf2_rank(rows: Sequence[Sequence[int]]) -> int:
    work = [[int(value) & 1 for value in row] for row in rows]
    if not work:
        return 0
    rank = 0
    for column in range(len(work[0])):
        pivot = next(
            (index for index in range(rank, len(work)) if work[index][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for index in range(len(work)):
            if index != rank and work[index][column]:
                work[index] = [
                    left ^ right for left, right in zip(work[index], work[rank])
                ]
        rank += 1
    return rank


def two_elementary_basis(
    gram: sp.Matrix, expected_length: int
) -> tuple[tuple[int, ...], ...]:
    """Basis bits v for the discriminant half-vectors v/2."""

    size = gram.rows
    work = [[int(gram[i, j]) & 1 for j in range(size)] for i in range(size)]
    pivots: list[int] = []
    rank = 0
    for column in range(size):
        pivot = next(
            (index for index in range(rank, size) if work[index][column]), None
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for index in range(size):
            if index != rank and work[index][column]:
                work[index] = [
                    left ^ right for left, right in zip(work[index], work[rank])
                ]
        pivots.append(column)
        rank += 1
    free = [column for column in range(size) if column not in pivots]
    basis = []
    for free_column in free:
        vector = [0] * size
        vector[free_column] = 1
        for row, pivot_column in reversed(list(enumerate(pivots))):
            vector[pivot_column] = sum(
                work[row][column] * vector[column] for column in free
            ) & 1
        basis.append(tuple(vector))
    if len(basis) != expected_length or gf2_rank(basis) != expected_length:
        raise AssertionError("unexpected two-elementary discriminant length")
    return tuple(basis)


def gf2_coordinates(
    vector: Sequence[int], basis: Sequence[Sequence[int]]
) -> tuple[int, ...]:
    for coefficients in itertools.product(range(2), repeat=len(basis)):
        candidate = tuple(
            sum(
                coefficients[index] * basis[index][column]
                for index in range(len(basis))
            )
            & 1
            for column in range(len(vector))
        )
        if candidate == tuple(int(value) & 1 for value in vector):
            return tuple(coefficients)
    raise AssertionError("vector is not in the displayed GF(2) span")


def discriminant_action(
    action: sp.Matrix, basis: Sequence[Sequence[int]]
) -> sp.Matrix:
    rows = []
    for vector in basis:
        image = tuple(
            sum(
                vector[index] * int(action[index, column])
                for index in range(action.rows)
            )
            & 1
            for column in range(action.cols)
        )
        rows.append(gf2_coordinates(image, basis))
    return sp.Matrix(rows)


def q_times_four_mod_eight(
    coordinates: Sequence[int],
    half_vector_basis: Sequence[Sequence[int]],
    gram: sp.Matrix,
) -> int:
    """Encode q in Q/2Z as 4q modulo eight."""

    bits = tuple(
        sum(
            coordinates[index] * half_vector_basis[index][column]
            for index in range(len(half_vector_basis))
        )
        & 1
        for column in range(gram.rows)
    )
    return (
        sum(
            bits[i] * int(gram[i, j]) * bits[j]
            for i in range(gram.rows)
            for j in range(gram.cols)
        )
        % 8
    )


def canonical_r14_basis(n_gram: sp.Matrix) -> sp.Matrix:
    """Return a unimodular U+E8(-1)+A1(-1)^4 basis of the source N."""

    e8_gram = n_gram[:8, :8]
    projection = sp.eye(14)
    projection[8:, :8] = -n_gram[8:, :8] * e8_gram.inv()
    split_gram = projection * n_gram * projection.T
    six_gram = split_gram[8:, 8:]

    u = sp.Matrix([[0, -1, 0, 0, 0, 0]])
    v = sp.Matrix([[0, -1, -1, 0, 0, 0]])
    complement_rows = [u, v]
    for column in (0, 3, 4, 5):
        vector = sp.zeros(1, 6)
        vector[0, column] = 1
        vector -= (vector * six_gram * v.T)[0] * u
        vector -= (vector * six_gram * u.T)[0] * v
        complement_rows.append(vector)
    complement = sp.Matrix.vstack(*complement_rows)
    roots = sp.eye(6)
    roots[2:, 2:] = sp.Matrix(
        [
            [1, -1, -1, -1],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ]
    )
    e8_u_a1_basis = sp.diag(sp.eye(8), roots * complement) * projection
    order = [8, 9, *range(8), 10, 13, 11, 12]
    canonical = e8_u_a1_basis[order, :]
    if abs(int(canonical.det())) != 1:
        raise AssertionError("the displayed R14 decomposition is not unimodular")
    return canonical


def glue_k3_lattice(
    n_gram: sp.Matrix,
    n_action: sp.Matrix,
    n_discriminant_basis: Sequence[Sequence[int]],
    t_gram: sp.Matrix,
    t_action: sp.Matrix,
    t_discriminant_basis: Sequence[Sequence[int]],
    identification: sp.Matrix,
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    target_bits = identification * sp.Matrix(t_discriminant_basis)
    glue_rows = sp.Matrix(
        [
            [
                *n_discriminant_basis[index],
                *map(int, target_bits.row(index)),
            ]
            for index in range(len(n_discriminant_basis))
        ]
    )
    doubled_rows = (2 * sp.eye(22)).col_join(glue_rows)
    k3_basis = hermite_normal_form(doubled_rows.T).T / 2
    k3_gram = sp.simplify(k3_basis * sp.diag(n_gram, t_gram) * k3_basis.T)
    k3_action = sp.simplify(
        k3_basis * sp.diag(n_action, t_action) * k3_basis.inv()
    )
    if not (
        all(sp.Rational(value).q == 1 for value in k3_gram)
        and all(k3_gram[index, index] % 2 == 0 for index in range(22))
        and k3_gram.det() == -1
        and exact_inertia(k3_gram) == (3, 19, 0)
        and all(sp.Rational(value).q == 1 for value in k3_action)
        and k3_action * k3_gram * k3_action.T == k3_gram
    ):
        raise AssertionError("the equivariant graph glue is not a K3 lattice action")
    return k3_basis, k3_gram, k3_action


def _parse_bh_row(line: str) -> dict[str, Any]:
    exact = re.sub(r"(?<![\w.])(\d+)/(\d+)", r"Fraction(\1,\2)", line)
    pairs = eval(exact, {"__builtins__": {}, "Fraction": Fraction})
    return dict(pairs)


def _rational(value: Any) -> sp.Rational:
    if isinstance(value, Fraction):
        return sp.Rational(value.numerator, value.denominator)
    return sp.Rational(value)


def replay_brandhorst_hofmann_order4(
    path: Path, canonical_gram: sp.Matrix, canonical_action: sp.Matrix
) -> dict[str, Any]:
    """Replay the complete table and the exact pair isometry for row 0.4.3.4."""

    raw = path.read_bytes()
    if sha256_bytes(raw) != BH_ORDER4_SHA256:
        raise AssertionError("the pinned Brandhorst--Hofmann order-four file changed")
    lines = raw.decode("utf-8").splitlines()
    if len(lines) != BH_ORDER4_ROW_COUNT:
        raise AssertionError("the Brandhorst--Hofmann order-four row count changed")

    hits: list[dict[str, Any]] = []
    target_data: dict[str, Any] | None = None
    for row_number, line in enumerate(lines, 1):
        data = _parse_bh_row(line)
        ambient_action = sp.Matrix(22, 22, data["distinguished_generator"])
        fixed_rank = 22 - (ambient_action - sp.eye(22)).rank()
        anti_rank = 22 - (ambient_action + sp.eye(22)).rank()
        gaussian_rank = (22 - fixed_rank - anti_rank) // 2
        if (fixed_rank, anti_rank, gaussian_rank) != (12, 2, 4):
            continue

        basis = sp.Matrix(22, 22, list(map(_rational, data["basis"])))
        ambient_gram = sp.Matrix(
            22, 22, list(map(_rational, data["inner_product_matrix"]))
        )
        gram = sp.simplify(basis * ambient_gram * basis.T)
        action = sp.simplify(basis * ambient_action * basis.inv())
        fixed_basis = saturated_row_kernel(action - sp.eye(22))
        fixed_gram = fixed_basis * gram * fixed_basis.T
        record = {
            "id": data["id"],
            "row_one_based": row_number,
            "row_sha256_with_newline": sha256_bytes((line + "\n").encode()),
            "fixed_determinant": int(fixed_gram.det()),
            "fixed_smith": smith_diagonal(fixed_gram),
        }
        hits.append(record)
        if data["id"] == "0.4.3.4":
            target_data = {
                "gram": gram,
                "action": action,
            }

    if tuple(hits) != BH_MATCHING_ROWS or target_data is None:
        raise AssertionError("the two classified (12,2,4) rows changed")

    gram = target_data["gram"]
    action = target_data["action"]
    if not (
        gram.det() == -1
        and exact_inertia(gram) == (3, 19, 0)
        and action**4 == sp.eye(22)
        and action * gram * action.T == gram
    ):
        raise AssertionError("the target Brandhorst--Hofmann K3 action changed")
    ns_basis = saturated_row_kernel(action**2 - sp.eye(22))
    ns_gram = ns_basis * gram * ns_basis.T
    ns_action = induced_row_action(ns_basis, action)

    rows = []
    for index in (0, 11, *range(3, 11)):
        vector = sp.zeros(1, 14)
        vector[0, index] = 1
        rows.append(vector)
    for positive, negative in ((13, None), (1, 13), (12, None), (2, 12)):
        vector = sp.zeros(1, 14)
        vector[0, positive] = 1
        if negative is not None:
            vector[0, negative] = -1
        rows.append(vector)
    raw_canonical = sp.Matrix.vstack(*rows)
    e8_isometry = sp.Matrix(
        [
            [1, 0, 0, 0, 0, 0, 0, 0],
            [-2, -1, 1, 2, -1, -3, -2, -2],
            [1, 0, -1, -2, 1, 2, 1, 1],
            [0, 0, 0, 1, -1, -1, 0, 0],
            [0, 0, 1, 1, 0, 0, 0, 0],
            [0, 1, -1, -1, 0, 1, 0, 1],
            [0, 0, 0, 0, 0, 0, 1, -1],
            [0, 0, 0, 0, 1, 0, 0, 1],
        ]
    )
    bh_canonical = sp.diag(sp.eye(2), e8_isometry, sp.eye(4)) * raw_canonical
    if not (
        abs(int(bh_canonical.det())) == 1
        and bh_canonical * ns_gram * bh_canonical.T == canonical_gram
        and bh_canonical * ns_action * bh_canonical.inv() == canonical_action
    ):
        raise AssertionError("the exact Brandhorst--Hofmann/R14 pair match failed")

    return {
        "member_sha256": BH_ORDER4_SHA256,
        "row_count": len(lines),
        "matching_character_rows": hits,
        "target_pair_exactly_integrally_conjugate_to_R14_C2D": True,
        "target_canonical_basis_sha256": matrix_digest(bh_canonical),
    }


def elliptic_normal_form_audit() -> dict[str, Any]:
    """Audit an explicit algebraic search chart carrying the order-four map."""

    t, x, y = sp.symbols("t x y")
    a1, a3, p, q = sp.symbols("a1 a3 p q")
    a = a1 * t + a3 * t**3
    b = (t**2 - p**2) * (t**2 - q**2)
    equation = y**2 - (x**3 + a * x**2 + b * x)
    transformed = equation.subs(
        {t: -t, x: -x, y: sp.I * y}, simultaneous=True
    )
    if sp.expand(transformed + equation) != 0:
        raise AssertionError("the elliptic order-four normal form is not invariant")

    example_a = t**3 + t
    example_b = (t**2 - 1) * (t**2 - 4)
    residual = sp.expand(example_a**2 - 4 * example_b)
    if not (
        sp.degree(residual, t) == 6
        and sp.gcd(residual, sp.diff(residual, t)) == 1
        and sp.gcd(residual, example_b) == 1
        and all(example_a.subs(t, value) != 0 for value in (1, -1, 2, -2))
    ):
        raise AssertionError("the displayed elliptic member is not generic enough")

    return {
        "family": (
            "y^2=x^3+(a1*t+a3*t^3)*x^2+"
            "(t^2-p^2)*(t^2-q^2)*x"
        ),
        "order_four_map": "sigma:(t,x,y)->(-t,-x,i*y)",
        "square": "sigma^2:(t,x,y)->(t,x,-y), elliptic inversion",
        "equation_pullback": "sigma^*(F)=-F",
        "two_torsion_section": "(x,y)=(0,0)",
        "discriminant_up_to_nonzero_constant": "B(t)^2*(A(t)^2-4*B(t))",
        "c4_up_to_nonzero_constant": "A(t)^2-3*B(t)",
        "explicit_member": (
            "y^2=x^3+(t^3+t)*x^2+(t^2-1)*(t^2-4)*x"
        ),
        "explicit_residual_factor": str(residual),
        "explicit_residual_discriminant": str(sp.discriminant(residual, t)),
        "generic_fiber_configuration": ["I4* at infinity", "4 I2", "6 I1"],
        "euler_number_check": "10+4*2+6=24",
        "lattice_chart_reason": (
            "the trivial lattice U+D8(-1)+A1(-1)^4 has determinant -64; "
            "the visible 2-torsion gives its index-two even overlattice of "
            "determinant -16, compatible with the certified R14 lattice"
        ),
        "scope_caution": (
            "the symbolic equation is an algebro-geometric normal-form chart; "
            "the exact NS=N and Torelli claims are certified independently by "
            "the K3-lattice period construction above, not by point counting "
            "on the displayed numerical member"
        ),
    }


def build_certificate(
    brandhorst_hofmann_order4: Path | None = None,
) -> dict[str, Any]:
    source_hash = sha256(SOURCE)
    if source_hash != EXPECTED_SOURCE_SHA256:
        raise AssertionError("the pinned R14/27 source certificate changed")
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    case = next(
        row for row in source["cases"] if row["id"] == "R14-UE8A1x4-C2D-H11-19"
    )
    vinberg = case["geometry"]["vinberg"]
    if not (
        source["all_checks_pass"]
        and case["all_checks_pass"]
        and all(case["checks"].values())
        and vinberg["classification"]
        == "finite_volume_noncompact_Coxeter_polytope"
        and vinberg["finite_volume"] is True
        and vinberg["compact"] is False
        and vinberg["dimension"] == 11
        and vinberg["facet_count"] == 19
        and vinberg["elliptic_edge_count"] == 795
        and vinberg["finite_vertex_count"] == 114
        and vinberg["ideal_vertex_count"] == 13
        and vinberg["bad_edge_count"] == 0
        and vinberg["bad_parabolic_neighborhood_count"] == 0
    ):
        raise AssertionError("the upstream R14 C2D Vinberg theorem changed")

    lane = source["source_reconstruction"]["authoritative_lane"]
    n_gram = sp.Matrix(lane["basis_gram"])
    wall_columns = sp.Matrix(lane["wall_coordinate_columns"])
    permutation = tuple(
        value - 1 for value in case["subgroup"]["elements_one_based"][1]
    )
    column_action = sp.Matrix.hstack(
        *(wall_columns[:, permutation[index]] for index in range(14))
    )
    n_action = column_action.T
    archived_fixed_basis = sp.Matrix(
        case["fixed_lattice"]["saturated_basis_columns_in_NS_basis"]
    ).T
    computed_fixed_basis = saturated_row_kernel(n_action - sp.eye(14))
    fixed_gram = archived_fixed_basis * n_gram * archived_fixed_basis.T
    wall_rows = wall_columns.T
    wall_lookup = {tuple(wall_rows.row(index)) for index in range(27)}
    if not (
        n_gram.det() == -16
        and exact_inertia(n_gram) == (1, 13, 0)
        and n_action**2 == sp.eye(14)
        and n_action * n_gram * n_action.T == n_gram
        and sp.factor(n_action.charpoly().as_expr())
        == (sp.Symbol("lambda") - 1) ** 12 * (sp.Symbol("lambda") + 1) ** 2
        and same_row_lattice(archived_fixed_basis, computed_fixed_basis)
        and fixed_gram.det() == -16
        and smith_diagonal(fixed_gram) == [*[1] * 10, 4, 4]
        and {tuple(wall_rows.row(index) * n_action) for index in range(27)}
        == wall_lookup
    ):
        raise AssertionError("the archived R14 C2D Picard action changed")

    canonical_basis = canonical_r14_basis(n_gram)
    canonical_gram = canonical_basis * n_gram * canonical_basis.T
    canonical_action = canonical_basis * n_action * canonical_basis.inv()
    u = sp.Matrix([[0, 1], [1, 0]])
    e8_gram = n_gram[:8, :8]
    a1_swap = sp.Matrix([[0, 1], [1, 0]])
    expected_canonical_gram = sp.diag(u, e8_gram, -2 * sp.eye(4))
    expected_canonical_action = sp.diag(
        sp.eye(10), a1_swap, a1_swap
    )
    if not (
        canonical_gram == expected_canonical_gram
        and canonical_action == expected_canonical_action
    ):
        raise AssertionError("the source action is not the canonical double transposition")

    anti_basis = saturated_row_kernel(n_action + sp.eye(14))
    canonical_anti = sp.Matrix.vstack(
        canonical_basis.row(10) - canonical_basis.row(11),
        canonical_basis.row(12) - canonical_basis.row(13),
    )
    anti_gram = anti_basis * n_gram * anti_basis.T
    if not (
        same_row_lattice(anti_basis, canonical_anti)
        and canonical_anti * n_gram * canonical_anti.T == -4 * sp.eye(2)
        and anti_gram.det() == 16
        and exact_inertia(anti_gram) == (0, 2, 0)
    ):
        raise AssertionError("the algebraic anti-invariant lattice is not <-4>^2")

    t_gram = sp.diag(u, u, -2 * sp.eye(4))
    rotation = sp.Matrix(
        [
            [0, 0, 0, -1],
            [0, 0, -1, 0],
            [0, 1, 0, 0],
            [1, 0, 0, 0],
        ]
    )
    t_action = sp.diag(rotation, rotation)
    if not (
        t_gram.det() == 16
        and exact_inertia(t_gram) == (2, 6, 0)
        and t_action * t_gram * t_action.T == t_gram
        and t_action**2 == -sp.eye(8)
    ):
        raise AssertionError("the transcendental order-four action changed")

    n_discriminant_basis = two_elementary_basis(n_gram, 4)
    t_discriminant_basis = two_elementary_basis(t_gram, 4)
    identification = sp.Matrix(
        [
            [1, 1, 1, 1],
            [1, 0, 1, 1],
            [1, 1, 0, 1],
            [1, 1, 1, 0],
        ]
    )
    n_discriminant_action = discriminant_action(n_action, n_discriminant_basis)
    t_discriminant_action = discriminant_action(t_action, t_discriminant_basis)
    if not (
        int(identification.det()) % 2 == 1
        and (
            n_discriminant_action * identification
            - identification * t_discriminant_action
        ).applyfunc(lambda value: int(value) % 2)
        == sp.zeros(4)
    ):
        raise AssertionError("the discriminant actions are not equivariant")

    q_table = []
    for coordinates in itertools.product(range(2), repeat=4):
        target = tuple(
            int(value) % 2
            for value in (sp.Matrix([coordinates]) * identification).row(0)
        )
        q_n = q_times_four_mod_eight(
            coordinates, n_discriminant_basis, n_gram
        )
        q_t = q_times_four_mod_eight(target, t_discriminant_basis, t_gram)
        if (q_n + q_t) % 8:
            raise AssertionError("the displayed discriminant graph is not isotropic")
        q_table.append(
            {
                "source": list(coordinates),
                "target": list(target),
                "4q_N_mod_8": q_n,
                "4q_T_mod_8": q_t,
            }
        )

    k3_basis, k3_gram, k3_action = glue_k3_lattice(
        n_gram,
        n_action,
        n_discriminant_basis,
        t_gram,
        t_action,
        t_discriminant_basis,
        identification,
    )
    lambda_symbol = sp.Symbol("lambda")
    expected_charpoly = (
        (lambda_symbol - 1) ** 12
        * (lambda_symbol + 1) ** 2
        * (lambda_symbol**2 + 1) ** 4
    )
    if not (
        k3_action**4 == sp.eye(22)
        and k3_action**2 != sp.eye(22)
        and sp.factor(k3_action.charpoly().as_expr()) == expected_charpoly
    ):
        raise AssertionError("the glued K3 action is not the required order-four action")

    n_embedding = sp.Matrix.hstack(sp.eye(14), sp.zeros(14, 8)) * k3_basis.inv()
    primitive_n = saturated_row_kernel(k3_basis[:, 14:22])
    square_fixed = saturated_row_kernel(k3_action**2 - sp.eye(22))
    global_fixed = saturated_row_kernel(k3_action - sp.eye(22))
    archived_fixed_in_k3 = archived_fixed_basis * n_embedding
    if not (
        same_row_lattice(primitive_n, n_embedding)
        and same_row_lattice(square_fixed, n_embedding)
        and same_row_lattice(global_fixed, archived_fixed_in_k3)
        and (global_fixed * k3_gram * global_fixed.T).det() == -16
        and smith_diagonal(global_fixed * k3_gram * global_fixed.T)
        == [*[1] * 10, 4, 4]
    ):
        raise AssertionError("the primitive Picard and fixed lattices changed after glue")

    eigenbasis = sp.Matrix(
        [
            [-sp.I, -sp.I, 1, 1, 0, 0, 0, 0],
            [sp.I, -sp.I, 1, -1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, -sp.I, 1, 0],
            [0, 0, 0, 0, -sp.I, 0, 0, 1],
        ]
    )
    hermitian_gram = sp.simplify(
        eigenbasis * t_gram * sp.conjugate(eigenbasis).T
    )
    if not (
        eigenbasis * t_action == sp.I * eigenbasis
        and eigenbasis * t_gram * eigenbasis.T == sp.zeros(4)
        and hermitian_gram == sp.diag(4, -4, -4, -4)
        and sp.Matrix.vstack(eigenbasis, sp.conjugate(eigenbasis)).rank() == 8
    ):
        raise AssertionError("the displayed positive i-eigenspace basis changed")
    global_eigenbasis = (
        sp.Matrix.hstack(sp.zeros(4, 14), eigenbasis) * k3_basis.inv()
    )
    if not (
        all(
            sp.Rational(sp.re(value)).q == 1
            and sp.Rational(sp.im(value)).q == 1
            for value in global_eigenbasis
        )
        and (
            global_eigenbasis * k3_action - sp.I * global_eigenbasis
        ).applyfunc(sp.simplify)
        == sp.zeros(4, 22)
        and global_eigenbasis * k3_gram * global_eigenbasis.T == sp.zeros(4)
        and sp.simplify(
            global_eigenbasis
            * k3_gram
            * sp.conjugate(global_eigenbasis).T
        )
        == hermitian_gram
    ):
        raise AssertionError("the i-eigenspace did not survive the K3 graph glue")

    alpha = sp.Symbol("alpha", real=True, positive=True)
    polynomial = sp.Poly(alpha**4 - 2, alpha, domain="QQ_I")
    if not polynomial.is_irreducible:
        raise AssertionError("alpha^4-2 is no longer irreducible over Q(i)")
    coefficients = sp.Matrix([[1, alpha / 10, alpha**2 / 100, alpha**3 / 1000]])
    period = coefficients * eigenbasis
    if not (
        (period * t_action - sp.I * period).applyfunc(sp.simplify)
        == sp.zeros(1, 8)
        and sp.simplify((period * t_gram * period.T)[0]) == 0
    ):
        raise AssertionError("the explicit algebraic period is not isotropic/equivariant")
    positive_lower_bound = Fraction(244949, 62500)

    d132 = sp.Matrix([2, 5, 9, 4, 10, 12, 15, 19, 9, 15, 12, 5, 5, 0])
    invariant_ample = d132 + column_action * d132
    ample_pairings = [
        sp.Rational(value) for value in list(invariant_ample.T * n_gram * wall_columns)
    ]
    global_ample = (
        sp.Matrix.hstack(invariant_ample.T, sp.zeros(1, 8)) * k3_basis.inv()
    )
    if not (
        column_action * invariant_ample == invariant_ample
        and (invariant_ample.T * n_gram * invariant_ample)[0] == 528
        and min(ample_pairings) == 2
        and all(value > 0 for value in ample_pairings)
        and all(sp.Rational(value).q == 1 for value in global_ample)
        and global_ample * k3_action == global_ample
        and (global_ample * k3_gram * global_ample.T)[0] == 528
    ):
        raise AssertionError("the invariant strict nef-chamber interior changed")

    bh_replay = None
    if brandhorst_hofmann_order4 is not None:
        bh_replay = replay_brandhorst_hofmann_order4(
            brandhorst_hofmann_order4,
            canonical_gram,
            canonical_action,
        )

    elliptic_audit = elliptic_normal_form_audit()

    certificate: dict[str, Any] = {
        "schema": "coxeter.roulleau_r14_c2d.special_k3_order4_lift.v1",
        "all_checks_pass": True,
        "source_file_sha256": source_hash,
        "case_id": "R14-UE8A1x4-C2D-H11-19",
        "theorem": (
            "There exists a complex projective K3 surface X with "
            "NS(X)=U+E8(-1)+A1(-1)^4 and a purely non-symplectic "
            "automorphism sigma of order four whose Picard action is the "
            "archived double transposition.  The fixed section of X's complete "
            "27-wall nef chamber is the finite-volume noncompact H11/19 "
            "Coxeter polytope R14-UE8A1x4-C2D-H11-19."
        ),
        "archived_picard_pair": {
            "lattice": "U+E8(-1)+A1(-1)^4",
            "rank_signature_determinant": [14, [1, 13], -16],
            "action_characteristic_polynomial": "(x-1)^12*(x+1)^2",
            "canonical_action": (
                "identity on U+E8(-1), and (12)(34) on the four A1(-1) summands"
            ),
            "canonical_basis_sha256": matrix_digest(canonical_basis),
            "canonical_gram_sha256": matrix_digest(canonical_gram),
            "fixed_lattice": "U+E8(-1)+<-4>+<-4>",
            "fixed_rank_signature_determinant": [12, [1, 11], -16],
            "fixed_smith": [*[1] * 10, 4, 4],
            "anti_invariant_lattice": "<-4>+<-4>",
            "anti_invariant_rootless": True,
        },
        "equivariant_k3_glue": {
            "transcendental_lattice": "U+U+A1(-1)^4",
            "transcendental_rank_signature_determinant": [8, [2, 6], 16],
            "transcendental_action_sha256": matrix_digest(t_action),
            "transcendental_action_square": "-identity",
            "anti_isometry_N_to_T_mod_two": matrix_payload(identification),
            "N_discriminant_action": [
                [int(value) % 2 for value in n_discriminant_action.row(index)]
                for index in range(4)
            ],
            "T_discriminant_action": [
                [int(value) % 2 for value in t_discriminant_action.row(index)]
                for index in range(4)
            ],
            "quadratic_form_table": q_table,
            "glued_basis_sha256": matrix_digest(k3_basis),
            "glued_gram_sha256": matrix_digest(k3_gram),
            "glued_action_sha256": matrix_digest(k3_action),
            "glued_lattice": "even unimodular II_(3,19)",
            "glued_action_characteristic_polynomial": (
                "(x-1)^12*(x+1)^2*(x^2+1)^4"
            ),
            "glued_action_order": 4,
            "square_fixed_lattice_is_exactly_N": True,
            "global_fixed_lattice_is_exactly_archived_fixed_lattice": True,
        },
        "explicit_period_and_picard": {
            "alpha": "positive real fourth root of 2",
            "minimal_polynomial_over_Q_i": "alpha^4-2",
            "minimal_polynomial_irreducible_over_Q_i": True,
            "i_eigenbasis_hermitian_gram": matrix_payload(hermitian_gram),
            "period_coefficients_in_i_eigenbasis": [
                "1",
                "alpha/10",
                "alpha^2/100",
                "alpha^3/1000",
            ],
            "global_i_eigenbasis_sha256": matrix_digest(global_eigenbasis),
            "period_square": "0",
            "period_eigenvalue": "i",
            "period_positive_norm_formula": (
                "4*(1-alpha^2/100-alpha^4/10000-alpha^6/1000000)"
            ),
            "period_positive_norm_strict_lower_bound": str(positive_lower_bound),
            "Q_i_linear_independence_reason": (
                "1,alpha,alpha^2,alpha^3 are independent because alpha^4-2 "
                "is irreducible over Q(i)"
            ),
            "no_nonzero_rational_T_orthogonal": True,
            "picard_lattice_is_exactly_N": True,
            "picard_rank": 14,
        },
        "kahler_and_torelli": {
            "source_chamber": "Roulleau's complete 27-wall K3 nef chamber",
            "all_27_walls_permuted": True,
            "invariant_ample_vector_in_source_basis": [
                int(value) for value in invariant_ample
            ],
            "invariant_ample_vector_in_glued_k3_basis_sha256": matrix_digest(
                global_ample
            ),
            "invariant_ample_square": 528,
            "minimum_wall_pairing": "2",
            "strictly_inside_source_chamber": True,
            "anti_invariant_NS_has_minus_two_vectors": False,
            "strong_torelli_conclusion": (
                "the order-four Hodge isometry preserving this Kahler class "
                "is induced by a genuine K3 automorphism"
            ),
        },
        "hyperbolic_fixed_section": {
            "dimension": 11,
            "facet_count": 19,
            "classification": "finite_volume_noncompact",
            "finite_volume": True,
            "compact": False,
            "elliptic_edge_count": 795,
            "finite_vertex_count": 114,
            "ideal_vertex_count": 13,
            "bad_edge_count": 0,
            "bad_parabolic_neighborhood_count": 0,
            "edge_endpoint_type_histogram": vinberg[
                "edge_endpoint_type_histogram"
            ],
            "finite_special_subset_counts_by_rank": vinberg[
                "finite_special_subset_counts_by_rank"
            ],
        },
        "brandhorst_hofmann_cross_identification": {
            "citation": (
                "Brandhorst--Hofmann, Finite subgroups of automorphisms of K3 surfaces"
            ),
            "doi": BH_DOI,
            "archive_version": BH_VERSION,
            "license": "CC-BY-4.0",
            "member": "K3Groups/results/purely_ns/order04.txt",
            "member_sha256": BH_ORDER4_SHA256,
            "complete_row_count": BH_ORDER4_ROW_COUNT,
            "cohomological_character_r_l_m": [12, 2, 4],
            "all_rows_with_that_character": list(BH_MATCHING_ROWS),
            "unique_row_with_fixed_determinant_minus_16": "0.4.3.4",
            "target_NS_pair_exactly_U_E8_A1x4_with_double_transposition": True,
            "external_full_replay_performed": bh_replay is not None,
            "external_replay": bh_replay,
            "role": (
                "independent exhaustive-classification cross-check; the positive "
                "existence theorem is already proved by the explicit graph glue"
            ),
        },
        "elliptic_algebro_geometric_chart": elliptic_audit,
        "novelty_scope": (
            "a new project-level genuine K3-automorphism realization of the "
            "already certified R14 H11/19 fixed-section polytope; not a new "
            "abstract Coxeter polytope, dimension, or arithmetic class"
        ),
    }
    canonical_payload = json.dumps(
        certificate, sort_keys=True, separators=(",", ":")
    ).encode()
    certificate["certificate_payload_sha256"] = sha256_bytes(canonical_payload)
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--brandhorst-hofmann-order4", type=Path)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate(args.brandhorst_hofmann_order4)
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != rendered:
            raise AssertionError("the archived R14 C2D order-four lift changed")
        print(f"archive OK: {args.output}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(args.output)
    print(certificate["certificate_payload_sha256"])


if __name__ == "__main__":
    main()
