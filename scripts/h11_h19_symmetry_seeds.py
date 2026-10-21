"""Exact infinite-order arithmetic symmetry seeds for H11 and H19.

These binary isometries are inputs for a chamber-stabilizer search.  They are
not, by themselves, certificates that either norm-two reflection subgroup has
infinite index: that requires an independently certified same-chamber witness.
"""

from __future__ import annotations

import json
from typing import Sequence

from coxeter_polytopes.exact.golden import Golden, ONE, PHI, ZERO
from coxeter_polytopes.exact.matrix import determinant


Matrix = Sequence[Sequence[Golden]]


def _matrix_product(left: Matrix, right: Matrix) -> list[list[Golden]]:
    return [
        [
            sum(
                (left[i][k] * right[k][j] for k in range(len(right))),
                ZERO,
            )
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def _transpose(matrix: Matrix) -> list[list[Golden]]:
    return [list(column) for column in zip(*matrix)]


def _encode_matrix(matrix: Matrix) -> list[list[str]]:
    return [[str(entry) for entry in row] for row in matrix]


def h11_binary_gram() -> list[list[Golden]]:
    return [[-2 * PHI, ONE], [ONE, Golden(2)]]


def h19_binary_gram() -> list[list[Golden]]:
    return [[-2 * PHI, PHI], [PHI, Golden(2)]]


def h11_seed() -> list[list[Golden]]:
    """A determinant-one loxodromic isometry of the H11 binary block."""

    return [[ONE, ONE - PHI], [-ONE, PHI]]


def h19_seed() -> list[list[Golden]]:
    """A determinant-one loxodromic isometry of the H19 binary block."""

    return [[ONE, -ONE], [-PHI, ONE + PHI]]


def _seed_record(
    name: str,
    gram: Matrix,
    seed: Matrix,
    pell_discriminant: Golden,
    pell_u: Golden,
    pell_t: Golden,
) -> dict[str, object]:
    pullback = _matrix_product(_transpose(seed), _matrix_product(gram, seed))
    trace = seed[0][0] + seed[1][1]
    pell_value = pell_t**2 - pell_discriminant * pell_u**2
    return {
        "lattice": name,
        "binary_gram": _encode_matrix(gram),
        "seed": _encode_matrix(seed),
        "preserves_binary_gram": pullback == gram,
        "determinant": str(determinant(seed)),
        "trace": str(trace),
        "trace_minus_two_identity_sign": (trace - 2).sign(),
        "pell_discriminant": str(pell_discriminant),
        "pell_u": str(pell_u),
        "pell_t": str(pell_t),
        "pell_equation_value": str(pell_value),
        "infinite_projective_order": (
            pullback == gram
            and determinant(seed) == ONE
            and (trace - 2).sign() > 0
        ),
        "same_chamber_certified": False,
    }


def analyze() -> dict[str, object]:
    """Return exact checks for the two search seeds.

    Trace greater than two at the Lorentzian embedding makes each
    determinant-one binary isometry loxodromic, hence of infinite projective
    order.  Extending it by the identity on E8 gives an integral isometry of
    the full rank-ten field lattice.
    """

    h11 = _seed_record(
        "H11",
        h11_binary_gram(),
        h11_seed(),
        ONE + 4 * PHI,
        PHI - ONE,
        ONE + PHI,
    )
    h19 = _seed_record(
        "H19",
        h19_binary_gram(),
        h19_seed(),
        ONE + 5 * PHI,
        ONE,
        Golden(2) + PHI,
    )
    return {
        "status": "search seeds only; no same-chamber certificate is claimed",
        "seeds": [h11, h19],
    }


def main() -> None:
    print(json.dumps(analyze(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
