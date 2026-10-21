"""Search golden rank-11 trace lattices realizable by supersingular K3s.

For an odd prime ``p``, Rudakov--Shafarevich characterize the Neron--Severi
lattice of a supersingular K3 surface of Artin invariant ``sigma`` as the
unique even hyperbolic p-elementary lattice of rank 22 and discriminant
``-p^(2*sigma)``.  This script searches ternary golden forms ``T`` for which

    T + E8

+has signature ``(10,1)`` at the identity embedding and is positive definite
at the conjugate embedding.  Restriction of scalars with the same codifferent
normalization used for the Bugaenko/K3 construction then has rank 22.  A
p-elementary trace lattice is therefore an exact supersingular-K3
realizability certificate and its negative golden line supplies an H^10
eigenspace.

The first search family fixes the diagonal ``(-2 phi, 2, 2)`` and enumerates
the three off-diagonal entries in a coefficient box in Z[phi].  Numerical
arithmetic is never used for acceptance.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from typing import Iterable

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.exact.matrix import determinant, inertia
from coxeter_polytopes.exact.trace_lattice import (
    INV_SQRT5,
    integral_gram,
    is_even_integral,
    trace_gram,
)


TRACE_SCALAR = -PHI * INV_SQRT5


@dataclass(frozen=True)
class Candidate:
    a: Golden
    b: Golden
    c: Golden
    x: Golden
    y: Golden
    z: Golden
    field_determinant: Golden
    artin_invariant: int
    prime: int
    type_i: bool
    trace_gram: tuple[tuple[int, ...], ...]
    smith_diagonal: tuple[int, ...]


def _ternary(
    a: Golden,
    b: Golden,
    c: Golden,
    x: Golden,
    y: Golden,
    z: Golden,
) -> list[list[Golden]]:
    return [
        [a, x, y],
        [x, b, z],
        [y, z, c],
    ]


def _ternary_determinant(
    a: Golden,
    b: Golden,
    c: Golden,
    x: Golden,
    y: Golden,
    z: Golden,
) -> Golden:
    return a * b * c + 2 * x * y * z - a * z * z - b * y * y - c * x * x


def _coefficient_box(bound: int) -> tuple[Golden, ...]:
    return tuple(
        Golden(a, b)
        for a in range(-bound, bound + 1)
        for b in range(-bound, bound + 1)
    )


def _smith_diagonal(gram: list[list[int]]) -> tuple[int, ...]:
    smith = smith_normal_form(Matrix(gram), domain=ZZ)
    return tuple(abs(int(smith[i, i])) for i in range(smith.rows))


def _is_type_i(gram: list[list[int]]) -> bool:
    """Test q(x) in Z/2Z on the discriminant group of a 2-elementary form."""

    inverse = Matrix(gram).inv()
    return all(inverse[i, i].q == 1 for i in range(inverse.rows))


def _artin_invariant_from_norm(norm: int, prime: int) -> int | None:
    if norm >= 0:
        return None
    absolute = -norm
    exponent = 0
    while absolute % prime == 0:
        absolute //= prime
        exponent += 1
    if absolute != 1 or exponent % 2 or not (2 <= exponent <= 20):
        return None
    return exponent // 2


def _accept(
    a: Golden,
    b: Golden,
    c: Golden,
    x: Golden,
    y: Golden,
    z: Golden,
    prime: int,
) -> Candidate | None:
    det = _ternary_determinant(a, b, c, x, y, z)
    form = _ternary(a, b, c, x, y, z)
    if det != determinant(form):
        raise ArithmeticError("closed determinant formula changed")
    norm = det.norm()
    if norm.denominator != 1:
        return None
    sigma = _artin_invariant_from_norm(norm.numerator, prime)
    if sigma is None:
        return None

    identity_signature = inertia(form)
    conjugate_signature = inertia(
        [[entry.conjugate() for entry in row] for row in form]
    )
    if (
        identity_signature.positive,
        identity_signature.negative,
        identity_signature.zero,
    ) != (2, 1, 0):
        return None
    if (
        conjugate_signature.positive,
        conjugate_signature.negative,
        conjugate_signature.zero,
    ) != (3, 0, 0):
        return None

    traced = trace_gram(form, TRACE_SCALAR)
    if not is_even_integral(traced):
        return None
    integer = integral_gram(traced)
    trace_det = determinant(integer)
    if trace_det.b or trace_det.a.denominator != 1 or trace_det.a.numerator != norm:
        raise ArithmeticError("trace determinant no longer equals the field norm")
    smith = _smith_diagonal(integer)
    if any(entry not in (1, prime) for entry in smith):
        return None
    if smith.count(prime) != 2 * sigma:
        return None
    type_i = _is_type_i(integer) if prime == 2 else False
    if prime == 2 and not type_i:
        return None
    trace_signature = inertia(integer)
    if (
        trace_signature.positive,
        trace_signature.negative,
        trace_signature.zero,
    ) != (1, 5, 0):
        raise ArithmeticError("accepted ternary trace form has wrong signature")
    return Candidate(
        a=a,
        b=b,
        c=c,
        x=x,
        y=y,
        z=z,
        field_determinant=det,
        artin_invariant=sigma,
        prime=prime,
        type_i=type_i,
        trace_gram=tuple(tuple(row) for row in integer),
        smith_diagonal=smith,
    )


def search(bound: int, prime: int = 5) -> Iterable[Candidate]:
    entries = _coefficient_box(bound)
    for x in entries:
        for y in entries:
            for z in entries:
                candidate = _accept(
                    -2 * PHI, Golden(2), Golden(2), x, y, z, prime
                )
                if candidate is not None:
                    yield candidate


def _even_diagonals(bound: int, sign_pair: tuple[int, int]) -> tuple[Golden, ...]:
    return tuple(
        2 * Golden(a, b)
        for a in range(-bound, bound + 1)
        for b in range(-bound, bound + 1)
        if Golden(a, b).sign() == sign_pair[0]
        and Golden(a, b).sign(1) == sign_pair[1]
    )


def search_solved_diagonal(
    diagonal_bound: int,
    off_diagonal_bound: int,
    derived_bound: int,
    unit_bound: int,
    prime: int = 5,
    sigma: int = 1,
) -> Iterable[Candidate]:
    """Solve the determinant equation for the third diagonal entry.

    Every target determinant of the form ``-p^sigma*phi^(2k+1)`` has norm
    ``-p^(2 sigma)`` and the required signs.  Solving linearly for ``c`` is
    much faster than scanning three diagonal entries independently.
    """

    negative_positive = _even_diagonals(diagonal_bound, (-1, 1))
    positive_positive = _even_diagonals(diagonal_bound, (1, 1))
    entries = _coefficient_box(off_diagonal_bound)
    targets = tuple(
        -(prime**sigma) * (PHI ** (2 * exponent + 1))
        for exponent in range(-unit_bound, unit_bound + 1)
    )
    seen: set[tuple[Golden, ...]] = set()
    for a in negative_positive:
        for b in positive_positive:
            for x in entries:
                coefficient = a * b - x * x
                if not coefficient:
                    continue
                coefficient_inverse = coefficient.inverse()
                scaled_targets = tuple(
                    target * coefficient_inverse for target in targets
                )
                for y in entries:
                    for z in entries:
                        constant = 2 * x * y * z - a * z * z - b * y * y
                        base = -constant * coefficient_inverse
                        for scaled_target in scaled_targets:
                            c = scaled_target + base
                            if (
                                not c.is_integral
                                or c.a.numerator % 2
                                or c.b.numerator % 2
                                or c.sign() <= 0
                                or c.sign(1) <= 0
                                or abs(c.a.numerator) > derived_bound
                                or abs(c.b.numerator) > derived_bound
                            ):
                                continue
                            key = (a, b, c, x, y, z)
                            if key in seen:
                                continue
                            seen.add(key)
                            candidate = _accept(a, b, c, x, y, z, prime)
                            if candidate is not None:
                                yield candidate


def _encode(candidate: Candidate) -> dict[str, object]:
    return {
        "diagonal": {
            "a": candidate.a.to_json(),
            "b": candidate.b.to_json(),
            "c": candidate.c.to_json(),
        },
        "off_diagonal": {
            "x": candidate.x.to_json(),
            "y": candidate.y.to_json(),
            "z": candidate.z.to_json(),
        },
        "ternary_field_gram": [
            [entry.to_json() for entry in row]
            for row in _ternary(
                candidate.a,
                candidate.b,
                candidate.c,
                candidate.x,
                candidate.y,
                candidate.z,
            )
        ],
        "field_determinant": candidate.field_determinant.to_json(),
        "field_determinant_norm": str(candidate.field_determinant.norm()),
        "field_signature_identity": {"positive": 2, "negative": 1, "zero": 0},
        "field_signature_conjugate": {"positive": 3, "negative": 0, "zero": 0},
        "full_field_lattice": "ternary + E8 over Z[phi]",
        "full_field_rank": 11,
        "hyperbolic_eigenspace_dimension": 10,
        "trace_rank": 22,
        "trace_signature": {"positive": 1, "negative": 21, "zero": 0},
        "trace_determinant": str(candidate.field_determinant.norm()),
        "trace_even_integral": True,
        "trace_smith_diagonal_ternary": list(candidate.smith_diagonal),
        "trace_smith_diagonal_full": [1] * 16 + list(candidate.smith_diagonal),
        "prime": candidate.prime,
        "artin_invariant": candidate.artin_invariant,
        "trace_2_elementary_type_I": candidate.type_i if candidate.prime == 2 else None,
        "supersingular_K3_NS_realizable": True,
        "realizability_theorem": (
            "Rudakov-Shafarevich uniqueness of even hyperbolic p-elementary "
            "rank-22 lattices of discriminant -p^(2 sigma), p odd"
        ),
        "scope": "lattice-and-K3-realizability certificate; chamber compactness not tested",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bound", type=int, default=2)
    parser.add_argument("--solve-diagonal", action="store_true")
    parser.add_argument("--diagonal-bound", type=int, default=2)
    parser.add_argument("--derived-bound", type=int, default=40)
    parser.add_argument("--unit-bound", type=int, default=4)
    parser.add_argument("--sigma", type=int, default=1)
    parser.add_argument("--prime", type=int, default=5)
    parser.add_argument("--limit", type=int, default=20)
    args = parser.parse_args()
    if args.bound < 0 or args.limit < 1:
        parser.error("bound must be nonnegative and limit positive")
    records = []
    candidates = (
        search_solved_diagonal(
            args.diagonal_bound,
            args.bound,
            args.derived_bound,
            args.unit_bound,
            args.prime,
            args.sigma,
        )
        if args.solve_diagonal
        else search(args.bound, args.prime)
    )
    for candidate in candidates:
        records.append(_encode(candidate))
        if len(records) >= args.limit:
            break
    print(
        json.dumps(
            {
                "search_family": (
                    "determinant-solved even diagonal ternary golden forms"
                    if args.solve_diagonal
                    else "diag(-2 phi,2,2), boxed golden off-diagonals"
                ),
                "coefficient_bound": args.bound,
                "diagonal_bound": args.diagonal_bound,
                "derived_diagonal_bound": args.derived_bound,
                "unit_bound": args.unit_bound,
                "target_artin_invariant": args.sigma,
                "prime": args.prime,
                "accepted": len(records),
                "candidates": records,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
